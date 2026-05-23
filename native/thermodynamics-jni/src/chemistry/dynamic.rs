use std::collections::BTreeMap;
use std::collections::BTreeSet;

use super::error::{ChemistryError, ChemistryResult};
use super::frowns::{parse_frowns, write_frowns};
use super::functional_group::{FunctionalGroup, FunctionalGroupType};
use super::molecule::MolecularStructure;
use super::organic;
use super::registry::{ChemistryRegistry, ChemistryRegistryBuilder};
use super::substance::{Substance, SubstanceId, SubstanceTagId};

const DEFAULT_DYNAMIC_DENSITY: f64 = 1000.0;
const DEFAULT_DYNAMIC_HEAT_CAPACITY: f64 = 100.0;
const DEFAULT_DYNAMIC_LATENT_HEAT: f64 = 20_000.0;
const DEFAULT_DYNAMIC_COLOR: u32 = 0x20FF_FFFF;
const MAX_DYNAMIC_ATOMS: usize = 100;
const MAX_DYNAMIC_WORK_ITEMS: usize = 1_000_000;
const MAX_DYNAMIC_QUEUE_ITEMS: usize = 100_000;

#[derive(Debug, Clone)]
pub struct DynamicChemistryRegistry {
    registry: ChemistryRegistry,
    canonical_to_id: BTreeMap<String, SubstanceId>,
    processed_generation_keys: BTreeSet<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DynamicGenerationReport {
    pub iterations: usize,
    pub added_substances: usize,
    pub added_reactions: usize,
    pub processed_work_items: usize,
    pub skipped_duplicates: usize,
    pub remaining_queue: usize,
    pub reached_fixed_point: bool,
    pub generator_errors: Vec<String>,
}

impl DynamicChemistryRegistry {
    pub fn from_destroy_catalog() -> ChemistryResult<Self> {
        Self::from_registry(super::destroy_registry_builder()?.build()?)
    }

    pub fn from_registry(registry: ChemistryRegistry) -> ChemistryResult<Self> {
        let mut result = Self {
            registry,
            canonical_to_id: BTreeMap::new(),
            processed_generation_keys: BTreeSet::new(),
        };
        result.rebuild_canonical_index()?;
        Ok(result)
    }

    pub fn registry(&self) -> &ChemistryRegistry {
        &self.registry
    }

    pub fn resolve_frowns(&mut self, code: &str) -> ChemistryResult<SubstanceId> {
        self.resolve_structure(parse_frowns(code)?)
    }

    pub fn resolve_structure(
        &mut self,
        structure: MolecularStructure,
    ) -> ChemistryResult<SubstanceId> {
        let canonical = write_frowns(&structure)?;
        if let Some(id) = self.canonical_to_id.get(&canonical) {
            return Ok(id.clone());
        }
        let substance = build_dynamic_substance(canonical.clone(), structure)?;
        let id = substance.id.clone();
        self.registry = ChemistryRegistryBuilder::from_registry(&self.registry)
            .substance(substance)
            .build()?;
        self.canonical_to_id.insert(canonical, id.clone());
        Ok(id)
    }

    pub fn generate_reactions_for(
        &mut self,
        substance_id: &SubstanceId,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        self.registry.substance(substance_id)?;
        let mut seeds = BTreeSet::from([substance_id.clone()]);
        self.generate_reactions_from_scope(&mut seeds, Some(max_iterations))
    }

    pub fn generate_reactions_for_to_fixed_point(
        &mut self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<DynamicGenerationReport> {
        self.registry.substance(substance_id)?;
        let mut seeds = BTreeSet::from([substance_id.clone()]);
        self.generate_reactions_from_scope(&mut seeds, None)
    }

    pub fn generate_reactions_for_substances(
        &mut self,
        substance_ids: impl IntoIterator<Item = SubstanceId>,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        let mut seeds = self.validated_substance_set(substance_ids)?;
        self.generate_reactions_from_scope(&mut seeds, Some(max_iterations))
    }

    pub fn generate_reactions_for_substances_to_fixed_point(
        &mut self,
        substance_ids: impl IntoIterator<Item = SubstanceId>,
    ) -> ChemistryResult<DynamicGenerationReport> {
        let mut seeds = self.validated_substance_set(substance_ids)?;
        self.generate_reactions_from_scope(&mut seeds, None)
    }

    pub fn generate_reactions(
        &mut self,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        let mut seeds = self
            .registry
            .substances()
            .map(|substance| substance.id.clone())
            .collect::<BTreeSet<_>>();
        self.generate_reactions_from_scope(&mut seeds, Some(max_iterations))
    }

    pub fn generate_to_fixed_point(&mut self) -> ChemistryResult<DynamicGenerationReport> {
        let mut seeds = self
            .registry
            .substances()
            .map(|substance| substance.id.clone())
            .collect::<BTreeSet<_>>();
        self.generate_reactions_from_scope(&mut seeds, None)
    }

    fn generate_reactions_from_scope(
        &mut self,
        seeds: &mut BTreeSet<SubstanceId>,
        max_iterations: Option<usize>,
    ) -> ChemistryResult<DynamicGenerationReport> {
        if max_iterations == Some(0) {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<dynamic-generation>".to_string(),
                reason: "max_iterations must be greater than zero".to_string(),
            });
        }

        let mut added_substances = 0usize;
        let mut added_reactions = 0usize;
        let mut processed_work_items = 0usize;
        let mut skipped_duplicates = 0usize;
        let mut queue = seeds.clone();
        let mut scope = seeds.clone();
        let mut iteration = 0usize;
        loop {
            if max_iterations.is_some_and(|max| iteration >= max) {
                return Ok(DynamicGenerationReport {
                    iterations: iteration,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: queue.len(),
                    reached_fixed_point: false,
                    generator_errors: Vec::new(),
                });
            }
            if queue.is_empty() {
                return Ok(DynamicGenerationReport {
                    iterations: iteration,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: 0,
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }
            let current_seeds = queue.clone();
            queue.clear();
            let mut unprocessed_seeds = BTreeSet::new();
            for seed in &current_seeds {
                let keys = self.generation_keys_for_substance(seed)?;
                if keys.is_empty() {
                    skipped_duplicates += 1;
                    continue;
                }
                let has_new_key = keys
                    .iter()
                    .any(|key| !self.processed_generation_keys.contains(key));
                if !has_new_key {
                    skipped_duplicates += 1;
                    continue;
                }
                self.processed_generation_keys.extend(keys);
                unprocessed_seeds.insert(seed.clone());
                processed_work_items += 1;
                if processed_work_items > MAX_DYNAMIC_WORK_ITEMS {
                    return Err(ChemistryError::GenerationInvariantViolation {
                        generator: "<dynamic-generation>".to_string(),
                        substance_id: seed.to_string(),
                        reason: format!(
                            "processed more than {MAX_DYNAMIC_WORK_ITEMS} work items without reaching a fixed point"
                        ),
                    });
                }
            }
            if unprocessed_seeds.is_empty() {
                return Ok(DynamicGenerationReport {
                    iterations: iteration + 1,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: 0,
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }

            let generated = organic::generate_organic_reactions_for_scope(
                &self.registry,
                &unprocessed_seeds,
                &scope,
            )?;
            let mut builder = ChemistryRegistryBuilder::from_registry(&self.registry);
            let mut changed = false;

            for substance in generated.substances {
                if self.registry.substance(&substance.id).is_ok() {
                    skipped_duplicates += 1;
                    continue;
                }
                let canonical = substance
                    .molecular_structure
                    .as_ref()
                    .map(write_frowns)
                    .transpose()?
                    .ok_or_else(|| ChemistryError::InvalidSubstance {
                        substance_id: substance.id.to_string(),
                        reason: "generated dynamic substance has no structure".to_string(),
                    })?;
                if self.canonical_to_id.contains_key(&canonical) {
                    skipped_duplicates += 1;
                    continue;
                }
                self.canonical_to_id.insert(canonical, substance.id.clone());
                queue.insert(substance.id.clone());
                scope.insert(substance.id.clone());
                if queue.len() > MAX_DYNAMIC_QUEUE_ITEMS {
                    return Err(ChemistryError::GenerationInvariantViolation {
                        generator: "<dynamic-generation>".to_string(),
                        substance_id: substance.id.to_string(),
                        reason: format!(
                            "dynamic generation queue exceeded {MAX_DYNAMIC_QUEUE_ITEMS} substances before reaching a fixed point"
                        ),
                    });
                }
                builder = builder.substance(substance);
                added_substances += 1;
                changed = true;
            }

            for reaction in generated.reactions {
                if self.registry.reaction(&reaction.id).is_ok() {
                    skipped_duplicates += 1;
                    continue;
                }
                builder = builder.reaction(reaction);
                added_reactions += 1;
                changed = true;
            }

            if !changed {
                return Ok(DynamicGenerationReport {
                    iterations: iteration + 1,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: queue.len(),
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }
            self.registry = builder.build()?;
            self.rebuild_canonical_index()?;
            *seeds = queue.clone();
            if queue.is_empty() {
                return Ok(DynamicGenerationReport {
                    iterations: iteration + 1,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: 0,
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }
            iteration += 1;
        }
    }

    fn rebuild_canonical_index(&mut self) -> ChemistryResult<()> {
        self.canonical_to_id.clear();
        for substance in self.registry.substances() {
            let Some(structure) = &substance.molecular_structure else {
                continue;
            };
            self.canonical_to_id
                .entry(write_frowns(structure)?)
                .or_insert_with(|| substance.id.clone());
        }
        Ok(())
    }

    fn validated_substance_set(
        &self,
        substance_ids: impl IntoIterator<Item = SubstanceId>,
    ) -> ChemistryResult<BTreeSet<SubstanceId>> {
        let mut result = BTreeSet::new();
        for substance_id in substance_ids {
            self.registry.substance(&substance_id)?;
            result.insert(substance_id);
        }
        Ok(result)
    }

    fn generation_keys_for_substance(
        &self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<BTreeSet<String>> {
        let substance = self.registry.substance(substance_id)?;
        let Some(structure) = &substance.molecular_structure else {
            return Ok(BTreeSet::new());
        };
        let canonical = write_frowns(structure)?;
        let mut keys = BTreeSet::new();
        for group in &substance.functional_groups {
            for generator in generators_for_group(&group.group_type) {
                keys.insert(format!(
                    "{canonical}|{generator}|{}",
                    functional_group_key(structure, group)
                ));
            }
        }
        Ok(keys)
    }
}

fn generators_for_group(group_type: &FunctionalGroupType) -> &'static [&'static str] {
    match group_type {
        FunctionalGroupType::Halide => &[
            "halide_hydroxide_substitution",
            "halide_ammonia_substitution",
            "halide_cyanide_substitution",
            "halide_amine_substitution",
        ],
        FunctionalGroupType::Alcohol => &[
            "alcohol_oxidation",
            "alcohol_dehydration",
            "thionyl_chloride_substitution",
            "carboxylic_acid_esterification",
            "acyl_chloride_esterification",
        ],
        FunctionalGroupType::Alkoxide => &["alkoxide_protonation"],
        FunctionalGroupType::Nitrile => &["nitrile_hydrolysis", "nitrile_hydrogenation"],
        FunctionalGroupType::Nitro => &["nitro_hydrogenation"],
        FunctionalGroupType::AcylChloride => {
            &["acyl_chloride_hydrolysis", "acyl_chloride_esterification"]
        }
        FunctionalGroupType::CarboxylicAcid => {
            &["acyl_chloride_formation", "carboxylic_acid_esterification"]
        }
        FunctionalGroupType::Carbonyl => &[
            "aldehyde_oxidation",
            "cyanide_nucleophilic_addition",
            "wolff_kishner_reduction",
        ],
        FunctionalGroupType::UnsubstitutedAmide => &["amide_hydrolysis"],
        FunctionalGroupType::PrimaryAmine => &["amine_phosgenation"],
        FunctionalGroupType::NonTertiaryAmine => {
            &["cyanamide_addition", "halide_amine_substitution"]
        }
        FunctionalGroupType::Isocyanate => &["isocyanate_hydrolysis"],
        FunctionalGroupType::Borane => &["borane_oxidation"],
        FunctionalGroupType::BorateEster => &["borate_ester_hydrolysis"],
        FunctionalGroupType::Alkene => &[
            "alkene_chlorination",
            "alkene_chlorohydrination",
            "alkene_hydrolysis",
            "alkene_hydroboration_with_borane",
            "alkene_hydrochlorination",
            "alkene_hydrogenation",
            "alkene_hydroiodination",
            "alkene_iodination",
        ],
        FunctionalGroupType::Alkyne => &[
            "alkyne_chlorination",
            "alkyne_chlorohydrination",
            "alkyne_hydrolysis",
            "alkyne_hydroboration_with_borane",
            "alkyne_hydrochlorination",
            "alkyne_hydrogenation",
            "alkyne_hydroiodination",
            "alkyne_iodination",
        ],
        _ => &[],
    }
}

fn functional_group_key(structure: &MolecularStructure, group: &FunctionalGroup) -> String {
    let mut atoms = group
        .atoms
        .iter()
        .map(|index| {
            let atom = &structure.atoms[*index];
            format!("{}:{}:{:.3}", index, atom.element, atom.charge)
        })
        .collect::<Vec<_>>();
    atoms.sort();
    format!("{:?}|{}", group.group_type, atoms.join(","))
}

fn build_dynamic_substance(
    canonical_frowns: String,
    structure: MolecularStructure,
) -> ChemistryResult<Substance> {
    if structure.atom_count() >= MAX_DYNAMIC_ATOMS {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: canonical_frowns,
            reason: "dynamic substance has too many atoms".to_string(),
        });
    }
    let summary = structure.summary()?;
    let tags = if structure.atoms.iter().any(|atom| atom.element == "R") {
        vec![SubstanceTagId::from("destroy:hypothetical")]
    } else {
        Vec::new()
    };
    Ok(Substance::new(
        SubstanceId::new(canonical_frowns.clone())?,
        summary.charge,
        summary.molar_mass_grams,
        DEFAULT_DYNAMIC_DENSITY,
        if summary.charge == 0 {
            estimate_dynamic_boiling_point(summary.molar_mass_grams)
        } else {
            f64::MAX
        },
        DEFAULT_DYNAMIC_HEAT_CAPACITY,
        DEFAULT_DYNAMIC_LATENT_HEAT,
    )
    .with_catalog_metadata(Some(canonical_frowns), None, DEFAULT_DYNAMIC_COLOR, tags)
    .with_molecular_structure(structure))
}

fn estimate_dynamic_boiling_point(molar_mass_grams: f64) -> f64 {
    2.042_598_921_281_41 * molar_mass_grams + 178.176_866_128_713
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::mixture::Mixture;

    #[test]
    fn resolves_known_substance_by_frowns() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acetone = registry.resolve_frowns("CC(=O)C").unwrap();
        assert_eq!(acetone.as_str(), "destroy:acetone");
    }

    #[test]
    fn creates_stable_dynamic_substance_without_duplicates() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let first = registry.resolve_frowns("CCCCCCCC").unwrap();
        let second = registry.resolve_frowns("CCCCCCCC").unwrap();
        assert_eq!(first, second);
        assert!(first.as_str().starts_with("destroy:linear:"));
    }

    #[test]
    fn hypothetical_dynamic_substance_cannot_enter_mixture() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let id = registry.resolve_frowns("R1C").unwrap();
        let substance = registry.registry().substance(&id).unwrap();
        assert!(substance
            .tags
            .iter()
            .any(|tag| tag.as_str() == "destroy:hypothetical"));

        let mut mixture = Mixture::empty();
        assert!(mixture.add_substance(registry.registry(), id, 1.0).is_err());
    }

    #[test]
    fn dynamic_products_are_available_for_later_generation() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let alkene = registry.resolve_frowns("CCCC=C").unwrap();
        let report = registry.generate_reactions_for(&alkene, 2).unwrap();
        assert!(report.added_substances > 0);
        assert!(report.added_reactions > 0);
        assert!(report.processed_work_items > 0);
        assert!(registry.registry().reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("alkene_hydrogenation/destroy_linear_CCCC_C/")
                || reaction.id.as_str().starts_with("alkene_hydrogenation/")
        }));
    }

    #[test]
    fn repeated_dynamic_generation_reaches_fixed_point_without_duplicates() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let alkene = registry.resolve_frowns("CCCC=C").unwrap();
        let first = registry.generate_reactions_for(&alkene, 1).unwrap();
        let second = registry.generate_reactions_for(&alkene, 1).unwrap();

        assert!(first.added_substances > 0);
        assert!(second.reached_fixed_point);
        assert_eq!(second.added_substances, 0);
        assert_eq!(second.added_reactions, 0);
        assert!(second.skipped_duplicates > 0 || second.processed_work_items > 0);
    }

    #[test]
    fn dynamic_generation_can_run_to_fixed_point_without_iteration_limit() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let methane = registry.resolve_frowns("C").unwrap();
        let first = registry
            .generate_reactions_for_to_fixed_point(&methane)
            .unwrap();
        let second = registry
            .generate_reactions_for_to_fixed_point(&methane)
            .unwrap();

        assert!(first.reached_fixed_point);
        assert_eq!(first.remaining_queue, 0);
        assert!(first.generator_errors.is_empty());
        assert!(second.reached_fixed_point);
        assert_eq!(second.added_substances, 0);
        assert_eq!(second.added_reactions, 0);
        assert_eq!(second.processed_work_items, 0);
    }

    #[test]
    fn pair_generators_use_generation_scope_not_whole_registry() {
        let mut acid_only = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acetic_acid = acid_only.resolve_frowns("CC(=O)O").unwrap();
        acid_only.generate_reactions_for(&acetic_acid, 1).unwrap();
        assert!(!acid_only.registry().reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("carboxylic_acid_esterification/destroy_acetic_acid/destroy_ethanol/")
        }));

        let mut acid_and_alcohol = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acetic_acid = acid_and_alcohol.resolve_frowns("CC(=O)O").unwrap();
        let ethanol = acid_and_alcohol.resolve_frowns("CCO").unwrap();
        acid_and_alcohol
            .generate_reactions_for_substances([acetic_acid, ethanol], 1)
            .unwrap();
        assert!(acid_and_alcohol.registry().reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("carboxylic_acid_esterification/destroy_acetic_acid/destroy_ethanol/")
        }));
    }
}
