use std::collections::BTreeMap;
use std::collections::BTreeSet;

use super::error::{ChemistryError, ChemistryResult};
use super::frowns::{parse_frowns, write_frowns};
use super::functional_group::{FunctionalGroup, FunctionalGroupType};
use super::molecule::MolecularStructure;
use super::organic;
use super::reaction::{Reaction, ReactionId};
use super::registry::ChemistryRegistry;
use super::substance::{Substance, SubstanceId, SubstanceTagId};

const DEFAULT_DYNAMIC_DENSITY: f64 = 1000.0;
const DEFAULT_DYNAMIC_HEAT_CAPACITY: f64 = 100.0;
const DEFAULT_DYNAMIC_LATENT_HEAT: f64 = 20_000.0;
const DEFAULT_DYNAMIC_COLOR: u32 = 0x20FF_FFFF;
const MAX_DYNAMIC_ATOMS: usize = 100;
const MAX_DYNAMIC_WORK_ITEMS: usize = 1_000_000;
const MAX_DYNAMIC_QUEUE_ITEMS: usize = 100_000;
const MASS_TOLERANCE_GRAMS_PER_MOL: f64 = 1.0e-6;

#[derive(Debug, Clone)]
pub struct DynamicChemistryRegistry {
    static_registry: ChemistryRegistry,
    dynamic_substances: BTreeMap<SubstanceId, Substance>,
    dynamic_reactions: BTreeMap<ReactionId, Reaction>,
    dynamic_reaction_index_by_substance: BTreeMap<SubstanceId, BTreeSet<ReactionId>>,
    dynamic_unindexed_reaction_ids: BTreeSet<ReactionId>,
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
            static_registry: registry,
            dynamic_substances: BTreeMap::new(),
            dynamic_reactions: BTreeMap::new(),
            dynamic_reaction_index_by_substance: BTreeMap::new(),
            dynamic_unindexed_reaction_ids: BTreeSet::new(),
            canonical_to_id: BTreeMap::new(),
            processed_generation_keys: BTreeSet::new(),
        };
        result.rebuild_canonical_index()?;
        Ok(result)
    }

    pub fn static_registry(&self) -> &ChemistryRegistry {
        &self.static_registry
    }

    pub fn substance(&self, id: &SubstanceId) -> ChemistryResult<&Substance> {
        self.dynamic_substances
            .get(id)
            .or_else(|| self.static_registry.substance(id).ok())
            .ok_or_else(|| ChemistryError::InvalidMixtureState(format!("unknown substance '{id}'")))
    }

    pub fn reaction(&self, id: &ReactionId) -> ChemistryResult<&Reaction> {
        self.dynamic_reactions
            .get(id)
            .or_else(|| self.static_registry.reaction(id).ok())
            .ok_or_else(|| ChemistryError::UnknownReaction(id.to_string()))
    }

    pub fn substances(&self) -> impl Iterator<Item = &Substance> {
        self.static_registry
            .substances()
            .chain(self.dynamic_substances.values())
    }

    pub fn reactions(&self) -> impl Iterator<Item = &Reaction> {
        self.static_registry
            .reactions()
            .chain(self.dynamic_reactions.values())
    }

    pub fn validate_substance_can_enter_mixture(
        &self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<()> {
        let substance = self.substance(substance_id)?;
        if substance
            .tags
            .iter()
            .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
        {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "hypothetical substance '{substance_id}' cannot be added to a mixture"
            )));
        }
        Ok(())
    }

    pub fn reaction_candidates_for_substances<'registry, 'substances, I>(
        &'registry self,
        substances: I,
    ) -> Vec<&'registry Reaction>
    where
        I: IntoIterator<Item = &'substances SubstanceId>,
    {
        let substance_ids = substances.into_iter().collect::<Vec<_>>();
        let mut dynamic_reaction_ids = self.dynamic_unindexed_reaction_ids.clone();
        for substance_id in &substance_ids {
            if let Some(indexed_reactions) = self.dynamic_reaction_index_by_substance.get(substance_id)
            {
                dynamic_reaction_ids.extend(indexed_reactions.iter().cloned());
            }
        }
        let mut result = self
            .static_registry
            .reaction_candidates_for_substances(substance_ids)
            .into_iter()
            .collect::<Vec<_>>();
        result.extend(
            dynamic_reaction_ids
                .into_iter()
                .filter_map(|reaction_id| self.dynamic_reactions.get(&reaction_id)),
        );
        result
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
        self.add_dynamic_substance(substance)?;
        self.canonical_to_id.insert(canonical, id.clone());
        Ok(id)
    }

    pub fn generate_reactions_for(
        &mut self,
        substance_id: &SubstanceId,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        self.substance(substance_id)?;
        let mut seeds = BTreeSet::from([substance_id.clone()]);
        self.generate_reactions_from_scope(&mut seeds, Some(max_iterations))
    }

    pub fn generate_reactions_for_to_fixed_point(
        &mut self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<DynamicGenerationReport> {
        self.substance(substance_id)?;
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
            .substances()
            .map(|substance| substance.id.clone())
            .collect::<BTreeSet<_>>();
        self.generate_reactions_from_scope(&mut seeds, Some(max_iterations))
    }

    pub fn generate_to_fixed_point(&mut self) -> ChemistryResult<DynamicGenerationReport> {
        let mut seeds = self
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

            let available_substances = self.substances().cloned().collect::<Vec<_>>();
            let generated = organic::generate_organic_reactions_for_substances(
                &available_substances,
                &unprocessed_seeds,
                &scope,
            )?;
            let mut changed = false;

            for substance in generated.substances {
                if self.substance(&substance.id).is_ok() {
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
                self.add_dynamic_substance(substance)?;
                added_substances += 1;
                changed = true;
            }

            for reaction in generated.reactions {
                if self.reaction(&reaction.id).is_ok() {
                    skipped_duplicates += 1;
                    continue;
                }
                self.add_dynamic_reaction(reaction)?;
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
        let canonical_entries = self
            .substances()
            .filter_map(|substance| {
                substance
                    .molecular_structure
                    .as_ref()
                    .map(|structure| (substance.id.clone(), structure))
            })
            .map(|(id, structure)| write_frowns(structure).map(|canonical| (canonical, id)))
            .collect::<ChemistryResult<Vec<_>>>()?;
        for (canonical, id) in canonical_entries {
            self.canonical_to_id.entry(canonical).or_insert(id);
        }
        Ok(())
    }

    fn add_dynamic_substance(&mut self, substance: Substance) -> ChemistryResult<()> {
        substance.validate()?;
        if self.static_registry.substance(&substance.id).is_ok()
            || self.dynamic_substances.contains_key(&substance.id)
        {
            return Err(ChemistryError::DuplicateSubstance(substance.id.to_string()));
        }
        for tag in &substance.tags {
            if !self.static_registry.has_substance_tag(tag) {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: substance.id.to_string(),
                    reason: format!("unknown substance tag '{tag}'"),
                });
            }
        }
        self.dynamic_substances.insert(substance.id.clone(), substance);
        Ok(())
    }

    fn add_dynamic_reaction(&mut self, reaction: Reaction) -> ChemistryResult<()> {
        if self.static_registry.reaction(&reaction.id).is_ok()
            || self.dynamic_reactions.contains_key(&reaction.id)
        {
            return Err(ChemistryError::DuplicateReaction(reaction.id.to_string()));
        }
        self.validate_dynamic_reaction(&reaction)?;
        index_dynamic_reaction(
            &mut self.dynamic_reaction_index_by_substance,
            &mut self.dynamic_unindexed_reaction_ids,
            &reaction,
        );
        self.dynamic_reactions.insert(reaction.id.clone(), reaction);
        Ok(())
    }

    fn validate_dynamic_reaction(&self, reaction: &Reaction) -> ChemistryResult<()> {
        reaction.validate_shape()?;
        for term in reaction.reactants.iter().chain(reaction.products.iter()) {
            self.substance(&term.substance_id)
                .map_err(|_| ChemistryError::UnknownSubstance {
                    reaction_id: reaction.id.to_string(),
                    substance_id: term.substance_id.to_string(),
                })?;
        }
        for ordered_substance in reaction.orders.keys() {
            self.substance(ordered_substance)
                .map_err(|_| ChemistryError::UnknownSubstance {
                    reaction_id: reaction.id.to_string(),
                    substance_id: ordered_substance.to_string(),
                })?;
        }
        let external_reactant_charge = reaction
            .external_reactants
            .iter()
            .filter_map(|requirement| {
                requirement
                    .charge
                    .map(|charge| charge * requirement.moles_per_reaction.round() as i32)
            })
            .sum::<i32>();
        let reactant_charge = reaction
            .reactants
            .iter()
            .map(|term| {
                self.substance(&term.substance_id)
                    .map(|substance| substance.charge * term.coefficient as i32)
            })
            .sum::<ChemistryResult<i32>>()?
            + external_reactant_charge;
        let product_charge = reaction
            .products
            .iter()
            .map(|term| {
                self.substance(&term.substance_id)
                    .map(|substance| substance.charge * term.coefficient as i32)
            })
            .sum::<ChemistryResult<i32>>()?;
        if reactant_charge != product_charge && !reaction.allow_charge_imbalance {
            return Err(ChemistryError::ChargeNotConserved {
                reaction_id: reaction.id.to_string(),
                reactants: reactant_charge,
                products: product_charge,
            });
        }

        let external_reactant_mass = reaction
            .external_reactants
            .iter()
            .filter_map(|requirement| {
                requirement
                    .molar_mass_grams
                    .map(|mass| mass * requirement.moles_per_reaction)
            })
            .sum::<f64>();
        let reactant_mass = reaction
            .reactants
            .iter()
            .map(|term| {
                self.substance(&term.substance_id).map(|substance| {
                    substance.molar_mass_grams * term.coefficient as f64
                })
            })
            .sum::<ChemistryResult<f64>>()?
            + external_reactant_mass;
        let product_mass = reaction
            .products
            .iter()
            .map(|term| {
                self.substance(&term.substance_id).map(|substance| {
                    substance.molar_mass_grams * term.coefficient as f64
                })
            })
            .sum::<ChemistryResult<f64>>()?;
        if (reactant_mass - product_mass).abs() > MASS_TOLERANCE_GRAMS_PER_MOL
            && !reaction.allow_mass_imbalance
        {
            return Err(ChemistryError::MassNotConserved {
                reaction_id: reaction.id.to_string(),
                reactants: reactant_mass,
                products: product_mass,
            });
        }
        Ok(())
    }

    fn validated_substance_set(
        &self,
        substance_ids: impl IntoIterator<Item = SubstanceId>,
    ) -> ChemistryResult<BTreeSet<SubstanceId>> {
        let mut result = BTreeSet::new();
        for substance_id in substance_ids {
            self.substance(&substance_id)?;
            result.insert(substance_id);
        }
        Ok(result)
    }

    fn generation_keys_for_substance(
        &self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<BTreeSet<String>> {
        let substance = self.substance(substance_id)?;
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

fn index_dynamic_reaction(
    by_substance: &mut BTreeMap<SubstanceId, BTreeSet<ReactionId>>,
    unindexed: &mut BTreeSet<ReactionId>,
    reaction: &Reaction,
) {
    let mut substances = BTreeSet::new();
    for reactant in &reaction.reactants {
        substances.insert(reactant.substance_id.clone());
    }
    for ordered_substance in reaction.orders.keys() {
        substances.insert(ordered_substance.clone());
    }

    if substances.is_empty() {
        unindexed.insert(reaction.id.clone());
        return;
    }

    for substance_id in substances {
        by_substance
            .entry(substance_id)
            .or_default()
            .insert(reaction.id.clone());
    }
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
        let substance = registry.substance(&id).unwrap();
        assert!(substance
            .tags
            .iter()
            .any(|tag| tag.as_str() == "destroy:hypothetical"));
        assert!(registry.validate_substance_can_enter_mixture(&id).is_err());
    }

    #[test]
    fn dynamic_products_are_available_for_later_generation() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let alkene = registry.resolve_frowns("CCCC=C").unwrap();
        let report = registry.generate_reactions_for(&alkene, 2).unwrap();
        assert!(report.added_substances > 0);
        assert!(report.added_reactions > 0);
        assert!(report.processed_work_items > 0);
        assert!(registry.reactions().any(|reaction| {
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
        assert!(!acid_only.reactions().any(|reaction| {
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
        assert!(acid_and_alcohol.reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("carboxylic_acid_esterification/destroy_acetic_acid/destroy_ethanol/")
        }));
    }
}
