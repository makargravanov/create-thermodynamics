use std::collections::BTreeMap;
use std::collections::BTreeSet;

use super::error::{ChemistryError, ChemistryResult};
use super::frowns::{parse_frowns, write_frowns};
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

#[derive(Debug, Clone)]
pub struct DynamicChemistryRegistry {
    registry: ChemistryRegistry,
    canonical_to_id: BTreeMap<String, SubstanceId>,
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
}

impl DynamicChemistryRegistry {
    pub fn from_destroy_catalog() -> ChemistryResult<Self> {
        Self::from_registry(super::destroy_registry_builder()?.build()?)
    }

    pub fn from_registry(registry: ChemistryRegistry) -> ChemistryResult<Self> {
        let mut result = Self {
            registry,
            canonical_to_id: BTreeMap::new(),
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
        self.generate_reactions_from_seeds(&mut seeds, max_iterations)
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
        self.generate_reactions_from_seeds(&mut seeds, max_iterations)
    }

    fn generate_reactions_from_seeds(
        &mut self,
        seeds: &mut BTreeSet<SubstanceId>,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        if max_iterations == 0 {
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
        let mut processed = BTreeSet::new();
        for iteration in 0..max_iterations {
            if queue.is_empty() {
                return Ok(DynamicGenerationReport {
                    iterations: iteration,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: 0,
                    reached_fixed_point: true,
                });
            }
            let current_seeds = queue.clone();
            queue.clear();
            let mut unprocessed_seeds = BTreeSet::new();
            for seed in &current_seeds {
                let canonical = self
                    .registry
                    .substance(seed)?
                    .molecular_structure
                    .as_ref()
                    .map(write_frowns)
                    .transpose()?
                    .unwrap_or_else(|| seed.to_string());
                if !processed.insert(canonical) {
                    skipped_duplicates += 1;
                    continue;
                }
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
                });
            }

            let generated =
                organic::generate_organic_reactions_for(&self.registry, &unprocessed_seeds)?;
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
                });
            }
        }
        Ok(DynamicGenerationReport {
            iterations: max_iterations,
            added_substances,
            added_reactions,
            processed_work_items,
            skipped_duplicates,
            remaining_queue: queue.len(),
            reached_fixed_point: false,
        })
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
}
