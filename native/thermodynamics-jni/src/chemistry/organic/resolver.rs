use std::collections::BTreeMap;

use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::MolecularStructure;
use crate::chemistry::substance::{Substance, SubstanceId};

const DEFAULT_DERIVED_DENSITY: f64 = 1000.0;
const DEFAULT_DERIVED_HEAT_CAPACITY: f64 = 100.0;
const DEFAULT_DERIVED_LATENT_HEAT: f64 = 20_000.0;

pub(crate) struct DerivedSubstanceResolver {
    canonical_to_id: BTreeMap<String, SubstanceId>,
    known_structures_by_id: BTreeMap<SubstanceId, MolecularStructure>,
    generated_id_to_canonical: BTreeMap<SubstanceId, String>,
    generation_context: Option<String>,
    pub(crate) substances: Vec<Substance>,
}

impl DerivedSubstanceResolver {
    #[cfg(test)]
    pub(crate) fn new_from_canonical_to_id(canonical_to_id: BTreeMap<String, SubstanceId>) -> Self {
        Self::new_from_known_structures(canonical_to_id, BTreeMap::new())
    }

    pub(crate) fn new_from_known_structures(
        canonical_to_id: BTreeMap<String, SubstanceId>,
        known_structures_by_id: BTreeMap<SubstanceId, MolecularStructure>,
    ) -> Self {
        Self {
            canonical_to_id,
            known_structures_by_id,
            generated_id_to_canonical: BTreeMap::new(),
            generation_context: None,
            substances: Vec::new(),
        }
    }

    pub(crate) fn set_generation_context(&mut self, generation_context: impl Into<String>) {
        self.generation_context = Some(generation_context.into());
    }

    pub(crate) fn known_structure(
        &self,
        substance_id: &SubstanceId,
    ) -> Option<&MolecularStructure> {
        self.known_structures_by_id.get(substance_id)
    }

    pub(crate) fn resolve(
        &mut self,
        structure: MolecularStructure,
    ) -> ChemistryResult<SubstanceId> {
        let canonical = crate::chemistry::frowns::write_frowns(&structure).map_err(|error| {
            ChemistryError::InvalidSubstance {
                substance_id: "generated".to_string(),
                reason: format!(
                    "{error}; context={}; atoms={:?}; bonds={:?}; stereo={:?}; source={}",
                    self.generation_context.as_deref().unwrap_or("<unknown>"),
                    structure.atoms,
                    structure.bonds,
                    structure.stereochemistry,
                    structure.source_code
                ),
            }
        })?;
        if let Some(id) = self.canonical_to_id.get(&canonical) {
            return Ok(id.clone());
        }
        let id = SubstanceId::new(canonical.clone())?;
        if let Some(existing) = self.generated_id_to_canonical.get(&id) {
            if existing != &canonical {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: id.to_string(),
                    reason: "derived substance id collision".to_string(),
                });
            }
        }
        let summary = structure
            .summary()
            .map_err(|error| ChemistryError::InvalidSubstance {
                substance_id: "generated".to_string(),
                reason: format!(
                    "{error}; context={}; atoms={:?}; bonds={:?}; stereo={:?}; source={}",
                    self.generation_context.as_deref().unwrap_or("<unknown>"),
                    structure.atoms,
                    structure.bonds,
                    structure.stereochemistry,
                    structure.source_code
                ),
            })?;
        let substance = Substance::new(
            id.clone(),
            summary.charge,
            summary.molar_mass_grams,
            DEFAULT_DERIVED_DENSITY,
            if summary.charge == 0 {
                1000.0
            } else {
                f64::MAX
            },
            DEFAULT_DERIVED_HEAT_CAPACITY,
            DEFAULT_DERIVED_LATENT_HEAT,
        )
        .with_catalog_metadata(Some(canonical.clone()), None, 0x20FF_FFFF, Vec::new())
        .with_molecular_structure(structure);
        self.canonical_to_id.insert(canonical.clone(), id.clone());
        self.generated_id_to_canonical.insert(id.clone(), canonical);
        self.substances.push(substance);
        Ok(id)
    }

    pub(crate) fn resolve_substance(
        &mut self,
        substance: Substance,
    ) -> ChemistryResult<SubstanceId> {
        substance
            .validate()
            .map_err(|error| ChemistryError::InvalidSubstance {
                substance_id: substance.id.to_string(),
                reason: match &substance.molecular_structure {
                    Some(structure) => format!(
                        "{error}; context={}; atoms={:?}; bonds={:?}; stereo={:?}; source={}",
                        self.generation_context.as_deref().unwrap_or("<unknown>"),
                        structure.atoms,
                        structure.bonds,
                        structure.stereochemistry,
                        structure.source_code
                    ),
                    None => error.to_string(),
                },
            })?;
        let id = substance.id.clone();
        let canonical = format!("material:{}", id.as_str());
        if let Some(existing) = self.canonical_to_id.get(&canonical) {
            return Ok(existing.clone());
        }
        if let Some(existing) = self.generated_id_to_canonical.get(&id) {
            if existing != &canonical {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: id.to_string(),
                    reason: "derived material id collision".to_string(),
                });
            }
            return Ok(id);
        }
        if self
            .canonical_to_id
            .values()
            .any(|existing| existing == &id)
        {
            return Ok(id);
        }
        self.canonical_to_id.insert(canonical.clone(), id.clone());
        self.generated_id_to_canonical.insert(id.clone(), canonical);
        self.substances.push(substance);
        Ok(id)
    }
}
