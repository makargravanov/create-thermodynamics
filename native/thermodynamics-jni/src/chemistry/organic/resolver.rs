use std::collections::BTreeMap;

use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::MolecularStructure;
use crate::chemistry::substance::{Substance, SubstanceId};

const DEFAULT_DERIVED_DENSITY: f64 = 1000.0;
const DEFAULT_DERIVED_HEAT_CAPACITY: f64 = 100.0;
const DEFAULT_DERIVED_LATENT_HEAT: f64 = 20_000.0;

pub(crate) struct DerivedSubstanceResolver {
    canonical_to_id: BTreeMap<String, SubstanceId>,
    generated_id_to_canonical: BTreeMap<SubstanceId, String>,
    pub(crate) substances: Vec<Substance>,
}

impl DerivedSubstanceResolver {
    pub(crate) fn new_from_canonical_to_id(canonical_to_id: BTreeMap<String, SubstanceId>) -> Self {
        Self {
            canonical_to_id,
            generated_id_to_canonical: BTreeMap::new(),
            substances: Vec::new(),
        }
    }

    pub(crate) fn resolve(
        &mut self,
        structure: MolecularStructure,
    ) -> ChemistryResult<SubstanceId> {
        let canonical = crate::chemistry::frowns::write_frowns(&structure)?;
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
        let summary = structure.summary()?;
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
}
