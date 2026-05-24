use crate::chemistry::{Element, ElementId, Species, SpeciesId};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SpeciesRegistryError {
    DuplicateElement(ElementId),
    DuplicateSpecies(SpeciesId),
    UnknownElement {
        species_id: SpeciesId,
        element_id: ElementId,
    },
    EmptyComposition(SpeciesId),
    MissingThermoData(SpeciesId),
}

#[derive(Debug, Clone)]
pub struct SpeciesRegistry {
    elements: BTreeMap<ElementId, Element>,
    species: BTreeMap<SpeciesId, Species>,
}

impl SpeciesRegistry {
    pub fn new(
        elements: Vec<Element>,
        species: Vec<Species>,
    ) -> Result<Self, SpeciesRegistryError> {
        let mut element_map = BTreeMap::new();
        for element in elements {
            let element_id = element.id;
            if element_map.insert(element_id, element).is_some() {
                return Err(SpeciesRegistryError::DuplicateElement(element_id));
            }
        }

        let element_ids: BTreeSet<ElementId> = element_map.keys().copied().collect();
        let mut species_map = BTreeMap::new();
        for species_record in species {
            let species_id = species_record.id;
            if species_record.composition.is_empty() {
                return Err(SpeciesRegistryError::EmptyComposition(species_id));
            }
            if !species_record
                .thermo
                .standard_gibbs_energy_joule_per_mol_298_15
                .is_finite()
            {
                return Err(SpeciesRegistryError::MissingThermoData(species_id));
            }
            for element_id in species_record.composition.keys() {
                if !element_ids.contains(element_id) {
                    return Err(SpeciesRegistryError::UnknownElement {
                        species_id,
                        element_id: *element_id,
                    });
                }
            }
            if species_map.insert(species_id, species_record).is_some() {
                return Err(SpeciesRegistryError::DuplicateSpecies(species_id));
            }
        }

        Ok(Self {
            elements: element_map,
            species: species_map,
        })
    }

    pub fn element(&self, id: ElementId) -> Option<&Element> {
        self.elements.get(&id)
    }

    pub fn species(&self, id: SpeciesId) -> Option<&Species> {
        self.species.get(&id)
    }

    pub fn elements(&self) -> impl Iterator<Item = &Element> {
        self.elements.values()
    }

    pub fn species_records(&self) -> impl Iterator<Item = &Species> {
        self.species.values()
    }
}
