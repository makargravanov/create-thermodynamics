use crate::chemistry::{ActivityModel, Element, ElementId, PhaseKind, Species, SpeciesId};
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
    InvalidTemperatureRange(SpeciesId),
    MissingDataSource(SpeciesId),
    IncompatibleActivityModel {
        species_id: SpeciesId,
        phase: PhaseKind,
        activity_model: ActivityModel,
    },
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
                .standard_gibbs_energy
                .value_joule_per_mol
                .is_finite()
                || !species_record
                    .thermo
                    .standard_gibbs_energy
                    .reference_temperature_kelvin
                    .is_finite()
            {
                return Err(SpeciesRegistryError::MissingThermoData(species_id));
            }
            if !species_record
                .thermo
                .standard_enthalpy_of_formation
                .value_joule_per_mol
                .is_finite()
                || !species_record
                    .thermo
                    .standard_enthalpy_of_formation
                    .reference_temperature_kelvin
                    .is_finite()
                || !species_record
                    .thermo
                    .constant_pressure_heat_capacity
                    .value_joule_per_mol_kelvin
                    .is_finite()
            {
                return Err(SpeciesRegistryError::MissingThermoData(species_id));
            }
            let valid_range = species_record.thermo.valid_temperature_range;
            if !valid_range.min_kelvin.is_finite()
                || !valid_range.max_kelvin.is_finite()
                || valid_range.min_kelvin <= 0.0
                || valid_range.min_kelvin > valid_range.max_kelvin
                || species_record
                    .thermo
                    .standard_gibbs_energy
                    .reference_temperature_kelvin
                    < valid_range.min_kelvin
                || species_record
                    .thermo
                    .standard_gibbs_energy
                    .reference_temperature_kelvin
                    > valid_range.max_kelvin
                || species_record
                    .thermo
                    .standard_enthalpy_of_formation
                    .reference_temperature_kelvin
                    < valid_range.min_kelvin
                || species_record
                    .thermo
                    .standard_enthalpy_of_formation
                    .reference_temperature_kelvin
                    > valid_range.max_kelvin
            {
                return Err(SpeciesRegistryError::InvalidTemperatureRange(species_id));
            }
            for source in [
                species_record.thermo.standard_gibbs_energy.source,
                species_record.thermo.standard_enthalpy_of_formation.source,
                species_record.thermo.constant_pressure_heat_capacity.source,
            ] {
                if source.citation.trim().is_empty() || source.note.trim().is_empty() {
                    return Err(SpeciesRegistryError::MissingDataSource(species_id));
                }
            }
            match (species_record.phase, species_record.activity_model) {
                (PhaseKind::Aqueous, ActivityModel::DaviesAqueous)
                | (PhaseKind::Aqueous, ActivityModel::IdealMolalityAqueous)
                | (PhaseKind::Aqueous, ActivityModel::UnitActivity)
                | (PhaseKind::Solid, ActivityModel::UnitActivity)
                | (PhaseKind::Gas, ActivityModel::IdealGas)
                | (PhaseKind::Gas, ActivityModel::UnitActivity) => {}
                (phase, activity_model) => {
                    return Err(SpeciesRegistryError::IncompatibleActivityModel {
                        species_id,
                        phase,
                        activity_model,
                    });
                }
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
