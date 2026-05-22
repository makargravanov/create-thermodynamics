use crate::chemistry::{ElementId, PhaseKind, SpeciesAmount, SpeciesId};
use crate::registry::SpeciesRegistry;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct CandidateSelectionRequest {
    pub temperature_kelvin: f64,
    pub initial_species_amounts_mol: Vec<SpeciesAmount>,
    pub phase_filter: CandidatePhaseFilter,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidatePhaseFilter {
    pub aqueous: bool,
    pub solid: bool,
    pub gas: bool,
}

impl CandidatePhaseFilter {
    pub const fn all_supported() -> Self {
        Self {
            aqueous: true,
            solid: true,
            gas: true,
        }
    }

    pub const fn allows(self, phase: PhaseKind) -> bool {
        match phase {
            PhaseKind::Aqueous => self.aqueous,
            PhaseKind::Solid => self.solid,
            PhaseKind::Gas => self.gas,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CandidateSelection {
    pub candidate_species: Vec<SpeciesId>,
    pub available_elements: Vec<ElementId>,
    pub initial_species: Vec<SpeciesAmount>,
    pub excluded_species: Vec<CandidateExclusion>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CandidateExclusion {
    pub species_id: SpeciesId,
    pub reason: CandidateExclusionReason,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CandidateExclusionReason {
    MissingElement(ElementId),
    PhaseDisabled(PhaseKind),
    TemperatureOutOfRange,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CandidateSelectionError {
    InvalidTemperatureKelvin(f64),
    NegativeAmount {
        species_id: SpeciesId,
        amount_mol: f64,
    },
    UnknownInputSpecies(SpeciesId),
    InputSpeciesPhaseDisabled {
        species_id: SpeciesId,
        phase: PhaseKind,
    },
    InputSpeciesTemperatureOutOfRange {
        species_id: SpeciesId,
        temperature_kelvin: f64,
        valid_min_temperature_kelvin: f64,
        valid_max_temperature_kelvin: f64,
    },
    NoPositiveInputAmounts,
}

pub fn select_candidate_species(
    registry: &SpeciesRegistry,
    request: &CandidateSelectionRequest,
) -> Result<CandidateSelection, CandidateSelectionError> {
    if !request.temperature_kelvin.is_finite() || request.temperature_kelvin <= 0.0 {
        return Err(CandidateSelectionError::InvalidTemperatureKelvin(
            request.temperature_kelvin,
        ));
    }

    let initial_amounts = normalized_positive_amounts(&request.initial_species_amounts_mol)?;
    if initial_amounts.is_empty() {
        return Err(CandidateSelectionError::NoPositiveInputAmounts);
    }

    let mut available_elements = BTreeSet::new();
    for amount in &initial_amounts {
        let species = registry.species(amount.species_id).ok_or(
            CandidateSelectionError::UnknownInputSpecies(amount.species_id),
        )?;
        if !request.phase_filter.allows(species.phase) {
            return Err(CandidateSelectionError::InputSpeciesPhaseDisabled {
                species_id: amount.species_id,
                phase: species.phase,
            });
        }
        let valid_range = species.thermo.valid_temperature_range;
        if request.temperature_kelvin < valid_range.min_kelvin
            || request.temperature_kelvin > valid_range.max_kelvin
        {
            return Err(CandidateSelectionError::InputSpeciesTemperatureOutOfRange {
                species_id: amount.species_id,
                temperature_kelvin: request.temperature_kelvin,
                valid_min_temperature_kelvin: valid_range.min_kelvin,
                valid_max_temperature_kelvin: valid_range.max_kelvin,
            });
        }
        available_elements.extend(species.composition.keys().copied());
    }

    let mut candidate_species = Vec::new();
    let mut excluded_species = Vec::new();
    for species in registry.species_records() {
        if let Some(reason) = candidate_exclusion_reason(
            species.phase,
            species.composition.keys().copied(),
            species.thermo.valid_temperature_range.min_kelvin,
            species.thermo.valid_temperature_range.max_kelvin,
            &available_elements,
            request.phase_filter,
            request.temperature_kelvin,
        ) {
            excluded_species.push(CandidateExclusion {
                species_id: species.id,
                reason,
            });
        } else {
            candidate_species.push(species.id);
        }
    }

    Ok(CandidateSelection {
        candidate_species,
        available_elements: available_elements.into_iter().collect(),
        initial_species: initial_amounts,
        excluded_species,
    })
}

fn normalized_positive_amounts(
    amounts: &[SpeciesAmount],
) -> Result<Vec<SpeciesAmount>, CandidateSelectionError> {
    let mut normalized = BTreeMap::<SpeciesId, f64>::new();
    for amount in amounts {
        if !amount.amount_mol.is_finite() || amount.amount_mol < 0.0 {
            return Err(CandidateSelectionError::NegativeAmount {
                species_id: amount.species_id,
                amount_mol: amount.amount_mol,
            });
        }
        if amount.amount_mol == 0.0 {
            continue;
        }
        *normalized.entry(amount.species_id).or_default() += amount.amount_mol;
    }

    Ok(normalized
        .into_iter()
        .map(|(species_id, amount_mol)| SpeciesAmount {
            species_id,
            amount_mol,
        })
        .collect())
}

fn candidate_exclusion_reason(
    phase: PhaseKind,
    composition: impl Iterator<Item = ElementId>,
    valid_min_temperature_kelvin: f64,
    valid_max_temperature_kelvin: f64,
    available_elements: &BTreeSet<ElementId>,
    phase_filter: CandidatePhaseFilter,
    temperature_kelvin: f64,
) -> Option<CandidateExclusionReason> {
    if !phase_filter.allows(phase) {
        return Some(CandidateExclusionReason::PhaseDisabled(phase));
    }
    if temperature_kelvin < valid_min_temperature_kelvin
        || temperature_kelvin > valid_max_temperature_kelvin
    {
        return Some(CandidateExclusionReason::TemperatureOutOfRange);
    }
    composition
        .filter(|element_id| !available_elements.contains(element_id))
        .min()
        .map(CandidateExclusionReason::MissingElement)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::{
        ActivityModel, ConstantPressureHeatCapacity, DataSource, Element, Species,
        StandardEnthalpyOfFormation, StandardGibbsEnergy, StandardThermo, TemperatureRange,
    };

    const H: ElementId = ElementId(1);
    const C: ElementId = ElementId(6);
    const O: ElementId = ElementId(8);
    const WATER: SpeciesId = SpeciesId(1);
    const WATER_GAS: SpeciesId = SpeciesId(2);
    const CARBON_DIOXIDE: SpeciesId = SpeciesId(3);
    const SOLID_CARBON: SpeciesId = SpeciesId(4);
    const HOT_ONLY_WATER: SpeciesId = SpeciesId(5);

    fn amount(species_id: SpeciesId, amount_mol: f64) -> SpeciesAmount {
        SpeciesAmount {
            species_id,
            amount_mol,
        }
    }

    fn registry() -> SpeciesRegistry {
        SpeciesRegistry::new(
            vec![
                Element {
                    id: H,
                    atomic_number: 1,
                    symbol: "H".to_owned(),
                },
                Element {
                    id: C,
                    atomic_number: 6,
                    symbol: "C".to_owned(),
                },
                Element {
                    id: O,
                    atomic_number: 8,
                    symbol: "O".to_owned(),
                },
            ],
            vec![
                species(
                    WATER,
                    "H2O(l)",
                    &[(H, 2), (O, 1)],
                    PhaseKind::Aqueous,
                    273.15,
                    373.15,
                ),
                species(
                    WATER_GAS,
                    "H2O(g)",
                    &[(H, 2), (O, 1)],
                    PhaseKind::Gas,
                    273.15,
                    373.15,
                ),
                species(
                    CARBON_DIOXIDE,
                    "CO2(g)",
                    &[(C, 1), (O, 2)],
                    PhaseKind::Gas,
                    273.15,
                    373.15,
                ),
                species(
                    SOLID_CARBON,
                    "C(s)",
                    &[(C, 1)],
                    PhaseKind::Solid,
                    273.15,
                    373.15,
                ),
                species(
                    HOT_ONLY_WATER,
                    "H2O(hot)",
                    &[(H, 2), (O, 1)],
                    PhaseKind::Aqueous,
                    500.0,
                    600.0,
                ),
            ],
        )
        .unwrap()
    }

    fn species(
        id: SpeciesId,
        symbol: &'static str,
        composition: &[(ElementId, u16)],
        phase: PhaseKind,
        min_kelvin: f64,
        max_kelvin: f64,
    ) -> Species {
        let source = DataSource {
            citation: "test".to_owned(),
            note: "test".to_owned(),
        };
        Species {
            id,
            symbol: symbol.to_owned(),
            composition: composition.iter().copied().collect(),
            charge_number: 0,
            phase,
            activity_model: match phase {
                PhaseKind::Aqueous => ActivityModel::UnitActivity,
                PhaseKind::Solid => ActivityModel::UnitActivity,
                PhaseKind::Gas => ActivityModel::IdealGas,
            },
            thermo: StandardThermo {
                standard_gibbs_energy: StandardGibbsEnergy {
                    value_joule_per_mol: 0.0,
                    reference_temperature_kelvin: min_kelvin,
                    source: source.clone(),
                },
                standard_enthalpy_of_formation: Some(StandardEnthalpyOfFormation {
                    value_joule_per_mol: 0.0,
                    reference_temperature_kelvin: min_kelvin,
                    source: source.clone(),
                }),
                constant_pressure_heat_capacity: Some(ConstantPressureHeatCapacity {
                    value_joule_per_mol_kelvin: 1.0,
                    source,
                }),
                valid_temperature_range: TemperatureRange {
                    min_kelvin,
                    max_kelvin,
                },
            },
        }
    }

    fn request(initial_species_amounts_mol: Vec<SpeciesAmount>) -> CandidateSelectionRequest {
        CandidateSelectionRequest {
            temperature_kelvin: 298.15,
            initial_species_amounts_mol,
            phase_filter: CandidatePhaseFilter::all_supported(),
        }
    }

    #[test]
    fn selects_only_species_composed_from_available_elements() {
        let selection =
            select_candidate_species(&registry(), &request(vec![amount(WATER, 1.0)])).unwrap();

        assert_eq!(selection.available_elements, vec![H, O]);
        assert!(selection.candidate_species.contains(&WATER));
        assert!(selection.candidate_species.contains(&WATER_GAS));
        assert!(!selection.candidate_species.contains(&CARBON_DIOXIDE));
        assert!(selection.excluded_species.iter().any(|excluded| {
            excluded.species_id == CARBON_DIOXIDE
                && excluded.reason == CandidateExclusionReason::MissingElement(C)
        }));
    }

    #[test]
    fn does_not_add_water_without_water_elements() {
        let selection =
            select_candidate_species(&registry(), &request(vec![amount(SOLID_CARBON, 1.0)]))
                .unwrap();

        assert_eq!(selection.available_elements, vec![C]);
        assert!(selection.candidate_species.contains(&SOLID_CARBON));
        assert!(!selection.candidate_species.contains(&WATER));
    }

    #[test]
    fn normalizes_and_preserves_positive_input_species() {
        let selection = select_candidate_species(
            &registry(),
            &request(vec![amount(WATER, 0.25), amount(WATER, 0.75)]),
        )
        .unwrap();

        assert_eq!(selection.initial_species, vec![amount(WATER, 1.0)]);
        assert!(selection.candidate_species.contains(&WATER));
    }

    #[test]
    fn excludes_disabled_phases() {
        let mut request = request(vec![amount(WATER, 1.0)]);
        request.phase_filter = CandidatePhaseFilter {
            aqueous: true,
            solid: true,
            gas: false,
        };

        let selection = select_candidate_species(&registry(), &request).unwrap();

        assert!(!selection.candidate_species.contains(&WATER_GAS));
        assert!(selection.excluded_species.iter().any(|excluded| {
            excluded.species_id == WATER_GAS
                && excluded.reason == CandidateExclusionReason::PhaseDisabled(PhaseKind::Gas)
        }));
    }

    #[test]
    fn rejects_input_species_disabled_by_phase_filter() {
        let mut request = request(vec![amount(WATER_GAS, 1.0)]);
        request.phase_filter = CandidatePhaseFilter {
            aqueous: true,
            solid: true,
            gas: false,
        };

        assert!(matches!(
            select_candidate_species(&registry(), &request),
            Err(CandidateSelectionError::InputSpeciesPhaseDisabled {
                species_id: WATER_GAS,
                phase: PhaseKind::Gas
            })
        ));
    }

    #[test]
    fn excludes_species_outside_temperature_range() {
        let selection =
            select_candidate_species(&registry(), &request(vec![amount(WATER, 1.0)])).unwrap();

        assert!(!selection.candidate_species.contains(&HOT_ONLY_WATER));
        assert!(selection.excluded_species.iter().any(|excluded| {
            excluded.species_id == HOT_ONLY_WATER
                && excluded.reason == CandidateExclusionReason::TemperatureOutOfRange
        }));
    }

    #[test]
    fn rejects_input_species_outside_temperature_range() {
        assert!(matches!(
            select_candidate_species(&registry(), &request(vec![amount(HOT_ONLY_WATER, 1.0)])),
            Err(CandidateSelectionError::InputSpeciesTemperatureOutOfRange {
                species_id: HOT_ONLY_WATER,
                temperature_kelvin: 298.15,
                valid_min_temperature_kelvin: 500.0,
                valid_max_temperature_kelvin: 600.0
            })
        ));
    }

    #[test]
    fn rejects_invalid_temperature() {
        let mut request = request(vec![amount(WATER, 1.0)]);
        request.temperature_kelvin = 0.0;

        assert!(matches!(
            select_candidate_species(&registry(), &request),
            Err(CandidateSelectionError::InvalidTemperatureKelvin(0.0))
        ));
    }

    #[test]
    fn rejects_negative_amounts() {
        assert!(matches!(
            select_candidate_species(&registry(), &request(vec![amount(WATER, -1.0)])),
            Err(CandidateSelectionError::NegativeAmount {
                species_id: WATER,
                amount_mol: -1.0
            })
        ));
    }

    #[test]
    fn rejects_unknown_input_species() {
        assert!(matches!(
            select_candidate_species(&registry(), &request(vec![amount(SpeciesId(99), 1.0)])),
            Err(CandidateSelectionError::UnknownInputSpecies(SpeciesId(99)))
        ));
    }

    #[test]
    fn rejects_empty_positive_input() {
        assert!(matches!(
            select_candidate_species(&registry(), &request(vec![amount(WATER, 0.0)])),
            Err(CandidateSelectionError::NoPositiveInputAmounts)
        ));
    }

    #[test]
    fn result_is_deterministic_for_input_order() {
        let first = select_candidate_species(
            &registry(),
            &request(vec![amount(WATER, 1.0), amount(SOLID_CARBON, 1.0)]),
        )
        .unwrap();
        let second = select_candidate_species(
            &registry(),
            &request(vec![amount(SOLID_CARBON, 1.0), amount(WATER, 1.0)]),
        )
        .unwrap();

        assert_eq!(first, second);
    }
}
