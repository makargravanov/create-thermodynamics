use thermodynamics_core::{
    select_candidate_species, CandidatePhaseFilter, CandidateSelectionError,
    CandidateSelectionRequest, Element, ElementId, EquilibriumProblem, Species, SpeciesAmount,
    SpeciesId, SpeciesRegistry,
};
use thermodynamics_datafile::{registry_from_database_file, DatafileError};

const EMBEDDED_DATABASE_BYTES: &[u8] =
    include_bytes!("../../generated-data/create_thermodynamics_db_v1.ctdb");
const SUPPORTED_RUNTIME_SYMBOLS: &[&str] = &[
    "H2O", "H+", "OH-", "Na+", "Cl-", "Ca+2", "CO3-2", "HCO3-", "CO2", "CaCO3(s)", "H2O(g)",
    "CO2(g)",
];

#[derive(Debug)]
pub enum DefaultDataError {
    Database(DatafileError),
    UnknownSpeciesSymbol(String),
}

#[derive(Debug)]
pub enum DefaultProblemError {
    Database(DatafileError),
    CandidateSelection(CandidateSelectionError),
}

pub fn default_registry() -> Result<SpeciesRegistry, DefaultDataError> {
    let full_registry =
        registry_from_database_file(EMBEDDED_DATABASE_BYTES).map_err(DefaultDataError::Database)?;
    subset_registry_by_symbols(&full_registry, SUPPORTED_RUNTIME_SYMBOLS)
}

pub fn default_species_ids() -> Result<Vec<SpeciesId>, DefaultDataError> {
    let registry = default_registry()?;
    Ok(registry
        .species_records()
        .map(|species| species.id)
        .collect())
}

pub fn species_id_by_symbol(
    registry: &SpeciesRegistry,
    symbol: &str,
) -> Result<SpeciesId, DefaultDataError> {
    registry
        .species_records()
        .find(|species| species.symbol == symbol)
        .map(|species| species.id)
        .ok_or_else(|| DefaultDataError::UnknownSpeciesSymbol(symbol.to_owned()))
}

pub fn source_registry() -> Result<SpeciesRegistry, DefaultDataError> {
    registry_from_database_file(EMBEDDED_DATABASE_BYTES).map_err(DefaultDataError::Database)
}

pub fn equilibrium_problem(
    temperature_kelvin: f64,
    pressure_pascal: f64,
    initial_species_amounts_mol: Vec<SpeciesAmount>,
) -> Result<EquilibriumProblem, DefaultProblemError> {
    selected_equilibrium_problem(
        temperature_kelvin,
        pressure_pascal,
        initial_species_amounts_mol,
    )
}

pub fn selected_equilibrium_problem(
    temperature_kelvin: f64,
    pressure_pascal: f64,
    initial_species_amounts_mol: Vec<SpeciesAmount>,
) -> Result<EquilibriumProblem, DefaultProblemError> {
    let registry = default_registry().map_err(|error| match error {
        DefaultDataError::Database(error) => DefaultProblemError::Database(error),
        DefaultDataError::UnknownSpeciesSymbol(symbol) => {
            unreachable!("embedded database lookup unexpectedly failed for symbol {symbol}")
        }
    })?;
    let selection = select_candidate_species(
        &registry,
        &CandidateSelectionRequest {
            temperature_kelvin,
            initial_species_amounts_mol: initial_species_amounts_mol.clone(),
            phase_filter: CandidatePhaseFilter::all_supported(),
        },
    )
    .map_err(DefaultProblemError::CandidateSelection)?;

    Ok(EquilibriumProblem {
        temperature_kelvin,
        pressure_pascal,
        initial_species_amounts_mol,
        candidate_species: selection.candidate_species,
    })
}

fn subset_registry_by_symbols(
    full_registry: &SpeciesRegistry,
    symbols: &[&str],
) -> Result<SpeciesRegistry, DefaultDataError> {
    let mut species = Vec::<Species>::new();
    let mut required_elements = Vec::<ElementId>::new();

    for symbol in symbols {
        let species_id = species_id_by_symbol(full_registry, symbol)?;
        let record = full_registry.species(species_id).unwrap().clone();
        required_elements.extend(record.composition.keys().copied());
        species.push(record);
    }

    required_elements.sort();
    required_elements.dedup();
    let elements = required_elements
        .into_iter()
        .map(|element_id| full_registry.element(element_id).unwrap().clone())
        .collect::<Vec<Element>>();

    SpeciesRegistry::new(elements, species)
        .map_err(|error| DefaultDataError::Database(DatafileError::Registry(error)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use thermodynamics_core::{
        analyze_aqueous_equilibrium, analyze_phase_equilibrium, mixture_enthalpy_joule,
        solve_closed_gas_equilibrium, solve_equilibrium, ClosedGasEquilibriumProblem, PhaseKind,
    };

    fn amount(species_id: SpeciesId, amount_mol: f64) -> SpeciesAmount {
        SpeciesAmount {
            species_id,
            amount_mol,
        }
    }

    fn result_amount(
        result: &thermodynamics_core::EquilibriumResult,
        species_id: SpeciesId,
    ) -> f64 {
        result
            .species_amounts_mol
            .iter()
            .find(|amount| amount.species_id == species_id)
            .map(|amount| amount.amount_mol)
            .unwrap_or_default()
    }

    fn aqueous_summary(
        registry: &SpeciesRegistry,
        result: &thermodynamics_core::EquilibriumResult,
    ) -> thermodynamics_core::AqueousEquilibriumSummary {
        analyze_aqueous_equilibrium(
            registry,
            result,
            species_id_by_symbol(registry, "H2O").unwrap(),
            species_id_by_symbol(registry, "H+").unwrap(),
        )
        .unwrap()
    }

    fn registry_subset(symbols: &[&str]) -> SpeciesRegistry {
        let registry = source_registry().unwrap();
        subset_registry_by_symbols(&registry, symbols).unwrap()
    }

    fn solvent_water_kg(water_mol: f64) -> f64 {
        water_mol * 0.018_015_28
    }

    #[test]
    fn embedded_database_loads_registry() {
        let registry = default_registry().unwrap();
        let source_registry = source_registry().unwrap();

        assert_eq!(
            registry.species_records().count(),
            SUPPORTED_RUNTIME_SYMBOLS.len()
        );
        assert!(source_registry.species_records().count() >= 50);
    }

    #[test]
    fn symbol_lookup_uses_embedded_database_records() {
        let registry = source_registry().unwrap();

        for symbol in [
            "H2O", "H2O(g)", "H+", "OH-", "CO2", "CO2(g)", "HCO3-", "CO3-2", "Ca+2", "CaCO3(s)",
            "O2", "O2(g)",
        ] {
            let species_id = species_id_by_symbol(&registry, symbol).unwrap();
            assert_eq!(registry.species(species_id).unwrap().symbol, symbol);
        }
    }

    #[test]
    fn selected_candidates_for_pure_water_stay_within_hydrogen_oxygen_species() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let water_gas = species_id_by_symbol(&registry, "H2O(g)").unwrap();
        let hydrogen = species_id_by_symbol(&registry, "H+").unwrap();
        let hydroxide = species_id_by_symbol(&registry, "OH-").unwrap();
        let carbon_dioxide = species_id_by_symbol(&registry, "CO2").unwrap();

        let problem =
            selected_equilibrium_problem(298.15, 101_325.0, vec![amount(water, 55.5)]).unwrap();

        assert!(problem.candidate_species.contains(&water));
        assert!(problem.candidate_species.contains(&water_gas));
        assert!(problem.candidate_species.contains(&hydrogen));
        assert!(problem.candidate_species.contains(&hydroxide));
        assert!(!problem.candidate_species.contains(&carbon_dioxide));
    }

    #[test]
    fn selected_candidates_for_water_and_carbon_dioxide_include_carbonate_family() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let carbon_dioxide_gas = species_id_by_symbol(&registry, "CO2(g)").unwrap();
        let carbon_dioxide_aq = species_id_by_symbol(&registry, "CO2").unwrap();
        let bicarbonate = species_id_by_symbol(&registry, "HCO3-").unwrap();
        let carbonate = species_id_by_symbol(&registry, "CO3-2").unwrap();
        let calcium = species_id_by_symbol(&registry, "Ca+2").unwrap();

        let problem = selected_equilibrium_problem(
            298.15,
            101_325.0,
            vec![amount(water, 55.5), amount(carbon_dioxide_gas, 1.0)],
        )
        .unwrap();

        assert!(problem.candidate_species.contains(&water));
        assert!(problem.candidate_species.contains(&carbon_dioxide_gas));
        assert!(problem.candidate_species.contains(&carbon_dioxide_aq));
        assert!(problem.candidate_species.contains(&bicarbonate));
        assert!(problem.candidate_species.contains(&carbonate));
        assert!(!problem.candidate_species.contains(&calcium));
    }

    #[test]
    fn selected_problem_rejects_zero_only_input() {
        let registry = default_registry().unwrap();
        let carbon_dioxide_gas = species_id_by_symbol(&registry, "CO2(g)").unwrap();

        assert!(matches!(
            selected_equilibrium_problem(298.15, 101_325.0, vec![amount(carbon_dioxide_gas, 0.0)]),
            Err(DefaultProblemError::CandidateSelection(
                CandidateSelectionError::NoPositiveInputAmounts
            ))
        ));
    }

    #[test]
    fn liquid_water_enthalpy_and_vaporization_come_from_embedded_data() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let water_gas = species_id_by_symbol(&registry, "H2O(g)").unwrap();
        let liquid = mixture_enthalpy_joule(&registry, &[amount(water, 1.0)], 298.15).unwrap();
        let vapor = mixture_enthalpy_joule(&registry, &[amount(water_gas, 1.0)], 298.15).unwrap();

        assert!(liquid.is_finite());
        assert!(vapor.is_finite());
        assert!((vapor - liquid) > 1.0e4);
    }

    #[test]
    fn pure_water_equilibrium_solves_with_runtime_subset() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let problem = equilibrium_problem(298.15, 101_325.0, vec![amount(water, 55.5)]).unwrap();
        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);

        assert!(summary.ph.is_finite());
        assert!(summary.solvent_water_mass_kg > 0.9);
    }

    #[test]
    fn carbon_dioxide_dissolution_uses_embedded_aqueous_and_gas_forms() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let carbon_dioxide_aq = species_id_by_symbol(&registry, "CO2").unwrap();
        let carbon_dioxide_gas = species_id_by_symbol(&registry, "CO2(g)").unwrap();
        let problem = EquilibriumProblem {
            temperature_kelvin: 298.15,
            pressure_pascal: 101_325.0,
            initial_species_amounts_mol: vec![amount(water, 55.5), amount(carbon_dioxide_gas, 1.0)],
            candidate_species: vec![water, carbon_dioxide_aq, carbon_dioxide_gas],
        };

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let dissolved_molality = result_amount(&result, carbon_dioxide_aq)
            / solvent_water_kg(result_amount(&result, water));

        assert!(dissolved_molality.is_finite());
        assert!(dissolved_molality >= 0.0);
        assert!(result_amount(&result, carbon_dioxide_gas) > 0.0);
    }

    #[test]
    fn oxygen_dissolution_uses_embedded_aqueous_and_gas_forms() {
        let registry = registry_subset(&["H2O", "O2", "O2(g)"]);
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let oxygen_aq = species_id_by_symbol(&registry, "O2").unwrap();
        let oxygen_gas = species_id_by_symbol(&registry, "O2(g)").unwrap();
        let problem = EquilibriumProblem {
            temperature_kelvin: 298.15,
            pressure_pascal: 101_325.0,
            initial_species_amounts_mol: vec![amount(water, 55.5), amount(oxygen_gas, 1.0)],
            candidate_species: vec![water, oxygen_aq, oxygen_gas],
        };

        let result = solve_equilibrium(&registry, &problem).unwrap();

        assert!(result_amount(&result, oxygen_aq).is_finite());
        assert!(result_amount(&result, oxygen_aq) >= 0.0);
        assert!(result_amount(&result, oxygen_gas) > 0.0);
    }

    #[test]
    fn phase_analysis_reports_vaporized_water_boundary() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let water_gas = species_id_by_symbol(&registry, "H2O(g)").unwrap();
        let problem = EquilibriumProblem {
            temperature_kelvin: 298.15,
            pressure_pascal: 1_500.0,
            initial_species_amounts_mol: vec![amount(water, 1.0)],
            candidate_species: vec![water, water_gas],
        };

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = analyze_phase_equilibrium(&registry, &result).unwrap();
        let gas_phase = summary
            .phases
            .iter()
            .find(|phase| phase.phase == PhaseKind::Gas)
            .unwrap();

        assert!(gas_phase.total_amount_mol > 0.99);
        assert!(summary
            .boundary_species
            .iter()
            .any(|amount| amount.species_id == water));
    }

    #[test]
    fn closed_gas_equilibrium_solves_pressure_from_embedded_data() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let carbon_dioxide_aq = species_id_by_symbol(&registry, "CO2").unwrap();
        let carbon_dioxide_gas = species_id_by_symbol(&registry, "CO2(g)").unwrap();
        let problem = ClosedGasEquilibriumProblem {
            temperature_kelvin: 298.15,
            gas_volume_cubic_meter: 0.024_465,
            initial_species_amounts_mol: vec![amount(water, 55.5), amount(carbon_dioxide_gas, 1.0)],
            candidate_species: vec![water, carbon_dioxide_aq, carbon_dioxide_gas],
        };

        let result = solve_closed_gas_equilibrium(&registry, &problem).unwrap();

        assert!(result.pressure_residual_pascal.abs() < 1.0e-3);
        assert!((90_000.0..110_000.0).contains(&result.pressure_pascal));
        assert!(result_amount(&result.equilibrium, carbon_dioxide_aq).is_finite());
        assert!(result_amount(&result.equilibrium, carbon_dioxide_aq) >= 0.0);
    }

    #[test]
    fn runtime_registry_contains_carbonate_runtime_slice() {
        let registry = default_registry().unwrap();
        let bicarbonate = species_id_by_symbol(&registry, "HCO3-").unwrap();
        let carbonate = species_id_by_symbol(&registry, "CO3-2").unwrap();
        let carbon_dioxide = species_id_by_symbol(&registry, "CO2").unwrap();

        assert!(registry.species(carbon_dioxide).is_some());
        assert!(registry.species(bicarbonate).is_some());
        assert!(registry.species(carbonate).is_some());
    }

    #[test]
    fn selected_candidates_include_calcium_carbonate_runtime_slice() {
        let registry = default_registry().unwrap();
        let water = species_id_by_symbol(&registry, "H2O").unwrap();
        let calcium = species_id_by_symbol(&registry, "Ca+2").unwrap();
        let carbonate = species_id_by_symbol(&registry, "CO3-2").unwrap();
        let calcite = species_id_by_symbol(&registry, "CaCO3(s)").unwrap();
        let problem = selected_equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(water, 55.5),
                amount(calcium, 2.0e-2),
                amount(carbonate, 2.0e-2),
            ],
        )
        .unwrap();

        assert!(problem.candidate_species.contains(&calcite));
    }

    #[test]
    fn supported_runtime_species_have_complete_thermal_data_when_expected() {
        let registry = default_registry().unwrap();
        for symbol in [
            "H2O", "OH-", "Na+", "Cl-", "Ca+2", "CO2", "CaCO3(s)", "H2O(g)", "CO2(g)",
        ] {
            let species = registry
                .species(species_id_by_symbol(&registry, symbol).unwrap())
                .unwrap();
            assert!(
                species.thermo.standard_enthalpy_of_formation.is_some(),
                "{symbol}"
            );
            assert!(
                species.thermo.constant_pressure_heat_capacity.is_some(),
                "{symbol}"
            );
        }
    }

    #[test]
    fn default_species_ids_match_registry_contents() {
        let registry = default_registry().unwrap();
        let ids = default_species_ids().unwrap();

        assert_eq!(ids.len(), registry.species_records().count());
        assert!(ids.windows(2).all(|pair| pair[0] < pair[1]));
    }
}
