use std::collections::BTreeMap;
use thermodynamics_core::{
    Element, ElementId, EquilibriumProblem, PhaseKind, Species, SpeciesAmount, SpeciesId,
    SpeciesRegistry, SpeciesRegistryError, StandardThermo,
};

pub const ELEMENT_H: ElementId = ElementId(1);
pub const ELEMENT_C: ElementId = ElementId(6);
pub const ELEMENT_O: ElementId = ElementId(8);
pub const ELEMENT_NA: ElementId = ElementId(11);
pub const ELEMENT_CL: ElementId = ElementId(17);
pub const ELEMENT_CA: ElementId = ElementId(20);

pub const H2O_L: SpeciesId = SpeciesId(1);
pub const H_PLUS: SpeciesId = SpeciesId(2);
pub const OH_MINUS: SpeciesId = SpeciesId(3);
pub const NA_PLUS: SpeciesId = SpeciesId(4);
pub const CL_MINUS: SpeciesId = SpeciesId(5);
pub const CA_2_PLUS: SpeciesId = SpeciesId(6);
pub const CO3_2_MINUS: SpeciesId = SpeciesId(7);
pub const HCO3_MINUS: SpeciesId = SpeciesId(8);
pub const CO2_AQ: SpeciesId = SpeciesId(9);
pub const CACO3_S: SpeciesId = SpeciesId(10);

const VALID_MIN_TEMPERATURE_KELVIN: f64 = 273.15;
const VALID_MAX_TEMPERATURE_KELVIN: f64 = 373.15;
const THERMO_PROVENANCE: &str =
    "Curated first-iteration constants near 298.15 K; values derived from common tabulated standard molar Gibbs energies of formation.";

pub fn curated_registry() -> Result<SpeciesRegistry, SpeciesRegistryError> {
    SpeciesRegistry::new(curated_elements(), curated_species())
}

pub fn curated_species_ids() -> Vec<SpeciesId> {
    vec![
        H2O_L,
        H_PLUS,
        OH_MINUS,
        NA_PLUS,
        CL_MINUS,
        CA_2_PLUS,
        CO3_2_MINUS,
        HCO3_MINUS,
        CO2_AQ,
        CACO3_S,
    ]
}

pub fn equilibrium_problem(
    temperature_kelvin: f64,
    pressure_pascal: f64,
    initial_species_amounts_mol: Vec<SpeciesAmount>,
) -> EquilibriumProblem {
    EquilibriumProblem {
        temperature_kelvin,
        pressure_pascal,
        initial_species_amounts_mol,
        candidate_species: curated_species_ids(),
    }
}

fn curated_elements() -> Vec<Element> {
    vec![
        Element {
            id: ELEMENT_H,
            atomic_number: 1,
            symbol: "H",
        },
        Element {
            id: ELEMENT_C,
            atomic_number: 6,
            symbol: "C",
        },
        Element {
            id: ELEMENT_O,
            atomic_number: 8,
            symbol: "O",
        },
        Element {
            id: ELEMENT_NA,
            atomic_number: 11,
            symbol: "Na",
        },
        Element {
            id: ELEMENT_CL,
            atomic_number: 17,
            symbol: "Cl",
        },
        Element {
            id: ELEMENT_CA,
            atomic_number: 20,
            symbol: "Ca",
        },
    ]
}

fn curated_species() -> Vec<Species> {
    vec![
        aqueous(
            H2O_L,
            "H2O(l)",
            &[(ELEMENT_H, 2), (ELEMENT_O, 1)],
            0,
            -237_130.0,
        ),
        aqueous(H_PLUS, "H+", &[(ELEMENT_H, 1)], 1, 0.0),
        aqueous(
            OH_MINUS,
            "OH-",
            &[(ELEMENT_O, 1), (ELEMENT_H, 1)],
            -1,
            -157_240.0,
        ),
        aqueous(NA_PLUS, "Na+", &[(ELEMENT_NA, 1)], 1, -261_900.0),
        aqueous(CL_MINUS, "Cl-", &[(ELEMENT_CL, 1)], -1, -131_200.0),
        aqueous(CA_2_PLUS, "Ca2+", &[(ELEMENT_CA, 1)], 2, -553_600.0),
        aqueous(
            CO3_2_MINUS,
            "CO3--",
            &[(ELEMENT_C, 1), (ELEMENT_O, 3)],
            -2,
            -527_900.0,
        ),
        aqueous(
            HCO3_MINUS,
            "HCO3-",
            &[(ELEMENT_H, 1), (ELEMENT_C, 1), (ELEMENT_O, 3)],
            -1,
            -586_900.0,
        ),
        aqueous(
            CO2_AQ,
            "CO2(aq)",
            &[(ELEMENT_C, 1), (ELEMENT_O, 2)],
            0,
            -386_000.0,
        ),
        solid(
            CACO3_S,
            "CaCO3(s)",
            &[(ELEMENT_CA, 1), (ELEMENT_C, 1), (ELEMENT_O, 3)],
            0,
            -1_128_800.0,
        ),
    ]
}

fn aqueous(
    id: SpeciesId,
    symbol: &'static str,
    composition: &[(ElementId, u16)],
    charge_number: i8,
    standard_gibbs_energy_joule_per_mol_298_15: f64,
) -> Species {
    species(
        id,
        symbol,
        composition,
        charge_number,
        PhaseKind::Aqueous,
        standard_gibbs_energy_joule_per_mol_298_15,
    )
}

fn solid(
    id: SpeciesId,
    symbol: &'static str,
    composition: &[(ElementId, u16)],
    charge_number: i8,
    standard_gibbs_energy_joule_per_mol_298_15: f64,
) -> Species {
    species(
        id,
        symbol,
        composition,
        charge_number,
        PhaseKind::Solid,
        standard_gibbs_energy_joule_per_mol_298_15,
    )
}

fn species(
    id: SpeciesId,
    symbol: &'static str,
    composition: &[(ElementId, u16)],
    charge_number: i8,
    phase: PhaseKind,
    standard_gibbs_energy_joule_per_mol_298_15: f64,
) -> Species {
    Species {
        id,
        symbol,
        composition: composition.iter().copied().collect::<BTreeMap<_, _>>(),
        charge_number,
        phase,
        thermo: StandardThermo {
            standard_gibbs_energy_joule_per_mol_298_15,
            valid_min_temperature_kelvin: VALID_MIN_TEMPERATURE_KELVIN,
            valid_max_temperature_kelvin: VALID_MAX_TEMPERATURE_KELVIN,
            provenance: THERMO_PROVENANCE,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use thermodynamics_core::{
        analyze_aqueous_equilibrium, solve_equilibrium, EquilibriumError, SpeciesAmount,
        DAVIES_MAX_IONIC_STRENGTH_MOLAL,
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
        analyze_aqueous_equilibrium(registry, result, H2O_L, H_PLUS).unwrap()
    }

    fn activity(
        summary: &thermodynamics_core::AqueousEquilibriumSummary,
        species_id: SpeciesId,
    ) -> f64 {
        summary
            .aqueous_species
            .iter()
            .find(|species| species.species_id == species_id)
            .map(|species| species.activity)
            .unwrap_or_default()
    }

    fn expected_equilibrium_constant(
        registry: &SpeciesRegistry,
        products: &[(SpeciesId, f64)],
        reactants: &[(SpeciesId, f64)],
    ) -> f64 {
        const GAS_CONSTANT_JOULE_PER_MOL_KELVIN: f64 = 8.314_462_618_153_24;
        const REFERENCE_TEMPERATURE_KELVIN: f64 = 298.15;

        let product_gibbs: f64 = products
            .iter()
            .map(|(species_id, coefficient)| {
                coefficient
                    * registry
                        .species(*species_id)
                        .unwrap()
                        .thermo
                        .standard_gibbs_energy_joule_per_mol_298_15
            })
            .sum();
        let reactant_gibbs: f64 = reactants
            .iter()
            .map(|(species_id, coefficient)| {
                coefficient
                    * registry
                        .species(*species_id)
                        .unwrap()
                        .thermo
                        .standard_gibbs_energy_joule_per_mol_298_15
            })
            .sum();

        let reaction_gibbs_joule_per_mol = product_gibbs - reactant_gibbs;
        (-reaction_gibbs_joule_per_mol
            / (GAS_CONSTANT_JOULE_PER_MOL_KELVIN * REFERENCE_TEMPERATURE_KELVIN))
            .exp()
    }

    fn assert_relative_error(actual: f64, expected: f64, tolerance: f64) {
        let relative_error = ((actual / expected) - 1.0).abs();
        assert!(
            relative_error <= tolerance,
            "actual={actual:e}, expected={expected:e}, relative_error={relative_error:e}, tolerance={tolerance:e}"
        );
    }

    #[test]
    fn curated_registry_has_complete_species_data() {
        let registry = curated_registry().unwrap();

        for species_id in curated_species_ids() {
            let species = registry.species(species_id).unwrap();
            assert!(!species.composition.is_empty());
            assert!(species
                .thermo
                .standard_gibbs_energy_joule_per_mol_298_15
                .is_finite());
            assert!(!species.thermo.provenance.is_empty());
        }
    }

    #[test]
    fn pure_water_equilibrium_is_near_neutral_at_298_kelvin() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(298.15, 101_325.0, vec![amount(H2O_L, 55.5)]);

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);
        let h_plus = summary
            .aqueous_species
            .iter()
            .find(|species| species.species_id == H_PLUS)
            .unwrap();
        let oh_minus = summary
            .aqueous_species
            .iter()
            .find(|species| species.species_id == OH_MINUS)
            .unwrap();

        assert!(
            (h_plus.molality_mol_per_kg_water / oh_minus.molality_mol_per_kg_water - 1.0).abs()
                < 0.05
        );
        assert!((6.0..8.0).contains(&summary.ph));
        assert!((1.0e-8..1.0e-6).contains(&h_plus.molality_mol_per_kg_water));
    }

    #[test]
    fn pure_water_matches_reference_ion_product_from_standard_gibbs_data() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(298.15, 101_325.0, vec![amount(H2O_L, 55.5)]);

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);
        let actual_ion_product = activity(&summary, H_PLUS) * activity(&summary, OH_MINUS);
        let expected_ion_product = expected_equilibrium_constant(
            &registry,
            &[(H_PLUS, 1.0), (OH_MINUS, 1.0)],
            &[(H2O_L, 1.0)],
        );

        assert_relative_error(actual_ion_product, expected_ion_product, 0.05);
    }

    #[test]
    fn hcl_naoh_input_preserves_charge_and_elements() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(H_PLUS, 1.0e-3),
                amount(CL_MINUS, 1.0e-3),
                amount(NA_PLUS, 1.0e-3),
                amount(OH_MINUS, 1.0e-3),
            ],
        );

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);

        assert!(result.residuals.max_element_balance_residual_mol < 1.0e-8);
        assert!(result.residuals.charge_balance_residual_mol.abs() < 1.0e-8);
        assert!(result_amount(&result, NA_PLUS) > 9.0e-4);
        assert!(result_amount(&result, CL_MINUS) > 9.0e-4);
        assert!(summary.ionic_strength_molal > 1.0e-3);
    }

    #[test]
    fn calcium_carbonate_precipitates_when_supersaturated() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(CA_2_PLUS, 1.0e-3),
                amount(CO3_2_MINUS, 1.0e-3),
            ],
        );

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);

        let precipitated_caco3_mol = summary
            .solid_species_amounts_mol
            .iter()
            .find(|amount| amount.species_id == CACO3_S)
            .map(|amount| amount.amount_mol)
            .unwrap_or_default();

        assert!(precipitated_caco3_mol > 5.0e-4);
    }

    #[test]
    fn calcite_solubility_product_matches_standard_gibbs_data_when_solid_is_present() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(CA_2_PLUS, 1.0e-3),
                amount(CO3_2_MINUS, 1.0e-3),
            ],
        );

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);
        let actual_solubility_product =
            activity(&summary, CA_2_PLUS) * activity(&summary, CO3_2_MINUS);
        let expected_solubility_product = expected_equilibrium_constant(
            &registry,
            &[(CA_2_PLUS, 1.0), (CO3_2_MINUS, 1.0)],
            &[(CACO3_S, 1.0)],
        );

        assert_relative_error(actual_solubility_product, expected_solubility_product, 0.20);
    }

    #[test]
    fn calcium_carbonate_does_not_precipitate_when_undersaturated() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(CA_2_PLUS, 1.0e-8),
                amount(CO3_2_MINUS, 1.0e-8),
            ],
        );

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);

        assert!(result_amount(&result, CACO3_S) < 1.0e-10);
        assert!(summary.ionic_strength_molal < 1.0e-5);
    }

    #[test]
    fn aqueous_summary_reports_ph_ionic_strength_and_activities() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(NA_PLUS, 1.0e-3),
                amount(CL_MINUS, 1.0e-3),
            ],
        );

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);
        let sodium = summary
            .aqueous_species
            .iter()
            .find(|species| species.species_id == NA_PLUS)
            .unwrap();

        assert!(summary.solvent_water_mass_kg > 0.99);
        assert!((1.0e-3..2.0e-3).contains(&summary.ionic_strength_molal));
        assert!((6.0..8.0).contains(&summary.ph));
        assert!(sodium.activity > 0.0);
        assert!(sodium.activity_coefficient < 1.0);
    }

    #[test]
    fn carbonate_acid_equilibria_match_standard_gibbs_data() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(CO2_AQ, 1.0e-3),
                amount(NA_PLUS, 1.0e-3),
                amount(HCO3_MINUS, 1.0e-3),
            ],
        );

        let result = solve_equilibrium(&registry, &problem).unwrap();
        let summary = aqueous_summary(&registry, &result);
        let first_dissociation = activity(&summary, H_PLUS) * activity(&summary, HCO3_MINUS)
            / activity(&summary, CO2_AQ);
        let expected_first_dissociation = expected_equilibrium_constant(
            &registry,
            &[(H_PLUS, 1.0), (HCO3_MINUS, 1.0)],
            &[(CO2_AQ, 1.0), (H2O_L, 1.0)],
        );
        let second_dissociation = activity(&summary, H_PLUS) * activity(&summary, CO3_2_MINUS)
            / activity(&summary, HCO3_MINUS);
        let expected_second_dissociation = expected_equilibrium_constant(
            &registry,
            &[(H_PLUS, 1.0), (CO3_2_MINUS, 1.0)],
            &[(HCO3_MINUS, 1.0)],
        );

        assert_relative_error(first_dissociation, expected_first_dissociation, 0.25);
        assert_relative_error(second_dissociation, expected_second_dissociation, 0.25);
    }

    #[test]
    fn non_neutral_input_is_rejected() {
        let registry = curated_registry().unwrap();
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![amount(H2O_L, 55.5), amount(H_PLUS, 1.0e-3)],
        );

        assert!(matches!(
            solve_equilibrium(&registry, &problem),
            Err(EquilibriumError::NonNeutralCharge { .. })
        ));
    }

    #[test]
    fn davies_range_is_enforced() {
        let registry = curated_registry().unwrap();
        let excessive_ca_mol = DAVIES_MAX_IONIC_STRENGTH_MOLAL * 55.5 * 0.018_015_28;
        let problem = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(CA_2_PLUS, excessive_ca_mol),
                amount(CO3_2_MINUS, excessive_ca_mol),
            ],
        );

        assert!(matches!(
            solve_equilibrium(&registry, &problem),
            Err(EquilibriumError::DaviesModelOutOfRange { .. })
        ));
    }

    #[test]
    fn deterministic_result_does_not_depend_on_input_order() {
        let registry = curated_registry().unwrap();
        let first = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(H2O_L, 55.5),
                amount(CA_2_PLUS, 1.0e-4),
                amount(CO3_2_MINUS, 1.0e-4),
            ],
        );
        let second = equilibrium_problem(
            298.15,
            101_325.0,
            vec![
                amount(CO3_2_MINUS, 1.0e-4),
                amount(H2O_L, 55.5),
                amount(CA_2_PLUS, 1.0e-4),
            ],
        );

        let first_result = solve_equilibrium(&registry, &first).unwrap();
        let second_result = solve_equilibrium(&registry, &second).unwrap();

        assert_eq!(
            first_result.species_amounts_mol.len(),
            second_result.species_amounts_mol.len()
        );
        for (left, right) in first_result
            .species_amounts_mol
            .iter()
            .zip(second_result.species_amounts_mol.iter())
        {
            assert_eq!(left.species_id, right.species_id);
            assert!((left.amount_mol - right.amount_mol).abs() < 1.0e-12);
        }
    }
}
