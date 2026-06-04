use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::molecule::MolecularStructure;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::substance::{Substance, SubstanceId, SubstanceTagId};

const OXYGEN_ID: &str = "destroy:oxygen";
const CARBON_DIOXIDE_ID: &str = "destroy:carbon_dioxide";
const WATER_ID: &str = "destroy:water";
const COMBUSTION_PRE_EXPONENTIAL: f64 = 2.5e11;
const COMBUSTION_ACTIVATION_ENERGY_KJ_PER_MOL: f64 = 115.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CombustionFormula {
    carbon: u32,
    hydrogen: u32,
    oxygen: u32,
    unsupported: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CombustionEquation {
    multiplier: u32,
    oxygen: u32,
    carbon_dioxide: u32,
    water: u32,
}

pub(crate) fn generate_complete_combustion(fuel: &Substance) -> ChemistryResult<Option<Reaction>> {
    if fuel.charge != 0 || is_hypothetical(fuel) {
        return Ok(None);
    }
    let Some(structure) = fuel.molecular_structure.as_ref() else {
        return Ok(None);
    };
    if structure.atoms.iter().any(|atom| atom.element == "R") {
        return Ok(None);
    }
    let formula = combustion_formula(&fuel.id, structure)?;
    let Some(equation) = complete_cho_combustion_equation(formula) else {
        return Ok(None);
    };
    let fuel_coefficient = equation.multiplier;
    let enthalpy = estimate_combustion_enthalpy_kj_per_reaction(&equation);
    let mut builder =
        Reaction::builder(format!("combustion/{}/complete", stable_id_part(&fuel.id)))
            .reactant(fuel.id.clone(), fuel_coefficient, 1)
            .reactant(OXYGEN_ID, equation.oxygen, 1)
            .product(CARBON_DIOXIDE_ID, equation.carbon_dioxide)
            .reactant_phase_access(fuel.id.clone(), [MixturePhase::Gas])
            .reactant_phase_access(OXYGEN_ID, [MixturePhase::Gas])
            .product_phase(CARBON_DIOXIDE_ID, MixturePhase::Gas)
            .pre_exponential_factor(COMBUSTION_PRE_EXPONENTIAL)
            .activation_energy_kj_per_mol(COMBUSTION_ACTIVATION_ENERGY_KJ_PER_MOL)
            .enthalpy_change_kj_per_mol(enthalpy)
            .reaction_result("estimated complete combustion enthalpy", 1.0);
    if equation.water > 0 {
        builder = builder
            .product(WATER_ID, equation.water)
            .product_phase(WATER_ID, MixturePhase::Gas);
    }
    Ok(Some(builder.build()))
}

fn combustion_formula(
    substance_id: &SubstanceId,
    structure: &MolecularStructure,
) -> ChemistryResult<CombustionFormula> {
    if structure
        .atoms
        .iter()
        .any(|atom| !matches!(atom.element.as_str(), "C" | "H" | "O"))
    {
        return Ok(CombustionFormula {
            carbon: 0,
            hydrogen: 0,
            oxygen: 0,
            unsupported: true,
        });
    }
    let mut formula = CombustionFormula {
        carbon: 0,
        hydrogen: 0,
        oxygen: 0,
        unsupported: false,
    };
    for atom in &structure.atoms {
        if atom.charge.abs() > 1.0e-9 {
            return Err(ChemistryError::GenerationInvariantViolation {
                generator: "combustion".to_string(),
                substance_id: substance_id.to_string(),
                reason: format!(
                    "charged atom '{}' cannot be used as neutral fuel",
                    atom.element
                ),
            });
        }
        match atom.element.as_str() {
            "C" => formula.carbon += 1,
            "H" => formula.hydrogen += 1,
            "O" => formula.oxygen += 1,
            _ => unreachable!("unsupported elements are checked before formula counting"),
        }
    }
    Ok(formula)
}

fn complete_cho_combustion_equation(formula: CombustionFormula) -> Option<CombustionEquation> {
    if formula.unsupported || formula.carbon == 0 {
        return None;
    }
    let oxygen_quarters =
        4_i64 * formula.carbon as i64 + formula.hydrogen as i64 - 2_i64 * formula.oxygen as i64;
    if oxygen_quarters <= 0 {
        return None;
    }
    let multiplier = 4_u32 / gcd_u32(oxygen_quarters as u32, 4);
    Some(CombustionEquation {
        multiplier,
        oxygen: (oxygen_quarters as u32 * multiplier) / 4,
        carbon_dioxide: formula.carbon * multiplier,
        water: (formula.hydrogen * multiplier) / 2,
    })
}

fn estimate_combustion_enthalpy_kj_per_reaction(equation: &CombustionEquation) -> f64 {
    let carbon_dioxide = equation.carbon_dioxide as f64;
    let water = equation.water as f64;
    -393.5 * carbon_dioxide - 241.8 * water
}

fn stable_id_part(id: &SubstanceId) -> String {
    id.as_str()
        .chars()
        .map(|character| match character {
            ':' => '_',
            '/' => '_',
            other => other,
        })
        .collect()
}

fn is_hypothetical(substance: &Substance) -> bool {
    substance
        .tags
        .iter()
        .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
}

fn gcd_u32(mut left: u32, mut right: u32) -> u32 {
    while right != 0 {
        let next = left % right;
        left = right;
        right = next;
    }
    left
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;
    use crate::chemistry::substance::Substance;

    fn fuel(id: &str, code: &str) -> Substance {
        let structure = parse_frowns(code).unwrap();
        let summary = structure.summary().unwrap();
        Substance::new(
            id,
            summary.charge,
            summary.molar_mass_grams,
            700.0,
            250.0,
            80.0,
            20_000.0,
        )
        .with_molecular_structure(structure)
    }

    #[test]
    fn methane_combustion_stoichiometry_is_balanced() {
        let reaction = generate_complete_combustion(&fuel("test:methane", "C"))
            .unwrap()
            .unwrap();
        assert_term(&reaction.reactants[0], "test:methane", 1);
        assert_term(&reaction.reactants[1], OXYGEN_ID, 2);
        assert_term(&reaction.products[0], CARBON_DIOXIDE_ID, 1);
        assert_term(&reaction.products[1], WATER_ID, 2);
    }

    #[test]
    fn ethanol_combustion_stoichiometry_is_balanced() {
        let reaction = generate_complete_combustion(&fuel("test:ethanol", "CCO"))
            .unwrap()
            .unwrap();
        assert_term(&reaction.reactants[0], "test:ethanol", 1);
        assert_term(&reaction.reactants[1], OXYGEN_ID, 3);
        assert_term(&reaction.products[0], CARBON_DIOXIDE_ID, 2);
        assert_term(&reaction.products[1], WATER_ID, 3);
    }

    #[test]
    fn oxidized_or_non_carbon_substance_does_not_generate_combustion() {
        assert!(generate_complete_combustion(&fuel("test:water", "O"))
            .unwrap()
            .is_none());
    }

    fn assert_term(
        term: &crate::chemistry::reaction::StoichiometricTerm,
        substance_id: &str,
        coefficient: u32,
    ) {
        assert_eq!(term.substance_id.as_str(), substance_id);
        assert_eq!(term.coefficient, coefficient);
    }
}
