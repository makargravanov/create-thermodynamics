use crate::activity::{davies_log10_gamma, DAVIES_MAX_IONIC_STRENGTH_MOLAL};
use crate::chemistry::{ElementId, PhaseKind, Species, SpeciesAmount, SpeciesId};
use crate::registry::SpeciesRegistry;
use std::collections::{BTreeMap, BTreeSet};

const GAS_CONSTANT_JOULE_PER_MOL_KELVIN: f64 = 8.314_462_618_153_24;
const WATER_MOLAR_MASS_KILOGRAM_PER_MOL: f64 = 0.018_015_28;
const MIN_ACTIVITY: f64 = 1.0e-300;
const MIN_AMOUNT_MOL: f64 = 1.0e-30;
const BALANCE_TOLERANCE_MOL: f64 = 1.0e-8;
const CHARGE_TOLERANCE_MOL: f64 = 1.0e-8;
const STEP_GRADIENT_TOLERANCE_JOULE_PER_MOL: f64 = 1.0;
const FINAL_GRADIENT_TOLERANCE_JOULE_PER_MOL: f64 = 100.0;
const MAX_ITERATIONS: usize = 4_000;

#[derive(Debug, Clone)]
pub struct EquilibriumProblem {
    pub temperature_kelvin: f64,
    pub pressure_pascal: f64,
    pub initial_species_amounts_mol: Vec<SpeciesAmount>,
    pub candidate_species: Vec<SpeciesId>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumResult {
    pub species_amounts_mol: Vec<SpeciesAmount>,
    pub diagnostics: Vec<EquilibriumDiagnostic>,
    pub residuals: EquilibriumResiduals,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumResiduals {
    pub max_element_balance_residual_mol: f64,
    pub charge_balance_residual_mol: f64,
    pub max_projected_gradient_joule_per_mol: f64,
    pub iterations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EquilibriumDiagnostic {
    DaviesActivityModel {
        ionic_strength_molal: f64,
        max_valid_ionic_strength_molal: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub enum EquilibriumError {
    InvalidTemperatureKelvin(f64),
    InvalidPressurePascal(f64),
    NegativeAmount {
        species_id: SpeciesId,
        amount_mol: f64,
    },
    NonNeutralCharge {
        charge_balance_mol: f64,
    },
    MissingSpeciesData(SpeciesId),
    UnsupportedPhase {
        species_id: SpeciesId,
        phase: PhaseKind,
    },
    UnsupportedTemperatureRange {
        species_id: SpeciesId,
        temperature_kelvin: f64,
        valid_min_temperature_kelvin: f64,
        valid_max_temperature_kelvin: f64,
    },
    DaviesModelOutOfRange {
        ionic_strength_molal: f64,
        max_valid_ionic_strength_molal: f64,
    },
    InfeasibleInitialState {
        max_element_balance_residual_mol: f64,
        charge_balance_residual_mol: f64,
    },
    NonConvergence {
        iterations: usize,
        max_projected_gradient_joule_per_mol: f64,
    },
}

pub fn solve_equilibrium(
    registry: &SpeciesRegistry,
    problem: &EquilibriumProblem,
) -> Result<EquilibriumResult, EquilibriumError> {
    if !problem.temperature_kelvin.is_finite() || problem.temperature_kelvin <= 0.0 {
        return Err(EquilibriumError::InvalidTemperatureKelvin(
            problem.temperature_kelvin,
        ));
    }
    if !problem.pressure_pascal.is_finite() || problem.pressure_pascal <= 0.0 {
        return Err(EquilibriumError::InvalidPressurePascal(
            problem.pressure_pascal,
        ));
    }

    let candidate_ids = sorted_unique_species(&problem.candidate_species);
    let species = load_species(registry, &candidate_ids, problem.temperature_kelvin)?;
    let initial_amount_by_species =
        normalized_initial_amounts(&problem.initial_species_amounts_mol)?;

    for species_id in initial_amount_by_species.keys() {
        if !candidate_ids.contains(species_id) {
            return Err(EquilibriumError::MissingSpeciesData(*species_id));
        }
    }

    let mut initial_amounts = vec![0.0; candidate_ids.len()];
    for (index, species_id) in candidate_ids.iter().enumerate() {
        initial_amounts[index] = *initial_amount_by_species.get(species_id).unwrap_or(&0.0);
    }

    let charge_balance = charge_balance_mol(&species, &initial_amounts);
    if charge_balance.abs() > CHARGE_TOLERANCE_MOL {
        return Err(EquilibriumError::NonNeutralCharge {
            charge_balance_mol: charge_balance,
        });
    }

    let element_ids = sorted_element_ids(&species);
    let constraint_matrix = build_constraint_matrix(&species, &element_ids);
    let conserved_totals = multiply_matrix_vector(&constraint_matrix, &initial_amounts);
    let search_directions = build_search_directions(&constraint_matrix);
    let mut amounts = initial_amounts;
    let mut diagnostics = Vec::new();
    let mut iterations = 0;

    for iteration in 0..MAX_ITERATIONS {
        iterations = iteration + 1;
        let potentials = chemical_potentials_joule_per_mol(
            &species,
            &amounts,
            problem.temperature_kelvin,
            &mut diagnostics,
        )?;
        let max_projected_gradient =
            feasible_descent_measure(&search_directions, &potentials, &amounts);

        if max_projected_gradient <= STEP_GRADIENT_TOLERANCE_JOULE_PER_MOL {
            break;
        }

        if !take_coordinate_root_step(
            &species,
            &mut amounts,
            &search_directions,
            problem.temperature_kelvin,
            &mut diagnostics,
        )? {
            break;
        }
    }

    let residuals = compute_residuals(
        &constraint_matrix,
        &conserved_totals,
        &species,
        &amounts,
        &search_directions,
        problem.temperature_kelvin,
    )?;

    if residuals.max_element_balance_residual_mol > BALANCE_TOLERANCE_MOL
        || residuals.charge_balance_residual_mol.abs() > CHARGE_TOLERANCE_MOL
    {
        return Err(EquilibriumError::InfeasibleInitialState {
            max_element_balance_residual_mol: residuals.max_element_balance_residual_mol,
            charge_balance_residual_mol: residuals.charge_balance_residual_mol,
        });
    }

    if residuals.max_projected_gradient_joule_per_mol > FINAL_GRADIENT_TOLERANCE_JOULE_PER_MOL {
        return Err(EquilibriumError::NonConvergence {
            iterations: residuals.iterations,
            max_projected_gradient_joule_per_mol: residuals.max_projected_gradient_joule_per_mol,
        });
    }

    Ok(EquilibriumResult {
        species_amounts_mol: candidate_ids
            .into_iter()
            .zip(amounts)
            .map(|(species_id, amount_mol)| SpeciesAmount {
                species_id,
                amount_mol: if amount_mol.abs() < MIN_AMOUNT_MOL {
                    0.0
                } else {
                    amount_mol
                },
            })
            .collect(),
        diagnostics,
        residuals: EquilibriumResiduals {
            iterations,
            ..residuals
        },
    })
}

fn sorted_unique_species(species: &[SpeciesId]) -> Vec<SpeciesId> {
    species
        .iter()
        .copied()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn load_species(
    registry: &SpeciesRegistry,
    ids: &[SpeciesId],
    temperature_kelvin: f64,
) -> Result<Vec<Species>, EquilibriumError> {
    ids.iter()
        .map(|species_id| {
            let species = registry
                .species(*species_id)
                .ok_or(EquilibriumError::MissingSpeciesData(*species_id))?;
            match species.phase {
                PhaseKind::Aqueous | PhaseKind::Solid => {}
                PhaseKind::Gas => {
                    return Err(EquilibriumError::UnsupportedPhase {
                        species_id: *species_id,
                        phase: species.phase,
                    });
                }
            }
            if temperature_kelvin < species.thermo.valid_min_temperature_kelvin
                || temperature_kelvin > species.thermo.valid_max_temperature_kelvin
            {
                return Err(EquilibriumError::UnsupportedTemperatureRange {
                    species_id: *species_id,
                    temperature_kelvin,
                    valid_min_temperature_kelvin: species.thermo.valid_min_temperature_kelvin,
                    valid_max_temperature_kelvin: species.thermo.valid_max_temperature_kelvin,
                });
            }
            Ok(species.clone())
        })
        .collect()
}

fn normalized_initial_amounts(
    amounts: &[SpeciesAmount],
) -> Result<BTreeMap<SpeciesId, f64>, EquilibriumError> {
    let mut by_species = BTreeMap::new();
    for amount in amounts {
        if !amount.amount_mol.is_finite() || amount.amount_mol < 0.0 {
            return Err(EquilibriumError::NegativeAmount {
                species_id: amount.species_id,
                amount_mol: amount.amount_mol,
            });
        }
        *by_species.entry(amount.species_id).or_insert(0.0) += amount.amount_mol;
    }
    Ok(by_species)
}

fn sorted_element_ids(species: &[Species]) -> Vec<ElementId> {
    species
        .iter()
        .flat_map(|species_record| species_record.composition.keys().copied())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn build_constraint_matrix(species: &[Species], element_ids: &[ElementId]) -> Vec<Vec<f64>> {
    let mut rows = Vec::new();
    for element_id in element_ids {
        rows.push(
            species
                .iter()
                .map(|species_record| species_record.element_count(*element_id))
                .collect(),
        );
    }
    rows.push(
        species
            .iter()
            .map(|species_record| f64::from(species_record.charge_number))
            .collect(),
    );
    rows
}

fn multiply_matrix_vector(matrix: &[Vec<f64>], vector: &[f64]) -> Vec<f64> {
    matrix
        .iter()
        .map(|row| row.iter().zip(vector).map(|(a, b)| a * b).sum())
        .collect()
}

fn charge_balance_mol(species: &[Species], amounts: &[f64]) -> f64 {
    species
        .iter()
        .zip(amounts)
        .map(|(species_record, amount)| f64::from(species_record.charge_number) * amount)
        .sum()
}

fn chemical_potentials_joule_per_mol(
    species: &[Species],
    amounts: &[f64],
    temperature_kelvin: f64,
    diagnostics: &mut Vec<EquilibriumDiagnostic>,
) -> Result<Vec<f64>, EquilibriumError> {
    let solvent_kg = solvent_water_kilograms(species, amounts);
    let ionic_strength = ionic_strength_molal(species, amounts, solvent_kg);
    if ionic_strength > DAVIES_MAX_IONIC_STRENGTH_MOLAL {
        return Err(EquilibriumError::DaviesModelOutOfRange {
            ionic_strength_molal: ionic_strength,
            max_valid_ionic_strength_molal: DAVIES_MAX_IONIC_STRENGTH_MOLAL,
        });
    }
    push_unique_davies_diagnostic(diagnostics, ionic_strength);

    Ok(species
        .iter()
        .zip(amounts)
        .map(|(species_record, amount_mol)| {
            let standard_gibbs = species_record
                .thermo
                .standard_gibbs_energy_joule_per_mol_298_15;
            match species_record.phase {
                PhaseKind::Solid => standard_gibbs,
                PhaseKind::Gas => standard_gibbs,
                PhaseKind::Aqueous if is_water_solvent(species_record) => standard_gibbs,
                PhaseKind::Aqueous => {
                    let molality =
                        (*amount_mol).max(MIN_AMOUNT_MOL) / solvent_kg.max(MIN_AMOUNT_MOL);
                    let log10_gamma =
                        davies_log10_gamma(species_record.charge_number, ionic_strength)
                            .expect("ionic strength range checked before gamma calculation");
                    let activity = (molality * 10.0_f64.powf(log10_gamma)).max(MIN_ACTIVITY);
                    standard_gibbs
                        + GAS_CONSTANT_JOULE_PER_MOL_KELVIN * temperature_kelvin * activity.ln()
                }
            }
        })
        .collect())
}

fn is_water_solvent(species: &Species) -> bool {
    species.symbol == "H2O(l)" && species.phase == PhaseKind::Aqueous
}

fn solvent_water_kilograms(species: &[Species], amounts: &[f64]) -> f64 {
    species
        .iter()
        .zip(amounts)
        .find(|(species_record, _)| is_water_solvent(species_record))
        .map(|(_, amount_mol)| amount_mol.max(MIN_AMOUNT_MOL) * WATER_MOLAR_MASS_KILOGRAM_PER_MOL)
        .unwrap_or(MIN_AMOUNT_MOL)
}

fn ionic_strength_molal(species: &[Species], amounts: &[f64], solvent_kg: f64) -> f64 {
    0.5 * species
        .iter()
        .zip(amounts)
        .filter(|(species_record, _)| species_record.phase == PhaseKind::Aqueous)
        .map(|(species_record, amount_mol)| {
            let charge = f64::from(species_record.charge_number);
            (amount_mol.max(0.0) / solvent_kg.max(MIN_AMOUNT_MOL)) * charge * charge
        })
        .sum::<f64>()
}

fn push_unique_davies_diagnostic(
    diagnostics: &mut Vec<EquilibriumDiagnostic>,
    ionic_strength_molal: f64,
) {
    if diagnostics.iter().any(|diagnostic| {
        matches!(
            diagnostic,
            EquilibriumDiagnostic::DaviesActivityModel { .. }
        )
    }) {
        return;
    }

    diagnostics.push(EquilibriumDiagnostic::DaviesActivityModel {
        ionic_strength_molal,
        max_valid_ionic_strength_molal: DAVIES_MAX_IONIC_STRENGTH_MOLAL,
    });
}

fn feasible_descent_measure(
    search_directions: &[Vec<f64>],
    potentials: &[f64],
    amounts: &[f64],
) -> f64 {
    let mut max_descent = 0.0;
    for basis in search_directions {
        for sign in [-1.0, 1.0] {
            let direction: Vec<f64> = basis.iter().map(|value| sign * value).collect();
            if !direction_is_feasible_at_boundary(amounts, &direction) {
                continue;
            }
            let directional_derivative: f64 = potentials
                .iter()
                .zip(direction.iter())
                .map(|(potential, delta)| potential * delta)
                .sum();
            if directional_derivative < 0.0 {
                let descent = -directional_derivative;
                if descent > max_descent {
                    max_descent = descent;
                }
            }
        }
    }
    max_descent
}

fn direction_is_feasible_at_boundary(amounts: &[f64], direction: &[f64]) -> bool {
    amounts
        .iter()
        .zip(direction)
        .all(|(amount, delta)| *amount > MIN_AMOUNT_MOL || *delta >= 0.0)
}

fn take_coordinate_root_step(
    species: &[Species],
    amounts: &mut [f64],
    search_directions: &[Vec<f64>],
    temperature_kelvin: f64,
    diagnostics: &mut Vec<EquilibriumDiagnostic>,
) -> Result<bool, EquilibriumError> {
    if search_directions.is_empty() {
        return Ok(false);
    }

    let current_potentials =
        chemical_potentials_joule_per_mol(species, amounts, temperature_kelvin, diagnostics)?;
    let mut coordinate_order: Vec<usize> = (0..search_directions.len()).collect();
    coordinate_order.sort_by(|left, right| {
        directional_derivative(&current_potentials, &search_directions[*right])
            .abs()
            .total_cmp(
                &directional_derivative(&current_potentials, &search_directions[*left]).abs(),
            )
    });

    for coordinate in coordinate_order {
        let direction = &search_directions[coordinate];
        if let Some(candidate_amounts) =
            solve_directional_root(species, amounts, direction, temperature_kelvin)?
        {
            let changed = amounts
                .iter()
                .zip(candidate_amounts.iter())
                .any(|(left, right)| (left - right).abs() > MIN_AMOUNT_MOL);
            if !changed {
                continue;
            }
            amounts.copy_from_slice(&candidate_amounts);
            let _ = chemical_potentials_joule_per_mol(
                species,
                amounts,
                temperature_kelvin,
                diagnostics,
            )?;
            return Ok(true);
        }
    }

    Ok(false)
}

fn solve_directional_root(
    species: &[Species],
    amounts: &[f64],
    direction: &[f64],
    temperature_kelvin: f64,
) -> Result<Option<Vec<f64>>, EquilibriumError> {
    let (lower_bound, upper_bound) = feasible_step_bounds(amounts, direction);
    if lower_bound >= upper_bound {
        return Ok(None);
    }
    let lower_bound =
        valid_step_towards(species, amounts, direction, lower_bound, temperature_kelvin)?;
    let upper_bound =
        valid_step_towards(species, amounts, direction, upper_bound, temperature_kelvin)?;
    if lower_bound >= upper_bound {
        return Ok(None);
    }

    let current_derivative =
        directional_derivative_at_step(species, amounts, direction, 0.0, temperature_kelvin)?;
    if current_derivative.abs() <= STEP_GRADIENT_TOLERANCE_JOULE_PER_MOL {
        return Ok(None);
    }

    let width = upper_bound - lower_bound;
    let epsilon = (width.abs() * 1.0e-12).max(1.0e-18);
    let lower = if lower_bound < 0.0 {
        (lower_bound + epsilon).min(0.0)
    } else {
        lower_bound
    };
    let upper = if upper_bound > 0.0 {
        (upper_bound - epsilon).max(0.0)
    } else {
        upper_bound
    };

    let lower_derivative =
        directional_derivative_at_step(species, amounts, direction, lower, temperature_kelvin)?;
    let upper_derivative =
        directional_derivative_at_step(species, amounts, direction, upper, temperature_kelvin)?;

    let target_step = if lower_derivative <= 0.0 && upper_derivative >= 0.0 {
        bisect_directional_derivative(
            species,
            amounts,
            direction,
            lower,
            upper,
            temperature_kelvin,
        )?
    } else if current_derivative > 0.0 && lower < 0.0 {
        lower
    } else if current_derivative < 0.0 && upper > 0.0 {
        upper
    } else {
        return Ok(None);
    };

    Ok(Some(amounts_at_step(amounts, direction, target_step)))
}

fn feasible_step_bounds(amounts: &[f64], direction: &[f64]) -> (f64, f64) {
    let mut lower_bound = f64::NEG_INFINITY;
    let mut upper_bound = f64::INFINITY;
    for (amount, delta) in amounts.iter().zip(direction) {
        if *delta < 0.0 {
            upper_bound = upper_bound.min(amount.max(0.0) / -delta);
        } else if *delta > 0.0 {
            lower_bound = lower_bound.max(-amount.max(0.0) / delta);
        }
    }

    if !lower_bound.is_finite() {
        lower_bound = -1.0;
    }
    if !upper_bound.is_finite() {
        upper_bound = 1.0;
    }

    (lower_bound, upper_bound)
}

fn valid_step_towards(
    species: &[Species],
    amounts: &[f64],
    direction: &[f64],
    target_step: f64,
    temperature_kelvin: f64,
) -> Result<f64, EquilibriumError> {
    if target_step == 0.0 {
        return Ok(0.0);
    }
    if directional_derivative_at_step(species, amounts, direction, target_step, temperature_kelvin)
        .is_ok()
    {
        return Ok(target_step);
    }

    let mut valid = 0.0;
    let mut invalid = target_step;
    for _ in 0..96 {
        let middle = (valid + invalid) * 0.5;
        match directional_derivative_at_step(
            species,
            amounts,
            direction,
            middle,
            temperature_kelvin,
        ) {
            Ok(_) => valid = middle,
            Err(EquilibriumError::DaviesModelOutOfRange { .. }) => invalid = middle,
            Err(error) => return Err(error),
        }
    }
    Ok(valid)
}

fn directional_derivative_at_step(
    species: &[Species],
    amounts: &[f64],
    direction: &[f64],
    step: f64,
    temperature_kelvin: f64,
) -> Result<f64, EquilibriumError> {
    let trial_amounts = amounts_at_step(amounts, direction, step);
    let mut diagnostics = Vec::new();
    let potentials = chemical_potentials_joule_per_mol(
        species,
        &trial_amounts,
        temperature_kelvin,
        &mut diagnostics,
    )?;
    Ok(directional_derivative(&potentials, direction))
}

fn directional_derivative(potentials: &[f64], direction: &[f64]) -> f64 {
    potentials
        .iter()
        .zip(direction.iter())
        .map(|(potential, delta)| potential * delta)
        .sum()
}

fn bisect_directional_derivative(
    species: &[Species],
    amounts: &[f64],
    direction: &[f64],
    mut lower: f64,
    mut upper: f64,
    temperature_kelvin: f64,
) -> Result<f64, EquilibriumError> {
    for _ in 0..96 {
        let middle = (lower + upper) * 0.5;
        let middle_derivative = directional_derivative_at_step(
            species,
            amounts,
            direction,
            middle,
            temperature_kelvin,
        )?;
        if middle_derivative.abs() <= STEP_GRADIENT_TOLERANCE_JOULE_PER_MOL {
            return Ok(middle);
        }
        if middle_derivative < 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    Ok((lower + upper) * 0.5)
}

fn amounts_at_step(amounts: &[f64], direction: &[f64], step: f64) -> Vec<f64> {
    amounts
        .iter()
        .zip(direction)
        .map(|(amount, delta)| {
            let value = amount + step * delta;
            if value.abs() < MIN_AMOUNT_MOL {
                0.0
            } else {
                value.max(0.0)
            }
        })
        .collect()
}

fn compute_residuals(
    constraint_matrix: &[Vec<f64>],
    conserved_totals: &[f64],
    species: &[Species],
    amounts: &[f64],
    search_directions: &[Vec<f64>],
    temperature_kelvin: f64,
) -> Result<EquilibriumResiduals, EquilibriumError> {
    let actual_totals = multiply_matrix_vector(constraint_matrix, amounts);
    let mut element_residual = 0.0;
    for (index, (actual, expected)) in actual_totals.iter().zip(conserved_totals).enumerate() {
        let residual = (actual - expected).abs();
        if index + 1 == actual_totals.len() {
            continue;
        }
        if residual > element_residual {
            element_residual = residual;
        }
    }

    let mut diagnostics = Vec::new();
    let potentials =
        chemical_potentials_joule_per_mol(species, amounts, temperature_kelvin, &mut diagnostics)?;
    Ok(EquilibriumResiduals {
        max_element_balance_residual_mol: element_residual,
        charge_balance_residual_mol: actual_totals.last().copied().unwrap_or_default(),
        max_projected_gradient_joule_per_mol: feasible_descent_measure(
            search_directions,
            &potentials,
            amounts,
        ),
        iterations: 0,
    })
}

fn build_search_directions(constraint_matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut directions = Vec::new();
    for basis_direction in nullspace(constraint_matrix) {
        push_unique_direction(&mut directions, basis_direction);
    }

    let species_count = constraint_matrix
        .first()
        .map(|row| row.len())
        .unwrap_or_default();
    for support_size in 2..=4 {
        let mut support = Vec::with_capacity(support_size);
        enumerate_supports(
            constraint_matrix,
            species_count,
            support_size,
            0,
            &mut support,
            &mut directions,
        );
    }

    directions
}

fn enumerate_supports(
    constraint_matrix: &[Vec<f64>],
    species_count: usize,
    support_size: usize,
    start_index: usize,
    support: &mut Vec<usize>,
    directions: &mut Vec<Vec<f64>>,
) {
    if support.len() == support_size {
        let mut coefficients = vec![0_i8; support_size];
        enumerate_coefficients(
            constraint_matrix,
            species_count,
            support,
            &mut coefficients,
            0,
            directions,
        );
        return;
    }

    for species_index in start_index..species_count {
        support.push(species_index);
        enumerate_supports(
            constraint_matrix,
            species_count,
            support_size,
            species_index + 1,
            support,
            directions,
        );
        support.pop();
    }
}

fn enumerate_coefficients(
    constraint_matrix: &[Vec<f64>],
    species_count: usize,
    support: &[usize],
    coefficients: &mut [i8],
    coefficient_index: usize,
    directions: &mut Vec<Vec<f64>>,
) {
    const COEFFICIENTS: [i8; 4] = [-2, -1, 1, 2];

    if coefficient_index == coefficients.len() {
        let mut direction = vec![0.0; species_count];
        for (species_index, coefficient) in support.iter().zip(coefficients.iter()) {
            direction[*species_index] = f64::from(*coefficient);
        }
        if constraint_matrix.iter().all(|row| {
            row.iter()
                .zip(direction.iter())
                .map(|(left, right)| left * right)
                .sum::<f64>()
                .abs()
                < 1.0e-12
        }) {
            push_unique_direction(directions, direction);
        }
        return;
    }

    for coefficient in COEFFICIENTS {
        coefficients[coefficient_index] = coefficient;
        enumerate_coefficients(
            constraint_matrix,
            species_count,
            support,
            coefficients,
            coefficient_index + 1,
            directions,
        );
    }
}

fn push_unique_direction(directions: &mut Vec<Vec<f64>>, mut direction: Vec<f64>) {
    if direction.iter().all(|value| value.abs() < 1.0e-12) {
        return;
    }

    let max_abs = direction
        .iter()
        .fold(0.0_f64, |acc, value| acc.max(value.abs()));
    if max_abs == 0.0 {
        return;
    }
    for value in &mut direction {
        *value /= max_abs;
    }

    if let Some(first_non_zero) = direction.iter().find(|value| value.abs() >= 1.0e-12) {
        if *first_non_zero < 0.0 {
            for value in &mut direction {
                *value = -*value;
            }
        }
    }

    if directions.iter().any(|existing| {
        existing
            .iter()
            .zip(direction.iter())
            .all(|(left, right)| (left - right).abs() < 1.0e-12)
    }) {
        return;
    }

    directions.push(direction);
}

fn nullspace(matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
    if matrix.is_empty() || matrix[0].is_empty() {
        return Vec::new();
    }

    let row_count = matrix.len();
    let column_count = matrix[0].len();
    let mut rref = matrix.to_vec();
    let mut pivot_columns = Vec::new();
    let mut pivot_row = 0;

    for column in 0..column_count {
        let Some(best_row) = (pivot_row..row_count)
            .max_by(|a, b| rref[*a][column].abs().total_cmp(&rref[*b][column].abs()))
        else {
            break;
        };

        if rref[best_row][column].abs() < 1.0e-12 {
            continue;
        }

        rref.swap(pivot_row, best_row);
        let pivot = rref[pivot_row][column];
        for value in &mut rref[pivot_row] {
            *value /= pivot;
        }

        for row in 0..row_count {
            if row == pivot_row {
                continue;
            }
            let factor = rref[row][column];
            if factor == 0.0 {
                continue;
            }
            for col in column..column_count {
                rref[row][col] -= factor * rref[pivot_row][col];
            }
        }

        pivot_columns.push(column);
        pivot_row += 1;
        if pivot_row == row_count {
            break;
        }
    }

    let pivot_set: BTreeSet<usize> = pivot_columns.iter().copied().collect();
    let free_columns: Vec<usize> = (0..column_count)
        .filter(|column| !pivot_set.contains(column))
        .collect();

    let mut basis = Vec::new();
    for free_column in free_columns {
        let mut vector = vec![0.0; column_count];
        vector[free_column] = 1.0;
        for (row, pivot_column) in pivot_columns.iter().enumerate() {
            vector[*pivot_column] = -rref[row][free_column];
        }
        let norm = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
        if norm > 0.0 {
            for value in &mut vector {
                *value /= norm;
            }
        }
        basis.push(vector);
    }
    basis
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::{PhaseKind, StandardThermo};
    use crate::registry::SpeciesRegistry;

    const H: ElementId = ElementId(1);
    const O: ElementId = ElementId(8);

    fn test_species(id: u16, symbol: &'static str, h: u16, o: u16, charge: i8) -> Species {
        let mut composition = BTreeMap::new();
        if h > 0 {
            composition.insert(H, h);
        }
        if o > 0 {
            composition.insert(O, o);
        }
        Species {
            id: SpeciesId(id),
            symbol,
            composition,
            charge_number: charge,
            phase: PhaseKind::Aqueous,
            thermo: StandardThermo {
                standard_gibbs_energy_joule_per_mol_298_15: 0.0,
                valid_min_temperature_kelvin: 273.15,
                valid_max_temperature_kelvin: 373.15,
                provenance: "test",
            },
        }
    }

    #[test]
    fn non_neutral_input_is_rejected() {
        let registry = SpeciesRegistry::new(
            vec![
                crate::chemistry::Element {
                    id: H,
                    atomic_number: 1,
                    symbol: "H",
                },
                crate::chemistry::Element {
                    id: O,
                    atomic_number: 8,
                    symbol: "O",
                },
            ],
            vec![
                test_species(1, "H+", 1, 0, 1),
                test_species(2, "OH-", 1, 1, -1),
            ],
        )
        .unwrap();

        let problem = EquilibriumProblem {
            temperature_kelvin: 298.15,
            pressure_pascal: 101_325.0,
            initial_species_amounts_mol: vec![SpeciesAmount {
                species_id: SpeciesId(1),
                amount_mol: 1.0e-6,
            }],
            candidate_species: vec![SpeciesId(1), SpeciesId(2)],
        };

        assert!(matches!(
            solve_equilibrium(&registry, &problem),
            Err(EquilibriumError::NonNeutralCharge { .. })
        ));
    }
}
