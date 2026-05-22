use crate::chemistry::{PhaseKind, SpeciesAmount, SpeciesId};
use crate::equilibrium::{
    solve_equilibrium, EquilibriumError, EquilibriumProblem, EquilibriumResult,
};
use crate::registry::SpeciesRegistry;

const GAS_CONSTANT_JOULE_PER_MOL_KELVIN: f64 = 8.314_462_618_153_24;
const MIN_PRESSURE_PASCAL: f64 = 1.0;
const MAX_PRESSURE_PASCAL: f64 = 100_000_000.0;
const PRESSURE_TOLERANCE_PASCAL: f64 = 1.0e-3;
const RELATIVE_PRESSURE_TOLERANCE: f64 = 1.0e-9;
const MIN_STABLE_GAS_AMOUNT_MOL: f64 = 1.0e-12;
const MAX_PRESSURE_ITERATIONS: usize = 96;

#[derive(Debug, Clone)]
pub struct ClosedGasEquilibriumProblem {
    pub temperature_kelvin: f64,
    pub gas_volume_cubic_meter: f64,
    pub initial_species_amounts_mol: Vec<SpeciesAmount>,
    pub candidate_species: Vec<SpeciesId>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ClosedGasEquilibriumResult {
    pub equilibrium: EquilibriumResult,
    pub pressure_pascal: f64,
    pub gas_amount_mol: f64,
    pub pressure_residual_pascal: f64,
    pub iterations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ClosedGasEquilibriumError {
    InvalidGasVolumeCubicMeter(f64),
    MissingGasCandidate,
    Equilibrium(EquilibriumError),
    PressureBracketFailure {
        lower_pressure_pascal: f64,
        lower_residual_pascal: f64,
        upper_pressure_pascal: f64,
        upper_residual_pascal: f64,
    },
    NoStableGasPhase {
        gas_amount_mol: f64,
    },
    NonConvergence {
        iterations: usize,
        pressure_residual_pascal: f64,
    },
}

pub fn solve_closed_gas_equilibrium(
    registry: &SpeciesRegistry,
    problem: &ClosedGasEquilibriumProblem,
) -> Result<ClosedGasEquilibriumResult, ClosedGasEquilibriumError> {
    if !problem.gas_volume_cubic_meter.is_finite() || problem.gas_volume_cubic_meter <= 0.0 {
        return Err(ClosedGasEquilibriumError::InvalidGasVolumeCubicMeter(
            problem.gas_volume_cubic_meter,
        ));
    }
    if !has_gas_candidate(registry, &problem.candidate_species) {
        return Err(ClosedGasEquilibriumError::MissingGasCandidate);
    }

    let mut lower = evaluate_pressure(
        registry,
        problem,
        MIN_PRESSURE_PASCAL,
        problem.gas_volume_cubic_meter,
    )?;
    let mut upper = evaluate_pressure(
        registry,
        problem,
        MAX_PRESSURE_PASCAL,
        problem.gas_volume_cubic_meter,
    )?;

    if pressure_converged(lower.pressure_pascal, lower.residual_pascal) {
        return lower.into_result(1);
    }
    if pressure_converged(upper.pressure_pascal, upper.residual_pascal) {
        return upper.into_result(1);
    }
    if lower.residual_pascal.signum() == upper.residual_pascal.signum() {
        return Err(ClosedGasEquilibriumError::PressureBracketFailure {
            lower_pressure_pascal: lower.pressure_pascal,
            lower_residual_pascal: lower.residual_pascal,
            upper_pressure_pascal: upper.pressure_pascal,
            upper_residual_pascal: upper.residual_pascal,
        });
    }

    let mut best = if lower.residual_pascal.abs() <= upper.residual_pascal.abs() {
        lower.clone()
    } else {
        upper.clone()
    };

    for iteration in 1..=MAX_PRESSURE_ITERATIONS {
        let middle_pressure = (lower.pressure_pascal + upper.pressure_pascal) * 0.5;
        let middle = evaluate_pressure(
            registry,
            problem,
            middle_pressure,
            problem.gas_volume_cubic_meter,
        )?;
        if middle.residual_pascal.abs() < best.residual_pascal.abs() {
            best = middle.clone();
        }
        if pressure_converged(middle.pressure_pascal, middle.residual_pascal) {
            return middle.into_result(iteration);
        }

        if middle.residual_pascal.signum() == lower.residual_pascal.signum() {
            lower = middle;
        } else {
            upper = middle;
        }
    }

    Err(ClosedGasEquilibriumError::NonConvergence {
        iterations: MAX_PRESSURE_ITERATIONS,
        pressure_residual_pascal: best.residual_pascal,
    })
}

#[derive(Clone)]
struct PressureEvaluation {
    equilibrium: EquilibriumResult,
    pressure_pascal: f64,
    gas_amount_mol: f64,
    residual_pascal: f64,
}

impl PressureEvaluation {
    fn into_result(
        self,
        iterations: usize,
    ) -> Result<ClosedGasEquilibriumResult, ClosedGasEquilibriumError> {
        if self.gas_amount_mol <= MIN_STABLE_GAS_AMOUNT_MOL {
            return Err(ClosedGasEquilibriumError::NoStableGasPhase {
                gas_amount_mol: self.gas_amount_mol,
            });
        }
        Ok(ClosedGasEquilibriumResult {
            equilibrium: self.equilibrium,
            pressure_pascal: self.pressure_pascal,
            gas_amount_mol: self.gas_amount_mol,
            pressure_residual_pascal: self.residual_pascal,
            iterations,
        })
    }
}

fn evaluate_pressure(
    registry: &SpeciesRegistry,
    problem: &ClosedGasEquilibriumProblem,
    pressure_pascal: f64,
    gas_volume_cubic_meter: f64,
) -> Result<PressureEvaluation, ClosedGasEquilibriumError> {
    let fixed_pressure_problem = EquilibriumProblem {
        temperature_kelvin: problem.temperature_kelvin,
        pressure_pascal,
        initial_species_amounts_mol: problem.initial_species_amounts_mol.clone(),
        candidate_species: problem.candidate_species.clone(),
    };
    let equilibrium = solve_equilibrium(registry, &fixed_pressure_problem)
        .map_err(ClosedGasEquilibriumError::Equilibrium)?;
    let gas_amount_mol = gas_amount_mol(registry, &equilibrium)?;
    let ideal_pressure_pascal =
        gas_amount_mol * GAS_CONSTANT_JOULE_PER_MOL_KELVIN * problem.temperature_kelvin
            / gas_volume_cubic_meter;
    Ok(PressureEvaluation {
        equilibrium,
        pressure_pascal,
        gas_amount_mol,
        residual_pascal: pressure_pascal - ideal_pressure_pascal,
    })
}

fn has_gas_candidate(registry: &SpeciesRegistry, candidate_species: &[SpeciesId]) -> bool {
    candidate_species.iter().any(|species_id| {
        registry
            .species(*species_id)
            .map(|species| species.phase == PhaseKind::Gas)
            .unwrap_or(false)
    })
}

fn gas_amount_mol(
    registry: &SpeciesRegistry,
    result: &EquilibriumResult,
) -> Result<f64, ClosedGasEquilibriumError> {
    let mut amount_mol = 0.0;
    for amount in &result.species_amounts_mol {
        let species =
            registry
                .species(amount.species_id)
                .ok_or(ClosedGasEquilibriumError::Equilibrium(
                    EquilibriumError::MissingSpeciesData(amount.species_id),
                ))?;
        if species.phase == PhaseKind::Gas {
            amount_mol += amount.amount_mol.max(0.0);
        }
    }
    Ok(amount_mol)
}

fn pressure_converged(pressure_pascal: f64, residual_pascal: f64) -> bool {
    residual_pascal.abs()
        <= PRESSURE_TOLERANCE_PASCAL.max(pressure_pascal.abs() * RELATIVE_PRESSURE_TOLERANCE)
}
