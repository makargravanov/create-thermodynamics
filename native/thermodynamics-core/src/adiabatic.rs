use crate::chemistry::{SpeciesAmount, SpeciesId};
use crate::equilibrium::{
    solve_equilibrium, EquilibriumError, EquilibriumProblem, EquilibriumResult,
};
use crate::registry::SpeciesRegistry;
use crate::thermal::{mixture_enthalpy_joule, MixtureThermalState, ThermalError};

const ENTHALPY_TOLERANCE_JOULE: f64 = 1.0e-6;
const TEMPERATURE_TOLERANCE_KELVIN: f64 = 1.0e-9;
const MAX_ITERATIONS: usize = 128;

#[derive(Debug, Clone)]
pub struct AdiabaticEquilibriumProblem {
    pub initial_temperature_kelvin: f64,
    pub pressure_pascal: f64,
    pub initial_species_amounts_mol: Vec<SpeciesAmount>,
    pub candidate_species: Vec<SpeciesId>,
    pub min_temperature_kelvin: f64,
    pub max_temperature_kelvin: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AdiabaticEquilibriumResult {
    pub equilibrium: EquilibriumResult,
    pub thermal_state: MixtureThermalState,
    pub initial_enthalpy_joule: f64,
    pub enthalpy_residual_joule: f64,
    pub iterations: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum AdiabaticEquilibriumError {
    InvalidTemperatureBounds {
        min_temperature_kelvin: f64,
        max_temperature_kelvin: f64,
    },
    Equilibrium(EquilibriumError),
    Thermal(ThermalError),
    TargetEnthalpyNotBracketed {
        target_enthalpy_joule: f64,
        min_temperature_enthalpy_joule: f64,
        max_temperature_enthalpy_joule: f64,
    },
    NonConvergence {
        iterations: usize,
        enthalpy_residual_joule: f64,
    },
}

pub fn solve_adiabatic_equilibrium(
    registry: &SpeciesRegistry,
    problem: &AdiabaticEquilibriumProblem,
) -> Result<AdiabaticEquilibriumResult, AdiabaticEquilibriumError> {
    if !problem.min_temperature_kelvin.is_finite()
        || !problem.max_temperature_kelvin.is_finite()
        || problem.min_temperature_kelvin <= 0.0
        || problem.max_temperature_kelvin < problem.min_temperature_kelvin
    {
        return Err(AdiabaticEquilibriumError::InvalidTemperatureBounds {
            min_temperature_kelvin: problem.min_temperature_kelvin,
            max_temperature_kelvin: problem.max_temperature_kelvin,
        });
    }

    let initial_enthalpy_joule = mixture_enthalpy_joule(
        registry,
        &problem.initial_species_amounts_mol,
        problem.initial_temperature_kelvin,
    )
    .map_err(AdiabaticEquilibriumError::Thermal)?;

    let mut lower_temperature = problem.min_temperature_kelvin;
    let mut upper_temperature = problem.max_temperature_kelvin;
    let lower_state = equilibrium_enthalpy_at_temperature(registry, problem, lower_temperature)?;
    let upper_state = equilibrium_enthalpy_at_temperature(registry, problem, upper_temperature)?;
    let mut lower_residual = lower_state.thermal_state.enthalpy_joule - initial_enthalpy_joule;
    let upper_residual = upper_state.thermal_state.enthalpy_joule - initial_enthalpy_joule;

    if lower_residual.abs() <= ENTHALPY_TOLERANCE_JOULE {
        return Ok(finish_result(
            lower_state,
            initial_enthalpy_joule,
            lower_residual,
            0,
        ));
    }
    if upper_residual.abs() <= ENTHALPY_TOLERANCE_JOULE {
        return Ok(finish_result(
            upper_state,
            initial_enthalpy_joule,
            upper_residual,
            0,
        ));
    }

    if lower_residual.signum() == upper_residual.signum() {
        return Err(AdiabaticEquilibriumError::TargetEnthalpyNotBracketed {
            target_enthalpy_joule: initial_enthalpy_joule,
            min_temperature_enthalpy_joule: lower_state.thermal_state.enthalpy_joule,
            max_temperature_enthalpy_joule: upper_state.thermal_state.enthalpy_joule,
        });
    }

    let mut best_residual = lower_residual;
    for iteration in 1..=MAX_ITERATIONS {
        let midpoint_temperature = (lower_temperature + upper_temperature) * 0.5;
        let midpoint_state =
            equilibrium_enthalpy_at_temperature(registry, problem, midpoint_temperature)?;
        let midpoint_residual =
            midpoint_state.thermal_state.enthalpy_joule - initial_enthalpy_joule;

        if midpoint_residual.abs() < best_residual.abs() {
            best_residual = midpoint_residual;
        }

        if midpoint_residual.abs() <= ENTHALPY_TOLERANCE_JOULE
            || (upper_temperature - lower_temperature).abs() <= TEMPERATURE_TOLERANCE_KELVIN
        {
            return Ok(finish_result(
                midpoint_state,
                initial_enthalpy_joule,
                midpoint_residual,
                iteration,
            ));
        }

        if midpoint_residual.signum() == lower_residual.signum() {
            lower_temperature = midpoint_temperature;
            lower_residual = midpoint_residual;
        } else {
            upper_temperature = midpoint_temperature;
        }
    }

    Err(AdiabaticEquilibriumError::NonConvergence {
        iterations: MAX_ITERATIONS,
        enthalpy_residual_joule: best_residual,
    })
}

#[derive(Debug, Clone)]
struct EquilibriumThermalState {
    equilibrium: EquilibriumResult,
    thermal_state: MixtureThermalState,
}

fn equilibrium_enthalpy_at_temperature(
    registry: &SpeciesRegistry,
    problem: &AdiabaticEquilibriumProblem,
    temperature_kelvin: f64,
) -> Result<EquilibriumThermalState, AdiabaticEquilibriumError> {
    let equilibrium_problem = EquilibriumProblem {
        temperature_kelvin,
        pressure_pascal: problem.pressure_pascal,
        initial_species_amounts_mol: problem.initial_species_amounts_mol.clone(),
        candidate_species: problem.candidate_species.clone(),
    };
    let equilibrium = solve_equilibrium(registry, &equilibrium_problem)
        .map_err(AdiabaticEquilibriumError::Equilibrium)?;
    let enthalpy_joule = mixture_enthalpy_joule(
        registry,
        &equilibrium.species_amounts_mol,
        temperature_kelvin,
    )
    .map_err(AdiabaticEquilibriumError::Thermal)?;
    let heat_capacity_joule_per_kelvin = crate::thermal::mixture_heat_capacity_joule_per_kelvin(
        registry,
        &equilibrium.species_amounts_mol,
    )
    .map_err(AdiabaticEquilibriumError::Thermal)?;

    Ok(EquilibriumThermalState {
        equilibrium,
        thermal_state: MixtureThermalState {
            temperature_kelvin,
            enthalpy_joule,
            heat_capacity_joule_per_kelvin,
        },
    })
}

fn finish_result(
    state: EquilibriumThermalState,
    initial_enthalpy_joule: f64,
    enthalpy_residual_joule: f64,
    iterations: usize,
) -> AdiabaticEquilibriumResult {
    AdiabaticEquilibriumResult {
        equilibrium: state.equilibrium,
        thermal_state: state.thermal_state,
        initial_enthalpy_joule,
        enthalpy_residual_joule,
        iterations,
    }
}
