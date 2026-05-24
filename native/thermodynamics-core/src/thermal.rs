use crate::chemistry::{SpeciesAmount, SpeciesId};
use crate::registry::SpeciesRegistry;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MixtureThermalState {
    pub temperature_kelvin: f64,
    pub enthalpy_joule: f64,
    pub heat_capacity_joule_per_kelvin: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ThermalError {
    InvalidTemperatureKelvin(f64),
    NegativeAmount {
        species_id: SpeciesId,
        amount_mol: f64,
    },
    MissingSpeciesData(SpeciesId),
    UnsupportedTemperatureRange {
        species_id: SpeciesId,
        temperature_kelvin: f64,
        valid_min_temperature_kelvin: f64,
        valid_max_temperature_kelvin: f64,
    },
    NonPositiveHeatCapacity {
        heat_capacity_joule_per_kelvin: f64,
    },
    TargetEnthalpyOutsideTemperatureRange {
        target_enthalpy_joule: f64,
        min_enthalpy_joule: f64,
        max_enthalpy_joule: f64,
    },
}

pub fn thermal_state_for_composition(
    registry: &SpeciesRegistry,
    amounts: &[SpeciesAmount],
    temperature_kelvin: f64,
) -> Result<MixtureThermalState, ThermalError> {
    Ok(MixtureThermalState {
        temperature_kelvin,
        enthalpy_joule: mixture_enthalpy_joule(registry, amounts, temperature_kelvin)?,
        heat_capacity_joule_per_kelvin: mixture_heat_capacity_joule_per_kelvin(registry, amounts)?,
    })
}

pub fn mixture_enthalpy_joule(
    registry: &SpeciesRegistry,
    amounts: &[SpeciesAmount],
    temperature_kelvin: f64,
) -> Result<f64, ThermalError> {
    if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
        return Err(ThermalError::InvalidTemperatureKelvin(temperature_kelvin));
    }

    let mut enthalpy_joule = 0.0;
    for amount in amounts {
        if !amount.amount_mol.is_finite() || amount.amount_mol < 0.0 {
            return Err(ThermalError::NegativeAmount {
                species_id: amount.species_id,
                amount_mol: amount.amount_mol,
            });
        }
        let species = registry
            .species(amount.species_id)
            .ok_or(ThermalError::MissingSpeciesData(amount.species_id))?;
        let valid_range = species.thermo.valid_temperature_range;
        if temperature_kelvin < valid_range.min_kelvin
            || temperature_kelvin > valid_range.max_kelvin
        {
            return Err(ThermalError::UnsupportedTemperatureRange {
                species_id: amount.species_id,
                temperature_kelvin,
                valid_min_temperature_kelvin: valid_range.min_kelvin,
                valid_max_temperature_kelvin: valid_range.max_kelvin,
            });
        }

        let standard_enthalpy = species.thermo.standard_enthalpy_of_formation;
        let heat_capacity = species.thermo.constant_pressure_heat_capacity;
        enthalpy_joule += amount.amount_mol
            * (standard_enthalpy.value_joule_per_mol
                + heat_capacity.value_joule_per_mol_kelvin
                    * (temperature_kelvin - standard_enthalpy.reference_temperature_kelvin));
    }

    Ok(enthalpy_joule)
}

pub fn mixture_heat_capacity_joule_per_kelvin(
    registry: &SpeciesRegistry,
    amounts: &[SpeciesAmount],
) -> Result<f64, ThermalError> {
    let mut heat_capacity_joule_per_kelvin = 0.0;
    for amount in amounts {
        if !amount.amount_mol.is_finite() || amount.amount_mol < 0.0 {
            return Err(ThermalError::NegativeAmount {
                species_id: amount.species_id,
                amount_mol: amount.amount_mol,
            });
        }
        let species = registry
            .species(amount.species_id)
            .ok_or(ThermalError::MissingSpeciesData(amount.species_id))?;
        heat_capacity_joule_per_kelvin += amount.amount_mol
            * species
                .thermo
                .constant_pressure_heat_capacity
                .value_joule_per_mol_kelvin;
    }

    if !heat_capacity_joule_per_kelvin.is_finite() || heat_capacity_joule_per_kelvin <= 0.0 {
        return Err(ThermalError::NonPositiveHeatCapacity {
            heat_capacity_joule_per_kelvin,
        });
    }

    Ok(heat_capacity_joule_per_kelvin)
}

pub fn solve_temperature_for_enthalpy(
    registry: &SpeciesRegistry,
    amounts: &[SpeciesAmount],
    target_enthalpy_joule: f64,
    min_temperature_kelvin: f64,
    max_temperature_kelvin: f64,
) -> Result<MixtureThermalState, ThermalError> {
    if !target_enthalpy_joule.is_finite() {
        return Err(ThermalError::TargetEnthalpyOutsideTemperatureRange {
            target_enthalpy_joule,
            min_enthalpy_joule: f64::NAN,
            max_enthalpy_joule: f64::NAN,
        });
    }
    if !min_temperature_kelvin.is_finite()
        || min_temperature_kelvin <= 0.0
        || !max_temperature_kelvin.is_finite()
        || max_temperature_kelvin < min_temperature_kelvin
    {
        return Err(ThermalError::InvalidTemperatureKelvin(
            min_temperature_kelvin.min(max_temperature_kelvin),
        ));
    }

    let min_enthalpy_joule = mixture_enthalpy_joule(registry, amounts, min_temperature_kelvin)?;
    let max_enthalpy_joule = mixture_enthalpy_joule(registry, amounts, max_temperature_kelvin)?;
    if target_enthalpy_joule < min_enthalpy_joule || target_enthalpy_joule > max_enthalpy_joule {
        return Err(ThermalError::TargetEnthalpyOutsideTemperatureRange {
            target_enthalpy_joule,
            min_enthalpy_joule,
            max_enthalpy_joule,
        });
    }

    let heat_capacity_joule_per_kelvin = mixture_heat_capacity_joule_per_kelvin(registry, amounts)?;
    let temperature_kelvin = min_temperature_kelvin
        + (target_enthalpy_joule - min_enthalpy_joule) / heat_capacity_joule_per_kelvin;

    Ok(MixtureThermalState {
        temperature_kelvin,
        enthalpy_joule: target_enthalpy_joule,
        heat_capacity_joule_per_kelvin,
    })
}
