use crate::chemistry::error::{ChemistryError, ChemistryResult};

use super::types::*;

pub(super) fn validate_phase_model(
    system: &MetallurgicalSystem,
    phase: &MetallurgicalPhaseModel,
) -> ChemistryResult<()> {
    if phase.id.trim().is_empty() {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "metallurgical system '{}' has a phase with empty id",
            system.id
        )));
    }
    for limit in &phase.component_limits {
        if !system.components.contains(&limit.component) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "phase '{}' references component '{}' outside metallurgical system '{}'",
                phase.id,
                limit.component.as_str(),
                system.id
            )));
        }
        validate_fraction(limit.min_fraction, "phase component minimum")?;
        validate_fraction(limit.max_fraction, "phase component maximum")?;
        if limit.min_fraction > limit.max_fraction {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "phase '{}' has component limit minimum above maximum",
                phase.id
            )));
        }
    }
    validate_phase_free_energy_model(&phase.free_energy_model)?;
    if let Some(hint) = &phase.fraction_hint {
        validate_fraction(hint.target_fraction, "phase fraction hint target")?;
        validate_non_negative_finite(hint.strength, "phase fraction hint strength")?;
        if hint.reason.trim().is_empty() {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "phase '{}' has a fraction hint without a reason",
                phase.id
            )));
        }
    }
    validate_property_model(&phase.property_model)?;
    validate_kinetic_model(&phase.kinetic_model)?;
    Ok(())
}

pub(super) fn validate_phase_boundary_point(
    system: &MetallurgicalSystem,
    point: &PhaseBoundaryPoint,
) -> ChemistryResult<()> {
    if point.composition.is_empty() {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "metallurgical system '{}' has an empty phase-boundary composition",
            system.id
        )));
    }
    let mut total = 0.0;
    for (component, fraction) in &point.composition {
        if !system.components.contains(component) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "phase-boundary point references component '{}' outside metallurgical system '{}'",
                component.as_str(),
                system.id
            )));
        }
        validate_fraction(*fraction, "phase-boundary component fraction")?;
        total += fraction;
    }
    if (total - 1.0).abs() > 1.0e-6 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "phase-boundary composition for system '{}' must sum to 1.0, got {total}",
            system.id
        )));
    }
    validate_phase_boundary_temperatures(point.solidus_kelvin, point.liquidus_kelvin)
}

pub(super) fn validate_phase_boundary_temperatures(
    solidus_kelvin: f64,
    liquidus_kelvin: f64,
) -> ChemistryResult<()> {
    validate_positive_finite(solidus_kelvin, "phase-boundary solidus")?;
    validate_positive_finite(liquidus_kelvin, "phase-boundary liquidus")?;
    if solidus_kelvin > liquidus_kelvin {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "phase-boundary solidus {solidus_kelvin} K is above liquidus {liquidus_kelvin} K"
        )));
    }
    Ok(())
}

pub(super) fn validate_phase_free_energy_model(
    model: &PhaseFreeEnergyModel,
) -> ChemistryResult<()> {
    validate_finite(
        model.reference_gibbs_j_per_mol,
        "phase reference Gibbs free energy",
    )?;
    validate_finite(model.entropy_j_per_mol_kelvin, "phase entropy")?;
    validate_non_negative_finite(model.low_temperature_kelvin, "phase low temperature")?;
    validate_positive_finite(model.high_temperature_kelvin, "phase high temperature")?;
    if model.low_temperature_kelvin >= model.high_temperature_kelvin {
        return Err(ChemistryError::InvalidMixtureState(
            "phase temperature window low bound must be below high bound".to_string(),
        ));
    }
    for term in &model.composition_terms {
        validate_fraction(term.center_fraction, "phase composition center")?;
        validate_positive_finite(term.width_fraction, "phase composition width")?;
        validate_non_negative_finite(term.penalty_j_per_mol, "phase composition penalty")?;
    }
    if let Some(threshold) = model.cooling_rate_stabilization_threshold_kelvin_per_second {
        validate_non_negative_finite(threshold, "phase cooling-rate threshold")?;
        validate_non_negative_finite(
            model.cooling_rate_stabilization_j_per_mol,
            "phase cooling-rate stabilization",
        )?;
    }
    Ok(())
}

pub(super) fn validate_property_model(
    model: &MetallurgicalPhasePropertyModel,
) -> ChemistryResult<()> {
    validate_non_negative_finite(model.hardness_hv, "phase hardness")?;
    validate_non_negative_finite(model.yield_strength_mpa, "phase yield strength")?;
    validate_fraction(model.ductility_fraction, "phase ductility")?;
    validate_non_negative_finite(
        model.electrical_resistivity_micro_ohm_meter,
        "phase electrical resistivity",
    )?;
    validate_non_negative_finite(
        model.thermal_conductivity_w_per_meter_kelvin,
        "phase thermal conductivity",
    )?;
    validate_fraction(
        model.corrosion_resistance_score,
        "phase corrosion resistance",
    )?;
    Ok(())
}

pub(super) fn validate_alloy_properties(properties: &AlloyPropertySnapshot) -> ChemistryResult<()> {
    validate_non_negative_finite(properties.hardness_hv, "alloy hardness")?;
    validate_non_negative_finite(properties.yield_strength_mpa, "alloy yield strength")?;
    validate_fraction(properties.ductility_fraction, "alloy ductility")?;
    validate_non_negative_finite(
        properties.electrical_resistivity_micro_ohm_meter,
        "alloy electrical resistivity",
    )?;
    validate_non_negative_finite(
        properties.thermal_conductivity_w_per_meter_kelvin,
        "alloy thermal conductivity",
    )?;
    validate_fraction(
        properties.corrosion_resistance_score,
        "alloy corrosion resistance",
    )?;
    Ok(())
}

pub(super) fn validate_service_properties(
    properties: &AlloyServicePropertySnapshot,
) -> ChemistryResult<()> {
    validate_non_negative_finite(
        properties.fracture_toughness_mpa_sqrt_meter,
        "alloy fracture toughness",
    )?;
    validate_fraction(properties.brittleness_score, "alloy brittleness")?;
    validate_fraction(properties.wear_resistance_score, "alloy wear resistance")?;
    validate_non_negative_finite(
        properties.electrical_conductivity_percent_iacs,
        "alloy electrical conductivity",
    )?;
    validate_fraction(
        properties.high_temperature_stability_score,
        "alloy high-temperature stability",
    )?;
    validate_non_negative_finite(
        properties.softening_temperature_kelvin,
        "alloy softening temperature",
    )?;
    Ok(())
}

pub(super) fn validate_kinetic_model(model: &PhaseKineticModel) -> ChemistryResult<()> {
    validate_non_negative_finite(
        model.diffusion_prefactor_square_meters_per_second,
        "phase diffusion prefactor",
    )?;
    validate_non_negative_finite(
        model.diffusion_activation_energy_j_per_mol,
        "phase diffusion activation energy",
    )?;
    validate_non_negative_finite(
        model.transformation_rate_per_second,
        "phase transformation rate",
    )?;
    validate_non_negative_finite(
        model.grain_growth_rate_micrometers_per_second,
        "phase grain-growth rate",
    )?;
    validate_non_negative_finite(model.recovery_rate_per_second, "phase recovery rate")?;
    validate_non_negative_finite(
        model.precipitation_rate_per_second,
        "phase precipitation rate",
    )?;
    Ok(())
}

pub(super) fn validate_fraction(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be finite and within 0.0..=1.0"
        )));
    }
    Ok(())
}

pub(super) fn validate_positive_finite(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be positive and finite"
        )));
    }
    Ok(())
}

pub(super) fn validate_non_negative_finite(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be non-negative and finite"
        )));
    }
    Ok(())
}

pub(super) fn validate_finite(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be finite"
        )));
    }
    Ok(())
}
