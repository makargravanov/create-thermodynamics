use std::collections::BTreeMap;

use crate::chemistry::error::ChemistryResult;

use super::super::validation::{validate_finite, validate_fraction, validate_positive_finite};
use super::GeneratedComponent;
use crate::chemistry::metallurgy::{
    CompositionEnergyTerm, MetallurgicalComponentId, MetallurgicalComposition,
    MetallurgicalPhasePropertyModel, MetallurgicalPropertyCalibration, PhaseFreeEnergyModel,
};

pub(super) fn weighted_properties(
    components: &[GeneratedComponent<'_>],
) -> ChemistryResult<MetallurgicalPhasePropertyModel> {
    let hardness_hv = weighted_average(components, |data| data.base_property_model.hardness_hv)?;
    let yield_strength_mpa = weighted_average(components, |data| {
        data.base_property_model.yield_strength_mpa
    })?;
    let ductility_fraction = weighted_average(components, |data| {
        data.base_property_model.ductility_fraction
    })?;
    let electrical_resistivity_micro_ohm_meter = weighted_average(components, |data| {
        data.base_property_model
            .electrical_resistivity_micro_ohm_meter
    })?;
    let thermal_conductivity_w_per_meter_kelvin = weighted_average(components, |data| {
        data.base_property_model
            .thermal_conductivity_w_per_meter_kelvin
    })?;
    let corrosion_resistance_score = weighted_average(components, |data| {
        data.base_property_model.corrosion_resistance_score
    })?;
    Ok(MetallurgicalPhasePropertyModel {
        hardness_hv,
        yield_strength_mpa,
        ductility_fraction,
        electrical_resistivity_micro_ohm_meter,
        thermal_conductivity_w_per_meter_kelvin,
        corrosion_resistance_score,
    })
}

pub(super) fn composition_energy(
    composition: &MetallurgicalComposition,
    reference_gibbs_j_per_mol: f64,
    penalty_j_per_mol: f64,
) -> PhaseFreeEnergyModel {
    composition_energy_for_map(
        &composition.components,
        reference_gibbs_j_per_mol,
        penalty_j_per_mol,
    )
}

pub(super) fn composition_energy_for_map(
    components: &BTreeMap<MetallurgicalComponentId, f64>,
    reference_gibbs_j_per_mol: f64,
    penalty_j_per_mol: f64,
) -> PhaseFreeEnergyModel {
    let mut model = PhaseFreeEnergyModel::new(reference_gibbs_j_per_mol, 6.0);
    for (component, fraction) in components {
        model = model.composition_term(CompositionEnergyTerm::new(
            component.clone(),
            *fraction,
            0.45,
            penalty_j_per_mol,
        ));
    }
    model
}

pub(super) fn weighted_average(
    components: &[GeneratedComponent<'_>],
    value: impl Fn(&crate::chemistry::metallurgy::MetallurgicalElementData) -> f64,
) -> ChemistryResult<f64> {
    let mut total = 0.0;
    for (_, fraction, data) in components {
        validate_fraction(*fraction, "generated metallurgical component fraction")?;
        let value = value(data);
        validate_finite(value, "generated metallurgical property")?;
        total += fraction * value;
    }
    validate_finite(total, "generated metallurgical weighted average")?;
    Ok(total)
}

pub(super) fn weighted_property_calibration(
    components: &[GeneratedComponent<'_>],
) -> ChemistryResult<MetallurgicalPropertyCalibration> {
    let yield_strength = weighted_average(components, |data| {
        data.base_property_model.yield_strength_mpa
    })?;
    let ductility = weighted_average(components, |data| {
        data.base_property_model.ductility_fraction
    })?;
    let resistivity = weighted_average(components, |data| {
        data.base_property_model
            .electrical_resistivity_micro_ohm_meter
    })?;
    let conductivity = weighted_average(components, |data| {
        data.base_property_model
            .thermal_conductivity_w_per_meter_kelvin
    })?;
    let strengthening = weighted_average(components, |data| {
        data.solid_solution_strengthening_mpa_per_fraction
    })?;
    let precipitation_response = (80.0 + strengthening * 0.16).clamp(100.0, 420.0);
    let cold_work_response = (yield_strength * (1.5 + ductility)).clamp(80.0, 520.0);
    let hall_petch = (35.0 + yield_strength.sqrt() * 5.0).clamp(35.0, 150.0);
    let dislocation = (yield_strength.sqrt() * 1.2).clamp(7.0, 35.0);
    let transport_sensitivity = if resistivity <= 0.03 || conductivity >= 180.0 {
        1.0
    } else if resistivity <= 0.12 || conductivity >= 80.0 {
        0.45
    } else {
        0.18
    };
    let calibration = MetallurgicalPropertyCalibration::neutral()
        .strength_response(
            hall_petch,
            dislocation,
            precipitation_response,
            cold_work_response,
        )
        .hardness_per_strength((0.22 + (1.0 - ductility).clamp(0.0, 1.0) * 0.16).clamp(0.20, 0.42))
        .ductility_penalties(
            (700.0 + (1.0 - ductility).clamp(0.0, 1.0) * 900.0).clamp(500.0, 1800.0),
            (0.20 + (1.0 - ductility).clamp(0.0, 1.0) * 0.40).clamp(0.15, 0.70),
            (0.30 + (1.0 - ductility).clamp(0.0, 1.0) * 0.40).clamp(0.25, 0.85),
        )
        .transport_penalties(
            resistivity * 0.45 * transport_sensitivity,
            resistivity * 0.30 * transport_sensitivity,
            0.35 * transport_sensitivity,
            0.25 * transport_sensitivity,
        );
    calibration.validate()?;
    Ok(calibration)
}

pub(super) fn radius_mismatch(
    components: &[GeneratedComponent<'_>],
    matrix_radius_pm: f64,
) -> ChemistryResult<f64> {
    validate_positive_finite(matrix_radius_pm, "generated metallurgical matrix radius")?;
    let mismatch = components
        .iter()
        .map(|(_, fraction, data)| {
            fraction * ((data.atomic_radius_pm - matrix_radius_pm).abs() / matrix_radius_pm)
        })
        .sum::<f64>();
    validate_finite(mismatch, "generated metallurgical radius mismatch")?;
    Ok(mismatch)
}
