use std::collections::{BTreeMap, BTreeSet};

use crate::chemistry::error::{ChemistryError, ChemistryResult};

use super::types::*;
use super::validation::{validate_finite, validate_fraction, validate_positive_finite};

pub(super) fn generated_system_for_composition(
    composition: &MetallurgicalComposition,
    element_data: &[MetallurgicalElementData],
) -> ChemistryResult<Option<MetallurgicalSystem>> {
    let data_by_component = element_data
        .iter()
        .map(|data| (data.component.clone(), data))
        .collect::<BTreeMap<_, _>>();
    let mut components = Vec::new();
    for (component, fraction) in &composition.components {
        let Some(data) = data_by_component.get(component) else {
            return Ok(None);
        };
        components.push((component.clone(), *fraction, *data));
    }
    if components.is_empty() {
        return Ok(None);
    }
    for (_, fraction, data) in &components {
        validate_fraction(*fraction, "generated metallurgical component fraction")?;
        data.validate()?;
    }

    let matrix = components
        .iter()
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .ok_or_else(|| {
            ChemistryError::InvalidMixtureState(
                "generated metallurgy requires at least one component".to_string(),
            )
        })?;
    let average_melting_point = weighted_average(&components, |data| data.melting_point_kelvin)?;
    let minimum_melting_point = components
        .iter()
        .map(|(_, _, data)| data.melting_point_kelvin)
        .fold(f64::INFINITY, f64::min);
    let maximum_melting_point = components
        .iter()
        .map(|(_, _, data)| data.melting_point_kelvin)
        .fold(0.0_f64, f64::max);
    let melting_spread = maximum_melting_point - minimum_melting_point;
    let radius_mismatch = radius_mismatch(&components, matrix.2.atomic_radius_pm)?;
    let intermetallic_tendency =
        weighted_average(&components, |data| data.intermetallic_forming_tendency)?;
    let phase_separation_tendency =
        weighted_average(&components, |data| data.phase_separation_tendency)?;

    let solidus = (average_melting_point - 0.30 * melting_spread - radius_mismatch * 600.0)
        .max(minimum_melting_point * 0.72)
        .max(250.0);
    let liquidus =
        (average_melting_point + 0.35 * melting_spread + phase_separation_tendency * 180.0)
            .max(solidus + 1.0);

    let system_id = generated_system_id(&components);
    let mut system = MetallurgicalSystem::new(
        system_id.clone(),
        components
            .iter()
            .map(|(component, _, _)| component.clone())
            .collect::<Vec<_>>(),
    )
    .phase_boundary(PhaseBoundaryPoint {
        composition: composition.components.clone(),
        solidus_kelvin: solidus,
        liquidus_kelvin: liquidus,
    });
    for (component, _, data) in &components {
        system = system.phase_boundary(PhaseBoundaryPoint::new(
            components
                .iter()
                .map(|(candidate, _, _)| {
                    (
                        candidate.clone(),
                        if candidate == component { 1.0 } else { 0.0 },
                    )
                })
                .collect::<Vec<_>>(),
            data.melting_point_kelvin,
            data.melting_point_kelvin,
        ));
    }

    system = system.phase_model(
        MetallurgicalPhaseModel::new(
            format!("{system_id}/liquid"),
            MetallurgicalPhaseKind::Liquid,
            liquid_properties(&components)?,
        )
        .free_energy_model(
            PhaseFreeEnergyModel::new(-20_000.0, 18.0)
                .temperature_window(solidus * 0.92, liquidus * 1.35),
        ),
    );
    system = system.phase_model(solid_solution_phase(
        &system_id,
        composition,
        &components,
        matrix.0.clone(),
        solidus,
        liquidus,
        radius_mismatch,
    )?);
    if intermetallic_tendency >= 0.35 || carbide_phase_required(&components) {
        system = system.phase_model(intermetallic_phase(
            &system_id,
            composition,
            &components,
            intermetallic_tendency,
            solidus,
        )?);
    }
    if components.len() > 1 && (phase_separation_tendency >= 0.45 || radius_mismatch >= 0.14) {
        system = system.phase_model(component_rich_phase(
            &system_id,
            &components,
            &matrix.0,
            phase_separation_tendency,
            solidus,
        )?);
    }
    system.validate()?;
    Ok(Some(system))
}

fn generated_system_id(
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
) -> String {
    let suffix = components
        .iter()
        .map(|(component, _, _)| {
            component
                .as_str()
                .chars()
                .map(|ch| {
                    if ch.is_ascii_alphanumeric() {
                        ch.to_ascii_lowercase()
                    } else {
                        '_'
                    }
                })
                .collect::<String>()
        })
        .collect::<Vec<_>>()
        .join("_");
    format!("metallurgy:generated/{suffix}")
}

fn solid_solution_phase(
    system_id: &str,
    composition: &MetallurgicalComposition,
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
    matrix: MetallurgicalComponentId,
    solidus: f64,
    liquidus: f64,
    radius_mismatch: f64,
) -> ChemistryResult<MetallurgicalPhaseModel> {
    let mut property_model = weighted_properties(components)?;
    let strengthening = components
        .iter()
        .filter(|(component, _, _)| component != &matrix)
        .map(|(_, fraction, data)| fraction * data.solid_solution_strengthening_mpa_per_fraction)
        .sum::<f64>();
    property_model.yield_strength_mpa += strengthening + radius_mismatch * 900.0;
    property_model.hardness_hv += (strengthening + radius_mismatch * 900.0) * 0.28;
    property_model.ductility_fraction =
        (property_model.ductility_fraction * (1.0 - radius_mismatch * 2.2)).clamp(0.02, 0.85);
    property_model.electrical_resistivity_micro_ohm_meter *= 1.0 + radius_mismatch * 4.0;

    let mut model = MetallurgicalPhaseModel::new(
        format!("{system_id}/generated_solid_solution"),
        MetallurgicalPhaseKind::SolidSolution,
        property_model,
    )
    .free_energy_model(
        composition_energy(composition, -14_000.0, 12_000.0)
            .temperature_window(0.0, liquidus.max(solidus + 1.0)),
    );
    for (component, _, _) in components {
        model = model.limit(ComponentLimit::new(component.clone(), 0.0, 1.0));
    }
    Ok(model)
}

fn intermetallic_phase(
    system_id: &str,
    composition: &MetallurgicalComposition,
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
    tendency: f64,
    solidus: f64,
) -> ChemistryResult<MetallurgicalPhaseModel> {
    let mut property_model = weighted_properties(components)?;
    property_model.hardness_hv =
        (property_model.hardness_hv + 420.0 + tendency * 320.0).clamp(250.0, 1150.0);
    property_model.yield_strength_mpa =
        (property_model.yield_strength_mpa + 650.0 + tendency * 550.0).clamp(500.0, 1900.0);
    property_model.ductility_fraction =
        (property_model.ductility_fraction * (0.20 - tendency * 0.10).max(0.04)).clamp(0.005, 0.18);
    property_model.electrical_resistivity_micro_ohm_meter *= 1.6 + tendency;
    property_model.thermal_conductivity_w_per_meter_kelvin *= 0.55;

    let mut model = MetallurgicalPhaseModel::new(
        format!("{system_id}/generated_intermetallic"),
        MetallurgicalPhaseKind::Intermetallic,
        property_model,
    )
    .free_energy_model(
        composition_energy(composition, -8_000.0 - tendency * 14_000.0, 18_000.0)
            .temperature_window(0.0, solidus * 1.08),
    );
    for (component, fraction, _) in components {
        model = model.limit(ComponentLimit::new(
            component.clone(),
            (fraction * 0.25).clamp(0.0, 1.0),
            (fraction + 0.50).clamp(0.0, 1.0),
        ));
    }
    Ok(model)
}

fn component_rich_phase(
    system_id: &str,
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
    matrix: &MetallurgicalComponentId,
    tendency: f64,
    solidus: f64,
) -> ChemistryResult<MetallurgicalPhaseModel> {
    let richest_secondary = components
        .iter()
        .filter(|(component, _, _)| component != matrix)
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .ok_or_else(|| {
            ChemistryError::InvalidMixtureState(
                "generated component-rich phase requires a secondary component".to_string(),
            )
        })?;
    let mut property_model = richest_secondary.2.base_property_model.clone();
    property_model.yield_strength_mpa =
        (property_model.yield_strength_mpa + tendency * 280.0).clamp(80.0, 1400.0);
    property_model.hardness_hv = (property_model.hardness_hv + tendency * 150.0).clamp(35.0, 900.0);
    property_model.ductility_fraction =
        (property_model.ductility_fraction * (1.0 - tendency * 0.55)).clamp(0.01, 0.80);

    let mut model = MetallurgicalPhaseModel::new(
        format!(
            "{system_id}/{}_rich",
            richest_secondary.0.as_str().to_ascii_lowercase()
        ),
        MetallurgicalPhaseKind::SolidSolution,
        property_model,
    )
    .free_energy_model(
        PhaseFreeEnergyModel::new(-6_000.0 - tendency * 9_000.0, 4.0)
            .composition_term(CompositionEnergyTerm::new(
                richest_secondary.0.clone(),
                0.82,
                0.30,
                16_000.0,
            ))
            .temperature_window(0.0, solidus * 1.05),
    );
    for (component, _, _) in components {
        let (min, max) = if component == &richest_secondary.0 {
            (0.08, 1.0)
        } else {
            (0.0, 0.92)
        };
        model = model.limit(ComponentLimit::new(component.clone(), min, max));
    }
    Ok(model)
}

fn liquid_properties(
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
) -> ChemistryResult<MetallurgicalPhasePropertyModel> {
    let mut properties = weighted_properties(components)?;
    properties.hardness_hv = (properties.hardness_hv * 0.35).max(20.0);
    properties.yield_strength_mpa = (properties.yield_strength_mpa * 0.25).max(30.0);
    properties.ductility_fraction = (properties.ductility_fraction + 0.30).clamp(0.25, 0.85);
    properties.electrical_resistivity_micro_ohm_meter *= 1.25;
    properties.thermal_conductivity_w_per_meter_kelvin *= 0.70;
    Ok(properties)
}

fn weighted_properties(
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
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

fn composition_energy(
    composition: &MetallurgicalComposition,
    reference_gibbs_j_per_mol: f64,
    penalty_j_per_mol: f64,
) -> PhaseFreeEnergyModel {
    let mut model = PhaseFreeEnergyModel::new(reference_gibbs_j_per_mol, 6.0);
    for (component, fraction) in &composition.components {
        model = model.composition_term(CompositionEnergyTerm::new(
            component.clone(),
            *fraction,
            0.45,
            penalty_j_per_mol,
        ));
    }
    model
}

fn weighted_average(
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
    value: impl Fn(&MetallurgicalElementData) -> f64,
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

fn radius_mismatch(
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
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

fn carbide_phase_required(
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
) -> bool {
    let has_carbon = components.iter().any(|(component, fraction, _)| {
        component.as_str() == "destroy:carbon" && *fraction > 0.005
    });
    let has_carbide_former = components
        .iter()
        .any(|(_, fraction, data)| *fraction > 0.01 && data.carbide_forming_tendency > 0.45);
    has_carbon && has_carbide_former
}

pub(crate) fn validate_element_data(
    element_data: Vec<MetallurgicalElementData>,
) -> ChemistryResult<Vec<MetallurgicalElementData>> {
    let mut seen = BTreeSet::new();
    for data in &element_data {
        data.validate()?;
        if !seen.insert(data.component.clone()) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "duplicate metallurgical element data '{}'",
                data.component.as_str()
            )));
        }
    }
    Ok(element_data)
}
