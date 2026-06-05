use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::metallurgy::{
    ComponentLimit, CompositionEnergyTerm, MetallurgicalComponentId, MetallurgicalComposition,
    MetallurgicalCompoundPhaseData, MetallurgicalPhaseKind, MetallurgicalPhaseModel,
    MetallurgicalPhasePropertyModel, PhaseFractionHint, PhaseFreeEnergyModel,
};

use super::super::constants::TRACE_COMPONENT_FRACTION;
use super::pair::{EutecticPairSummary, PairInteractionSummary};
use super::properties::{composition_energy, composition_energy_for_map, weighted_properties};
use super::GeneratedComponent;

pub(super) fn solid_solution_phase(
    system_id: &str,
    composition: &MetallurgicalComposition,
    components: &[GeneratedComponent<'_>],
    matrix: MetallurgicalComponentId,
    solidus: f64,
    liquidus: f64,
    radius_mismatch: f64,
    pair_summary: &PairInteractionSummary,
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
    property_model.yield_strength_mpa += pair_summary.strengthening_mpa;
    property_model.hardness_hv += pair_summary.strengthening_mpa * 0.28;
    property_model.ductility_fraction = (property_model.ductility_fraction
        * (1.0 - pair_summary.ductility_penalty))
        .clamp(0.01, 0.85);
    property_model.electrical_resistivity_micro_ohm_meter += pair_summary.resistivity_penalty;

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
        let max_fraction = if component == &matrix {
            1.0
        } else {
            pair_summary
                .solid_solution_limit_for(component)
                .unwrap_or(1.0)
        };
        model = model.limit(ComponentLimit::new(component.clone(), 0.0, max_fraction));
    }
    Ok(model)
}

pub(super) fn compound_phase_model(
    system_id: &str,
    composition: &MetallurgicalComposition,
    components: &[GeneratedComponent<'_>],
    compound: &MetallurgicalCompoundPhaseData,
) -> ChemistryResult<MetallurgicalPhaseModel> {
    let mut model = MetallurgicalPhaseModel::new(
        format!("{system_id}/{}", compound.id.replace(':', "_")),
        compound.kind,
        compound.property_model.clone(),
    )
    .free_energy_model(
        composition_energy_for_map(
            &compound.components,
            compound.formation_energy_j_per_mol,
            22_000.0,
        )
        .temperature_window(
            compound.low_temperature_kelvin,
            compound.high_temperature_kelvin,
        ),
    )
    .fraction_hint(PhaseFractionHint::new(
        stoichiometric_phase_fraction(&compound.components, composition).clamp(0.0, 1.0),
        5.0,
        format!(
            "stoichiometric lever rule from available components for compound '{}'",
            compound.id
        ),
    ));
    if let Some(kinetic_model) = &compound.kinetic_model {
        model = model.kinetic_model(kinetic_model.clone());
    }
    for (component, _, _) in components {
        model = model.limit(ComponentLimit::new(component.clone(), 0.0, 1.0));
    }
    for component in compound.components.keys() {
        if !composition.components.contains_key(component) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "generated compound phase '{}' references absent component '{}'",
                compound.id,
                component.as_str()
            )));
        }
    }
    Ok(model)
}

pub(super) fn intermetallic_phase(
    system_id: &str,
    composition: &MetallurgicalComposition,
    components: &[GeneratedComponent<'_>],
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

pub(super) fn component_rich_phase(
    system_id: &str,
    components: &[GeneratedComponent<'_>],
    matrix: &MetallurgicalComponentId,
    tendency: f64,
    solidus: f64,
    eutectic_pair: Option<&EutecticPairSummary>,
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
    )
    .fraction_hint(component_rich_fraction_hint(
        &richest_secondary.0,
        richest_secondary.1,
        eutectic_pair,
    ));
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

pub(super) fn liquid_properties(
    components: &[GeneratedComponent<'_>],
) -> ChemistryResult<MetallurgicalPhasePropertyModel> {
    let mut properties = weighted_properties(components)?;
    properties.hardness_hv = (properties.hardness_hv * 0.35).max(20.0);
    properties.yield_strength_mpa = (properties.yield_strength_mpa * 0.25).max(30.0);
    properties.ductility_fraction = (properties.ductility_fraction + 0.30).clamp(0.25, 0.85);
    properties.electrical_resistivity_micro_ohm_meter *= 1.25;
    properties.thermal_conductivity_w_per_meter_kelvin *= 0.70;
    Ok(properties)
}

pub(super) fn compound_phases_for_composition<'a>(
    compound_phases: &'a [MetallurgicalCompoundPhaseData],
    composition: &MetallurgicalComposition,
) -> Vec<&'a MetallurgicalCompoundPhaseData> {
    compound_phases
        .iter()
        .filter(|phase| {
            phase
                .components
                .keys()
                .all(|component| composition.components.contains_key(component))
                && stoichiometric_phase_fraction(&phase.components, composition)
                    > TRACE_COMPONENT_FRACTION
        })
        .collect()
}

pub(super) fn considered_compound_phase_ids(
    compound_phases: &[MetallurgicalCompoundPhaseData],
    composition: &MetallurgicalComposition,
) -> Vec<String> {
    compound_phases
        .iter()
        .filter(|phase| {
            phase
                .components
                .keys()
                .all(|component| composition.components.contains_key(component))
        })
        .map(|phase| phase.id.clone())
        .collect()
}

fn stoichiometric_phase_fraction(
    phase_composition: &std::collections::BTreeMap<MetallurgicalComponentId, f64>,
    alloy_composition: &MetallurgicalComposition,
) -> f64 {
    phase_composition
        .iter()
        .filter(|(_, target)| **target > 0.0)
        .map(|(component, target)| alloy_composition.fraction_of(component) / target)
        .fold(1.0, f64::min)
        .clamp(0.0, 1.0)
}

fn component_rich_fraction_hint(
    rich_component: &MetallurgicalComponentId,
    rich_component_fraction: f64,
    eutectic_pair: Option<&EutecticPairSummary>,
) -> PhaseFractionHint {
    let rich_phase_component_fraction = eutectic_pair
        .and_then(|pair| {
            if &pair.second == rich_component {
                Some(pair.second_fraction.max(0.82))
            } else if &pair.first == rich_component {
                Some((1.0 - pair.second_fraction).max(0.82))
            } else {
                None
            }
        })
        .unwrap_or(0.82)
        .clamp(0.08, 1.0);
    let target_fraction = (rich_component_fraction / rich_phase_component_fraction).clamp(0.0, 1.0);
    PhaseFractionHint::new(
        target_fraction,
        4.0,
        format!(
            "component-rich lever rule: component '{}' fraction {:.6} divided by phase-side fraction {:.6}",
            rich_component.as_str(),
            rich_component_fraction,
            rich_phase_component_fraction
        ),
    )
}

pub(super) fn carbide_phase_required(components: &[GeneratedComponent<'_>]) -> bool {
    let has_carbon = components.iter().any(|(component, fraction, _)| {
        component.as_str() == "destroy:carbon" && *fraction > 0.005
    });
    let has_carbide_former = components
        .iter()
        .any(|(_, fraction, data)| *fraction > 0.01 && data.carbide_forming_tendency > 0.45);
    has_carbon && has_carbide_former
}
