use std::collections::{BTreeMap, BTreeSet};

use crate::chemistry::error::{ChemistryError, ChemistryResult};

use super::types::*;
use super::validation::{validate_finite, validate_fraction, validate_positive_finite};

pub(super) fn generated_system_for_composition(
    composition: &MetallurgicalComposition,
    element_data: &[MetallurgicalElementData],
    pair_interactions: &[MetallurgicalPairInteractionData],
    compound_phases: &[MetallurgicalCompoundPhaseData],
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
    let pair_summary = pair_summary(&components, pair_interactions);
    let element_intermetallic_tendency =
        weighted_average(&components, |data| data.intermetallic_forming_tendency)?;
    let intermetallic_tendency =
        element_intermetallic_tendency.max(pair_summary.intermetallic_tendency);
    let element_phase_separation_tendency =
        weighted_average(&components, |data| data.phase_separation_tendency)?;
    let phase_separation_tendency =
        element_phase_separation_tendency.max(pair_summary.phase_separation_tendency);

    let solidus = (average_melting_point - 0.30 * melting_spread - radius_mismatch * 600.0)
        .min(
            pair_summary
                .eutectic_temperature_kelvin
                .unwrap_or(f64::INFINITY),
        )
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
        &pair_summary,
    )?);
    let mut added_specific_compound = false;
    for compound in compound_phases_for_composition(compound_phases, composition) {
        system = system.phase_model(compound_phase_model(
            &system_id,
            composition,
            &components,
            compound,
        )?);
        added_specific_compound = true;
    }
    if !added_specific_compound
        && (intermetallic_tendency >= 0.35 || carbide_phase_required(&components))
    {
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

fn compound_phase_model(
    system_id: &str,
    composition: &MetallurgicalComposition,
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
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
    );
    if let Some(kinetic_model) = &compound.kinetic_model {
        model = model.kinetic_model(kinetic_model.clone());
    }
    for (component, _, _) in components {
        let target = compound.components.get(component).copied().unwrap_or(0.0);
        let min = (target - compound.composition_tolerance_fraction).clamp(0.0, 1.0);
        let max = (target + compound.composition_tolerance_fraction).clamp(0.0, 1.0);
        model = model.limit(ComponentLimit::new(component.clone(), min, max.max(min)));
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
    composition_energy_for_map(
        &composition.components,
        reference_gibbs_j_per_mol,
        penalty_j_per_mol,
    )
}

fn composition_energy_for_map(
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

#[derive(Debug, Clone)]
struct PairInteractionSummary {
    phase_separation_tendency: f64,
    intermetallic_tendency: f64,
    resistivity_penalty: f64,
    ductility_penalty: f64,
    strengthening_mpa: f64,
    eutectic_temperature_kelvin: Option<f64>,
    solid_solution_limits: BTreeMap<MetallurgicalComponentId, f64>,
}

impl PairInteractionSummary {
    fn solid_solution_limit_for(&self, component: &MetallurgicalComponentId) -> Option<f64> {
        self.solid_solution_limits.get(component).copied()
    }
}

fn pair_summary(
    components: &[(MetallurgicalComponentId, f64, &MetallurgicalElementData)],
    pair_interactions: &[MetallurgicalPairInteractionData],
) -> PairInteractionSummary {
    let mut summary = PairInteractionSummary {
        phase_separation_tendency: 0.0,
        intermetallic_tendency: 0.0,
        resistivity_penalty: 0.0,
        ductility_penalty: 0.0,
        strengthening_mpa: 0.0,
        eutectic_temperature_kelvin: None,
        solid_solution_limits: BTreeMap::new(),
    };
    for (left_index, (left, left_fraction, _)) in components.iter().enumerate() {
        for (right, right_fraction, _) in components.iter().skip(left_index + 1) {
            let pair_weight = 4.0 * left_fraction * right_fraction;
            let interaction = pair_interactions
                .iter()
                .find(|interaction| interaction.contains_pair(left, right));
            let solid_miscibility = interaction
                .map(|interaction| interaction.solid_miscibility)
                .unwrap_or(SolidMiscibility::Limited);
            let separation = solid_miscibility_separation_score(solid_miscibility);
            summary.phase_separation_tendency = summary
                .phase_separation_tendency
                .max(separation * pair_weight);
            summary.intermetallic_tendency = summary.intermetallic_tendency.max(
                interaction
                    .map(|interaction| interaction.interaction_strength_j_per_mol.abs() / 30_000.0)
                    .unwrap_or(0.0)
                    .clamp(0.0, 1.0)
                    * pair_weight,
            );
            summary.resistivity_penalty += interaction
                .map(|interaction| interaction.resistivity_penalty_per_fraction)
                .unwrap_or(0.04)
                * pair_weight;
            summary.ductility_penalty += interaction
                .map(|interaction| interaction.ductility_penalty_per_fraction)
                .unwrap_or(0.08)
                * pair_weight;
            summary.strengthening_mpa += interaction
                .map(|interaction| interaction.strengthening_mpa_per_fraction)
                .unwrap_or(120.0)
                * pair_weight;
            if let Some(interaction) = interaction {
                if let Some(eutectic_temperature) = interaction.eutectic_temperature_kelvin {
                    summary.eutectic_temperature_kelvin = Some(
                        summary
                            .eutectic_temperature_kelvin
                            .map(|current| current.min(eutectic_temperature))
                            .unwrap_or(eutectic_temperature),
                    );
                }
            }
            let limit = solid_miscibility_limit(solid_miscibility);
            merge_solution_limit(&mut summary.solid_solution_limits, left.clone(), limit);
            merge_solution_limit(&mut summary.solid_solution_limits, right.clone(), limit);
        }
    }
    summary.ductility_penalty = summary.ductility_penalty.clamp(0.0, 0.85);
    summary
}

fn solid_miscibility_separation_score(miscibility: SolidMiscibility) -> f64 {
    match miscibility {
        SolidMiscibility::Complete => 0.0,
        SolidMiscibility::High => 0.10,
        SolidMiscibility::Limited => 0.35,
        SolidMiscibility::VeryLimited => 0.70,
        SolidMiscibility::Immiscible => 1.0,
    }
}

fn solid_miscibility_limit(miscibility: SolidMiscibility) -> f64 {
    match miscibility {
        SolidMiscibility::Complete => 1.0,
        SolidMiscibility::High => 0.70,
        SolidMiscibility::Limited => 0.35,
        SolidMiscibility::VeryLimited => 0.12,
        SolidMiscibility::Immiscible => 0.04,
    }
}

fn merge_solution_limit(
    limits: &mut BTreeMap<MetallurgicalComponentId, f64>,
    component: MetallurgicalComponentId,
    limit: f64,
) {
    limits
        .entry(component)
        .and_modify(|current| *current = current.min(limit))
        .or_insert(limit);
}

fn compound_phases_for_composition<'a>(
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
                && compound_distance(phase, composition) <= phase.composition_tolerance_fraction
        })
        .collect()
}

fn compound_distance(
    phase: &MetallurgicalCompoundPhaseData,
    composition: &MetallurgicalComposition,
) -> f64 {
    phase
        .components
        .iter()
        .map(|(component, target)| (composition.fraction_of(component) - target).abs())
        .sum::<f64>()
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

pub(crate) fn validate_pair_interactions(
    pair_interactions: Vec<MetallurgicalPairInteractionData>,
) -> ChemistryResult<Vec<MetallurgicalPairInteractionData>> {
    let mut seen = BTreeSet::new();
    for interaction in &pair_interactions {
        interaction.validate()?;
        let key = ordered_pair_key(&interaction.first, &interaction.second);
        if !seen.insert(key.clone()) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "duplicate metallurgical pair interaction '{}:{}'",
                key.0.as_str(),
                key.1.as_str()
            )));
        }
    }
    Ok(pair_interactions)
}

pub(crate) fn validate_compound_phases(
    compound_phases: Vec<MetallurgicalCompoundPhaseData>,
) -> ChemistryResult<Vec<MetallurgicalCompoundPhaseData>> {
    let mut seen = BTreeSet::new();
    for phase in &compound_phases {
        phase.validate()?;
        if !seen.insert(phase.id.clone()) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "duplicate metallurgical compound phase '{}'",
                phase.id
            )));
        }
    }
    Ok(compound_phases)
}

fn ordered_pair_key(
    left: &MetallurgicalComponentId,
    right: &MetallurgicalComponentId,
) -> (MetallurgicalComponentId, MetallurgicalComponentId) {
    if left <= right {
        (left.clone(), right.clone())
    } else {
        (right.clone(), left.clone())
    }
}
