mod pair;
mod phases;
mod properties;
mod validation;

use std::collections::BTreeMap;

use crate::chemistry::error::{ChemistryError, ChemistryResult};

use super::types::*;
use super::validation::validate_fraction;
use pair::pair_summary;
use phases::{
    carbide_phase_required, component_rich_phase, compound_phase_model,
    compound_phases_for_composition, considered_compound_phase_ids, intermetallic_phase,
    liquid_properties, solid_solution_phase,
};
use properties::{radius_mismatch, weighted_average};
pub(crate) use validation::{
    validate_compound_phases, validate_element_data, validate_pair_interactions,
};

pub(super) type GeneratedComponent<'a> =
    (MetallurgicalComponentId, f64, &'a MetallurgicalElementData);

pub(super) struct GeneratedMetallurgicalSystem {
    pub system: MetallurgicalSystem,
    pub diagnostic: GeneratedMetallurgyDiagnostic,
}

pub(super) fn generated_system_for_composition(
    composition: &MetallurgicalComposition,
    element_data: &[MetallurgicalElementData],
    pair_interactions: &[MetallurgicalPairInteractionData],
    compound_phases: &[MetallurgicalCompoundPhaseData],
) -> ChemistryResult<Option<GeneratedMetallurgicalSystem>> {
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

    let considered_compound_phases = considered_compound_phase_ids(compound_phases, composition);
    let mut selected_compound_phases = Vec::new();
    for compound in compound_phases_for_composition(compound_phases, composition) {
        system = system.phase_model(compound_phase_model(
            &system_id,
            composition,
            &components,
            compound,
        )?);
        selected_compound_phases.push(compound.id.clone());
    }
    let used_generic_intermetallic = selected_compound_phases.is_empty()
        && (intermetallic_tendency >= 0.35 || carbide_phase_required(&components));
    if used_generic_intermetallic {
        system = system.phase_model(intermetallic_phase(
            &system_id,
            composition,
            &components,
            intermetallic_tendency,
            solidus,
        )?);
    }
    let used_component_rich_phase =
        components.len() > 1 && (phase_separation_tendency >= 0.45 || radius_mismatch >= 0.14);
    if used_component_rich_phase {
        system = system.phase_model(component_rich_phase(
            &system_id,
            &components,
            &matrix.0,
            phase_separation_tendency,
            solidus,
            pair_summary.eutectic_pair.as_ref(),
        )?);
    }
    system.validate()?;

    let diagnostic = GeneratedMetallurgyDiagnostic {
        system_id: system_id.clone(),
        matrix_component: matrix.0.clone(),
        generated_components: components
            .iter()
            .map(|(component, _, _)| component.clone())
            .collect(),
        missing_element_data: Vec::new(),
        used_pair_interactions: pair_summary.used_pair_interactions,
        missing_pair_interactions: pair_summary.missing_pair_interactions,
        considered_compound_phases,
        selected_compound_phases,
        used_generic_intermetallic,
        used_component_rich_phase,
        radius_mismatch,
        phase_separation_tendency,
        intermetallic_tendency,
        eutectic_temperature_kelvin: pair_summary.eutectic_temperature_kelvin,
        solidus_kelvin: solidus,
        liquidus_kelvin: liquidus,
        reason: "generated from metallurgical element, pair-interaction, and compound-phase data"
            .to_string(),
    };
    Ok(Some(GeneratedMetallurgicalSystem { system, diagnostic }))
}

fn generated_system_id(components: &[GeneratedComponent<'_>]) -> String {
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
