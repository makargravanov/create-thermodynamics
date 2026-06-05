use std::collections::BTreeMap;

use super::GeneratedComponent;
use crate::chemistry::metallurgy::{
    MetallurgicalComponentId, MetallurgicalPairInteractionData, SolidMiscibility,
};

#[derive(Debug, Clone)]
pub(super) struct PairInteractionSummary {
    pub phase_separation_tendency: f64,
    pub intermetallic_tendency: f64,
    pub resistivity_penalty: f64,
    pub ductility_penalty: f64,
    pub strengthening_mpa: f64,
    pub eutectic_temperature_kelvin: Option<f64>,
    pub eutectic_pair: Option<EutecticPairSummary>,
    pub solid_solution_limits: BTreeMap<MetallurgicalComponentId, f64>,
    pub used_pair_interactions: Vec<String>,
    pub missing_pair_interactions: Vec<String>,
}

#[derive(Debug, Clone)]
pub(super) struct EutecticPairSummary {
    pub first: MetallurgicalComponentId,
    pub second: MetallurgicalComponentId,
    pub second_fraction: f64,
    pub temperature_kelvin: f64,
}

impl PairInteractionSummary {
    pub(super) fn solid_solution_limit_for(
        &self,
        component: &MetallurgicalComponentId,
    ) -> Option<f64> {
        self.solid_solution_limits.get(component).copied()
    }
}

pub(super) fn pair_summary(
    components: &[GeneratedComponent<'_>],
    pair_interactions: &[MetallurgicalPairInteractionData],
) -> PairInteractionSummary {
    let mut summary = PairInteractionSummary {
        phase_separation_tendency: 0.0,
        intermetallic_tendency: 0.0,
        resistivity_penalty: 0.0,
        ductility_penalty: 0.0,
        strengthening_mpa: 0.0,
        eutectic_temperature_kelvin: None,
        eutectic_pair: None,
        solid_solution_limits: BTreeMap::new(),
        used_pair_interactions: Vec::new(),
        missing_pair_interactions: Vec::new(),
    };
    for (left_index, (left, left_fraction, _)) in components.iter().enumerate() {
        for (right, right_fraction, _) in components.iter().skip(left_index + 1) {
            let pair_weight = 4.0 * left_fraction * right_fraction;
            let interaction = pair_interactions
                .iter()
                .find(|interaction| interaction.contains_pair(left, right));
            let pair_id = format!("{}:{}", left.as_str(), right.as_str());
            if interaction.is_some() {
                summary.used_pair_interactions.push(pair_id);
            } else {
                summary.missing_pair_interactions.push(pair_id);
            }
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
                    let second_fraction = interaction
                        .eutectic_second_fraction
                        .expect("validated pair interaction has eutectic composition");
                    if summary
                        .eutectic_pair
                        .as_ref()
                        .is_none_or(|current| eutectic_temperature < current.temperature_kelvin)
                    {
                        summary.eutectic_pair = Some(EutecticPairSummary {
                            first: interaction.first.clone(),
                            second: interaction.second.clone(),
                            second_fraction,
                            temperature_kelvin: eutectic_temperature,
                        });
                    }
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
