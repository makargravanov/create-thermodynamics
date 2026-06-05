use std::collections::{BTreeMap, BTreeSet};

use crate::chemistry::alloy::AlloyPhaseSnapshot;
use crate::chemistry::error::{ChemistryError, ChemistryResult};

use super::constants::{
    DEFAULT_GRAIN_SIZE_MICROMETERS, DEFAULT_HOMOGENIZATION_LENGTH_MICROMETERS,
    GAS_CONSTANT_J_PER_MOL_KELVIN, TRACE_COMPONENT_FRACTION,
};
use super::generation::generated_system_for_composition;
use super::types::*;
use super::validation::*;

pub fn apply_mechanical_working(
    state: &MetallurgicalState,
    process: MechanicalWorkingProcess,
) -> ChemistryResult<MetallurgicalState> {
    validate_mechanical_working_process(&process)?;
    if matches!(state.kind, MetallurgicalStateKind::Unmodeled { .. }) || state.phases.is_empty() {
        return Err(ChemistryError::InvalidMixtureState(
            "mechanical working requires a modeled metallurgical state with phases".to_string(),
        ));
    }
    if state
        .phases
        .iter()
        .any(|phase| phase.kind == MetallurgicalPhaseKind::Liquid && phase.fraction > 0.05)
    {
        return Err(ChemistryError::InvalidMixtureState(
            "mechanical working is not defined for mostly molten alloy phases".to_string(),
        ));
    }

    let mut worked = state.clone();
    let work_factor = mechanical_work_factor(process.mode);
    let hot_work_softening = hot_work_softening_fraction(&process, state.phase_boundaries);
    let effective_strain = process.true_strain * work_factor * (1.0 - hot_work_softening);
    worked.mechanical_history.accumulated_true_strain =
        (worked.mechanical_history.accumulated_true_strain + effective_strain).clamp(0.0, 12.0);
    worked.mechanical_history.recent_true_strain = effective_strain;
    worked.mechanical_history.strain_rate_per_second = process.strain_rate_per_second;
    worked.mechanical_history.deformation_temperature_kelvin = process.temperature_kelvin;
    worked.mechanical_history.elapsed_work_seconds += process.duration_seconds;

    let cold_work_increment = 1.0 - (-effective_strain).exp();
    worked.defect_state.cold_work_fraction = (worked.defect_state.cold_work_fraction
        + (1.0 - worked.defect_state.cold_work_fraction) * cold_work_increment)
        .clamp(0.0, 1.0);
    let dislocation_multiplier =
        1.0 + effective_strain * (process.strain_rate_per_second + 1.0).log10().max(0.0) * 80.0;
    worked.defect_state.dislocation_density_per_square_meter =
        (worked.defect_state.dislocation_density_per_square_meter * dislocation_multiplier)
            .clamp(1.0e8, 1.0e16);

    let grain_refinement = 1.0
        + effective_strain
            * match process.mode {
                MechanicalWorkingMode::Forging | MechanicalWorkingMode::Rolling => 1.8,
                MechanicalWorkingMode::Drawing | MechanicalWorkingMode::Extrusion => 1.3,
                MechanicalWorkingMode::Machining => 0.35,
            };
    worked.grain_structure.average_grain_size_micrometers =
        (worked.grain_structure.average_grain_size_micrometers / grain_refinement.max(1.0))
            .clamp(0.5, 5000.0);
    worked.grain_structure.distribution_width =
        (worked.grain_structure.distribution_width + effective_strain * 0.12).clamp(0.1, 1.0);

    let recrystallized = hot_work_softening * cold_work_increment;
    worked.mechanical_history.recrystallized_fraction =
        (worked.mechanical_history.recrystallized_fraction
            + (1.0 - worked.mechanical_history.recrystallized_fraction) * recrystallized)
            .clamp(0.0, 1.0);
    if recrystallized > 0.0 {
        worked.defect_state.cold_work_fraction *= 1.0 - recrystallized;
        worked.defect_state.dislocation_density_per_square_meter =
            (worked.defect_state.dislocation_density_per_square_meter * (1.0 - recrystallized)
                + 1.0e10 * recrystallized)
                .clamp(1.0e8, 1.0e16);
    }

    worked.properties = estimate_properties(
        &worked.phases,
        &worked.grain_structure,
        &worked.defect_state,
        &worked.diffusion_state,
        &worked.property_calibration,
    )?;
    worked.service_properties = estimate_service_properties(
        &worked.properties,
        &worked.phases,
        &worked.grain_structure,
        &worked.defect_state,
        &worked.diffusion_state,
        worked.phase_boundaries,
    )?;
    worked.use_profile = estimate_use_profile(&worked.properties, &worked.service_properties)?;
    worked.diagnostics.phase_reasons.push(MetallurgicalPhaseDiagnostic {
        phase_id: "mechanical_working".to_string(),
        kind: MetallurgicalPhaseKind::SolidSolution,
        selected: true,
        fraction: 0.0,
        gibbs_j_per_mol: None,
        energy_above_minimum_j_per_mol: None,
        reason: format!(
            "mechanical working {:?} applied true strain {:.6}, effective strain {:.6}, hot-work softening {:.6}",
            process.mode, process.true_strain, effective_strain, hot_work_softening
        ),
    });
    Ok(worked)
}

pub fn metallurgical_state_from_alloy_phase(
    alloy: &AlloyPhaseSnapshot,
    systems: &[MetallurgicalSystem],
    element_data: &[MetallurgicalElementData],
    pair_interactions: &[MetallurgicalPairInteractionData],
    compound_phases: &[MetallurgicalCompoundPhaseData],
    previous: Option<&MetallurgicalState>,
    delta_seconds: f64,
) -> ChemistryResult<MetallurgicalState> {
    validate_positive_finite(delta_seconds, "metallurgical tick duration")?;
    let composition = MetallurgicalComposition::from_alloy_phase(alloy)?;
    let thermal_treatment = match previous {
        Some(previous) => previous
            .thermal_treatment
            .advance(alloy.temperature_kelvin, delta_seconds)?,
        None => ThermalTreatmentState::initial(alloy.temperature_kelvin)?,
    };
    let mut selected_system = None;
    let mut selected_distance = f64::INFINITY;
    let mut considered_systems = Vec::new();
    for system in systems {
        system.validate()?;
        let missing_components = missing_system_components(system, &composition);
        let covers_composition = missing_components.is_empty();
        if covers_composition {
            let distance = system_distance_to_composition(system, &composition)?;
            considered_systems.push(MetallurgicalSystemSelectionDiagnostic {
                system_id: system.id.clone(),
                covers_composition,
                missing_components,
                composition_distance: Some(distance),
            });
            if distance < selected_distance {
                selected_distance = distance;
                selected_system = Some(system);
            }
        } else {
            considered_systems.push(MetallurgicalSystemSelectionDiagnostic {
                system_id: system.id.clone(),
                covers_composition,
                missing_components,
                composition_distance: None,
            });
        }
    }
    if let Some(system) = selected_system {
        return modeled_state(
            alloy,
            composition,
            system,
            thermal_treatment,
            previous,
            delta_seconds,
            considered_systems,
            None,
        );
    }
    if let Some(generated_system) = generated_system_for_composition(
        &composition,
        element_data,
        pair_interactions,
        compound_phases,
    )? {
        let generated_distance =
            system_distance_to_composition(&generated_system.system, &composition)?;
        considered_systems.push(MetallurgicalSystemSelectionDiagnostic {
            system_id: generated_system.system.id.clone(),
            covers_composition: true,
            missing_components: Vec::new(),
            composition_distance: Some(generated_distance),
        });
        return modeled_state(
            alloy,
            composition,
            &generated_system.system,
            thermal_treatment,
            previous,
            delta_seconds,
            considered_systems,
            Some(generated_system.diagnostic),
        );
    }
    let reason = "no registered metallurgical system covers all components";
    Ok(unmodeled_state(
        alloy,
        composition,
        thermal_treatment,
        considered_systems,
        delta_seconds,
        reason,
    ))
}

fn modeled_state(
    alloy: &AlloyPhaseSnapshot,
    composition: MetallurgicalComposition,
    system: &MetallurgicalSystem,
    thermal_treatment: ThermalTreatmentState,
    previous: Option<&MetallurgicalState>,
    delta_seconds: f64,
    considered_systems: Vec<MetallurgicalSystemSelectionDiagnostic>,
    generated_system: Option<GeneratedMetallurgyDiagnostic>,
) -> ChemistryResult<MetallurgicalState> {
    let phase_boundaries = phase_boundaries_for_composition(system, &composition)?;
    let mut phase_energies = Vec::new();
    let mut phase_diagnostics = Vec::new();
    for model in &system.phase_models {
        if !phase_limits_match(model, &composition) {
            phase_diagnostics.push(MetallurgicalPhaseDiagnostic {
                phase_id: model.id.clone(),
                kind: model.kind,
                selected: false,
                fraction: 0.0,
                gibbs_j_per_mol: None,
                energy_above_minimum_j_per_mol: None,
                reason: phase_limit_mismatch_reason(model, &composition),
            });
            continue;
        }
        let energy = phase_free_energy(
            model,
            &composition,
            alloy.temperature_kelvin,
            &thermal_treatment,
            Some(phase_boundaries),
        )?;
        if !energy.is_finite() {
            phase_diagnostics.push(MetallurgicalPhaseDiagnostic {
                phase_id: model.id.clone(),
                kind: model.kind,
                selected: false,
                fraction: 0.0,
                gibbs_j_per_mol: None,
                energy_above_minimum_j_per_mol: None,
                reason: phase_energy_exclusion_reason(
                    model,
                    alloy.temperature_kelvin,
                    phase_boundaries,
                ),
            });
        }
        phase_energies.push((model, energy));
    }
    if phase_energies.is_empty() {
        let reason = "registered system has no stable phase for this composition and temperature";
        return Ok(unmodeled_state(
            alloy,
            composition,
            thermal_treatment,
            considered_systems,
            delta_seconds,
            reason,
        ));
    }
    let minimum_finite_energy = phase_energies
        .iter()
        .filter_map(|(_, energy)| energy.is_finite().then_some(*energy))
        .fold(f64::INFINITY, f64::min);
    let phase_fractions = equilibrium_phase_fractions(&phase_energies, alloy.temperature_kelvin)?;
    let phase_fractions = apply_thermal_treatment_bias(
        phase_fractions,
        &thermal_treatment,
        &system.thermal_treatment_profile,
    )?;
    let phase_fractions = relax_phase_fractions(
        phase_fractions,
        &thermal_treatment,
        &system.thermal_treatment_profile,
        previous,
        delta_seconds,
    )?;
    for (model, fraction) in &phase_fractions {
        let energy = phase_energies
            .iter()
            .find(|(candidate, _)| candidate.id == model.id)
            .map(|(_, energy)| *energy)
            .filter(|energy| energy.is_finite());
        phase_diagnostics.push(MetallurgicalPhaseDiagnostic {
            phase_id: model.id.clone(),
            kind: model.kind,
            selected: *fraction > TRACE_COMPONENT_FRACTION,
            fraction: *fraction,
            gibbs_j_per_mol: energy,
            energy_above_minimum_j_per_mol: energy.map(|value| value - minimum_finite_energy),
            reason: selected_phase_reason(
                model,
                *fraction,
                energy,
                &thermal_treatment,
                &system.thermal_treatment_profile,
                phase_boundaries,
            ),
        });
    }
    let phase_compositions = distribute_components_between_phases(&composition, &phase_fractions)?;
    let phases = phase_fractions
        .into_iter()
        .zip(phase_compositions)
        .filter(|((_, fraction), _)| *fraction > TRACE_COMPONENT_FRACTION)
        .map(
            |((model, fraction), phase_composition)| MetallurgicalPhaseAmount {
                phase_id: model.id.clone(),
                kind: model.kind,
                fraction,
                composition: phase_composition,
                property_model: model.property_model.clone(),
                kinetic_model: model.kinetic_model.clone(),
            },
        )
        .collect::<Vec<_>>();
    let diffusion_state = estimate_diffusion_state(
        &phases,
        alloy.temperature_kelvin,
        &thermal_treatment,
        &system.thermal_treatment_profile,
        previous,
        delta_seconds,
    )?;
    let grain_structure = estimate_grain_structure(
        &thermal_treatment,
        &diffusion_state,
        &system.thermal_treatment_profile,
        previous,
    );
    let mechanical_history = estimate_mechanical_history(
        &thermal_treatment,
        &diffusion_state,
        &system.thermal_treatment_profile,
        previous,
        alloy.temperature_kelvin,
        delta_seconds,
    );
    let defect_state = estimate_defect_state(
        &thermal_treatment,
        &diffusion_state,
        &system.thermal_treatment_profile,
        &mechanical_history,
        previous,
    );
    let properties = estimate_properties(
        &phases,
        &grain_structure,
        &defect_state,
        &diffusion_state,
        &system.property_calibration,
    )?;
    let service_properties = estimate_service_properties(
        &properties,
        &phases,
        &grain_structure,
        &defect_state,
        &diffusion_state,
        Some(phase_boundaries),
    )?;
    let use_profile = estimate_use_profile(&properties, &service_properties)?;
    let diagnostics = MetallurgicalDiagnosticReport {
        selected_system_id: Some(system.id.clone()),
        considered_systems,
        generated_system,
        phase_boundaries: Some(phase_boundaries),
        phase_reasons: phase_diagnostics,
        thermal_reason: thermal_diagnostic(
            alloy,
            &thermal_treatment,
            &system.thermal_treatment_profile,
            delta_seconds,
        ),
        unmodeled_reason: None,
    };
    Ok(MetallurgicalState {
        kind: MetallurgicalStateKind::Modeled {
            system_id: system.id.clone(),
        },
        composition,
        temperature_kelvin: alloy.temperature_kelvin,
        phase_boundaries: Some(phase_boundaries),
        phases,
        grain_structure,
        defect_state,
        mechanical_history,
        diffusion_state,
        thermal_treatment,
        property_calibration: system.property_calibration.clone(),
        properties,
        service_properties,
        use_profile,
        diagnostics,
    })
}

fn unmodeled_state(
    alloy: &AlloyPhaseSnapshot,
    composition: MetallurgicalComposition,
    thermal_treatment: ThermalTreatmentState,
    considered_systems: Vec<MetallurgicalSystemSelectionDiagnostic>,
    delta_seconds: f64,
    reason: impl Into<String>,
) -> MetallurgicalState {
    let reason = reason.into();
    let diagnostics = MetallurgicalDiagnosticReport {
        selected_system_id: None,
        considered_systems,
        generated_system: None,
        phase_boundaries: None,
        phase_reasons: Vec::new(),
        thermal_reason: thermal_diagnostic(
            alloy,
            &thermal_treatment,
            &ThermalTreatmentProfile::neutral(),
            delta_seconds,
        ),
        unmodeled_reason: Some(reason.clone()),
    };
    MetallurgicalState {
        kind: MetallurgicalStateKind::Unmodeled {
            reason: reason.clone(),
        },
        composition,
        temperature_kelvin: alloy.temperature_kelvin,
        phase_boundaries: None,
        phases: Vec::new(),
        grain_structure: GrainStructure {
            average_grain_size_micrometers: DEFAULT_GRAIN_SIZE_MICROMETERS,
            distribution_width: 0.3,
        },
        defect_state: DefectState {
            vacancy_fraction: 0.0,
            dislocation_density_per_square_meter: 1.0e10,
            cold_work_fraction: 0.0,
        },
        mechanical_history: MechanicalHistoryState::initial(alloy.temperature_kelvin),
        diffusion_state: DiffusionState {
            effective_diffusivity_square_meters_per_second: 0.0,
            diffusion_length_micrometers: 0.0,
            homogenization_fraction: 0.0,
            precipitation_fraction: 0.0,
            aging_fraction: 0.0,
        },
        thermal_treatment,
        property_calibration: MetallurgicalPropertyCalibration::neutral(),
        properties: AlloyPropertySnapshot {
            hardness_hv: 0.0,
            yield_strength_mpa: 0.0,
            ductility_fraction: 0.0,
            electrical_resistivity_micro_ohm_meter: 0.0,
            thermal_conductivity_w_per_meter_kelvin: 0.0,
            corrosion_resistance_score: 0.0,
        },
        service_properties: AlloyServicePropertySnapshot {
            fracture_toughness_mpa_sqrt_meter: 0.0,
            brittleness_score: 1.0,
            wear_resistance_score: 0.0,
            electrical_conductivity_percent_iacs: 0.0,
            high_temperature_stability_score: 0.0,
            softening_temperature_kelvin: 0.0,
        },
        use_profile: AlloyUseProfile {
            suitability: Vec::new(),
        },
        diagnostics,
    }
}

fn missing_system_components(
    system: &MetallurgicalSystem,
    composition: &MetallurgicalComposition,
) -> Vec<MetallurgicalComponentId> {
    composition
        .components
        .keys()
        .filter(|component| !system.components.contains(*component))
        .cloned()
        .collect()
}

fn distribute_components_between_phases(
    composition: &MetallurgicalComposition,
    phase_fractions: &[(&MetallurgicalPhaseModel, f64)],
) -> ChemistryResult<Vec<BTreeMap<MetallurgicalComponentId, f64>>> {
    validate_phase_fraction_sum(phase_fractions)?;
    composition.validate()?;

    let components = composition.components.keys().cloned().collect::<Vec<_>>();
    if components.is_empty() || phase_fractions.is_empty() {
        return Err(ChemistryError::InvalidMixtureState(
            "phase component distribution requires phases and components".to_string(),
        ));
    }

    let mut matrix = Vec::with_capacity(phase_fractions.len());
    for (phase, phase_fraction) in phase_fractions {
        validate_fraction(*phase_fraction, "metallurgical phase fraction")?;
        let mut row = Vec::with_capacity(components.len());
        for component in &components {
            let component_fraction = composition.fraction_of(component);
            let affinity = phase_component_affinity(phase, component);
            row.push((phase_fraction * component_fraction * affinity).max(1.0e-30));
        }
        matrix.push(row);
    }

    for _ in 0..96 {
        for component_index in 0..components.len() {
            let target = composition.fraction_of(&components[component_index]);
            let current = matrix.iter().map(|row| row[component_index]).sum::<f64>();
            if current <= 0.0 || !current.is_finite() {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "component '{}' cannot be distributed between metallurgical phases",
                    components[component_index].as_str()
                )));
            }
            let scale = target / current;
            for row in &mut matrix {
                row[component_index] *= scale;
            }
        }
        for (phase_index, (_, target)) in phase_fractions.iter().enumerate() {
            let current = matrix[phase_index].iter().sum::<f64>();
            if current <= 0.0 || !current.is_finite() {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "metallurgical phase '{}' cannot receive a valid component distribution",
                    phase_fractions[phase_index].0.id
                )));
            }
            let scale = target / current;
            for value in &mut matrix[phase_index] {
                *value *= scale;
            }
        }
    }

    validate_phase_component_distribution(composition, phase_fractions, &components, &matrix)?;

    Ok(matrix
        .into_iter()
        .zip(phase_fractions.iter())
        .map(|(row, (_, phase_fraction))| {
            if *phase_fraction <= TRACE_COMPONENT_FRACTION {
                return BTreeMap::new();
            }
            components
                .iter()
                .cloned()
                .zip(row)
                .filter_map(|(component, amount)| {
                    let fraction_in_phase = amount / phase_fraction;
                    (fraction_in_phase > TRACE_COMPONENT_FRACTION)
                        .then_some((component, fraction_in_phase))
                })
                .collect()
        })
        .collect())
}

fn phase_component_affinity(
    phase: &MetallurgicalPhaseModel,
    component: &MetallurgicalComponentId,
) -> f64 {
    let mut affinity: f64 = match phase.kind {
        MetallurgicalPhaseKind::Liquid
        | MetallurgicalPhaseKind::SolidSolution
        | MetallurgicalPhaseKind::Ferrite
        | MetallurgicalPhaseKind::Austenite
        | MetallurgicalPhaseKind::Martensite
        | MetallurgicalPhaseKind::Pearlite
        | MetallurgicalPhaseKind::Bainite
        | MetallurgicalPhaseKind::TemperedMartensite => 1.0,
        MetallurgicalPhaseKind::Intermetallic
        | MetallurgicalPhaseKind::Cementite
        | MetallurgicalPhaseKind::Graphite => 0.08,
    };

    for term in &phase.free_energy_model.composition_terms {
        if &term.component == component {
            affinity = affinity.max(1.0 + 80.0 * term.center_fraction.max(0.01));
        }
    }
    for limit in &phase.component_limits {
        if &limit.component == component {
            let center = ((limit.min_fraction + limit.max_fraction) * 0.5).max(0.01);
            affinity = affinity.max(1.0 + 40.0 * center);
        }
    }
    if matches!(phase.kind, MetallurgicalPhaseKind::Graphite)
        && component.as_str() == "destroy:carbon"
    {
        affinity = affinity.max(100.0);
    }
    affinity
}

fn validate_phase_component_distribution(
    composition: &MetallurgicalComposition,
    phase_fractions: &[(&MetallurgicalPhaseModel, f64)],
    components: &[MetallurgicalComponentId],
    matrix: &[Vec<f64>],
) -> ChemistryResult<()> {
    for (phase_index, row) in matrix.iter().enumerate() {
        for value in row {
            if !value.is_finite() || *value < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "metallurgical phase '{}' received invalid component amount {value}",
                    phase_fractions[phase_index].0.id
                )));
            }
        }
        let row_total = row.iter().sum::<f64>();
        let target = phase_fractions[phase_index].1;
        if (row_total - target).abs() > 1.0e-8 {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "metallurgical phase '{}' component total does not match phase fraction: {row_total} vs {target}",
                phase_fractions[phase_index].0.id
            )));
        }
    }
    for (component_index, component) in components.iter().enumerate() {
        let column_total = matrix.iter().map(|row| row[component_index]).sum::<f64>();
        let target = composition.fraction_of(component);
        if (column_total - target).abs() > 1.0e-8 {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "metallurgical component '{}' is not conserved across phases: {column_total} vs {target}",
                component.as_str()
            )));
        }
    }
    Ok(())
}

fn system_distance_to_composition(
    system: &MetallurgicalSystem,
    composition: &MetallurgicalComposition,
) -> ChemistryResult<f64> {
    let distance = system
        .phase_boundaries
        .iter()
        .map(|point| composition_distance(&composition.components, &point.composition))
        .fold(f64::INFINITY, f64::min);
    validate_finite(distance, "metallurgical system composition distance")?;
    Ok(distance)
}

fn phase_limits_match(
    phase: &MetallurgicalPhaseModel,
    composition: &MetallurgicalComposition,
) -> bool {
    phase.component_limits.iter().all(|limit| {
        let fraction = composition.fraction_of(&limit.component);
        fraction >= limit.min_fraction && fraction <= limit.max_fraction
    })
}

fn phase_limit_mismatch_reason(
    phase: &MetallurgicalPhaseModel,
    composition: &MetallurgicalComposition,
) -> String {
    let mismatches = phase
        .component_limits
        .iter()
        .filter_map(|limit| {
            let fraction = composition.fraction_of(&limit.component);
            (fraction < limit.min_fraction || fraction > limit.max_fraction).then(|| {
                format!(
                    "{} fraction {:.6} outside {:.6}..={:.6}",
                    limit.component.as_str(),
                    fraction,
                    limit.min_fraction,
                    limit.max_fraction
                )
            })
        })
        .collect::<Vec<_>>();
    if mismatches.is_empty() {
        format!("phase '{}' was excluded by component limits", phase.id)
    } else {
        format!(
            "phase '{}' component limits do not match: {}",
            phase.id,
            mismatches.join(", ")
        )
    }
}

fn phase_energy_exclusion_reason(
    phase: &MetallurgicalPhaseModel,
    temperature_kelvin: f64,
    boundaries: PhaseBoundarySnapshot,
) -> String {
    match phase.kind {
        MetallurgicalPhaseKind::Liquid if temperature_kelvin < boundaries.solidus_kelvin => {
            format!(
                "liquid phase '{}' excluded because temperature {:.3} K is below solidus {:.3} K",
                phase.id, temperature_kelvin, boundaries.solidus_kelvin
            )
        }
        MetallurgicalPhaseKind::Liquid
            if temperature_kelvin > phase.free_energy_model.high_temperature_kelvin =>
        {
            format!(
                "liquid phase '{}' excluded because temperature {:.3} K is above phase range {:.3} K",
                phase.id, temperature_kelvin, phase.free_energy_model.high_temperature_kelvin
            )
        }
        _ if temperature_kelvin > boundaries.liquidus_kelvin => {
            format!(
                "solid phase '{}' excluded because temperature {:.3} K is above liquidus {:.3} K",
                phase.id, temperature_kelvin, boundaries.liquidus_kelvin
            )
        }
        _ => format!(
            "phase '{}' excluded because temperature {:.3} K is outside phase range {:.3}..={:.3} K",
            phase.id,
            temperature_kelvin,
            phase.free_energy_model.low_temperature_kelvin,
            phase.free_energy_model.high_temperature_kelvin
        ),
    }
}

fn selected_phase_reason(
    phase: &MetallurgicalPhaseModel,
    fraction: f64,
    energy: Option<f64>,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
    boundaries: PhaseBoundarySnapshot,
) -> String {
    let energy_text = energy
        .map(|value| format!("G={value:.3} J/mol"))
        .unwrap_or_else(|| "finite energy unavailable".to_string());
    let boundary_text = match phase.kind {
        MetallurgicalPhaseKind::Liquid => format!(
            "liquid is allowed above solidus {:.3} K and constrained by liquidus {:.3} K",
            boundaries.solidus_kelvin, boundaries.liquidus_kelvin
        ),
        _ => format!(
            "solid phase is allowed below liquidus {:.3} K",
            boundaries.liquidus_kelvin
        ),
    };
    let cooling_text = phase
        .free_energy_model
        .cooling_rate_stabilization_threshold_kelvin_per_second
        .map(|threshold| {
            if thermal_treatment.cooling_rate_kelvin_per_second >= threshold {
                format!(
                    "cooling rate {:.3} K/s stabilizes this phase over threshold {:.3} K/s",
                    thermal_treatment.cooling_rate_kelvin_per_second, threshold
                )
            } else {
                format!(
                    "cooling rate {:.3} K/s is below stabilization threshold {:.3} K/s",
                    thermal_treatment.cooling_rate_kelvin_per_second, threshold
                )
            }
        })
        .unwrap_or_else(|| "no cooling-rate stabilization term".to_string());
    let fraction_hint_text = phase
        .fraction_hint
        .as_ref()
        .map(|hint| {
            format!(
                "fraction hint target {:.6} with strength {:.3}: {}",
                hint.target_fraction, hint.strength, hint.reason
            )
        })
        .unwrap_or_else(|| "no phase-fraction hint".to_string());
    let treatment_text = phase_treatment_reason(phase.kind, thermal_treatment, treatment_profile)
        .unwrap_or_else(|| "thermal-treatment profile adds no phase-specific bias".to_string());
    format!(
        "phase '{}' selected fraction {:.6}; {}; {}; {}; {}; {}",
        phase.id,
        fraction,
        energy_text,
        boundary_text,
        cooling_text,
        fraction_hint_text,
        treatment_text
    )
}

fn thermal_diagnostic(
    alloy: &AlloyPhaseSnapshot,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
    delta_seconds: f64,
) -> MetallurgicalThermalDiagnostic {
    MetallurgicalThermalDiagnostic {
        previous_temperature_kelvin: thermal_treatment.previous_temperature_kelvin,
        current_temperature_kelvin: alloy.temperature_kelvin,
        cooling_rate_kelvin_per_second: thermal_treatment.cooling_rate_kelvin_per_second,
        hold_time_seconds: thermal_treatment.hold_time_seconds,
        delta_seconds,
        treatment_profile_id: treatment_profile.id.clone(),
        treatment_events: treatment_events(
            alloy.temperature_kelvin,
            thermal_treatment,
            treatment_profile,
        ),
    }
}

fn treatment_events(
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> Vec<String> {
    let mut events = Vec::new();
    if austenitized(thermal_treatment, treatment_profile) {
        events.push("austenitized before cooling".to_string());
    }
    if martensite_bias(temperature_kelvin, thermal_treatment, treatment_profile) > 0.0 {
        events.push("cooling path favors martensite".to_string());
    }
    if bainite_bias(temperature_kelvin, thermal_treatment, treatment_profile) > 0.0 {
        events.push("cooling path favors bainite".to_string());
    }
    if tempering_bias(temperature_kelvin, thermal_treatment, treatment_profile) > 0.0 {
        events.push("temperature and hold favor tempering".to_string());
    }
    if solution_treated(thermal_treatment, treatment_profile) {
        events.push("previous peak temperature reached solution-treatment range".to_string());
    }
    if aging_bias(temperature_kelvin, thermal_treatment, treatment_profile) > 0.0 {
        events.push("temperature and hold favor precipitation aging".to_string());
    }
    events
}

fn phase_free_energy(
    phase: &MetallurgicalPhaseModel,
    composition: &MetallurgicalComposition,
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
    phase_boundaries: Option<PhaseBoundarySnapshot>,
) -> ChemistryResult<f64> {
    validate_non_negative_finite(temperature_kelvin, "metallurgical temperature")?;
    let model = &phase.free_energy_model;
    validate_phase_free_energy_model(model)?;
    if let Some(boundaries) = phase_boundaries {
        match phase.kind {
            MetallurgicalPhaseKind::Liquid if temperature_kelvin < boundaries.solidus_kelvin => {
                return Ok(f64::INFINITY);
            }
            MetallurgicalPhaseKind::Liquid => {
                if temperature_kelvin > model.high_temperature_kelvin {
                    return Ok(f64::INFINITY);
                }
            }
            _ if temperature_kelvin > boundaries.liquidus_kelvin => return Ok(f64::INFINITY),
            _ => {
                if temperature_kelvin < model.low_temperature_kelvin
                    || temperature_kelvin > model.high_temperature_kelvin
                {
                    return Ok(f64::INFINITY);
                }
            }
        }
    } else if temperature_kelvin < model.low_temperature_kelvin
        || temperature_kelvin > model.high_temperature_kelvin
    {
        return Ok(f64::INFINITY);
    }
    let mut energy =
        model.reference_gibbs_j_per_mol - temperature_kelvin * model.entropy_j_per_mol_kelvin;
    for term in &model.composition_terms {
        let fraction = composition.fraction_of(&term.component);
        let normalized_distance = (fraction - term.center_fraction) / term.width_fraction;
        energy += term.penalty_j_per_mol * normalized_distance * normalized_distance;
    }
    if let Some(threshold) = model.cooling_rate_stabilization_threshold_kelvin_per_second {
        if thermal_treatment.cooling_rate_kelvin_per_second >= threshold {
            energy -= model.cooling_rate_stabilization_j_per_mol;
        }
    }
    Ok(energy)
}

fn phase_boundaries_for_composition(
    system: &MetallurgicalSystem,
    composition: &MetallurgicalComposition,
) -> ChemistryResult<PhaseBoundarySnapshot> {
    let mut weighted_solidus = 0.0;
    let mut weighted_liquidus = 0.0;
    let mut total_weight = 0.0;

    for point in &system.phase_boundaries {
        let distance = composition_distance(&composition.components, &point.composition);
        let weight = 1.0 / (distance * distance + 1.0e-12);
        weighted_solidus += weight * point.solidus_kelvin;
        weighted_liquidus += weight * point.liquidus_kelvin;
        total_weight += weight;
    }
    validate_positive_finite(total_weight, "phase-boundary interpolation weight")?;
    let snapshot = PhaseBoundarySnapshot {
        solidus_kelvin: weighted_solidus / total_weight,
        liquidus_kelvin: weighted_liquidus / total_weight,
    };
    validate_phase_boundary_temperatures(snapshot.solidus_kelvin, snapshot.liquidus_kelvin)?;
    Ok(snapshot)
}

fn composition_distance(
    left: &BTreeMap<MetallurgicalComponentId, f64>,
    right: &BTreeMap<MetallurgicalComponentId, f64>,
) -> f64 {
    let mut components = BTreeSet::new();
    components.extend(left.keys().cloned());
    components.extend(right.keys().cloned());
    components
        .into_iter()
        .map(|component| {
            let delta = left.get(&component).copied().unwrap_or(0.0)
                - right.get(&component).copied().unwrap_or(0.0);
            delta * delta
        })
        .sum::<f64>()
        .sqrt()
}

fn equilibrium_phase_fractions<'a>(
    phase_energies: &[(&'a MetallurgicalPhaseModel, f64)],
    temperature_kelvin: f64,
) -> ChemistryResult<Vec<(&'a MetallurgicalPhaseModel, f64)>> {
    validate_positive_finite(temperature_kelvin, "metallurgical temperature")?;
    let finite = phase_energies
        .iter()
        .copied()
        .filter(|(_, energy)| energy.is_finite())
        .collect::<Vec<_>>();
    if finite.is_empty() {
        return Err(ChemistryError::InvalidMixtureState(
            "metallurgical energy minimization has no finite candidate phases".to_string(),
        ));
    }
    let minimum_energy = finite
        .iter()
        .map(|(_, energy)| *energy)
        .fold(f64::INFINITY, f64::min);
    let thermal_energy = GAS_CONSTANT_J_PER_MOL_KELVIN * temperature_kelvin;
    let mut weights = Vec::new();
    for (phase, energy) in finite {
        let exponent = -((energy - minimum_energy) / thermal_energy).clamp(0.0, 80.0);
        weights.push((phase, exponent.exp()));
    }
    let total = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
    validate_positive_finite(total, "metallurgical equilibrium partition function")?;
    apply_phase_fraction_hints(
        weights
            .into_iter()
            .map(|(phase, weight)| (phase, weight / total))
            .collect(),
    )
}

fn apply_phase_fraction_hints<'a>(
    phase_fractions: Vec<(&'a MetallurgicalPhaseModel, f64)>,
) -> ChemistryResult<Vec<(&'a MetallurgicalPhaseModel, f64)>> {
    if phase_fractions
        .iter()
        .all(|(phase, _)| phase.fraction_hint.is_none())
    {
        validate_phase_fraction_sum(&phase_fractions)?;
        return Ok(phase_fractions);
    }

    let mut adjusted = Vec::with_capacity(phase_fractions.len());
    for (phase, base_fraction) in phase_fractions {
        let fraction = if let Some(hint) = &phase.fraction_hint {
            validate_fraction(hint.target_fraction, "phase fraction hint target")?;
            validate_non_negative_finite(hint.strength, "phase fraction hint strength")?;
            let blend = hint.strength / (1.0 + hint.strength);
            base_fraction * (1.0 - blend) + hint.target_fraction * blend
        } else {
            base_fraction
        };
        adjusted.push((phase, fraction.clamp(0.0, 1.0)));
    }
    normalize_phase_fractions(adjusted)
}

fn apply_thermal_treatment_bias<'a>(
    phase_fractions: Vec<(&'a MetallurgicalPhaseModel, f64)>,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> ChemistryResult<Vec<(&'a MetallurgicalPhaseModel, f64)>> {
    let mut biased = Vec::with_capacity(phase_fractions.len());
    for (phase, fraction) in phase_fractions {
        let multiplier =
            phase_treatment_multiplier(phase.kind, thermal_treatment, treatment_profile);
        biased.push((phase, fraction * multiplier));
    }
    normalize_phase_fractions(biased)
}

fn phase_treatment_multiplier(
    phase_kind: MetallurgicalPhaseKind,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> f64 {
    let temperature_kelvin = thermal_treatment.previous_temperature_kelvin;
    match phase_kind {
        MetallurgicalPhaseKind::Martensite => {
            1.0 + 6.0 * martensite_bias(temperature_kelvin, thermal_treatment, treatment_profile)
        }
        MetallurgicalPhaseKind::Bainite => {
            1.0 + 4.0 * bainite_bias(temperature_kelvin, thermal_treatment, treatment_profile)
        }
        MetallurgicalPhaseKind::TemperedMartensite => {
            1.0 + 5.0 * tempering_bias(temperature_kelvin, thermal_treatment, treatment_profile)
        }
        MetallurgicalPhaseKind::Intermetallic => {
            1.0 + 2.5 * aging_bias(temperature_kelvin, thermal_treatment, treatment_profile)
        }
        _ => 1.0,
    }
}

fn phase_treatment_reason(
    phase_kind: MetallurgicalPhaseKind,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> Option<String> {
    let multiplier = phase_treatment_multiplier(phase_kind, thermal_treatment, treatment_profile);
    (multiplier > 1.0 + 1.0e-9).then(|| {
        format!(
            "thermal-treatment profile '{}' applies phase multiplier {:.6}",
            treatment_profile.id, multiplier
        )
    })
}

fn relax_phase_fractions<'a>(
    equilibrium: Vec<(&'a MetallurgicalPhaseModel, f64)>,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
    previous: Option<&MetallurgicalState>,
    delta_seconds: f64,
) -> ChemistryResult<Vec<(&'a MetallurgicalPhaseModel, f64)>> {
    validate_positive_finite(delta_seconds, "metallurgical phase relaxation duration")?;
    if previous.is_none()
        || equilibrium
            .iter()
            .any(|(phase, _)| phase.kind == MetallurgicalPhaseKind::Liquid)
        || thermal_treatment.cooling_rate_kelvin_per_second > 50.0
    {
        validate_phase_fraction_sum(&equilibrium)?;
        return Ok(equilibrium);
    }

    let previous = previous.expect("checked above");
    let mut relaxed = Vec::with_capacity(equilibrium.len());
    for (phase, equilibrium_fraction) in equilibrium {
        validate_kinetic_model(&phase.kinetic_model)?;
        let previous_fraction = previous
            .phases
            .iter()
            .find(|previous_phase| previous_phase.phase_id == phase.id)
            .map(|previous_phase| previous_phase.fraction)
            .unwrap_or(0.0);
        let treatment_multiplier =
            phase_treatment_multiplier(phase.kind, thermal_treatment, treatment_profile);
        let relaxation = 1.0
            - (-(phase.kinetic_model.transformation_rate_per_second * treatment_multiplier)
                * delta_seconds)
                .exp();
        let fraction = previous_fraction
            + (equilibrium_fraction - previous_fraction) * relaxation.clamp(0.0, 1.0);
        relaxed.push((phase, fraction.max(0.0)));
    }
    normalize_phase_fractions(relaxed)
}

fn normalize_phase_fractions<'a>(
    phase_fractions: Vec<(&'a MetallurgicalPhaseModel, f64)>,
) -> ChemistryResult<Vec<(&'a MetallurgicalPhaseModel, f64)>> {
    let total = phase_fractions
        .iter()
        .map(|(_, fraction)| *fraction)
        .sum::<f64>();
    validate_positive_finite(total, "metallurgical relaxed phase fraction total")?;
    Ok(phase_fractions
        .into_iter()
        .map(|(phase, fraction)| (phase, fraction / total))
        .collect())
}

fn validate_phase_fraction_sum(
    phase_fractions: &[(&MetallurgicalPhaseModel, f64)],
) -> ChemistryResult<()> {
    let total = phase_fractions
        .iter()
        .map(|(_, fraction)| *fraction)
        .sum::<f64>();
    if !total.is_finite() || (total - 1.0).abs() > 1.0e-6 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "metallurgical phase fractions must sum to 1.0, got {total}"
        )));
    }
    Ok(())
}

fn austenitized(
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> bool {
    treatment_profile
        .austenitizing_temperature_kelvin
        .is_some_and(|temperature| thermal_treatment.peak_temperature_kelvin >= temperature)
}

fn solution_treated(
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> bool {
    treatment_profile
        .solution_temperature_kelvin
        .is_some_and(|temperature| thermal_treatment.peak_temperature_kelvin >= temperature)
}

fn martensite_bias(
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> f64 {
    let Some(start_temperature) = treatment_profile.martensite_start_kelvin else {
        return 0.0;
    };
    let Some(rate_threshold) = treatment_profile.martensite_cooling_rate_kelvin_per_second else {
        return 0.0;
    };
    if !austenitized(thermal_treatment, treatment_profile) {
        return 0.0;
    }
    let temperature_factor = if temperature_kelvin <= start_temperature {
        1.0
    } else {
        (1.0 - (temperature_kelvin - start_temperature) / start_temperature.max(1.0))
            .clamp(0.0, 1.0)
    };
    let rate_factor = if rate_threshold <= 0.0 {
        1.0
    } else {
        (thermal_treatment.cooling_rate_kelvin_per_second / rate_threshold).clamp(0.0, 2.0) / 2.0
    };
    (temperature_factor * rate_factor).clamp(0.0, 1.0)
}

fn bainite_bias(
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> f64 {
    let Some((low_temperature, high_temperature)) =
        treatment_profile.bainite_temperature_window_kelvin
    else {
        return 0.0;
    };
    let Some((low_rate, high_rate)) =
        treatment_profile.bainite_cooling_rate_window_kelvin_per_second
    else {
        return 0.0;
    };
    if !austenitized(thermal_treatment, treatment_profile)
        || temperature_kelvin < low_temperature
        || temperature_kelvin > high_temperature
    {
        return 0.0;
    }
    let rate = thermal_treatment.cooling_rate_kelvin_per_second;
    if rate < low_rate || rate > high_rate {
        return 0.0;
    }
    let temperature_center = 0.5 * (low_temperature + high_temperature);
    let temperature_half_width = 0.5 * (high_temperature - low_temperature).max(1.0);
    let temperature_factor = (1.0
        - ((temperature_kelvin - temperature_center).abs() / temperature_half_width))
        .clamp(0.0, 1.0);
    let rate_center = 0.5 * (low_rate + high_rate);
    let rate_half_width = 0.5 * (high_rate - low_rate).max(1.0e-9);
    let rate_factor = (1.0 - ((rate - rate_center).abs() / rate_half_width)).clamp(0.0, 1.0);
    (temperature_factor * rate_factor).clamp(0.0, 1.0)
}

fn tempering_bias(
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> f64 {
    let Some((low_temperature, high_temperature)) =
        treatment_profile.tempering_temperature_window_kelvin
    else {
        return 0.0;
    };
    if temperature_kelvin < low_temperature || temperature_kelvin > high_temperature {
        return 0.0;
    }
    let temperature_factor = ((temperature_kelvin - low_temperature)
        / (high_temperature - low_temperature).max(1.0))
    .clamp(0.0, 1.0);
    let time_factor = 1.0 - (-(thermal_treatment.hold_time_seconds / 1800.0)).exp();
    (temperature_factor * time_factor).clamp(0.0, 1.0)
}

fn aging_bias(
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
) -> f64 {
    let Some((low_temperature, high_temperature)) =
        treatment_profile.aging_temperature_window_kelvin
    else {
        return 0.0;
    };
    if !solution_treated(thermal_treatment, treatment_profile)
        || temperature_kelvin < low_temperature
        || temperature_kelvin > high_temperature
    {
        return 0.0;
    }
    let center = 0.5 * (low_temperature + high_temperature);
    let half_width = 0.5 * (high_temperature - low_temperature).max(1.0);
    let temperature_factor =
        (1.0 - ((temperature_kelvin - center).abs() / half_width)).clamp(0.0, 1.0);
    let time_factor = 1.0 - (-(thermal_treatment.hold_time_seconds / 7200.0)).exp();
    (temperature_factor * time_factor).clamp(0.0, 1.0)
}

fn estimate_diffusion_state(
    phases: &[MetallurgicalPhaseAmount],
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
    treatment_profile: &ThermalTreatmentProfile,
    previous: Option<&MetallurgicalState>,
    delta_seconds: f64,
) -> ChemistryResult<DiffusionState> {
    validate_positive_finite(temperature_kelvin, "diffusion temperature")?;
    validate_positive_finite(delta_seconds, "diffusion duration")?;
    let mut effective_diffusivity = 0.0;
    let mut precipitation_rate = 0.0;
    let mut recovery_rate = 0.0;
    let has_liquid = phases
        .iter()
        .any(|phase| phase.kind == MetallurgicalPhaseKind::Liquid && phase.fraction > 0.0);
    for phase in phases {
        validate_kinetic_model(&phase.kinetic_model)?;
        effective_diffusivity += phase.fraction
            * phase
                .kinetic_model
                .diffusivity_at_kelvin(temperature_kelvin)?;
        precipitation_rate += phase.fraction * phase.kinetic_model.precipitation_rate_per_second;
        recovery_rate += phase.fraction * phase.kinetic_model.recovery_rate_per_second;
    }
    validate_non_negative_finite(effective_diffusivity, "effective metallurgical diffusivity")?;
    let previous_length = previous
        .map(|state| state.diffusion_state.diffusion_length_micrometers)
        .unwrap_or(0.0);
    let added_length_square_micrometers = 2.0 * effective_diffusivity * delta_seconds * 1.0e12;
    validate_non_negative_finite(
        added_length_square_micrometers,
        "metallurgical diffusion distance increment",
    )?;
    let diffusion_length =
        (previous_length * previous_length + added_length_square_micrometers).sqrt();
    let homogenization_fraction =
        1.0 - (-diffusion_length / DEFAULT_HOMOGENIZATION_LENGTH_MICROMETERS).exp();

    let previous_precipitation = previous
        .map(|state| state.diffusion_state.precipitation_fraction)
        .unwrap_or(0.0);
    let previous_aging = previous
        .map(|state| state.diffusion_state.aging_fraction)
        .unwrap_or(0.0);
    let aging_bias = aging_bias(temperature_kelvin, thermal_treatment, treatment_profile);
    precipitation_rate *=
        1.0 + aging_bias * treatment_profile.precipitation_strength_multiplier.max(0.0);
    recovery_rate *= treatment_profile.recovery_multiplier.max(0.0);

    let precipitation_fraction = if has_liquid {
        previous_precipitation * (-recovery_rate.max(1.0) * delta_seconds).exp()
    } else {
        previous_precipitation
            + (1.0 - previous_precipitation) * (1.0 - (-precipitation_rate * delta_seconds).exp())
    }
    .clamp(0.0, 1.0);
    let aging_fraction = if has_liquid {
        0.0
    } else {
        previous_aging
            + (1.0 - previous_aging)
                * homogenization_fraction
                * precipitation_fraction
                * (1.0 + aging_bias)
    }
    .clamp(0.0, 1.0);

    Ok(DiffusionState {
        effective_diffusivity_square_meters_per_second: effective_diffusivity,
        diffusion_length_micrometers: diffusion_length,
        homogenization_fraction: homogenization_fraction.clamp(0.0, 1.0),
        precipitation_fraction,
        aging_fraction,
    })
}

fn estimate_grain_structure(
    thermal_treatment: &ThermalTreatmentState,
    diffusion_state: &DiffusionState,
    treatment_profile: &ThermalTreatmentProfile,
    previous: Option<&MetallurgicalState>,
) -> GrainStructure {
    let previous_size = previous
        .map(|state| state.grain_structure.average_grain_size_micrometers)
        .unwrap_or(DEFAULT_GRAIN_SIZE_MICROMETERS);
    let cooling_refinement = 1.0 / (1.0 + thermal_treatment.cooling_rate_kelvin_per_second / 50.0);
    let hold_growth = 1.0
        + thermal_treatment.hold_time_seconds.min(3600.0) / 3600.0
            * diffusion_state.homogenization_fraction
            * treatment_profile.grain_growth_multiplier.max(0.0);
    let temperature_growth = if thermal_treatment.previous_temperature_kelvin > 1000.0 {
        thermal_treatment.previous_temperature_kelvin / 1000.0
    } else {
        1.0
    };
    GrainStructure {
        average_grain_size_micrometers: (previous_size
            * hold_growth
            * temperature_growth
            * cooling_refinement.max(0.15))
        .clamp(1.0, 5000.0),
        distribution_width: (0.2 + thermal_treatment.cooling_rate_kelvin_per_second / 500.0)
            .clamp(0.1, 1.0),
    }
}

fn estimate_defect_state(
    thermal_treatment: &ThermalTreatmentState,
    diffusion_state: &DiffusionState,
    treatment_profile: &ThermalTreatmentProfile,
    mechanical_history: &MechanicalHistoryState,
    previous: Option<&MetallurgicalState>,
) -> DefectState {
    let previous_dislocation = previous
        .map(|state| state.defect_state.dislocation_density_per_square_meter)
        .unwrap_or(1.0e10);
    let quench_factor = ((thermal_treatment.cooling_rate_kelvin_per_second / 200.0)
        * treatment_profile.quench_vacancy_multiplier.max(0.0))
    .clamp(0.0, 1.0);
    let recovery_factor = (1.0
        - diffusion_state.homogenization_fraction
            * treatment_profile.recovery_multiplier.max(0.0)
            * thermal_treatment.hold_time_seconds
            / 7200.0)
        .clamp(0.2, 1.0);
    let previous_cold_work = previous
        .map(|state| state.defect_state.cold_work_fraction)
        .unwrap_or(0.0);
    let cold_work_recovery = mechanical_history.recrystallized_fraction
        + diffusion_state.homogenization_fraction
            * treatment_profile.recovery_multiplier.max(0.0)
            * thermal_treatment.hold_time_seconds
            / 7200.0;
    let cold_work_fraction =
        (previous_cold_work * (1.0 - cold_work_recovery.clamp(0.0, 1.0))).clamp(0.0, 1.0);
    DefectState {
        vacancy_fraction: (quench_factor * 1.0e-4).clamp(0.0, 1.0e-3),
        dislocation_density_per_square_meter: (previous_dislocation * recovery_factor
            + quench_factor * 5.0e13)
            .clamp(1.0e8, 1.0e16),
        cold_work_fraction,
    }
}

fn estimate_properties(
    phases: &[MetallurgicalPhaseAmount],
    grain_structure: &GrainStructure,
    defect_state: &DefectState,
    diffusion_state: &DiffusionState,
    calibration: &MetallurgicalPropertyCalibration,
) -> ChemistryResult<AlloyPropertySnapshot> {
    calibration.validate()?;
    let total = phases.iter().map(|phase| phase.fraction).sum::<f64>();
    if !total.is_finite() || (total - 1.0).abs() > 1.0e-6 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "metallurgical phase fractions must sum to 1.0, got {total}"
        )));
    }
    let mut properties = AlloyPropertySnapshot {
        hardness_hv: 0.0,
        yield_strength_mpa: 0.0,
        ductility_fraction: 0.0,
        electrical_resistivity_micro_ohm_meter: 0.0,
        thermal_conductivity_w_per_meter_kelvin: 0.0,
        corrosion_resistance_score: 0.0,
    };
    for phase in phases {
        let f = phase.fraction;
        let model = &phase.property_model;
        properties.hardness_hv += f * model.hardness_hv;
        properties.yield_strength_mpa += f * model.yield_strength_mpa;
        properties.ductility_fraction += f * model.ductility_fraction;
        properties.electrical_resistivity_micro_ohm_meter +=
            f * model.electrical_resistivity_micro_ohm_meter;
        properties.thermal_conductivity_w_per_meter_kelvin +=
            f * model.thermal_conductivity_w_per_meter_kelvin;
        properties.corrosion_resistance_score += f * model.corrosion_resistance_score;
    }
    let grain_strengthening = calibration.hall_petch_mpa_sqrt_micrometer
        / grain_structure.average_grain_size_micrometers.sqrt();
    let dislocation_strengthening = (defect_state.dislocation_density_per_square_meter / 1.0e12)
        .sqrt()
        * calibration.dislocation_strengthening_mpa_at_1e12;
    let precipitation_strengthening =
        calibration.precipitation_strengthening_mpa * diffusion_state.aging_fraction;
    let cold_work_strengthening =
        calibration.cold_work_strengthening_mpa * defect_state.cold_work_fraction;
    properties.yield_strength_mpa += grain_strengthening
        + dislocation_strengthening
        + precipitation_strengthening
        + cold_work_strengthening;
    properties.hardness_hv += (grain_strengthening
        + dislocation_strengthening
        + precipitation_strengthening
        + cold_work_strengthening)
        * calibration.hardness_per_strength_mpa;
    properties.ductility_fraction = (properties.ductility_fraction
        * (1.0
            - defect_state.vacancy_fraction * calibration.vacancy_ductility_penalty_per_fraction))
        * (1.0
            - calibration.precipitation_ductility_penalty * diffusion_state.precipitation_fraction)
            .clamp(0.0, 1.0)
        * (1.0 - calibration.cold_work_ductility_penalty * defect_state.cold_work_fraction)
            .clamp(0.0, 1.0);
    properties.electrical_resistivity_micro_ohm_meter += calibration
        .resistivity_precipitation_penalty_micro_ohm_meter
        * diffusion_state.precipitation_fraction
        + calibration.resistivity_cold_work_penalty_micro_ohm_meter
            * defect_state.cold_work_fraction;
    let thermal_penalty = 1.0
        + calibration.thermal_conductivity_precipitation_penalty
            * diffusion_state.precipitation_fraction
        + calibration.thermal_conductivity_defect_penalty
            * (defect_state.cold_work_fraction
                + (defect_state.dislocation_density_per_square_meter / 1.0e14)
                    .sqrt()
                    .clamp(0.0, 1.0));
    properties.thermal_conductivity_w_per_meter_kelvin /= thermal_penalty.max(1.0);
    validate_alloy_properties(&properties)?;
    Ok(properties)
}

fn estimate_service_properties(
    properties: &AlloyPropertySnapshot,
    phases: &[MetallurgicalPhaseAmount],
    grain_structure: &GrainStructure,
    defect_state: &DefectState,
    diffusion_state: &DiffusionState,
    phase_boundaries: Option<PhaseBoundarySnapshot>,
) -> ChemistryResult<AlloyServicePropertySnapshot> {
    validate_alloy_properties(properties)?;
    validate_positive_finite(
        grain_structure.average_grain_size_micrometers,
        "average metallurgical grain size",
    )?;
    validate_non_negative_finite(
        defect_state.dislocation_density_per_square_meter,
        "dislocation density",
    )?;

    let hard_phase_fraction = phases
        .iter()
        .filter(|phase| {
            matches!(
                phase.kind,
                MetallurgicalPhaseKind::Cementite
                    | MetallurgicalPhaseKind::Intermetallic
                    | MetallurgicalPhaseKind::Martensite
                    | MetallurgicalPhaseKind::TemperedMartensite
                    | MetallurgicalPhaseKind::Bainite
            )
        })
        .map(|phase| phase.fraction)
        .sum::<f64>()
        .clamp(0.0, 1.0);
    let liquid_fraction = phases
        .iter()
        .filter(|phase| phase.kind == MetallurgicalPhaseKind::Liquid)
        .map(|phase| phase.fraction)
        .sum::<f64>()
        .clamp(0.0, 1.0);
    let brittleness_score = (score_lower(properties.ductility_fraction, 0.28) * 0.42
        + score_higher(properties.hardness_hv, 850.0) * 0.18
        + hard_phase_fraction * 0.25
        + score_higher(defect_state.vacancy_fraction, 4.0e-4) * 0.15)
        .clamp(0.0, 1.0);
    let fracture_toughness_mpa_sqrt_meter = (12.0
        + properties.ductility_fraction.clamp(0.0, 1.0) * 180.0
        + score_higher(properties.yield_strength_mpa, 1200.0) * 35.0
        - brittleness_score * 70.0
        - hard_phase_fraction * 25.0)
        .clamp(2.0, 250.0);
    let wear_resistance_score = (score_higher(properties.hardness_hv, 900.0) * 0.55
        + hard_phase_fraction * 0.25
        + score_higher(properties.yield_strength_mpa, 1500.0) * 0.20
        - properties.ductility_fraction.clamp(0.0, 1.0) * 0.08)
        .clamp(0.0, 1.0);
    let electrical_conductivity_percent_iacs =
        if properties.electrical_resistivity_micro_ohm_meter > 0.0 {
            (100.0 * 0.017241 / properties.electrical_resistivity_micro_ohm_meter).clamp(0.0, 120.0)
        } else {
            0.0
        };
    let softening_temperature_kelvin = phase_boundaries
        .map(|boundaries| boundaries.solidus_kelvin * (0.58 + 0.16 * hard_phase_fraction))
        .unwrap_or(0.0)
        .max(0.0);
    let high_temperature_stability_score = (score_higher(softening_temperature_kelvin, 1500.0)
        * 0.42
        + score_higher(properties.yield_strength_mpa, 1200.0) * 0.22
        + hard_phase_fraction * 0.18
        + diffusion_state.precipitation_fraction.clamp(0.0, 1.0) * 0.10
        + properties.corrosion_resistance_score.clamp(0.0, 1.0) * 0.08)
        * (1.0 - liquid_fraction);

    let snapshot = AlloyServicePropertySnapshot {
        fracture_toughness_mpa_sqrt_meter,
        brittleness_score,
        wear_resistance_score,
        electrical_conductivity_percent_iacs,
        high_temperature_stability_score: high_temperature_stability_score.clamp(0.0, 1.0),
        softening_temperature_kelvin,
    };
    validate_service_properties(&snapshot)?;
    Ok(snapshot)
}

fn estimate_use_profile(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> ChemistryResult<AlloyUseProfile> {
    validate_alloy_properties(properties)?;
    validate_service_properties(service)?;

    let mut suitability = vec![
        suitability(
            MetallurgicalUseKind::Structural,
            structural_score(properties, service),
            structural_limiting_factor(properties, service),
        ),
        suitability(
            MetallurgicalUseKind::CuttingTool,
            cutting_tool_score(properties, service),
            cutting_tool_limiting_factor(properties, service),
        ),
        suitability(
            MetallurgicalUseKind::Spring,
            spring_score(properties, service),
            spring_limiting_factor(properties, service),
        ),
        suitability(
            MetallurgicalUseKind::ElectricalConductor,
            electrical_conductor_score(properties, service),
            electrical_limiting_factor(properties, service),
        ),
        suitability(
            MetallurgicalUseKind::ThermalConductor,
            thermal_conductor_score(properties, service),
            thermal_limiting_factor(properties, service),
        ),
        suitability(
            MetallurgicalUseKind::CorrosionResistant,
            corrosion_resistant_score(properties, service),
            corrosion_limiting_factor(properties, service),
        ),
        suitability(
            MetallurgicalUseKind::HighTemperature,
            high_temperature_score(properties, service),
            high_temperature_limiting_factor(properties, service),
        ),
        suitability(
            MetallurgicalUseKind::WearResistant,
            wear_resistant_score(properties, service),
            wear_limiting_factor(properties, service),
        ),
    ];
    suitability.sort_by(|left, right| {
        right
            .score
            .total_cmp(&left.score)
            .then_with(|| left.kind.cmp(&right.kind))
    });
    Ok(AlloyUseProfile { suitability })
}

fn suitability(
    kind: MetallurgicalUseKind,
    score: f64,
    limiting_factor: &'static str,
) -> MetallurgicalUseSuitability {
    MetallurgicalUseSuitability {
        kind,
        score: score.clamp(0.0, 1.0),
        limiting_factor: limiting_factor.to_string(),
    }
}

fn structural_score(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> f64 {
    weighted_score(&[
        (score_higher(properties.yield_strength_mpa, 900.0), 0.35),
        (
            score_higher(service.fracture_toughness_mpa_sqrt_meter, 90.0),
            0.25,
        ),
        (score_higher(properties.ductility_fraction, 0.22), 0.20),
        (1.0 - service.brittleness_score, 0.20),
    ])
}

fn cutting_tool_score(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> f64 {
    weighted_score(&[
        (score_higher(properties.hardness_hv, 900.0), 0.38),
        (service.wear_resistance_score, 0.30),
        (service.high_temperature_stability_score, 0.22),
        (
            1.0 - score_higher(properties.ductility_fraction, 0.35),
            0.10,
        ),
    ])
}

fn spring_score(properties: &AlloyPropertySnapshot, service: &AlloyServicePropertySnapshot) -> f64 {
    weighted_score(&[
        (score_higher(properties.yield_strength_mpa, 1100.0), 0.35),
        (score_higher(properties.ductility_fraction, 0.18), 0.25),
        (1.0 - service.brittleness_score, 0.25),
        (
            score_higher(service.fracture_toughness_mpa_sqrt_meter, 75.0),
            0.15,
        ),
    ])
}

fn electrical_conductor_score(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> f64 {
    weighted_score(&[
        (
            score_higher(service.electrical_conductivity_percent_iacs, 75.0),
            0.65,
        ),
        (score_higher(properties.ductility_fraction, 0.18), 0.15),
        (properties.corrosion_resistance_score.clamp(0.0, 1.0), 0.12),
        (1.0 - service.brittleness_score, 0.08),
    ])
}

fn thermal_conductor_score(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> f64 {
    weighted_score(&[
        (
            score_higher(properties.thermal_conductivity_w_per_meter_kelvin, 260.0),
            0.68,
        ),
        (service.high_temperature_stability_score, 0.15),
        (properties.corrosion_resistance_score.clamp(0.0, 1.0), 0.10),
        (1.0 - service.brittleness_score, 0.07),
    ])
}

fn corrosion_resistant_score(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> f64 {
    weighted_score(&[
        (properties.corrosion_resistance_score.clamp(0.0, 1.0), 0.60),
        (1.0 - service.brittleness_score, 0.15),
        (
            score_higher(service.fracture_toughness_mpa_sqrt_meter, 55.0),
            0.12,
        ),
        (score_higher(properties.yield_strength_mpa, 550.0), 0.13),
    ])
}

fn high_temperature_score(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> f64 {
    weighted_score(&[
        (service.high_temperature_stability_score, 0.45),
        (
            score_higher(service.softening_temperature_kelvin, 1300.0),
            0.22,
        ),
        (score_higher(properties.yield_strength_mpa, 1000.0), 0.18),
        (properties.corrosion_resistance_score.clamp(0.0, 1.0), 0.15),
    ])
}

fn wear_resistant_score(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> f64 {
    weighted_score(&[
        (service.wear_resistance_score, 0.50),
        (score_higher(properties.hardness_hv, 750.0), 0.25),
        (1.0 - service.brittleness_score, 0.15),
        (score_higher(properties.yield_strength_mpa, 900.0), 0.10),
    ])
}

fn weighted_score(parts: &[(f64, f64)]) -> f64 {
    let weight = parts.iter().map(|(_, weight)| *weight).sum::<f64>();
    if weight <= 0.0 || !weight.is_finite() {
        return 0.0;
    }
    (parts
        .iter()
        .map(|(score, weight)| score.clamp(0.0, 1.0) * weight)
        .sum::<f64>()
        / weight)
        .clamp(0.0, 1.0)
}

fn score_higher(value: f64, excellent: f64) -> f64 {
    if !value.is_finite() || !excellent.is_finite() || excellent <= 0.0 {
        0.0
    } else {
        (value / excellent).clamp(0.0, 1.0)
    }
}

fn score_lower(value: f64, poor: f64) -> f64 {
    if !value.is_finite() || !poor.is_finite() || poor <= 0.0 {
        1.0
    } else {
        (1.0 - value / poor).clamp(0.0, 1.0)
    }
}

fn structural_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if service.brittleness_score > 0.55 {
        "brittleness"
    } else if service.fracture_toughness_mpa_sqrt_meter < 45.0 {
        "fracture_toughness"
    } else if properties.yield_strength_mpa < 350.0 {
        "yield_strength"
    } else {
        "none"
    }
}

fn cutting_tool_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if properties.hardness_hv < 450.0 {
        "hardness"
    } else if service.high_temperature_stability_score < 0.35 {
        "hot_softening"
    } else if service.wear_resistance_score < 0.45 {
        "wear_resistance"
    } else {
        "none"
    }
}

fn spring_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if service.brittleness_score > 0.45 {
        "brittleness"
    } else if properties.yield_strength_mpa < 550.0 {
        "yield_strength"
    } else if properties.ductility_fraction < 0.08 {
        "ductility"
    } else {
        "none"
    }
}

fn electrical_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if service.electrical_conductivity_percent_iacs < 25.0 {
        "electrical_resistivity"
    } else if service.brittleness_score > 0.60 {
        "brittleness"
    } else if properties.corrosion_resistance_score < 0.25 {
        "corrosion"
    } else {
        "none"
    }
}

fn thermal_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if properties.thermal_conductivity_w_per_meter_kelvin < 45.0 {
        "thermal_conductivity"
    } else if service.high_temperature_stability_score < 0.20 {
        "hot_softening"
    } else {
        "none"
    }
}

fn corrosion_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if properties.corrosion_resistance_score < 0.45 {
        "corrosion"
    } else if service.brittleness_score > 0.65 {
        "brittleness"
    } else {
        "none"
    }
}

fn high_temperature_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if service.softening_temperature_kelvin < 900.0 {
        "softening_temperature"
    } else if properties.yield_strength_mpa < 450.0 {
        "yield_strength"
    } else if properties.corrosion_resistance_score < 0.35 {
        "corrosion"
    } else {
        "none"
    }
}

fn wear_limiting_factor(
    properties: &AlloyPropertySnapshot,
    service: &AlloyServicePropertySnapshot,
) -> &'static str {
    if properties.hardness_hv < 300.0 {
        "hardness"
    } else if service.brittleness_score > 0.75 {
        "brittleness"
    } else if service.wear_resistance_score < 0.35 {
        "wear_resistance"
    } else {
        "none"
    }
}

fn estimate_mechanical_history(
    thermal_treatment: &ThermalTreatmentState,
    diffusion_state: &DiffusionState,
    treatment_profile: &ThermalTreatmentProfile,
    previous: Option<&MetallurgicalState>,
    temperature_kelvin: f64,
    delta_seconds: f64,
) -> MechanicalHistoryState {
    let mut history = previous
        .map(|state| state.mechanical_history.clone())
        .unwrap_or_else(|| MechanicalHistoryState::initial(temperature_kelvin));
    history.recent_true_strain = 0.0;
    history.strain_rate_per_second = 0.0;
    let effective_recovery_time = thermal_treatment.hold_time_seconds.max(delta_seconds);
    let recovery = (diffusion_state.homogenization_fraction
        * treatment_profile.recovery_multiplier.max(0.0)
        * effective_recovery_time
        / 3600.0)
        .clamp(0.0, 1.0);
    history.accumulated_true_strain *= 1.0 - recovery;
    history.recrystallized_fraction = (history.recrystallized_fraction
        + (1.0 - history.recrystallized_fraction) * recovery)
        .clamp(0.0, 1.0);
    history.deformation_temperature_kelvin = temperature_kelvin;
    history
}

fn validate_mechanical_working_process(process: &MechanicalWorkingProcess) -> ChemistryResult<()> {
    validate_positive_finite(process.true_strain, "mechanical true strain")?;
    validate_non_negative_finite(process.strain_rate_per_second, "mechanical strain rate")?;
    validate_positive_finite(process.temperature_kelvin, "mechanical working temperature")?;
    validate_positive_finite(process.duration_seconds, "mechanical working duration")?;
    Ok(())
}

fn mechanical_work_factor(mode: MechanicalWorkingMode) -> f64 {
    match mode {
        MechanicalWorkingMode::Forging => 1.15,
        MechanicalWorkingMode::Rolling => 1.0,
        MechanicalWorkingMode::Drawing => 1.25,
        MechanicalWorkingMode::Extrusion => 1.1,
        MechanicalWorkingMode::Machining => 0.35,
    }
}

fn hot_work_softening_fraction(
    process: &MechanicalWorkingProcess,
    phase_boundaries: Option<PhaseBoundarySnapshot>,
) -> f64 {
    let Some(boundaries) = phase_boundaries else {
        return 0.0;
    };
    let recrystallization_start = 0.45 * boundaries.solidus_kelvin;
    if process.temperature_kelvin <= recrystallization_start {
        return 0.0;
    }
    let thermal_fraction = ((process.temperature_kelvin - recrystallization_start)
        / (boundaries.solidus_kelvin - recrystallization_start).max(1.0))
    .clamp(0.0, 1.0);
    let time_fraction = 1.0 - (-process.duration_seconds / 600.0).exp();
    (thermal_fraction * time_fraction).clamp(0.0, 0.95)
}
