use super::*;
use crate::chemistry::alloy::{alloy_phase_snapshots, AlloyPhaseSnapshot};
use crate::chemistry::metallurgy_data::{
    default_metallurgical_compound_phases, default_metallurgical_element_data,
    default_metallurgical_pair_interactions, default_metallurgical_systems,
};
use crate::chemistry::mixture::Mixture;
use crate::chemistry::registry::ChemistryRegistryBuilder;
use crate::chemistry::substance::{
    LiquidPhasePreference, SolventRole, Substance, SubstancePhaseProperties,
    SubstanceRepresentation,
};

#[test]
fn default_metallurgical_systems_are_valid() {
    for system in default_metallurgical_systems() {
        system.validate().unwrap();
    }
}

#[test]
fn default_metallurgical_element_data_are_valid() {
    for data in default_metallurgical_element_data() {
        data.validate().unwrap();
    }
}

#[test]
fn default_generated_metallurgy_data_are_valid() {
    for interaction in default_metallurgical_pair_interactions() {
        interaction.validate().unwrap();
    }
    for phase in default_metallurgical_compound_phases() {
        phase.validate().unwrap();
    }
}

#[test]
fn iron_carbon_melt_gets_modeled_liquid_state() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1900.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_iron", 0.98)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:carbon", 0.02)
        .unwrap();

    let alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Modeled { ref system_id } if system_id == "metallurgy:fe_c"
    ));
    assert!(state
        .phases
        .iter()
        .any(|phase| phase.kind == MetallurgicalPhaseKind::Liquid && phase.fraction > 0.5));
}

#[test]
fn unknown_gold_silver_alloy_uses_generated_metallurgical_system() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1400.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_gold", 0.55)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:test_silver", 0.45)
        .unwrap();

    let mut alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    alloy.temperature_kelvin = 1000.0;
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Modeled { ref system_id }
            if system_id == "metallurgy:generated/ag_au"
    ));
    assert!(state
        .phases
        .iter()
        .any(|phase| phase.kind == MetallurgicalPhaseKind::SolidSolution));
    assert!(state.properties.electrical_resistivity_micro_ohm_meter > 0.0);
    let generated = state
        .diagnostics
        .generated_system
        .as_ref()
        .expect("generated Au-Ag alloy must carry generator diagnostics");
    assert_eq!(generated.system_id, "metallurgy:generated/ag_au");
    assert!(generated
        .used_pair_interactions
        .iter()
        .any(|pair| pair == "Ag:Au"));
    assert!(generated.missing_pair_interactions.is_empty());
}

#[test]
fn exact_registered_system_takes_priority_over_generated_system() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1500.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_copper", 0.70)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:test_zinc", 0.30)
        .unwrap();

    let alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Modeled { ref system_id } if system_id == "metallurgy:cu_zn"
    ));
}

#[test]
fn generated_copper_beryllium_system_uses_specific_cube_phase() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1800.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_copper", 0.50)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:test_beryllium", 0.50)
        .unwrap();
    let mut alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    alloy.temperature_kelvin = 1200.0;
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Modeled { ref system_id }
            if system_id == "metallurgy:generated/be_cu"
    ));
    assert!(
        state
            .diagnostics
            .phase_reasons
            .iter()
            .any(|phase| phase.phase_id.contains("compound/cube")
                || phase.phase_id.contains("compound_cube")),
        "phase diagnostics: {:?}",
        state.diagnostics.phase_reasons
    );
}

#[test]
fn generated_tin_lead_system_uses_eutectic_solidus() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(700.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_tin", 0.62)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:test_lead", 0.38)
        .unwrap();
    let alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Modeled { ref system_id }
            if system_id == "metallurgy:generated/pb_sn"
    ));
    assert!(
        state
            .phase_boundaries
            .is_some_and(|boundaries| boundaries.solidus_kelvin <= 457.0),
        "phase boundaries: {:?}",
        state.phase_boundaries
    );
    let generated = state
        .diagnostics
        .generated_system
        .as_ref()
        .expect("generated Sn-Pb alloy must carry generator diagnostics");
    assert_eq!(generated.eutectic_temperature_kelvin, Some(456.0));
    assert!(
        generated.solidus_kelvin <= 457.0,
        "generated diagnostic solidus {}",
        generated.solidus_kelvin
    );
}

#[test]
fn fast_cooling_promotes_martensite_and_harder_state() {
    let registry = test_registry().build().unwrap();
    let mut hot = Mixture::new(1900.0).unwrap();
    hot.add_substance(&registry, "destroy:test_iron", 0.97)
        .unwrap();
    hot.add_substance(&registry, "destroy:carbon", 0.03)
        .unwrap();
    let hot_alloy = alloy_phase_snapshots(&registry, &hot).unwrap().remove(0);
    let hot_state = registry
        .metallurgical_state_from_alloy_phase(&hot_alloy, None, 1.0)
        .unwrap();

    let mut slow_alloy = hot_alloy.clone();
    slow_alloy.temperature_kelvin = 500.0;
    let slow_state = registry
        .metallurgical_state_from_alloy_phase(&slow_alloy, Some(&hot_state), 50.0)
        .unwrap();

    let fast_state = registry
        .metallurgical_state_from_alloy_phase(&slow_alloy, Some(&hot_state), 1.0)
        .unwrap();

    let fast_martensite = phase_fraction(&fast_state, MetallurgicalPhaseKind::Martensite);
    let slow_martensite = phase_fraction(&slow_state, MetallurgicalPhaseKind::Martensite);
    assert!(
        fast_martensite > slow_martensite,
        "fast martensite {fast_martensite}, slow martensite {slow_martensite}"
    );
    assert!(
        fast_state.properties.hardness_hv > slow_state.properties.hardness_hv,
        "fast hardness {}, slow hardness {}",
        fast_state.properties.hardness_hv,
        slow_state.properties.hardness_hv
    );
}

#[test]
fn slow_cooled_hypoeutectoid_steel_contains_ferrite_and_pearlite() {
    let registry = test_registry().build().unwrap();
    let hot_state = steel_state_from_temperature(&registry, 0.97, 0.03, 1900.0, None, 1.0);
    let cold_alloy = steel_alloy(&registry, 0.97, 0.03, 700.0);
    let cold_state = registry
        .metallurgical_state_from_alloy_phase(&cold_alloy, Some(&hot_state), 7200.0)
        .unwrap();

    let ferrite = summed_phase_fraction(&cold_state, MetallurgicalPhaseKind::Ferrite);
    let pearlite = summed_phase_fraction(&cold_state, MetallurgicalPhaseKind::Pearlite);

    assert!(ferrite > 0.0, "ferrite fraction {ferrite}");
    assert!(pearlite > 0.0, "pearlite fraction {pearlite}");
    assert_eq!(
        summed_phase_fraction(&cold_state, MetallurgicalPhaseKind::Martensite),
        0.0
    );
}

#[test]
fn intermediate_cooling_can_form_bainite_without_martensitic_quench() {
    let registry = test_registry().build().unwrap();
    let hot_state = steel_state_from_temperature(&registry, 0.97, 0.03, 1900.0, None, 1.0);
    let bainitic_alloy = steel_alloy(&registry, 0.97, 0.03, 720.0);
    let bainitic_state = registry
        .metallurgical_state_from_alloy_phase(&bainitic_alloy, Some(&hot_state), 120.0)
        .unwrap();

    let bainite = summed_phase_fraction(&bainitic_state, MetallurgicalPhaseKind::Bainite);
    let martensite = summed_phase_fraction(&bainitic_state, MetallurgicalPhaseKind::Martensite);

    assert!(bainite > 0.0, "bainite fraction {bainite}");
    assert!(
        bainite > martensite,
        "bainite {bainite}, martensite {martensite}"
    );
}

#[test]
fn tempering_quenched_steel_replaces_martensite_with_tempered_martensite() {
    let registry = test_registry().build().unwrap();
    let hot_state = steel_state_from_temperature(&registry, 0.97, 0.03, 1900.0, None, 1.0);
    let quenched_alloy = steel_alloy(&registry, 0.97, 0.03, 500.0);
    let quenched_state = registry
        .metallurgical_state_from_alloy_phase(&quenched_alloy, Some(&hot_state), 1.0)
        .unwrap();
    let tempered_alloy = steel_alloy(&registry, 0.97, 0.03, 780.0);
    let tempered_state = registry
        .metallurgical_state_from_alloy_phase(&tempered_alloy, Some(&quenched_state), 3600.0)
        .unwrap();

    let tempered =
        summed_phase_fraction(&tempered_state, MetallurgicalPhaseKind::TemperedMartensite);
    let quenched_martensite =
        summed_phase_fraction(&quenched_state, MetallurgicalPhaseKind::Martensite);
    let tempered_martensite =
        summed_phase_fraction(&tempered_state, MetallurgicalPhaseKind::Martensite);

    assert!(
        quenched_martensite > 0.0,
        "quenched martensite {quenched_martensite}"
    );
    assert!(tempered > 0.0, "tempered martensite {tempered}");
    assert!(
        tempered_martensite < quenched_martensite,
        "tempered martensite {tempered_martensite}, quenched {quenched_martensite}"
    );
    assert!(
        tempered_state.properties.ductility_fraction > quenched_state.properties.ductility_fraction,
        "tempered ductility {}, quenched {}",
        tempered_state.properties.ductility_fraction,
        quenched_state.properties.ductility_fraction
    );
}

#[test]
fn phase_boundaries_are_interpolated_from_system_data() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1900.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_iron", 0.83)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:carbon", 0.17)
        .unwrap();

    let alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();
    let boundaries = state
        .phase_boundaries
        .expect("modeled alloy must carry phase-boundary data");

    assert!((boundaries.solidus_kelvin - 1420.0).abs() < 1.0);
    assert!((boundaries.liquidus_kelvin - 1420.0).abs() < 1.0);
}

#[test]
fn phase_boundaries_constrain_liquid_and_solid_phases() {
    let registry = test_registry().build().unwrap();
    let mut hot = Mixture::new(1900.0).unwrap();
    hot.add_substance(&registry, "destroy:test_iron", 0.98)
        .unwrap();
    hot.add_substance(&registry, "destroy:carbon", 0.02)
        .unwrap();
    let hot_alloy = alloy_phase_snapshots(&registry, &hot).unwrap().remove(0);
    let hot_state = registry
        .metallurgical_state_from_alloy_phase(&hot_alloy, None, 1.0)
        .unwrap();
    assert_eq!(hot_state.phases.len(), 1);
    assert_eq!(hot_state.phases[0].kind, MetallurgicalPhaseKind::Liquid);

    let mut cold_alloy = hot_alloy;
    cold_alloy.temperature_kelvin = 500.0;
    let cold_state = registry
        .metallurgical_state_from_alloy_phase(&cold_alloy, Some(&hot_state), 1.0)
        .unwrap();
    assert_eq!(
        phase_fraction(&cold_state, MetallurgicalPhaseKind::Liquid),
        0.0
    );
}

#[test]
fn holding_at_temperature_accumulates_diffusion_and_aging() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1900.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_iron", 0.97)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:carbon", 0.03)
        .unwrap();
    let mut alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    alloy.temperature_kelvin = 950.0;

    let initial = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();
    let held = registry
        .metallurgical_state_from_alloy_phase(&alloy, Some(&initial), 3600.0)
        .unwrap();

    assert!(
        held.diffusion_state.diffusion_length_micrometers
            > initial.diffusion_state.diffusion_length_micrometers,
        "held diffusion length {}, initial {}",
        held.diffusion_state.diffusion_length_micrometers,
        initial.diffusion_state.diffusion_length_micrometers
    );
    assert!(
        held.diffusion_state.homogenization_fraction
            > initial.diffusion_state.homogenization_fraction,
        "held homogenization {}, initial {}",
        held.diffusion_state.homogenization_fraction,
        initial.diffusion_state.homogenization_fraction
    );
    assert!(
        held.diffusion_state.aging_fraction >= initial.diffusion_state.aging_fraction,
        "held aging {}, initial {}",
        held.diffusion_state.aging_fraction,
        initial.diffusion_state.aging_fraction
    );
}

#[test]
fn cold_mechanical_work_accumulates_dislocations_and_strengthens_alloy() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1900.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_iron", 0.97)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:carbon", 0.03)
        .unwrap();
    let mut alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    alloy.temperature_kelvin = 700.0;
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    let worked = apply_mechanical_working(
        &state,
        MechanicalWorkingProcess::new(
            MechanicalWorkingMode::Rolling,
            0.45,
            2.0,
            alloy.temperature_kelvin,
            30.0,
        ),
    )
    .unwrap();

    assert!(
        worked.defect_state.dislocation_density_per_square_meter
            > state.defect_state.dislocation_density_per_square_meter,
        "worked dislocation {}, initial {}",
        worked.defect_state.dislocation_density_per_square_meter,
        state.defect_state.dislocation_density_per_square_meter
    );
    assert!(
        worked.properties.yield_strength_mpa > state.properties.yield_strength_mpa,
        "worked strength {}, initial {}",
        worked.properties.yield_strength_mpa,
        state.properties.yield_strength_mpa
    );
    assert!(
        worked.properties.ductility_fraction < state.properties.ductility_fraction,
        "worked ductility {}, initial {}",
        worked.properties.ductility_fraction,
        state.properties.ductility_fraction
    );
}

#[test]
fn service_profile_is_derived_from_modeled_state_properties() {
    let registry = test_registry().build().unwrap();
    let state = nonferrous_state(
        &registry,
        [("destroy:test_copper", 0.72), ("destroy:test_nickel", 0.28)],
        1300.0,
        None,
        1.0,
    );

    assert!(
        state.service_properties.fracture_toughness_mpa_sqrt_meter > 0.0,
        "fracture toughness {}",
        state.service_properties.fracture_toughness_mpa_sqrt_meter
    );
    assert!(
        state
            .service_properties
            .electrical_conductivity_percent_iacs
            > 0.0,
        "electrical conductivity {}",
        state
            .service_properties
            .electrical_conductivity_percent_iacs
    );
    assert_eq!(state.use_profile.suitability.len(), 8);
    assert!(state
        .use_profile
        .suitability
        .iter()
        .all(|entry| entry.score.is_finite() && (0.0..=1.0).contains(&entry.score)));
}

#[test]
fn mechanical_work_updates_service_properties_and_use_profile() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1900.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_iron", 0.97)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:carbon", 0.03)
        .unwrap();
    let mut alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    alloy.temperature_kelvin = 700.0;
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    let worked = apply_mechanical_working(
        &state,
        MechanicalWorkingProcess::new(
            MechanicalWorkingMode::Rolling,
            0.45,
            2.0,
            alloy.temperature_kelvin,
            30.0,
        ),
    )
    .unwrap();

    assert!(
        worked.service_properties.wear_resistance_score
            >= state.service_properties.wear_resistance_score,
        "worked wear {}, initial wear {}",
        worked.service_properties.wear_resistance_score,
        state.service_properties.wear_resistance_score
    );
    assert!(
        use_score(&worked, MetallurgicalUseKind::WearResistant)
            >= use_score(&state, MetallurgicalUseKind::WearResistant),
        "worked use {:?}, initial use {:?}",
        worked.use_profile.suitability,
        state.use_profile.suitability
    );
}

#[test]
fn heat_holding_recovers_cold_work_without_manual_reset() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(1900.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_iron", 0.97)
        .unwrap();
    mixture
        .add_substance(&registry, "destroy:carbon", 0.03)
        .unwrap();
    let mut alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    alloy.temperature_kelvin = 700.0;
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();
    let worked = apply_mechanical_working(
        &state,
        MechanicalWorkingProcess::new(
            MechanicalWorkingMode::Forging,
            0.6,
            1.0,
            alloy.temperature_kelvin,
            20.0,
        ),
    )
    .unwrap();

    alloy.temperature_kelvin = 950.0;
    let recovered = registry
        .metallurgical_state_from_alloy_phase(&alloy, Some(&worked), 7200.0)
        .unwrap();

    assert!(
        recovered.defect_state.cold_work_fraction < worked.defect_state.cold_work_fraction,
        "recovered cold work {}, worked {}",
        recovered.defect_state.cold_work_fraction,
        worked.defect_state.cold_work_fraction
    );
    assert!(
        recovered.mechanical_history.recrystallized_fraction
            > worked.mechanical_history.recrystallized_fraction,
        "recovered recrystallized {}, worked {}",
        recovered.mechanical_history.recrystallized_fraction,
        worked.mechanical_history.recrystallized_fraction
    );
}

#[test]
fn aluminum_copper_magnesium_aging_forms_strengthening_precipitates() {
    let registry = test_registry().build().unwrap();
    let composition = [
        ("destroy:test_aluminum", 0.91),
        ("destroy:test_copper", 0.06),
        ("destroy:test_magnesium", 0.03),
    ];
    let solution = nonferrous_state(&registry, composition, 880.0, None, 1.0);
    let aged = nonferrous_state(&registry, composition, 460.0, Some(&solution), 18_000.0);

    assert!(matches!(
        aged.kind,
        MetallurgicalStateKind::Modeled { ref system_id } if system_id == "metallurgy:al_cu_mg"
    ));
    assert!(
        summed_phase_fraction(&aged, MetallurgicalPhaseKind::Intermetallic) > 0.0,
        "aged phases: {:?}",
        aged.phases
    );
    assert!(
        aged.properties.yield_strength_mpa > solution.properties.yield_strength_mpa,
        "aged strength {}, solution strength {}",
        aged.properties.yield_strength_mpa,
        solution.properties.yield_strength_mpa
    );
    assert!(
        aged.diffusion_state.aging_fraction > solution.diffusion_state.aging_fraction,
        "aged fraction {}, solution fraction {}",
        aged.diffusion_state.aging_fraction,
        solution.diffusion_state.aging_fraction
    );
}

#[test]
fn nickel_chromium_aluminum_superalloy_forms_gamma_prime() {
    let registry = test_registry().build().unwrap();
    let state = nonferrous_state(
        &registry,
        [
            ("destroy:test_nickel", 0.72),
            ("destroy:test_chromium", 0.16),
            ("destroy:test_aluminum", 0.12),
        ],
        1100.0,
        None,
        3_600.0,
    );

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Modeled { ref system_id } if system_id == "metallurgy:ni_cr_al"
    ));
    assert!(
        state.phases.iter().any(|phase| {
            phase.phase_id == "metallurgy:ni_cr_al/gamma_prime_ni3al" && phase.fraction > 0.0
        }),
        "superalloy phases: {:?}",
        state.phases
    );
    assert!(
        state.properties.corrosion_resistance_score > 0.85,
        "corrosion score {}",
        state.properties.corrosion_resistance_score
    );
}

#[test]
fn copper_nickel_uses_continuous_solid_solution() {
    let registry = test_registry().build().unwrap();
    let state = nonferrous_state(
        &registry,
        [("destroy:test_copper", 0.70), ("destroy:test_nickel", 0.30)],
        900.0,
        None,
        1.0,
    );

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Modeled { ref system_id } if system_id == "metallurgy:cu_ni"
    ));
    assert!(
        state.phases.iter().any(|phase| {
            phase.phase_id == "metallurgy:cu_ni/continuous_solid_solution" && phase.fraction > 0.5
        }),
        "cu-ni phases: {:?}",
        state.phases
    );
}

#[test]
fn magnesium_aluminum_zinc_alloy_can_age_by_intermetallic_precipitation() {
    let registry = test_registry().build().unwrap();
    let composition = [
        ("destroy:test_magnesium", 0.86),
        ("destroy:test_aluminum", 0.10),
        ("destroy:test_zinc", 0.04),
    ];
    let solution = nonferrous_state(&registry, composition, 760.0, None, 1.0);
    let aged = nonferrous_state(&registry, composition, 430.0, Some(&solution), 14_400.0);

    assert!(matches!(
        aged.kind,
        MetallurgicalStateKind::Modeled { ref system_id } if system_id == "metallurgy:mg_al_zn"
    ));
    assert!(
        summed_phase_fraction(&aged, MetallurgicalPhaseKind::Intermetallic) > 0.0,
        "aged magnesium alloy phases: {:?}",
        aged.phases
    );
    assert!(
        aged.properties.hardness_hv > solution.properties.hardness_hv,
        "aged hardness {}, solution hardness {}",
        aged.properties.hardness_hv,
        solution.properties.hardness_hv
    );
}

#[test]
fn modeled_state_carries_metallurgical_diagnostics() {
    let registry = test_registry().build().unwrap();
    let state = steel_state_from_temperature(&registry, 0.97, 0.03, 950.0, None, 1.0);

    assert_eq!(
        state.diagnostics.selected_system_id.as_deref(),
        Some("metallurgy:fe_c")
    );
    assert!(state.diagnostics.phase_boundaries.is_some());
    assert!(state
        .diagnostics
        .considered_systems
        .iter()
        .any(|system| system.system_id == "metallurgy:fe_c"
            && system.covers_composition
            && system.composition_distance.is_some()));
    assert!(
        state.diagnostics.phase_reasons.iter().any(|phase| {
            phase.selected && phase.fraction > 0.0 && phase.gibbs_j_per_mol.is_some()
        }),
        "phase diagnostics: {:?}",
        state.diagnostics.phase_reasons
    );
    assert_eq!(state.diagnostics.unmodeled_reason, None);
}

#[test]
fn unmodeled_state_reports_missing_metallurgical_components() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(2000.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_unknown_metal", 1.0)
        .unwrap();
    let alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    assert!(state.diagnostics.selected_system_id.is_none());
    assert!(state.diagnostics.unmodeled_reason.is_some());
    assert!(state.diagnostics.considered_systems.iter().all(|system| {
        !system.covers_composition
            && system
                .missing_components
                .iter()
                .any(|component| component.as_str() == "Xx")
    }));
}

#[test]
fn invalid_phase_kinetic_data_is_rejected() {
    let mut system = default_metallurgical_systems()
        .into_iter()
        .find(|system| system.id == "metallurgy:fe_c")
        .unwrap();
    system.phase_models[0]
        .kinetic_model
        .diffusion_prefactor_square_meters_per_second = f64::NAN;

    assert!(system.validate().is_err());
}

#[test]
fn unknown_metal_system_is_explicitly_unmodeled() {
    let registry = test_registry().build().unwrap();
    let mut mixture = Mixture::new(2000.0).unwrap();
    mixture
        .add_substance(&registry, "destroy:test_unknown_metal", 1.0)
        .unwrap();

    let alloy = alloy_phase_snapshots(&registry, &mixture)
        .unwrap()
        .remove(0);
    let state = registry
        .metallurgical_state_from_alloy_phase(&alloy, None, 1.0)
        .unwrap();

    assert!(matches!(
        state.kind,
        MetallurgicalStateKind::Unmodeled { .. }
    ));
    assert!(state.phases.is_empty());
}

fn phase_fraction(state: &MetallurgicalState, kind: MetallurgicalPhaseKind) -> f64 {
    state
        .phases
        .iter()
        .find(|phase| phase.kind == kind)
        .map(|phase| phase.fraction)
        .unwrap_or(0.0)
}

fn summed_phase_fraction(state: &MetallurgicalState, kind: MetallurgicalPhaseKind) -> f64 {
    state
        .phases
        .iter()
        .filter(|phase| phase.kind == kind)
        .map(|phase| phase.fraction)
        .sum()
}

fn use_score(state: &MetallurgicalState, kind: MetallurgicalUseKind) -> f64 {
    state
        .use_profile
        .suitability
        .iter()
        .find(|entry| entry.kind == kind)
        .map(|entry| entry.score)
        .unwrap_or(0.0)
}

fn steel_state_from_temperature(
    registry: &crate::chemistry::registry::ChemistryRegistry,
    iron_fraction: f64,
    carbon_fraction: f64,
    temperature_kelvin: f64,
    previous: Option<&MetallurgicalState>,
    delta_seconds: f64,
) -> MetallurgicalState {
    let alloy = steel_alloy(registry, iron_fraction, carbon_fraction, temperature_kelvin);
    registry
        .metallurgical_state_from_alloy_phase(&alloy, previous, delta_seconds)
        .unwrap()
}

fn steel_alloy(
    registry: &crate::chemistry::registry::ChemistryRegistry,
    iron_fraction: f64,
    carbon_fraction: f64,
    temperature_kelvin: f64,
) -> AlloyPhaseSnapshot {
    let mut mixture = Mixture::new(1900.0).unwrap();
    mixture
        .add_substance(registry, "destroy:test_iron", iron_fraction)
        .unwrap();
    mixture
        .add_substance(registry, "destroy:carbon", carbon_fraction)
        .unwrap();
    let mut alloy = alloy_phase_snapshots(registry, &mixture).unwrap().remove(0);
    alloy.temperature_kelvin = temperature_kelvin;
    alloy
}

fn nonferrous_state<const N: usize>(
    registry: &crate::chemistry::registry::ChemistryRegistry,
    components: [(&'static str, f64); N],
    temperature_kelvin: f64,
    previous: Option<&MetallurgicalState>,
    delta_seconds: f64,
) -> MetallurgicalState {
    let mut mixture = Mixture::new(1900.0).unwrap();
    for (substance_id, fraction) in components {
        mixture
            .add_substance(registry, substance_id, fraction)
            .unwrap();
    }
    let mut alloy = alloy_phase_snapshots(registry, &mixture).unwrap().remove(0);
    alloy.temperature_kelvin = temperature_kelvin;
    registry
        .metallurgical_state_from_alloy_phase(&alloy, previous, delta_seconds)
        .unwrap()
}

fn test_registry() -> ChemistryRegistryBuilder {
    ChemistryRegistryBuilder::new()
        .metallurgical_systems(default_metallurgical_systems())
        .metallurgical_elements(default_metallurgical_element_data())
        .metallurgical_pair_interactions(default_metallurgical_pair_interactions())
        .metallurgical_compound_phases(default_metallurgical_compound_phases())
        .substance(test_metal("destroy:test_iron", "Fe", 55.845, 1811.0))
        .substance(test_metal("destroy:test_lead", "Pb", 207.2, 600.61))
        .substance(test_metal("destroy:test_aluminum", "Al", 26.982, 933.0))
        .substance(test_metal("destroy:test_copper", "Cu", 63.546, 1358.0))
        .substance(test_metal("destroy:test_magnesium", "Mg", 24.305, 923.0))
        .substance(test_metal("destroy:test_zinc", "Zn", 65.38, 692.7))
        .substance(test_metal("destroy:test_nickel", "Ni", 58.693, 1728.0))
        .substance(test_metal("destroy:test_chromium", "Cr", 51.996, 2180.0))
        .substance(test_metal("destroy:test_gold", "Au", 196.967, 1337.0))
        .substance(test_metal("destroy:test_silver", "Ag", 107.868, 1235.0))
        .substance(test_metal("destroy:test_tin", "Sn", 118.71, 505.0))
        .substance(test_metal("destroy:test_beryllium", "Be", 9.0122, 1560.0))
        .substance(test_metal(
            "destroy:test_unknown_metal",
            "Xx",
            100.0,
            1200.0,
        ))
        .substance(
            Substance::new("destroy:carbon", 0, 12.011, 2_200.0, 4300.0, 8.5, 0.0)
                .with_melting_point_kelvin(1000.0)
                .with_phase_properties(molten_metal_phase_properties())
                .with_representation(SubstanceRepresentation::UnspecifiedMaterial {
                    reason: "test metallurgical carbon component".to_string(),
                }),
        )
}

fn test_metal(
    id: &'static str,
    element: &'static str,
    molar_mass: f64,
    melting_point: f64,
) -> Substance {
    Substance::new(id, 0, molar_mass, 7_800.0, 3300.0, 25.0, 0.0)
        .with_solid_density_grams_per_bucket(7_800.0)
        .with_melting_point_kelvin(melting_point)
        .with_phase_properties(molten_metal_phase_properties())
        .with_representation(SubstanceRepresentation::Metal {
            element_symbol: element.to_string(),
        })
}

fn molten_metal_phase_properties() -> SubstancePhaseProperties {
    SubstancePhaseProperties {
        preferred_liquid_phase: LiquidPhasePreference::MoltenMetal,
        aqueous_solubility_mol_per_bucket: Some(0.0),
        organic_solubility_mol_per_bucket: Some(0.0),
        can_precipitate: true,
        can_form_liquid_phase: true,
        solvent_role: SolventRole::NotSolvent,
    }
}
