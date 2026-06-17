use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::reactor::peripheral::{
    ControlMode, ElectrodeState, Peripheral, SmartHeaterState,
};
use crate::chemistry::reactor::reactor::{
    PhaseEntry, Reactor, SubstanceEntry, TransitionMode, ZoneId, ZoneTransition,
};
use crate::chemistry::reactor::zone::ReactorZone;
use crate::chemistry::substance::SubstanceId;
use crate::chemistry::ChemistryResult;

pub fn batch_reactor(volume_m3: f64, target_kelvin: f64) -> ChemistryResult<Reactor> {
    let mut reactor = Reactor::new();
    let zone = ReactorZone::new(volume_m3)?.with_peripheral(Peripheral::SmartHeater(
        SmartHeaterState::new(
            target_kelvin,
            10.0,
            293.15,
            0.0004,
            1673.0,
            230.0,
            ControlMode::PID,
        ),
    ));
    reactor.add_zone(zone);
    Ok(reactor)
}

pub fn cstr(volume_m3: f64, target_kelvin: f64) -> ChemistryResult<Reactor> {
    let mut reactor = Reactor::new();

    let input = reactor.add_zone(ReactorZone::new(volume_m3)?);
    let reaction = reactor.add_zone(ReactorZone::new(volume_m3)?.with_peripheral(
        Peripheral::SmartHeater(SmartHeaterState::new(
            target_kelvin,
            10.0,
            293.15,
            0.0004,
            1673.0,
            230.0,
            ControlMode::PID,
        )),
    ));
    let output = reactor.add_zone(ReactorZone::new(volume_m3)?);

    reactor.add_transition(ZoneTransition::new(
        input,
        reaction,
        TransitionMode::All {
            rate_mol_per_second: 0.5,
        },
    ));
    reactor.add_transition(ZoneTransition::new(
        reaction,
        output,
        TransitionMode::All {
            rate_mol_per_second: 0.5,
        },
    ));

    Ok(reactor)
}

pub fn electrolysis_cell(volume_m3: f64) -> ChemistryResult<Reactor> {
    let mut reactor = Reactor::new();

    let electrolyte = reactor.add_zone(
        ReactorZone::new(volume_m3)?
            .with_peripheral(Peripheral::SmartHeater(SmartHeaterState::new(
                353.0,
                5.0,
                293.15,
                0.0004,
                1673.0,
                12.0,
                ControlMode::P,
            )))
            .with_peripheral(Peripheral::Electrode(ElectrodeState::new(12.0, 50.0, 1.23))),
    );
    let product = reactor.add_zone(ReactorZone::new(volume_m3)?);

    reactor.add_transition(ZoneTransition::new(
        electrolyte,
        product,
        TransitionMode::Substances {
            entries: vec![
                SubstanceEntry {
                    id: SubstanceId::from("destroy:hydrogen"),
                    rate_mol_per_second: 0.1,
                },
                SubstanceEntry {
                    id: SubstanceId::from("destroy:oxygen"),
                    rate_mol_per_second: 0.1,
                },
            ],
        },
    ));

    Ok(reactor)
}

pub fn distillation_column(stages: usize, volume_per_stage_m3: f64) -> ChemistryResult<Reactor> {
    if stages == 0 {
        return Err(crate::chemistry::ChemistryError::InvalidMixtureState(
            "distillation column must have at least one stage".to_string(),
        ));
    }
    let mut reactor = Reactor::new();
    reactor.set_vle_iterations(2);
    reactor.set_vle_relaxation(1.0);

    let mut zone_ids: Vec<ZoneId> = Vec::new();
    for _ in 0..stages {
        let zone =
            ReactorZone::new(volume_per_stage_m3)?.with_peripheral(Peripheral::HeatExchanger {
                coolant_temperature_kelvin: 373.0,
                u_kw_per_kelvin: 0.5,
            });
        zone_ids.push(reactor.add_zone(zone));
    }

    for i in 0..stages - 1 {
        reactor.add_transition(ZoneTransition::new(
            zone_ids[i],
            zone_ids[i + 1],
            TransitionMode::Phases {
                entries: vec![PhaseEntry {
                    phase: MixturePhase::Gas,
                    rate_mol_per_second: 0.2,
                }],
            },
        ));
        reactor.add_transition(ZoneTransition::new(
            zone_ids[i + 1],
            zone_ids[i],
            TransitionMode::Phases {
                entries: vec![PhaseEntry {
                    phase: MixturePhase::Aqueous,
                    rate_mol_per_second: 0.1,
                }],
            },
        ));
    }

    Ok(reactor)
}

pub fn arc_furnace(volume_m3: f64) -> ChemistryResult<Reactor> {
    let mut reactor = Reactor::new();

    let chamber = reactor.add_zone(
        ReactorZone::new(volume_m3)?
            .with_peripheral(Peripheral::Electrode(ElectrodeState::new(80.0, 500.0, 0.0)))
            .with_peripheral(Peripheral::HeatExchanger {
                coolant_temperature_kelvin: 293.0,
                u_kw_per_kelvin: 2.0,
            }),
    );

    let output = reactor.add_zone(ReactorZone::new(volume_m3)?);

    reactor.add_transition(ZoneTransition::new(
        chamber,
        output,
        TransitionMode::Phases {
            entries: vec![PhaseEntry {
                phase: MixturePhase::MoltenMetal,
                rate_mol_per_second: 1.0,
            }],
        },
    ));

    Ok(reactor)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::destroy_registry_builder;
    use crate::chemistry::reactor::io;
    use crate::chemistry::reactor::zone::ReactorVolumeMode;
    use crate::chemistry::registry::ChemistryRegistry;

    fn liquid_amount_for_condensed_volume(
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
        condensed_volume_cubic_meters: f64,
    ) -> f64 {
        let substance = registry.substance_index(substance_id).unwrap();
        let properties = registry.substance_properties();
        condensed_volume_cubic_meters / crate::chemistry::mixture::DEFAULT_GAS_VOLUME_CUBIC_METERS
            * properties.liquid_density_grams_per_bucket[substance.as_usize()]
            / properties.molar_mass_grams[substance.as_usize()]
    }

    #[test]
    fn distillation_column_separates_ethanol_from_water() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let ethanol_id = SubstanceId::from("destroy:ethanol");

        let temps = [298.0, 353.0, 373.0, 400.0];

        for &temp in &temps {
            let mut reactor = distillation_column(5, 0.001).unwrap();

            // Fill zone 0 with 60% water, 40% ethanol without overfilling the stage.
            {
                let zone = reactor.zone_mut(&ZoneId(0)).unwrap();
                zone.add_substance_checked(&registry, water_id.clone(), 3.0)
                    .unwrap();
                zone.add_substance_checked(&registry, ethanol_id.clone(), 2.0)
                    .unwrap();
            }

            // Set coolant temperature for all zones
            for i in 0..5 {
                let zone = reactor.zone_mut(&ZoneId(i)).unwrap();
                zone.peripherals_mut().clear();
                zone.add_peripheral(Peripheral::HeatExchanger {
                    coolant_temperature_kelvin: temp,
                    u_kw_per_kelvin: 0.5,
                });
            }

            // Run for 100 ticks
            for _ in 0..100 {
                reactor.tick(&registry, 1.0).unwrap();
            }

            // Collect results
            let mut results: Vec<(String, f64, f64, f64)> = Vec::new();
            for i in 0..5 {
                let snapshot = io::mixture_snapshot(reactor.zone(&ZoneId(i)).unwrap());
                let w = snapshot
                    .substances
                    .iter()
                    .find(|s| s.id == water_id)
                    .map(|s| s.total_mol_per_bucket)
                    .unwrap_or(0.0);
                let e = snapshot
                    .substances
                    .iter()
                    .find(|s| s.id == ethanol_id)
                    .map(|s| s.total_mol_per_bucket)
                    .unwrap_or(0.0);
                results.push((format!("zone_{}", i), w, e, snapshot.temperature_kelvin));
            }

            // Print table
            println!("\n=== Distillation at coolant T = {} K ===", temp);
            println!(
                "{:<10} {:>12} {:>12} {:>12} {:>10}",
                "Zone", "Water", "Ethanol", "T (K)", "EtOH frac"
            );
            println!("{}", "-".repeat(62));
            for (name, w, e, t) in &results {
                let total = w + e;
                let frac = if total > 0.001 { e / total } else { 0.0 };
                println!(
                    "{:<10} {:>12.4} {:>12.4} {:>12.2} {:>10.4}",
                    name, w, e, t, frac
                );
            }

            // Ethanol fraction should increase from bottom (zone_0) to top (zone_4)
            let frac_bottom = {
                let total = results[0].1 + results[0].2;
                if total > 0.001 {
                    results[0].2 / total
                } else {
                    0.0
                }
            };
            let frac_top = {
                let total = results[4].1 + results[4].2;
                if total > 0.001 {
                    results[4].2 / total
                } else {
                    0.0
                }
            };

            println!(
                "\nEthanol fraction: bottom={:.4}, top={:.4}",
                frac_bottom, frac_top
            );

            if temp >= 353.0 {
                // Ethanol should have moved upward: top zone either has
                // more ethanol fraction, or ethanol left zone 0 entirely
                let ethanol_in_top_zones: f64 = results[1..].iter().map(|r| r.2).sum();
                let ethanol_in_bottom = results[0].2;

                println!(
                    "Ethanol in bottom: {:.4}, in top zones: {:.4}",
                    ethanol_in_bottom, ethanol_in_top_zones
                );

                assert!(
                    ethanol_in_top_zones > 0.0,
                    "At T={}: ethanol should have moved to upper zones",
                    temp
                );

                if temp < 400.0 {
                    let water_in_top_zones: f64 = results[1..].iter().map(|r| r.1).sum();
                    let water_in_bottom = results[0].1;
                    assert!(
                        water_in_bottom > water_in_top_zones,
                        "At T={}: water should stay mostly in bottom zone",
                        temp
                    );
                }
            }
        }
    }

    #[test]
    fn input_and_output_process_independently() {
        use crate::chemistry::reactor::reactor::{Input, Output};

        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");

        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        let input_idx = reactor.add_input(Input::new(
            zone,
            TransitionMode::Substances {
                entries: vec![SubstanceEntry {
                    id: water_id.clone(),
                    rate_mol_per_second: 1.0,
                }],
            },
        ));

        let output_idx = reactor.add_output(Output::new(
            zone,
            TransitionMode::Substances {
                entries: vec![SubstanceEntry {
                    id: water_id.clone(),
                    rate_mol_per_second: 0.5,
                }],
            },
        ));

        // Apply input for 2 seconds в†’ 2 mol water
        reactor.apply_input(&registry, input_idx, 2.0).unwrap();
        assert!(
            (reactor.zone(&zone).unwrap().concentration_of(&water_id) - 2.0).abs() < 1.0e-9,
            "after input: {}",
            reactor.zone(&zone).unwrap().concentration_of(&water_id)
        );

        // Apply output for 1 second в†’ 0.5 mol removed
        let extracted = reactor.apply_output(&registry, output_idx, 1.0).unwrap();
        assert_eq!(extracted.len(), 1);
        assert!(
            (extracted[0].1 - 0.5).abs() < 1.0e-9,
            "extracted: {}",
            extracted[0].1
        );
        assert!(
            (reactor.zone(&zone).unwrap().concentration_of(&water_id) - 1.5).abs() < 1.0e-9,
            "after output: {}",
            reactor.zone(&zone).unwrap().concentration_of(&water_id)
        );
    }

    #[test]
    fn disabled_input_does_nothing() {
        use crate::chemistry::reactor::reactor::Input;

        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");

        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        let input_idx = reactor.add_input(
            Input::new(
                zone,
                TransitionMode::Substances {
                    entries: vec![SubstanceEntry {
                        id: water_id.clone(),
                        rate_mol_per_second: 1.0,
                    }],
                },
            )
            .with_enabled(false),
        );

        reactor.apply_input(&registry, input_idx, 5.0).unwrap();
        assert!(
            reactor.zone(&zone).unwrap().concentration_of(&water_id) < 1.0e-9,
            "disabled input should not add substances"
        );
    }

    #[test]
    fn per_substance_rates_in_transition() {
        use crate::chemistry::reactor::reactor::SubstanceEntry;

        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let ethanol_id = SubstanceId::from("destroy:ethanol");

        let mut reactor = Reactor::new();
        let from_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());
        let to_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let mixture = reactor.zone_mut(&from_zone).unwrap().mixture_mut();
            let _ = mixture.add_substance(&registry, water_id.clone(), 10.0);
            let _ = mixture.add_substance(&registry, ethanol_id.clone(), 10.0);
        }

        reactor.add_transition(ZoneTransition::new(
            from_zone,
            to_zone,
            TransitionMode::Substances {
                entries: vec![
                    SubstanceEntry {
                        id: water_id.clone(),
                        rate_mol_per_second: 1.0,
                    },
                    SubstanceEntry {
                        id: ethanol_id.clone(),
                        rate_mol_per_second: 0.5,
                    },
                ],
            },
        ));

        reactor.tick(&registry, 1.0).unwrap();

        let water_to = reactor.zone(&to_zone).unwrap().concentration_of(&water_id);
        let ethanol_to = reactor
            .zone(&to_zone)
            .unwrap()
            .concentration_of(&ethanol_id);

        assert!(
            (water_to - 1.0).abs() < 1.0e-6,
            "water should transfer at 1.0 mol/s, got {water_to}"
        );
        assert!(
            (ethanol_to - 0.5).abs() < 1.0e-6,
            "ethanol should transfer at 0.5 mol/s, got {ethanol_to}"
        );
    }

    #[test]
    fn per_phase_rates_in_transition() {
        use crate::chemistry::reactor::reactor::PhaseEntry;

        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let oxygen_id = SubstanceId::from("destroy:oxygen");

        let mut reactor = Reactor::new();
        let from_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());
        let to_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let mixture = reactor.zone_mut(&from_zone).unwrap().mixture_mut();
            let _ = mixture.add_substance(&registry, water_id.clone(), 5.0);
            let _ = mixture.add_substance(&registry, oxygen_id.clone(), 5.0);
        }

        reactor.add_transition(ZoneTransition::new(
            from_zone,
            to_zone,
            TransitionMode::Phases {
                entries: vec![
                    PhaseEntry {
                        phase: MixturePhase::Aqueous,
                        rate_mol_per_second: 1.0,
                    },
                    PhaseEntry {
                        phase: MixturePhase::Gas,
                        rate_mol_per_second: 0.5,
                    },
                ],
            },
        ));

        reactor.tick(&registry, 1.0).unwrap();

        let water_to = reactor.zone(&to_zone).unwrap().concentration_of(&water_id);
        let oxygen_to = reactor.zone(&to_zone).unwrap().concentration_of(&oxygen_id);

        assert!(
            water_to > 0.0,
            "aqueous phase should transfer, got water={water_to}"
        );
        assert!(
            oxygen_to > 0.0,
            "gas phase should transfer, got oxygen={oxygen_to}"
        );
    }

    #[test]
    fn disabled_transition_is_skipped() {
        use crate::chemistry::reactor::reactor::SubstanceEntry;

        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");

        let mut reactor = Reactor::new();
        let from_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());
        let to_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let mixture = reactor.zone_mut(&from_zone).unwrap().mixture_mut();
            let _ = mixture.add_substance(&registry, water_id.clone(), 10.0);
        }

        reactor.add_transition(
            ZoneTransition::new(
                from_zone,
                to_zone,
                TransitionMode::Substances {
                    entries: vec![SubstanceEntry {
                        id: water_id.clone(),
                        rate_mol_per_second: 1.0,
                    }],
                },
            )
            .with_enabled(false),
        );

        reactor.tick(&registry, 1.0).unwrap();

        let from_water = reactor
            .zone(&from_zone)
            .unwrap()
            .concentration_of(&water_id);
        let to_water = reactor.zone(&to_zone).unwrap().concentration_of(&water_id);

        assert!(
            (from_water - 10.0).abs() < 1.0e-6,
            "from zone should retain all water"
        );
        assert!(
            to_water < 1.0e-9,
            "to zone should receive nothing from disabled transition"
        );
    }

    #[test]
    fn reactor_zone_rejects_invalid_volume() {
        assert!(ReactorZone::new(0.0).is_err());
        assert!(ReactorZone::new(f64::NAN).is_err());
    }

    #[test]
    fn reactor_zone_tracks_condensed_volume_as_headspace() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let mut zone = ReactorZone::new(0.001).unwrap();

        let initial_headspace = zone.headspace_volume_cubic_meters(&registry).unwrap();
        zone.add_substance_checked(&registry, water_id, 10.0)
            .unwrap();

        let condensed = zone.condensed_volume_cubic_meters(&registry).unwrap();
        let headspace = zone.headspace_volume_cubic_meters(&registry).unwrap();
        assert!(condensed > 0.0);
        assert!(headspace < initial_headspace);
        assert!((zone.mixture().gas_volume_cubic_meters() - headspace).abs() < 1.0e-12);
    }

    #[test]
    fn reactor_zone_rejects_volume_shrink_that_would_overfill() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let mut zone = ReactorZone::new(0.001).unwrap();
        zone.add_substance_checked(&registry, water_id, 10.0)
            .unwrap();
        let original_volume = zone.volume_cubic_meters();
        let original_gas_volume = zone.mixture().gas_volume_cubic_meters();
        let condensed = zone.condensed_volume_cubic_meters(&registry).unwrap();

        let error = zone
            .set_volume_cubic_meters(&registry, condensed * 0.5)
            .unwrap_err();

        assert!(matches!(
            error,
            crate::chemistry::ChemistryError::InvalidMixtureState(_)
        ));
        assert_eq!(zone.volume_cubic_meters(), original_volume);
        assert_eq!(
            zone.mixture().gas_volume_cubic_meters(),
            original_gas_volume
        );
    }

    #[test]
    fn reactor_zone_rejects_liquid_overfill_without_mutating_mixture() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let mut zone = ReactorZone::new(0.00001).unwrap();

        let error = zone
            .add_substance_checked(&registry, water_id.clone(), 1.0)
            .unwrap_err();

        assert!(matches!(
            error,
            crate::chemistry::ChemistryError::InvalidMixtureState(_)
        ));
        assert!(zone.concentration_of(&water_id) < 1.0e-12);
    }

    #[test]
    fn liquid_filled_zone_allows_no_headspace_without_gas_phase() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let volume = 0.0005;
        let target_condensed_volume = volume - 0.5e-9;
        let amount =
            liquid_amount_for_condensed_volume(&registry, &water_id, target_condensed_volume);
        let mut zone = ReactorZone::new(volume)
            .unwrap()
            .with_volume_mode(ReactorVolumeMode::LiquidFilled);

        zone.add_substance_checked(&registry, water_id.clone(), amount)
            .unwrap();

        let headspace = zone.headspace_volume_cubic_meters(&registry).unwrap();
        let condensed = zone.condensed_volume_cubic_meters(&registry).unwrap();
        assert!(
            headspace <= 1.0e-9,
            "headspace={headspace}, condensed={condensed}, volume={volume}"
        );
        assert!(zone.pressure_pascal() <= 1.0e-9);
        assert!((zone.concentration_of(&water_id) - amount).abs() < 1.0e-9);
    }

    #[test]
    fn headspace_required_zone_rejects_exact_liquid_fill() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let volume = 0.0005;
        let target_condensed_volume = volume - 0.5e-9;
        let amount =
            liquid_amount_for_condensed_volume(&registry, &water_id, target_condensed_volume);
        let mut zone = ReactorZone::new(volume).unwrap();

        let error = zone
            .add_substance_checked(&registry, water_id.clone(), amount)
            .unwrap_err();

        assert!(matches!(
            error,
            crate::chemistry::ChemistryError::InvalidMixtureState(_)
        ));
        assert!(zone.concentration_of(&water_id) < 1.0e-12);
    }

    #[test]
    fn liquid_filled_zone_rejects_gas_without_headspace() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let oxygen_id = SubstanceId::from("destroy:oxygen");
        let volume = 0.0005;
        let target_condensed_volume = volume - 0.5e-9;
        let amount =
            liquid_amount_for_condensed_volume(&registry, &water_id, target_condensed_volume);
        let mut zone = ReactorZone::new(volume)
            .unwrap()
            .with_volume_mode(ReactorVolumeMode::LiquidFilled);
        zone.add_substance_checked(&registry, water_id, amount)
            .unwrap();

        let error = zone
            .add_substance_checked(&registry, oxygen_id.clone(), 0.1)
            .unwrap_err();

        assert!(matches!(
            error,
            crate::chemistry::ChemistryError::InvalidMixtureState(_)
        ));
        assert!(zone.concentration_of(&oxygen_id) < 1.0e-12);
    }

    #[test]
    fn transition_transfers_thermal_energy_with_substance() {
        use crate::chemistry::reactor::reactor::SubstanceEntry;

        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");

        let mut reactor = Reactor::new();
        let from_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());
        let to_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let mixture = reactor.zone_mut(&from_zone).unwrap().mixture_mut();
            mixture
                .add_substance(&registry, water_id.clone(), 5.0)
                .unwrap();
            mixture.heat(&registry, 50_000.0).unwrap();
        }
        {
            let mixture = reactor.zone_mut(&to_zone).unwrap().mixture_mut();
            mixture
                .add_substance(&registry, water_id.clone(), 1.0)
                .unwrap();
        }

        let source_initial_temperature = reactor.zone(&from_zone).unwrap().temperature_kelvin();
        let target_initial_temperature = reactor.zone(&to_zone).unwrap().temperature_kelvin();
        assert!(source_initial_temperature > target_initial_temperature);

        reactor.add_transition(ZoneTransition::new(
            from_zone,
            to_zone,
            TransitionMode::Substances {
                entries: vec![SubstanceEntry {
                    id: water_id,
                    rate_mol_per_second: 1.0,
                }],
            },
        ));

        reactor.tick(&registry, 1.0).unwrap();

        let target_final_temperature = reactor.zone(&to_zone).unwrap().temperature_kelvin();
        assert!(
            target_final_temperature > target_initial_temperature,
            "hot transferred stream should heat receiver: {target_initial_temperature} -> {target_final_temperature}"
        );
    }

    #[test]
    fn electrical_interface_on_smart_heater() {
        let water_id: SubstanceId = SubstanceId::from("destroy:water");
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut heater =
            SmartHeaterState::new(373.0, 10.0, 293.15, 0.0004, 1673.0, 120.0, ControlMode::PID);

        assert!((heater.voltage() - 120.0).abs() < 1.0e-9);
        heater.set_voltage(240.0);
        assert!((heater.voltage() - 240.0).abs() < 1.0e-9);

        let r = heater.resistance_at(300.0);
        assert!(r > 0.0);
        let expected_power = 240.0 * 240.0 / r;
        assert!(
            (heater.power_at(300.0) - expected_power).abs() < 1.0e-6,
            "power_at should use current voltage"
        );

        let mut zone = ReactorZone::new(0.001).unwrap();
        let _ = zone.mixture_mut().add_substance(&registry, water_id, 1.0);
        let mut peripheral = Peripheral::SmartHeater(heater);
        peripheral.apply(&mut zone, &registry, 1.0).unwrap();

        let state = match &peripheral {
            Peripheral::SmartHeater(s) => s,
            _ => unreachable!(),
        };
        assert!(
            state.electrical_draw_w() >= 0.0,
            "electrical_draw_w should be non-negative"
        );
        assert!(
            state.heating_power_w() >= 0.0,
            "heating_power_w should be non-negative"
        );
        assert!(state.last_resistance_ohm() > 0.0);
    }

    #[test]
    fn electrical_interface_on_electrode() {
        let water_id: SubstanceId = SubstanceId::from("destroy:water");
        let mut electrode = ElectrodeState::new(12.0, 50.0, 1.23);

        assert!((electrode.voltage() - 12.0).abs() < 1.0e-9);
        electrode.set_voltage(24.0);
        assert!((electrode.voltage() - 24.0).abs() < 1.0e-9);

        assert!((electrode.current_a() - 0.0).abs() < 1.0e-9);
        electrode.set_current(10.0);
        assert!((electrode.current_a() - 10.0).abs() < 1.0e-9);
        assert!((electrode.max_current_a() - 50.0).abs() < 1.0e-9);

        electrode.set_current(100.0);
        assert!(
            (electrode.current_a() - 50.0).abs() < 1.0e-9,
            "should clamp to max_current_a"
        );

        let draw = electrode.electrical_draw_w();
        assert!(draw >= 0.0);
        let energy = electrode.energy_delivered_j();
        assert!(energy >= 0.0);

        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut zone = ReactorZone::new(0.001).unwrap();
        let _ = zone.mixture_mut().add_substance(&registry, water_id, 1.0);
        let mut peripheral = Peripheral::Electrode(electrode);
        peripheral.apply(&mut zone, &registry, 1.0).unwrap();

        match &peripheral {
            Peripheral::Electrode(e) => {
                assert!(e.electrical_draw_w() >= 0.0);
                assert!(e.energy_delivered_j() >= 0.0);
            }
            _ => unreachable!(),
        }
    }

    #[test]
    fn reactor_total_electrical_draw() {
        let water_id: SubstanceId = SubstanceId::from("destroy:water");
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(
            ReactorZone::new(0.001)
                .unwrap()
                .with_peripheral(Peripheral::SmartHeater(SmartHeaterState::new(
                    373.0,
                    10.0,
                    293.15,
                    0.0004,
                    1673.0,
                    120.0,
                    ControlMode::PID,
                )))
                .with_peripheral(Peripheral::Electrode(
                    ElectrodeState::new(12.0, 50.0, 1.23).with_current(5.0),
                )),
        );

        {
            let z = reactor.zone_mut(&zone).unwrap();
            let _ = z.mixture_mut().add_substance(&registry, water_id, 1.0);
        }

        reactor.tick(&registry, 1.0).unwrap();

        let total = reactor.total_electrical_draw_w();
        assert!(
            total >= 0.0,
            "total electrical draw should be non-negative, got {total}"
        );
    }

    #[test]
    fn ambient_heat_exchange_cools_hot_zone() {
        let water_id: SubstanceId = SubstanceId::from("destroy:water");
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let z = reactor.zone_mut(&zone).unwrap();
            let _ = z.mixture_mut().add_substance(&registry, water_id, 1.0);
        }

        reactor.set_ambient_temperature(293.0);
        reactor.set_heat_transfer_coefficient(0.1);

        let initial_temp = reactor.zone(&zone).unwrap().temperature_kelvin();
        assert!(
            initial_temp > 293.0,
            "initial temp should be above ambient, got {initial_temp}"
        );

        for _ in 0..10 {
            reactor.tick(&registry, 1.0).unwrap();
        }

        let final_temp = reactor.zone(&zone).unwrap().temperature_kelvin();
        assert!(
            final_temp < initial_temp,
            "zone should cool toward ambient: {initial_temp} -> {final_temp}"
        );

        let energy_exchanged = reactor.last_ambient_energy_exchange_j();
        assert!(
            energy_exchanged < 0.0,
            "energy should be negative (lost to ambient), got {energy_exchanged}"
        );
    }

    #[test]
    fn ambient_heat_exchange_heats_cold_zone() {
        let water_id: SubstanceId = SubstanceId::from("destroy:water");
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let z = reactor.zone_mut(&zone).unwrap();
            let _ = z.mixture_mut().add_substance(&registry, water_id, 1.0);
        }

        reactor.set_ambient_temperature(400.0);
        reactor.set_heat_transfer_coefficient(0.1);

        let initial_temp = reactor.zone(&zone).unwrap().temperature_kelvin();
        assert!(
            initial_temp < 400.0,
            "initial temp should be below ambient, got {initial_temp}"
        );

        for _ in 0..10 {
            reactor.tick(&registry, 1.0).unwrap();
        }

        let final_temp = reactor.zone(&zone).unwrap().temperature_kelvin();
        assert!(
            final_temp > initial_temp,
            "zone should heat toward ambient: {initial_temp} -> {final_temp}"
        );

        let energy_exchanged = reactor.last_ambient_energy_exchange_j();
        assert!(
            energy_exchanged > 0.0,
            "energy should be positive (gained from ambient), got {energy_exchanged}"
        );
    }

    #[test]
    fn ambient_heat_exchange_independent_zones() {
        let water_id: SubstanceId = SubstanceId::from("destroy:water");
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut reactor = Reactor::new();

        let hot_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());
        let cold_zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let z = reactor.zone_mut(&hot_zone).unwrap();
            let _ = z
                .mixture_mut()
                .add_substance(&registry, water_id.clone(), 5.0);
        }
        {
            let z = reactor.zone_mut(&cold_zone).unwrap();
            let _ = z.mixture_mut().add_substance(&registry, water_id, 0.5);
        }

        reactor
            .zone_mut(&hot_zone)
            .unwrap()
            .mixture_mut()
            .heat(&registry, 50_000.0)
            .unwrap();
        reactor
            .zone_mut(&cold_zone)
            .unwrap()
            .mixture_mut()
            .heat(&registry, 50_000.0)
            .unwrap();

        let hot_initial = reactor.zone(&hot_zone).unwrap().temperature_kelvin();
        let cold_initial = reactor.zone(&cold_zone).unwrap().temperature_kelvin();

        reactor.set_ambient_temperature(293.0);
        reactor.set_heat_transfer_coefficient(0.01);

        for _ in 0..100 {
            reactor.tick(&registry, 1.0).unwrap();
        }

        let hot_final = reactor.zone(&hot_zone).unwrap().temperature_kelvin();
        let cold_final = reactor.zone(&cold_zone).unwrap().temperature_kelvin();

        assert!(
            hot_final < hot_initial,
            "hot zone should cool: {hot_initial} -> {hot_final}"
        );
        assert!(
            cold_final < cold_initial,
            "cold zone should also cool: {cold_initial} -> {cold_final}"
        );
        assert!(
            (hot_initial - cold_initial).abs() > (hot_final - cold_final).abs(),
            "temperature difference should decrease: О”T_initial={} vs О”T_final={}",
            hot_initial - cold_initial,
            hot_final - cold_final
        );
    }

    #[test]
    fn no_ambient_exchange_when_not_configured() {
        let water_id: SubstanceId = SubstanceId::from("destroy:water");
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());

        {
            let z = reactor.zone_mut(&zone).unwrap();
            let _ = z.mixture_mut().add_substance(&registry, water_id, 1.0);
        }

        let initial_temp = reactor.zone(&zone).unwrap().temperature_kelvin();

        for _ in 0..10 {
            reactor.tick(&registry, 1.0).unwrap();
        }

        let final_temp = reactor.zone(&zone).unwrap().temperature_kelvin();
        assert!(
            (final_temp - initial_temp).abs() < 1.0,
            "without ambient config, temperature should not change: {initial_temp} -> {final_temp}"
        );
        assert!(
            (reactor.last_ambient_energy_exchange_j()).abs() < 1.0e-12,
            "energy exchange should be zero"
        );
    }
}
