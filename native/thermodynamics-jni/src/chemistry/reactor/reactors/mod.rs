use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::reactor::peripheral::{ControlMode, ElectrodeState, Peripheral, SmartHeaterState};
use crate::chemistry::reactor::reactor::{Reactor, TransitionMode, ZoneId, ZoneTransition};
use crate::chemistry::reactor::zone::ReactorZone;

pub fn batch_reactor(volume_m3: f64, target_kelvin: f64) -> Reactor {
    let mut reactor = Reactor::new();
    let zone = ReactorZone::new(volume_m3)
        .with_peripheral(Peripheral::SmartHeater(SmartHeaterState::new(
            target_kelvin,
            10.0,
            293.15,
            0.0004,
            1673.0,
            230.0,
            ControlMode::PID,
        )));
    reactor.add_zone(zone);
    reactor
}

pub fn cstr(volume_m3: f64, target_kelvin: f64) -> Reactor {
    let mut reactor = Reactor::new();

    let input = reactor.add_zone(ReactorZone::new(volume_m3));
    let reaction = reactor.add_zone(
        ReactorZone::new(volume_m3)
            .with_peripheral(Peripheral::SmartHeater(SmartHeaterState::new(
                target_kelvin,
                10.0,
                293.15,
                0.0004,
                1673.0,
                230.0,
                ControlMode::PID,
            ))),
    );
    let output = reactor.add_zone(ReactorZone::new(volume_m3));

    reactor.add_transition(ZoneTransition {
        from: input,
        to: reaction,
        mode: TransitionMode::All {
            rate_mol_per_second: 0.5,
        },
    });
    reactor.add_transition(ZoneTransition {
        from: reaction,
        to: output,
        mode: TransitionMode::All {
            rate_mol_per_second: 0.5,
        },
    });

    reactor
}

pub fn electrolysis_cell(volume_m3: f64) -> Reactor {
    let mut reactor = Reactor::new();

    let electrolyte = reactor.add_zone(
        ReactorZone::new(volume_m3)
            .with_peripheral(Peripheral::SmartHeater(SmartHeaterState::new(
                353.0,
                5.0,
                293.15,
                0.0004,
                1673.0,
                12.0,
                ControlMode::P,
            )))
            .with_peripheral(Peripheral::Electrode(ElectrodeState::new(
                12.0, 50.0, 1.23,
            ))),
    );
    let product = reactor.add_zone(ReactorZone::new(volume_m3));

    reactor.add_transition(ZoneTransition {
        from: electrolyte,
        to: product,
        mode: TransitionMode::Substances {
            ids: vec!["destroy:hydrogen".into(), "destroy:oxygen".into()],
            rate_mol_per_second: 0.1,
        },
    });

    reactor
}

pub fn distillation_column(stages: usize, volume_per_stage_m3: f64) -> Reactor {
    let mut reactor = Reactor::new();

    let mut zone_ids: Vec<ZoneId> = Vec::new();
    for _ in 0..stages {
        let zone = ReactorZone::new(volume_per_stage_m3)
            .with_peripheral(Peripheral::HeatExchanger {
                coolant_temperature_kelvin: 373.0,
                u_kw_per_kelvin: 0.5,
            });
        zone_ids.push(reactor.add_zone(zone));
    }

    for i in 0..stages - 1 {
        reactor.add_transition(ZoneTransition {
            from: zone_ids[i],
            to: zone_ids[i + 1],
            mode: TransitionMode::Phases {
                phases: vec![MixturePhase::Gas],
                rate_mol_per_second: 0.2,
            },
        });
        reactor.add_transition(ZoneTransition {
            from: zone_ids[i + 1],
            to: zone_ids[i],
            mode: TransitionMode::Phases {
                phases: vec![MixturePhase::Aqueous],
                rate_mol_per_second: 0.1,
            },
        });
    }

    reactor
}

pub fn arc_furnace(volume_m3: f64) -> Reactor {
    let mut reactor = Reactor::new();

    let chamber = reactor.add_zone(
        ReactorZone::new(volume_m3)
            .with_peripheral(Peripheral::Electrode(ElectrodeState::new(
                80.0, 500.0, 0.0,
            )))
            .with_peripheral(Peripheral::HeatExchanger {
                coolant_temperature_kelvin: 293.0,
                u_kw_per_kelvin: 2.0,
            }),
    );

    let output = reactor.add_zone(ReactorZone::new(volume_m3));

    reactor.add_transition(ZoneTransition {
        from: chamber,
        to: output,
        mode: TransitionMode::Phases {
            phases: vec![MixturePhase::MoltenMetal],
            rate_mol_per_second: 1.0,
        },
    });

    reactor
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::destroy_registry_builder;
    use crate::chemistry::reactor::io;
    use crate::chemistry::substance::SubstanceId;
    use crate::chemistry::mixture::MixturePhase;

    #[test]
    fn distillation_column_separates_ethanol_from_water() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let water_id = SubstanceId::from("destroy:water");
        let ethanol_id = SubstanceId::from("destroy:ethanol");

        let temps = [298.0, 353.0, 373.0, 400.0];

        for &temp in &temps {
            let mut reactor = distillation_column(5, 0.001);

            // Fill zone 0 with 60% water, 40% ethanol
            {
                let zone = reactor.zone_mut(&ZoneId(0)).unwrap();
                let mixture = zone.mixture_mut();
                let _ = mixture.add_substance(&registry, water_id.clone(), 30.0);
                let _ = mixture.add_substance(&registry, ethanol_id.clone(), 20.0);
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
            println!("{:<10} {:>12} {:>12} {:>12} {:>10}", "Zone", "Water", "Ethanol", "T (K)", "EtOH frac");
            println!("{}", "-".repeat(62));
            for (name, w, e, t) in &results {
                let total = w + e;
                let frac = if total > 0.001 { e / total } else { 0.0 };
                println!("{:<10} {:>12.4} {:>12.4} {:>12.2} {:>10.4}", name, w, e, t, frac);
            }

            // Ethanol fraction should increase from bottom (zone_0) to top (zone_4)
            let frac_bottom = {
                let total = results[0].1 + results[0].2;
                if total > 0.001 { results[0].2 / total } else { 0.0 }
            };
            let frac_top = {
                let total = results[4].1 + results[4].2;
                if total > 0.001 { results[4].2 / total } else { 0.0 }
            };

            println!("\nEthanol fraction: bottom={:.4}, top={:.4}", frac_bottom, frac_top);

            if temp >= 353.0 {
                // Ethanol should have moved upward — top zone either has
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

                // Water should stay mostly in bottom
                let water_in_top_zones: f64 = results[1..].iter().map(|r| r.1).sum();
                let water_in_bottom = results[0].1;
                assert!(
                    water_in_bottom > water_in_top_zones,
                    "At T={}: water should stay in bottom zone",
                    temp
                );
            }
        }
    }
}
