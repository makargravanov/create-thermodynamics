use crate::chemistry::metallurgy::{
    ComponentLimit, MetallurgicalPhaseKind, MetallurgicalPhaseModel,
    MetallurgicalPhasePropertyModel, MetallurgicalSystem, PhaseWeightRule,
};

pub fn default_metallurgical_systems() -> Vec<MetallurgicalSystem> {
    vec![
        iron_carbon_system(),
        copper_zinc_system(),
        copper_tin_system(),
        aluminum_silicon_system(),
    ]
}

fn iron_carbon_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:fe_c", ["Fe", "destroy:carbon"])
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(60.0, 120.0, 0.55, 0.95, 35.0, 0.45),
            )
            .weight_rule(PhaseWeightRule::constant(6.0).temperature_window(1700.0, 3400.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/ferrite",
                MetallurgicalPhaseKind::Ferrite,
                properties(95.0, 220.0, 0.35, 0.10, 80.0, 0.45),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.0, 0.025))
            .weight_rule(
                PhaseWeightRule::composition_centered("destroy:carbon", 0.004, 0.03, 4.0)
                    .temperature_window(0.0, 1185.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/austenite",
                MetallurgicalPhaseKind::Austenite,
                properties(150.0, 310.0, 0.28, 0.12, 18.0, 0.40),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.0, 0.09))
            .weight_rule(
                PhaseWeightRule::composition_centered("destroy:carbon", 0.035, 0.08, 5.0)
                    .temperature_window(900.0, 1800.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/cementite",
                MetallurgicalPhaseKind::Cementite,
                properties(800.0, 1200.0, 0.02, 0.55, 11.0, 0.20),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.01, 0.25))
            .weight_rule(
                PhaseWeightRule::composition_centered("destroy:carbon", 0.25, 0.25, 2.0)
                    .temperature_window(0.0, 1500.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/pearlite",
                MetallurgicalPhaseKind::Pearlite,
                properties(250.0, 520.0, 0.18, 0.18, 45.0, 0.35),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.02, 0.08))
            .weight_rule(
                PhaseWeightRule::composition_centered("destroy:carbon", 0.035, 0.04, 3.5)
                    .temperature_window(0.0, 1000.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/martensite",
                MetallurgicalPhaseKind::Martensite,
                properties(650.0, 1500.0, 0.04, 0.30, 28.0, 0.25),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.01, 0.12))
            .weight_rule(
                PhaseWeightRule::composition_centered("destroy:carbon", 0.04, 0.08, 0.0)
                    .temperature_window(0.0, 650.0)
                    .cooling_rate_bonus(80.0, 8.0),
            ),
        )
}

fn copper_zinc_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:cu_zn", ["Cu", "Zn"])
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_zn/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(50.0, 100.0, 0.50, 0.20, 55.0, 0.55),
            )
            .weight_rule(PhaseWeightRule::constant(5.0).temperature_window(1150.0, 3000.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_zn/alpha_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(100.0, 260.0, 0.45, 0.08, 85.0, 0.55),
            )
            .limit(ComponentLimit::new("Zn", 0.0, 0.38))
            .weight_rule(PhaseWeightRule::composition_centered("Zn", 0.18, 0.38, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_zn/beta_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(170.0, 480.0, 0.25, 0.12, 45.0, 0.45),
            )
            .limit(ComponentLimit::new("Zn", 0.30, 0.55))
            .weight_rule(PhaseWeightRule::composition_centered("Zn", 0.45, 0.25, 3.0)),
        )
}

fn copper_tin_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:cu_sn", ["Cu", "Sn"])
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_sn/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(45.0, 90.0, 0.45, 0.18, 50.0, 0.50),
            )
            .weight_rule(PhaseWeightRule::constant(5.0).temperature_window(1050.0, 2900.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_sn/alpha_bronze",
                MetallurgicalPhaseKind::SolidSolution,
                properties(130.0, 340.0, 0.35, 0.10, 70.0, 0.55),
            )
            .limit(ComponentLimit::new("Sn", 0.0, 0.16))
            .weight_rule(PhaseWeightRule::composition_centered("Sn", 0.08, 0.16, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_sn/intermetallic",
                MetallurgicalPhaseKind::Intermetallic,
                properties(360.0, 800.0, 0.06, 0.22, 25.0, 0.35),
            )
            .limit(ComponentLimit::new("Sn", 0.12, 0.45))
            .weight_rule(PhaseWeightRule::composition_centered("Sn", 0.25, 0.22, 2.5)),
        )
}

fn aluminum_silicon_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_si", ["Al", "destroy:silicon"])
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_si/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(35.0, 70.0, 0.50, 0.05, 95.0, 0.60),
            )
            .weight_rule(PhaseWeightRule::constant(5.0).temperature_window(850.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_si/aluminum_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(55.0, 170.0, 0.35, 0.04, 120.0, 0.60),
            )
            .limit(ComponentLimit::new("destroy:silicon", 0.0, 0.15))
            .weight_rule(PhaseWeightRule::composition_centered(
                "destroy:silicon",
                0.07,
                0.15,
                4.0,
            )),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_si/silicon_rich",
                MetallurgicalPhaseKind::Intermetallic,
                properties(450.0, 700.0, 0.03, 0.25, 30.0, 0.35),
            )
            .limit(ComponentLimit::new("destroy:silicon", 0.12, 0.60))
            .weight_rule(PhaseWeightRule::composition_centered(
                "destroy:silicon",
                0.25,
                0.28,
                2.0,
            )),
        )
}

fn properties(
    hardness_hv: f64,
    yield_strength_mpa: f64,
    ductility_fraction: f64,
    electrical_resistivity_micro_ohm_meter: f64,
    thermal_conductivity_w_per_meter_kelvin: f64,
    corrosion_resistance_score: f64,
) -> MetallurgicalPhasePropertyModel {
    MetallurgicalPhasePropertyModel {
        hardness_hv,
        yield_strength_mpa,
        ductility_fraction,
        electrical_resistivity_micro_ohm_meter,
        thermal_conductivity_w_per_meter_kelvin,
        corrosion_resistance_score,
    }
}
