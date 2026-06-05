use crate::chemistry::metallurgy::{
    ComponentLimit, CompositionEnergyTerm, CrystalStructure, LiquidMiscibility,
    MetallurgicalCompoundPhaseData, MetallurgicalElementData, MetallurgicalPairInteractionData,
    MetallurgicalPhaseKind, MetallurgicalPhaseModel, MetallurgicalPhasePropertyModel,
    MetallurgicalPropertyCalibration, MetallurgicalSystem, PhaseBoundaryPoint,
    PhaseFreeEnergyModel, PhaseKineticModel, SolidMiscibility, ThermalTreatmentProfile,
};

pub fn default_metallurgical_element_data() -> Vec<MetallurgicalElementData> {
    vec![
        element(
            "Fe",
            1811.0,
            126.0,
            CrystalStructure::BodyCenteredCubic,
            properties(90.0, 210.0, 0.35, 0.10, 80.0, 0.45),
        )
        .solid_solution_strengthening(850.0),
        element(
            "destroy:carbon",
            3915.0,
            70.0,
            CrystalStructure::Complex,
            properties(900.0, 1200.0, 0.01, 1.50, 8.0, 0.15),
        )
        .solid_solution_strengthening(2400.0)
        .intermetallic_forming_tendency(0.70),
        element(
            "Cr",
            2180.0,
            128.0,
            CrystalStructure::BodyCenteredCubic,
            properties(150.0, 360.0, 0.26, 0.13, 94.0, 0.90),
        )
        .solid_solution_strengthening(950.0)
        .carbide_forming_tendency(0.80),
        element(
            "Ni",
            1728.0,
            124.0,
            CrystalStructure::FaceCenteredCubic,
            properties(90.0, 170.0, 0.45, 0.07, 91.0, 0.72),
        )
        .solid_solution_strengthening(700.0),
        element(
            "Mn",
            1519.0,
            127.0,
            CrystalStructure::Complex,
            properties(120.0, 260.0, 0.25, 1.45, 7.8, 0.30),
        )
        .solid_solution_strengthening(900.0),
        element(
            "destroy:silicon",
            1687.0,
            111.0,
            CrystalStructure::DiamondCubic,
            properties(1150.0, 700.0, 0.01, 2300.0, 150.0, 0.58),
        )
        .solid_solution_strengthening(1200.0)
        .intermetallic_forming_tendency(0.35),
        element(
            "Mo",
            2896.0,
            139.0,
            CrystalStructure::BodyCenteredCubic,
            properties(250.0, 550.0, 0.18, 0.053, 139.0, 0.78),
        )
        .solid_solution_strengthening(1200.0)
        .carbide_forming_tendency(0.90),
        element(
            "V",
            2183.0,
            134.0,
            CrystalStructure::BodyCenteredCubic,
            properties(260.0, 650.0, 0.16, 0.20, 31.0, 0.70),
        )
        .solid_solution_strengthening(1300.0)
        .carbide_forming_tendency(0.95),
        element(
            "Al",
            933.0,
            143.0,
            CrystalStructure::FaceCenteredCubic,
            properties(25.0, 35.0, 0.42, 0.028, 237.0, 0.76),
        )
        .solid_solution_strengthening(480.0),
        element(
            "Cu",
            1358.0,
            128.0,
            CrystalStructure::FaceCenteredCubic,
            properties(45.0, 70.0, 0.50, 0.017, 401.0, 0.62),
        )
        .solid_solution_strengthening(520.0),
        element(
            "Mg",
            923.0,
            160.0,
            CrystalStructure::HexagonalClosePacked,
            properties(35.0, 80.0, 0.28, 0.044, 156.0, 0.32),
        )
        .solid_solution_strengthening(430.0),
        element(
            "Zn",
            693.0,
            134.0,
            CrystalStructure::HexagonalClosePacked,
            properties(35.0, 90.0, 0.22, 0.059, 116.0, 0.50),
        )
        .solid_solution_strengthening(420.0),
        element(
            "Sn",
            505.0,
            145.0,
            CrystalStructure::Tetragonal,
            properties(12.0, 18.0, 0.35, 0.115, 67.0, 0.54),
        )
        .solid_solution_strengthening(260.0)
        .phase_separation_tendency(0.35),
        element(
            "Pb",
            601.0,
            175.0,
            CrystalStructure::FaceCenteredCubic,
            properties(5.0, 12.0, 0.40, 0.208, 35.0, 0.42),
        )
        .solid_solution_strengthening(180.0)
        .phase_separation_tendency(0.70),
        element(
            "Ag",
            1235.0,
            144.0,
            CrystalStructure::FaceCenteredCubic,
            properties(25.0, 55.0, 0.55, 0.016, 429.0, 0.80),
        )
        .solid_solution_strengthening(330.0),
        element(
            "Au",
            1337.0,
            144.0,
            CrystalStructure::FaceCenteredCubic,
            properties(25.0, 65.0, 0.55, 0.024, 318.0, 0.95),
        )
        .solid_solution_strengthening(350.0),
        element(
            "Bi",
            545.0,
            160.0,
            CrystalStructure::Rhombohedral,
            properties(9.0, 20.0, 0.03, 1.29, 8.0, 0.60),
        )
        .phase_separation_tendency(0.75),
        element(
            "Ti",
            1941.0,
            147.0,
            CrystalStructure::HexagonalClosePacked,
            properties(160.0, 330.0, 0.24, 0.42, 22.0, 0.82),
        )
        .solid_solution_strengthening(900.0),
        element(
            "Co",
            1768.0,
            125.0,
            CrystalStructure::HexagonalClosePacked,
            properties(125.0, 260.0, 0.28, 0.062, 100.0, 0.70),
        )
        .solid_solution_strengthening(780.0),
        element(
            "Be",
            1560.0,
            112.0,
            CrystalStructure::HexagonalClosePacked,
            properties(100.0, 240.0, 0.10, 0.040, 200.0, 0.68),
        )
        .solid_solution_strengthening(1700.0)
        .intermetallic_forming_tendency(0.55),
    ]
}

pub fn default_metallurgical_systems() -> Vec<MetallurgicalSystem> {
    vec![
        iron_carbon_system(),
        iron_carbon_chromium_nickel_system(),
        iron_carbon_manganese_silicon_system(),
        copper_zinc_system(),
        copper_tin_system(),
        tin_lead_system(),
        tin_silver_copper_system(),
        bismuth_tin_system(),
        aluminum_silicon_system(),
        aluminum_copper_system(),
        aluminum_copper_magnesium_system(),
        aluminum_magnesium_system(),
        aluminum_magnesium_silicon_system(),
        aluminum_zinc_magnesium_system(),
        aluminum_zinc_magnesium_copper_system(),
        nickel_chromium_system(),
        nickel_chromium_aluminum_system(),
        nickel_chromium_cobalt_molybdenum_system(),
        copper_nickel_system(),
        copper_beryllium_system(),
        copper_aluminum_system(),
        magnesium_aluminum_zinc_system(),
        titanium_aluminum_vanadium_system(),
        titanium_aluminum_molybdenum_system(),
        iron_carbon_chromium_molybdenum_vanadium_system(),
    ]
}

pub fn default_metallurgical_pair_interactions() -> Vec<MetallurgicalPairInteractionData> {
    use LiquidMiscibility::Complete as LiquidComplete;
    use SolidMiscibility::{Complete as SolidComplete, High, Immiscible, Limited, VeryLimited};
    vec![
        pair("Au", "Ag", SolidComplete, LiquidComplete)
            .strengthening(160.0)
            .resistivity_penalty(0.018),
        pair("Cu", "Ni", SolidComplete, LiquidComplete)
            .strengthening(260.0)
            .resistivity_penalty(0.060),
        pair("Cu", "Zn", High, LiquidComplete)
            .strengthening(360.0)
            .resistivity_penalty(0.075),
        pair("Cu", "Sn", Limited, LiquidComplete)
            .strengthening(620.0)
            .ductility_penalty(0.20)
            .interaction_strength(-18_000.0),
        pair("Al", "Cu", Limited, LiquidComplete)
            .strengthening(760.0)
            .ductility_penalty(0.22)
            .interaction_strength(-22_000.0),
        pair("Al", "Mg", Limited, LiquidComplete)
            .strengthening(520.0)
            .ductility_penalty(0.16)
            .interaction_strength(-16_000.0),
        pair("Al", "Zn", High, LiquidComplete)
            .eutectic(655.0, 0.95)
            .strengthening(480.0)
            .ductility_penalty(0.12),
        pair("Al", "destroy:silicon", VeryLimited, LiquidComplete)
            .eutectic(850.0, 0.12)
            .strengthening(420.0)
            .ductility_penalty(0.24),
        pair("Sn", "Pb", Limited, LiquidComplete)
            .eutectic(456.0, 0.38)
            .ductility_penalty(0.12),
        pair("Bi", "Sn", Limited, LiquidComplete)
            .eutectic(412.0, 0.43)
            .ductility_penalty(0.28),
        pair("Sn", "Ag", Limited, LiquidComplete)
            .strengthening(480.0)
            .interaction_strength(-20_000.0),
        pair("Ni", "Al", Limited, LiquidComplete)
            .strengthening(900.0)
            .interaction_strength(-32_000.0)
            .ductility_penalty(0.20),
        pair("Cu", "Be", Limited, LiquidComplete)
            .strengthening(1100.0)
            .interaction_strength(-28_000.0)
            .ductility_penalty(0.28),
        pair("Fe", "destroy:carbon", VeryLimited, LiquidComplete)
            .strengthening(1800.0)
            .interaction_strength(-24_000.0)
            .ductility_penalty(0.35),
        pair("Cr", "destroy:carbon", VeryLimited, LiquidComplete)
            .strengthening(1600.0)
            .interaction_strength(-30_000.0)
            .ductility_penalty(0.38),
        pair("Mo", "destroy:carbon", VeryLimited, LiquidComplete)
            .strengthening(1700.0)
            .interaction_strength(-34_000.0)
            .ductility_penalty(0.36),
        pair("V", "destroy:carbon", VeryLimited, LiquidComplete)
            .strengthening(1750.0)
            .interaction_strength(-36_000.0)
            .ductility_penalty(0.36),
        pair("Fe", "Cr", High, LiquidComplete)
            .strengthening(520.0)
            .resistivity_penalty(0.11),
        pair("Fe", "Ni", High, LiquidComplete)
            .strengthening(430.0)
            .resistivity_penalty(0.08),
        pair("Fe", "Mn", High, LiquidComplete)
            .strengthening(520.0)
            .resistivity_penalty(0.10),
        pair("Fe", "destroy:silicon", Limited, LiquidComplete)
            .strengthening(620.0)
            .resistivity_penalty(0.16),
        pair("Ag", "Cu", Limited, LiquidComplete)
            .eutectic(1052.0, 0.40)
            .strengthening(300.0)
            .ductility_penalty(0.10),
        pair("Pb", "Cu", Immiscible, LiquidComplete)
            .ductility_penalty(0.40)
            .resistivity_penalty(0.08),
        pair("Pb", "Ag", VeryLimited, LiquidComplete)
            .eutectic(578.0, 0.97)
            .ductility_penalty(0.18),
        pair("Bi", "Pb", Limited, LiquidComplete)
            .eutectic(398.0, 0.56)
            .ductility_penalty(0.24),
        pair("Mg", "destroy:silicon", VeryLimited, LiquidComplete)
            .strengthening(620.0)
            .interaction_strength(-18_000.0)
            .ductility_penalty(0.20),
        pair("Mg", "Zn", Limited, LiquidComplete)
            .strengthening(460.0)
            .interaction_strength(-14_000.0)
            .ductility_penalty(0.14),
        pair("Cu", "Mg", VeryLimited, LiquidComplete)
            .strengthening(540.0)
            .interaction_strength(-12_000.0)
            .ductility_penalty(0.20),
        pair("Fe", "Mo", High, LiquidComplete)
            .strengthening(820.0)
            .interaction_strength(-14_000.0)
            .resistivity_penalty(0.08),
        pair("Fe", "V", High, LiquidComplete)
            .strengthening(900.0)
            .interaction_strength(-16_000.0)
            .resistivity_penalty(0.09),
        pair("Cr", "Mo", High, LiquidComplete)
            .strengthening(780.0)
            .interaction_strength(-18_000.0)
            .resistivity_penalty(0.07),
        pair("Cr", "V", High, LiquidComplete)
            .strengthening(860.0)
            .interaction_strength(-20_000.0)
            .resistivity_penalty(0.08),
        pair("Mo", "V", High, LiquidComplete)
            .strengthening(920.0)
            .interaction_strength(-22_000.0)
            .resistivity_penalty(0.06),
        pair("Ni", "Co", SolidComplete, LiquidComplete)
            .strengthening(420.0)
            .resistivity_penalty(0.05),
        pair("Ni", "Mo", High, LiquidComplete)
            .strengthening(760.0)
            .interaction_strength(-18_000.0)
            .resistivity_penalty(0.07),
        pair("Cr", "Co", High, LiquidComplete)
            .strengthening(620.0)
            .interaction_strength(-12_000.0)
            .resistivity_penalty(0.06),
        pair("Co", "Mo", High, LiquidComplete)
            .strengthening(780.0)
            .interaction_strength(-18_000.0)
            .resistivity_penalty(0.05),
        pair("Ti", "Al", Limited, LiquidComplete)
            .strengthening(900.0)
            .interaction_strength(-24_000.0)
            .ductility_penalty(0.18),
        pair("Ti", "Mo", High, LiquidComplete)
            .strengthening(820.0)
            .interaction_strength(-14_000.0)
            .resistivity_penalty(0.08),
        pair("Al", "Mo", VeryLimited, LiquidComplete)
            .strengthening(700.0)
            .interaction_strength(-12_000.0)
            .ductility_penalty(0.20),
    ]
}

pub fn default_metallurgical_compound_phases() -> Vec<MetallurgicalCompoundPhaseData> {
    vec![
        compound(
            "metallurgy:compound/ni3al",
            [("Ni", 0.75), ("Al", 0.25)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(420.0, 1050.0, 0.10, 0.65, 28.0, 0.82),
            -34_000.0,
        )
        .composition_tolerance(0.18)
        .temperature_window(300.0, 1700.0),
        compound(
            "metallurgy:compound/al2cu",
            [("Al", 0.667), ("Cu", 0.333)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(520.0, 850.0, 0.04, 0.42, 38.0, 0.58),
            -28_000.0,
        )
        .composition_tolerance(0.20),
        compound(
            "metallurgy:compound/mg17al12",
            [("Mg", 0.586), ("Al", 0.414)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(300.0, 520.0, 0.05, 0.22, 52.0, 0.42),
            -19_000.0,
        )
        .composition_tolerance(0.18),
        compound(
            "metallurgy:compound/fe3c",
            [("Fe", 0.75), ("destroy:carbon", 0.25)],
            MetallurgicalPhaseKind::Cementite,
            properties(820.0, 1250.0, 0.015, 0.55, 11.0, 0.18),
            -30_000.0,
        )
        .composition_tolerance(0.22),
        compound(
            "metallurgy:compound/cr_carbide",
            [("Cr", 0.75), ("destroy:carbon", 0.25)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(1050.0, 1700.0, 0.01, 0.72, 8.0, 0.50),
            -38_000.0,
        )
        .composition_tolerance(0.24),
        compound(
            "metallurgy:compound/mo_carbide",
            [("Mo", 0.50), ("destroy:carbon", 0.50)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(1150.0, 1800.0, 0.008, 0.62, 15.0, 0.62),
            -40_000.0,
        )
        .composition_tolerance(0.24),
        compound(
            "metallurgy:compound/v_carbide",
            [("V", 0.50), ("destroy:carbon", 0.50)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(1100.0, 1750.0, 0.008, 0.58, 18.0, 0.58),
            -42_000.0,
        )
        .composition_tolerance(0.24),
        compound(
            "metallurgy:compound/cube",
            [("Cu", 0.50), ("Be", 0.50)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(560.0, 1200.0, 0.03, 0.18, 70.0, 0.70),
            -32_000.0,
        )
        .composition_tolerance(0.20),
        compound(
            "metallurgy:compound/cu6sn5",
            [("Cu", 0.545), ("Sn", 0.455)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(430.0, 760.0, 0.025, 0.34, 36.0, 0.48),
            -24_000.0,
        )
        .composition_tolerance(0.18),
        compound(
            "metallurgy:compound/cu3sn",
            [("Cu", 0.75), ("Sn", 0.25)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(480.0, 820.0, 0.02, 0.30, 42.0, 0.50),
            -25_000.0,
        )
        .composition_tolerance(0.18),
        compound(
            "metallurgy:compound/ag3sn",
            [("Ag", 0.75), ("Sn", 0.25)],
            MetallurgicalPhaseKind::Intermetallic,
            properties(260.0, 430.0, 0.06, 0.11, 80.0, 0.66),
            -21_000.0,
        )
        .composition_tolerance(0.18),
    ]
}

fn iron_carbon_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:fe_c", ["Fe", "destroy:carbon"])
        .phase_boundary(binary_boundary("Fe", "destroy:carbon", 0.0, 1811.0, 1811.0))
        .phase_boundary(binary_boundary(
            "Fe",
            "destroy:carbon",
            0.035,
            1420.0,
            1720.0,
        ))
        .phase_boundary(binary_boundary(
            "Fe",
            "destroy:carbon",
            0.17,
            1420.0,
            1420.0,
        ))
        .phase_boundary(binary_boundary(
            "Fe",
            "destroy:carbon",
            0.25,
            1500.0,
            1620.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(60.0, 120.0, 0.55, 0.95, 35.0, 0.45),
            )
            .free_energy_model(phase_energy(6.0).temperature_window(1700.0, 3400.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/hypoeutectoid_ferrite",
                MetallurgicalPhaseKind::Ferrite,
                properties(95.0, 220.0, 0.35, 0.10, 80.0, 0.45),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.0, 0.055))
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.004, 0.025, 4.2)
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
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.035, 0.08, 5.0)
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
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.25, 0.25, 2.0)
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
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.035, 0.04, 3.5)
                    .temperature_window(0.0, 1000.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/hypereutectoid_cementite_network",
                MetallurgicalPhaseKind::Cementite,
                properties(900.0, 1250.0, 0.01, 0.62, 9.0, 0.18),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.055, 0.25))
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.12, 0.16, 3.2)
                    .temperature_window(0.0, 1100.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/bainite",
                MetallurgicalPhaseKind::Bainite,
                properties(430.0, 1050.0, 0.11, 0.24, 32.0, 0.32),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.012, 0.10))
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.035, 0.07, 1.7)
                    .temperature_window(520.0, 850.0)
                    .cooling_rate_stabilization(8.0, 9_000.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/martensite",
                MetallurgicalPhaseKind::Martensite,
                properties(650.0, 1500.0, 0.04, 0.30, 28.0, 0.25),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.01, 0.12))
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.04, 0.08, 0.0)
                    .temperature_window(0.0, 650.0)
                    .cooling_rate_stabilization(80.0, 40_000.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:fe_c/tempered_martensite",
                MetallurgicalPhaseKind::TemperedMartensite,
                properties(430.0, 1150.0, 0.12, 0.26, 31.0, 0.30),
            )
            .limit(ComponentLimit::new("destroy:carbon", 0.01, 0.12))
            .free_energy_model(
                composition_phase_energy("destroy:carbon", 0.04, 0.08, 4.6)
                    .temperature_window(650.0, 950.0),
            ),
        )
        .thermal_treatment_profile(plain_carbon_steel_treatment())
        .property_calibration(steel_property_calibration())
}

fn iron_carbon_chromium_nickel_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new(
        "metallurgy:fe_c_cr_ni",
        ["Fe", "destroy:carbon", "Cr", "Ni"],
    )
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.70),
            ("destroy:carbon", 0.01),
            ("Cr", 0.18),
            ("Ni", 0.11),
        ],
        1670.0,
        1770.0,
    ))
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.58),
            ("destroy:carbon", 0.02),
            ("Cr", 0.28),
            ("Ni", 0.12),
        ],
        1620.0,
        1740.0,
    ))
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.73),
            ("destroy:carbon", 0.01),
            ("Cr", 0.10),
            ("Ni", 0.16),
        ],
        1660.0,
        1760.0,
    ))
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_ni/liquid",
            MetallurgicalPhaseKind::Liquid,
            properties(65.0, 135.0, 0.50, 0.92, 30.0, 0.65),
        )
        .free_energy_model(phase_energy(6.0).temperature_window(1650.0, 3400.0)),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_ni/austenitic_matrix",
            MetallurgicalPhaseKind::Austenite,
            properties(190.0, 420.0, 0.42, 0.72, 18.0, 0.90),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.0, 0.04))
        .limit(ComponentLimit::new("Cr", 0.08, 0.30))
        .limit(ComponentLimit::new("Ni", 0.04, 0.30))
        .free_energy_model(
            composition_phase_energy("Ni", 0.12, 0.24, 5.0).temperature_window(300.0, 1750.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_ni/ferritic_matrix",
            MetallurgicalPhaseKind::Ferrite,
            properties(170.0, 360.0, 0.32, 0.66, 22.0, 0.86),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.0, 0.035))
        .limit(ComponentLimit::new("Cr", 0.10, 0.35))
        .limit(ComponentLimit::new("Ni", 0.0, 0.08))
        .free_energy_model(
            composition_phase_energy("Cr", 0.18, 0.25, 4.0).temperature_window(250.0, 1600.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_ni/chromium_carbides",
            MetallurgicalPhaseKind::Intermetallic,
            properties(950.0, 1300.0, 0.03, 0.80, 12.0, 0.35),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.006, 0.08))
        .limit(ComponentLimit::new("Cr", 0.08, 0.45))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.025, 0.06, 1.8)
                .temperature_window(300.0, 1250.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_ni/martensite",
            MetallurgicalPhaseKind::Martensite,
            properties(700.0, 1600.0, 0.05, 0.70, 20.0, 0.55),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.01, 0.10))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.035, 0.08, 0.0)
                .temperature_window(0.0, 650.0)
                .cooling_rate_stabilization(100.0, 25_000.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_ni/tempered_martensite",
            MetallurgicalPhaseKind::TemperedMartensite,
            properties(480.0, 1250.0, 0.16, 0.66, 22.0, 0.66),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.01, 0.10))
        .limit(ComponentLimit::new("Cr", 0.08, 0.30))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.035, 0.08, 4.2)
                .temperature_window(650.0, 980.0),
        ),
    )
    .thermal_treatment_profile(stainless_steel_treatment())
    .property_calibration(stainless_property_calibration())
}

fn iron_carbon_manganese_silicon_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new(
        "metallurgy:fe_c_mn_si",
        ["Fe", "destroy:carbon", "Mn", "destroy:silicon"],
    )
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.94),
            ("destroy:carbon", 0.02),
            ("Mn", 0.03),
            ("destroy:silicon", 0.01),
        ],
        1680.0,
        1780.0,
    ))
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.86),
            ("destroy:carbon", 0.04),
            ("Mn", 0.08),
            ("destroy:silicon", 0.02),
        ],
        1600.0,
        1740.0,
    ))
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.80),
            ("destroy:carbon", 0.08),
            ("Mn", 0.10),
            ("destroy:silicon", 0.02),
        ],
        1500.0,
        1680.0,
    ))
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_mn_si/liquid",
            MetallurgicalPhaseKind::Liquid,
            properties(62.0, 125.0, 0.52, 0.95, 33.0, 0.46),
        )
        .free_energy_model(phase_energy(6.0).temperature_window(1650.0, 3400.0)),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_mn_si/ferrite",
            MetallurgicalPhaseKind::Ferrite,
            properties(125.0, 300.0, 0.33, 0.18, 70.0, 0.48),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.0, 0.035))
        .limit(ComponentLimit::new("Mn", 0.0, 0.12))
        .limit(ComponentLimit::new("destroy:silicon", 0.0, 0.08))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.006, 0.04, 4.0)
                .temperature_window(250.0, 1185.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_mn_si/pearlite",
            MetallurgicalPhaseKind::Pearlite,
            properties(290.0, 620.0, 0.16, 0.25, 40.0, 0.38),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.018, 0.09))
        .limit(ComponentLimit::new("Mn", 0.0, 0.16))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.035, 0.05, 4.0)
                .temperature_window(250.0, 1000.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_mn_si/bainite",
            MetallurgicalPhaseKind::Bainite,
            properties(470.0, 1150.0, 0.10, 0.29, 30.0, 0.34),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.012, 0.10))
        .limit(ComponentLimit::new("Mn", 0.0, 0.18))
        .limit(ComponentLimit::new("destroy:silicon", 0.0, 0.08))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.04, 0.07, 2.2)
                .composition_term(CompositionEnergyTerm::new("Mn", 0.06, 0.12, 8_000.0))
                .temperature_window(520.0, 880.0)
                .cooling_rate_stabilization(6.0, 12_000.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_mn_si/austenite",
            MetallurgicalPhaseKind::Austenite,
            properties(165.0, 360.0, 0.30, 0.16, 16.0, 0.42),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.0, 0.10))
        .limit(ComponentLimit::new("Mn", 0.0, 0.20))
        .free_energy_model(
            composition_phase_energy("Mn", 0.05, 0.16, 3.0).temperature_window(850.0, 1800.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_mn_si/martensite",
            MetallurgicalPhaseKind::Martensite,
            properties(720.0, 1700.0, 0.035, 0.34, 25.0, 0.28),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.012, 0.12))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.045, 0.09, 0.0)
                .temperature_window(0.0, 650.0)
                .cooling_rate_stabilization(70.0, 35_000.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_mn_si/tempered_martensite",
            MetallurgicalPhaseKind::TemperedMartensite,
            properties(460.0, 1200.0, 0.13, 0.30, 29.0, 0.34),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.012, 0.12))
        .limit(ComponentLimit::new("Mn", 0.0, 0.20))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.045, 0.09, 4.4)
                .temperature_window(650.0, 980.0),
        ),
    )
    .thermal_treatment_profile(manganese_steel_treatment())
    .property_calibration(steel_property_calibration())
}

fn copper_zinc_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:cu_zn", ["Cu", "Zn"])
        .phase_boundary(binary_boundary("Cu", "Zn", 0.0, 1358.0, 1358.0))
        .phase_boundary(binary_boundary("Cu", "Zn", 0.30, 1180.0, 1250.0))
        .phase_boundary(binary_boundary("Cu", "Zn", 0.45, 1120.0, 1210.0))
        .phase_boundary(binary_boundary("Cu", "Zn", 1.0, 692.7, 692.7))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_zn/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(50.0, 100.0, 0.50, 0.20, 55.0, 0.55),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1150.0, 3000.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_zn/alpha_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(100.0, 260.0, 0.45, 0.08, 85.0, 0.55),
            )
            .limit(ComponentLimit::new("Zn", 0.0, 0.38))
            .free_energy_model(composition_phase_energy("Zn", 0.18, 0.38, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_zn/beta_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(170.0, 480.0, 0.25, 0.12, 45.0, 0.45),
            )
            .limit(ComponentLimit::new("Zn", 0.30, 0.55))
            .free_energy_model(composition_phase_energy("Zn", 0.45, 0.25, 3.0)),
        )
        .property_calibration(copper_property_calibration())
}

fn copper_tin_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:cu_sn", ["Cu", "Sn"])
        .phase_boundary(binary_boundary("Cu", "Sn", 0.0, 1358.0, 1358.0))
        .phase_boundary(binary_boundary("Cu", "Sn", 0.10, 1120.0, 1250.0))
        .phase_boundary(binary_boundary("Cu", "Sn", 0.22, 1070.0, 1120.0))
        .phase_boundary(binary_boundary("Cu", "Sn", 1.0, 505.0, 505.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_sn/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(45.0, 90.0, 0.45, 0.18, 50.0, 0.50),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1050.0, 2900.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_sn/alpha_bronze",
                MetallurgicalPhaseKind::SolidSolution,
                properties(130.0, 340.0, 0.35, 0.10, 70.0, 0.55),
            )
            .limit(ComponentLimit::new("Sn", 0.0, 0.16))
            .free_energy_model(composition_phase_energy("Sn", 0.08, 0.16, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_sn/intermetallic",
                MetallurgicalPhaseKind::Intermetallic,
                properties(360.0, 800.0, 0.06, 0.22, 25.0, 0.35),
            )
            .limit(ComponentLimit::new("Sn", 0.12, 0.45))
            .free_energy_model(composition_phase_energy("Sn", 0.25, 0.22, 2.5)),
        )
        .property_calibration(copper_property_calibration())
}

fn tin_lead_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:sn_pb", ["Sn", "Pb"])
        .phase_boundary(binary_boundary("Sn", "Pb", 0.0, 505.0, 505.0))
        .phase_boundary(binary_boundary("Sn", "Pb", 0.38, 456.0, 456.0))
        .phase_boundary(binary_boundary("Sn", "Pb", 0.85, 470.0, 585.0))
        .phase_boundary(binary_boundary("Sn", "Pb", 1.0, 601.0, 601.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_pb/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(8.0, 16.0, 0.55, 0.13, 48.0, 0.45),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(455.0, 1400.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_pb/tin_rich_solid",
                MetallurgicalPhaseKind::SolidSolution,
                properties(14.0, 24.0, 0.36, 0.13, 62.0, 0.52),
            )
            .limit(ComponentLimit::new("Pb", 0.0, 0.18))
            .free_energy_model(composition_phase_energy("Pb", 0.04, 0.12, 4.2)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_pb/lead_rich_solid",
                MetallurgicalPhaseKind::SolidSolution,
                properties(7.0, 16.0, 0.42, 0.19, 38.0, 0.42),
            )
            .limit(ComponentLimit::new("Pb", 0.60, 1.0))
            .free_energy_model(composition_phase_energy("Pb", 0.90, 0.22, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_pb/eutectic",
                MetallurgicalPhaseKind::Intermetallic,
                properties(18.0, 34.0, 0.30, 0.15, 48.0, 0.46),
            )
            .limit(ComponentLimit::new("Pb", 0.25, 0.55))
            .free_energy_model(
                composition_phase_energy("Pb", 0.38, 0.14, 2.0).temperature_window(0.0, 520.0),
            ),
        )
        .property_calibration(solder_property_calibration())
}

fn tin_silver_copper_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:sn_ag_cu", ["Sn", "Ag", "Cu"])
        .phase_boundary(multi_boundary(
            [("Sn", 0.965), ("Ag", 0.030), ("Cu", 0.005)],
            490.0,
            510.0,
        ))
        .phase_boundary(multi_boundary(
            [("Sn", 0.930), ("Ag", 0.050), ("Cu", 0.020)],
            500.0,
            535.0,
        ))
        .phase_boundary(multi_boundary(
            [("Sn", 0.880), ("Ag", 0.080), ("Cu", 0.040)],
            525.0,
            610.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_ag_cu/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(12.0, 24.0, 0.48, 0.12, 58.0, 0.56),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(490.0, 1500.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_ag_cu/tin_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(18.0, 42.0, 0.32, 0.12, 62.0, 0.56),
            )
            .limit(ComponentLimit::new("Ag", 0.0, 0.08))
            .limit(ComponentLimit::new("Cu", 0.0, 0.04))
            .free_energy_model(composition_phase_energy("Ag", 0.03, 0.06, 4.2)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_ag_cu/ag3sn",
                MetallurgicalPhaseKind::Intermetallic,
                properties(260.0, 430.0, 0.06, 0.11, 80.0, 0.66),
            )
            .limit(ComponentLimit::new("Ag", 0.015, 0.20))
            .free_energy_model(
                composition_phase_energy("Ag", 0.06, 0.10, 2.1).temperature_window(300.0, 900.0),
            )
            .kinetic_model(precipitation_kinetics(1.3e-4)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:sn_ag_cu/cu6sn5",
                MetallurgicalPhaseKind::Intermetallic,
                properties(430.0, 760.0, 0.025, 0.34, 36.0, 0.48),
            )
            .limit(ComponentLimit::new("Cu", 0.002, 0.12))
            .free_energy_model(
                composition_phase_energy("Cu", 0.015, 0.08, 1.9).temperature_window(300.0, 950.0),
            )
            .kinetic_model(precipitation_kinetics(1.7e-4)),
        )
        .property_calibration(solder_property_calibration())
}

fn bismuth_tin_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:bi_sn", ["Bi", "Sn"])
        .phase_boundary(binary_boundary("Bi", "Sn", 0.0, 545.0, 545.0))
        .phase_boundary(binary_boundary("Bi", "Sn", 0.43, 412.0, 412.0))
        .phase_boundary(binary_boundary("Bi", "Sn", 1.0, 505.0, 505.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:bi_sn/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(9.0, 18.0, 0.40, 0.55, 34.0, 0.55),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(410.0, 1300.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:bi_sn/bismuth_rich_solid",
                MetallurgicalPhaseKind::SolidSolution,
                properties(12.0, 28.0, 0.10, 0.95, 14.0, 0.58),
            )
            .limit(ComponentLimit::new("Sn", 0.0, 0.25))
            .free_energy_model(composition_phase_energy("Sn", 0.10, 0.18, 3.8)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:bi_sn/tin_rich_solid",
                MetallurgicalPhaseKind::SolidSolution,
                properties(15.0, 30.0, 0.30, 0.18, 58.0, 0.54),
            )
            .limit(ComponentLimit::new("Sn", 0.35, 1.0))
            .free_energy_model(composition_phase_energy("Sn", 0.80, 0.25, 3.8)),
        )
        .property_calibration(solder_property_calibration())
}

fn aluminum_silicon_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_si", ["Al", "destroy:silicon"])
        .phase_boundary(binary_boundary("Al", "destroy:silicon", 0.0, 933.0, 933.0))
        .phase_boundary(binary_boundary("Al", "destroy:silicon", 0.12, 850.0, 850.0))
        .phase_boundary(binary_boundary(
            "Al",
            "destroy:silicon",
            0.30,
            900.0,
            1200.0,
        ))
        .phase_boundary(binary_boundary(
            "Al",
            "destroy:silicon",
            1.0,
            1687.0,
            1687.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_si/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(35.0, 70.0, 0.50, 0.05, 95.0, 0.60),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(850.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_si/aluminum_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(55.0, 170.0, 0.35, 0.04, 120.0, 0.60),
            )
            .limit(ComponentLimit::new("destroy:silicon", 0.0, 0.15))
            .free_energy_model(composition_phase_energy(
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
            .free_energy_model(composition_phase_energy(
                "destroy:silicon",
                0.25,
                0.28,
                2.0,
            )),
        )
        .property_calibration(aluminum_property_calibration())
}

fn aluminum_copper_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_cu", ["Al", "Cu"])
        .phase_boundary(binary_boundary("Al", "Cu", 0.0, 933.0, 933.0))
        .phase_boundary(binary_boundary("Al", "Cu", 0.17, 821.0, 821.0))
        .phase_boundary(binary_boundary("Al", "Cu", 0.33, 820.0, 900.0))
        .phase_boundary(binary_boundary("Al", "Cu", 1.0, 1358.0, 1358.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_cu/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(35.0, 75.0, 0.52, 0.06, 85.0, 0.45),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(850.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_cu/aluminum_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(75.0, 240.0, 0.30, 0.055, 110.0, 0.48),
            )
            .limit(ComponentLimit::new("Cu", 0.0, 0.08))
            .free_energy_model(composition_phase_energy("Cu", 0.025, 0.08, 4.5)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_cu/theta_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(420.0, 900.0, 0.04, 0.16, 45.0, 0.30),
            )
            .limit(ComponentLimit::new("Cu", 0.02, 0.35))
            .free_energy_model(
                composition_phase_energy("Cu", 0.10, 0.18, 2.5).temperature_window(250.0, 825.0),
            ),
        )
        .thermal_treatment_profile(aluminum_precipitation_treatment())
        .property_calibration(aluminum_property_calibration())
}

fn aluminum_copper_magnesium_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_cu_mg", ["Al", "Cu", "Mg"])
        .phase_boundary(multi_boundary(
            [("Al", 0.93), ("Cu", 0.05), ("Mg", 0.02)],
            780.0,
            910.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.88), ("Cu", 0.08), ("Mg", 0.04)],
            760.0,
            890.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.80), ("Cu", 0.14), ("Mg", 0.06)],
            735.0,
            870.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_cu_mg/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(34.0, 72.0, 0.50, 0.065, 84.0, 0.44),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(850.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_cu_mg/aluminum_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(105.0, 320.0, 0.24, 0.075, 92.0, 0.46),
            )
            .limit(ComponentLimit::new("Cu", 0.0, 0.12))
            .limit(ComponentLimit::new("Mg", 0.0, 0.08))
            .free_energy_model(
                composition_phase_energy("Cu", 0.04, 0.10, 4.2)
                    .composition_term(CompositionEnergyTerm::new("Mg", 0.02, 0.06, 14_000.0)),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_cu_mg/theta_al2cu",
                MetallurgicalPhaseKind::Intermetallic,
                properties(460.0, 980.0, 0.035, 0.17, 42.0, 0.30),
            )
            .limit(ComponentLimit::new("Cu", 0.025, 0.28))
            .free_energy_model(
                composition_phase_energy("Cu", 0.13, 0.16, 2.4).temperature_window(280.0, 820.0),
            )
            .kinetic_model(precipitation_kinetics(1.6e-4)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_cu_mg/s_phase_al2cumg",
                MetallurgicalPhaseKind::Intermetallic,
                properties(540.0, 1120.0, 0.028, 0.20, 34.0, 0.28),
            )
            .limit(ComponentLimit::new("Cu", 0.035, 0.24))
            .limit(ComponentLimit::new("Mg", 0.012, 0.16))
            .free_energy_model(
                composition_phase_energy("Cu", 0.09, 0.14, 2.1)
                    .composition_term(CompositionEnergyTerm::new("Mg", 0.04, 0.10, 12_000.0))
                    .temperature_window(280.0, 780.0),
            )
            .kinetic_model(precipitation_kinetics(2.2e-4)),
        )
        .thermal_treatment_profile(aluminum_precipitation_treatment())
        .property_calibration(aluminum_property_calibration())
}

fn aluminum_magnesium_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_mg", ["Al", "Mg"])
        .phase_boundary(binary_boundary("Al", "Mg", 0.0, 933.0, 933.0))
        .phase_boundary(binary_boundary("Al", "Mg", 0.18, 720.0, 880.0))
        .phase_boundary(binary_boundary("Al", "Mg", 0.35, 710.0, 820.0))
        .phase_boundary(binary_boundary("Al", "Mg", 1.0, 923.0, 923.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_mg/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(32.0, 68.0, 0.55, 0.055, 95.0, 0.50),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(850.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_mg/solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(85.0, 260.0, 0.34, 0.07, 95.0, 0.56),
            )
            .limit(ComponentLimit::new("Mg", 0.0, 0.18))
            .free_energy_model(composition_phase_energy("Mg", 0.06, 0.16, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_mg/beta_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(300.0, 650.0, 0.06, 0.18, 40.0, 0.35),
            )
            .limit(ComponentLimit::new("Mg", 0.08, 0.45))
            .free_energy_model(
                composition_phase_energy("Mg", 0.22, 0.24, 2.0).temperature_window(250.0, 800.0),
            ),
        )
        .thermal_treatment_profile(aluminum_precipitation_treatment())
        .property_calibration(aluminum_property_calibration())
}

fn aluminum_magnesium_silicon_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_mg_si", ["Al", "Mg", "destroy:silicon"])
        .phase_boundary(multi_boundary(
            [("Al", 0.970), ("Mg", 0.010), ("destroy:silicon", 0.020)],
            840.0,
            925.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.930), ("Mg", 0.030), ("destroy:silicon", 0.040)],
            820.0,
            900.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.860), ("Mg", 0.080), ("destroy:silicon", 0.060)],
            790.0,
            880.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_mg_si/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(34.0, 70.0, 0.52, 0.06, 95.0, 0.58),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(820.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_mg_si/aluminum_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(90.0, 260.0, 0.30, 0.055, 110.0, 0.58),
            )
            .limit(ComponentLimit::new("Mg", 0.0, 0.08))
            .limit(ComponentLimit::new("destroy:silicon", 0.0, 0.08))
            .free_energy_model(
                composition_phase_energy("Mg", 0.025, 0.065, 4.1).composition_term(
                    CompositionEnergyTerm::new("destroy:silicon", 0.025, 0.060, 16_000.0),
                ),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_mg_si/mg2si_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(360.0, 760.0, 0.05, 0.18, 55.0, 0.48),
            )
            .limit(ComponentLimit::new("Mg", 0.008, 0.18))
            .limit(ComponentLimit::new("destroy:silicon", 0.008, 0.18))
            .free_energy_model(
                composition_phase_energy("Mg", 0.045, 0.08, 2.0)
                    .composition_term(CompositionEnergyTerm::new(
                        "destroy:silicon",
                        0.035,
                        0.08,
                        10_000.0,
                    ))
                    .temperature_window(300.0, 820.0),
            )
            .kinetic_model(precipitation_kinetics(2.0e-4)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_mg_si/silicon_rich",
                MetallurgicalPhaseKind::Intermetallic,
                properties(450.0, 700.0, 0.03, 0.25, 30.0, 0.35),
            )
            .limit(ComponentLimit::new("destroy:silicon", 0.04, 0.35))
            .free_energy_model(composition_phase_energy(
                "destroy:silicon",
                0.12,
                0.20,
                1.8,
            )),
        )
        .thermal_treatment_profile(aluminum_precipitation_treatment())
        .property_calibration(aluminum_property_calibration())
}

fn aluminum_zinc_magnesium_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_zn_mg", ["Al", "Zn", "Mg"])
        .phase_boundary(multi_boundary(
            [("Al", 0.90), ("Zn", 0.07), ("Mg", 0.03)],
            750.0,
            890.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.82), ("Zn", 0.12), ("Mg", 0.06)],
            720.0,
            850.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.70), ("Zn", 0.20), ("Mg", 0.10)],
            700.0,
            830.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_zn_mg/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(34.0, 70.0, 0.50, 0.065, 82.0, 0.44),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(820.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_zn_mg/aluminum_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(115.0, 360.0, 0.22, 0.085, 78.0, 0.45),
            )
            .limit(ComponentLimit::new("Zn", 0.0, 0.16))
            .limit(ComponentLimit::new("Mg", 0.0, 0.10))
            .free_energy_model(composition_phase_energy("Zn", 0.06, 0.14, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_zn_mg/eta_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(520.0, 1050.0, 0.035, 0.20, 35.0, 0.25),
            )
            .limit(ComponentLimit::new("Zn", 0.03, 0.35))
            .limit(ComponentLimit::new("Mg", 0.01, 0.20))
            .free_energy_model(
                composition_phase_energy("Zn", 0.14, 0.22, 2.6).temperature_window(250.0, 760.0),
            ),
        )
        .thermal_treatment_profile(aluminum_precipitation_treatment())
        .property_calibration(aluminum_property_calibration())
}

fn aluminum_zinc_magnesium_copper_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:al_zn_mg_cu", ["Al", "Zn", "Mg", "Cu"])
        .phase_boundary(multi_boundary(
            [("Al", 0.900), ("Zn", 0.055), ("Mg", 0.025), ("Cu", 0.020)],
            735.0,
            880.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.820), ("Zn", 0.105), ("Mg", 0.040), ("Cu", 0.035)],
            710.0,
            850.0,
        ))
        .phase_boundary(multi_boundary(
            [("Al", 0.720), ("Zn", 0.160), ("Mg", 0.060), ("Cu", 0.060)],
            690.0,
            830.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_zn_mg_cu/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(34.0, 72.0, 0.48, 0.07, 78.0, 0.45),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(800.0, 2800.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_zn_mg_cu/aluminum_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(140.0, 430.0, 0.18, 0.095, 72.0, 0.45),
            )
            .limit(ComponentLimit::new("Zn", 0.0, 0.18))
            .limit(ComponentLimit::new("Mg", 0.0, 0.10))
            .limit(ComponentLimit::new("Cu", 0.0, 0.08))
            .free_energy_model(
                composition_phase_energy("Zn", 0.075, 0.16, 4.0)
                    .composition_term(CompositionEnergyTerm::new("Mg", 0.035, 0.08, 12_000.0))
                    .composition_term(CompositionEnergyTerm::new("Cu", 0.025, 0.07, 10_000.0)),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_zn_mg_cu/eta_mgzn2",
                MetallurgicalPhaseKind::Intermetallic,
                properties(560.0, 1180.0, 0.030, 0.22, 32.0, 0.25),
            )
            .limit(ComponentLimit::new("Zn", 0.04, 0.35))
            .limit(ComponentLimit::new("Mg", 0.015, 0.18))
            .free_energy_model(
                composition_phase_energy("Zn", 0.16, 0.20, 2.1)
                    .composition_term(CompositionEnergyTerm::new("Mg", 0.055, 0.12, 11_000.0))
                    .temperature_window(300.0, 760.0),
            )
            .kinetic_model(precipitation_kinetics(2.3e-4)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:al_zn_mg_cu/al2cu_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(500.0, 950.0, 0.035, 0.18, 40.0, 0.30),
            )
            .limit(ComponentLimit::new("Cu", 0.015, 0.18))
            .free_energy_model(
                composition_phase_energy("Cu", 0.07, 0.12, 2.0).temperature_window(300.0, 820.0),
            )
            .kinetic_model(precipitation_kinetics(1.8e-4)),
        )
        .thermal_treatment_profile(aluminum_precipitation_treatment())
        .property_calibration(aluminum_property_calibration())
}

fn nickel_chromium_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:ni_cr", ["Ni", "Cr"])
        .phase_boundary(binary_boundary("Ni", "Cr", 0.0, 1728.0, 1728.0))
        .phase_boundary(binary_boundary("Ni", "Cr", 0.20, 1680.0, 1780.0))
        .phase_boundary(binary_boundary("Ni", "Cr", 0.50, 1750.0, 1900.0))
        .phase_boundary(binary_boundary("Ni", "Cr", 1.0, 2180.0, 2180.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(70.0, 150.0, 0.45, 1.05, 20.0, 0.85),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1650.0, 3500.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr/nickel_chromium_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(210.0, 520.0, 0.28, 1.10, 18.0, 0.92),
            )
            .limit(ComponentLimit::new("Cr", 0.0, 0.45))
            .free_energy_model(composition_phase_energy("Cr", 0.18, 0.35, 4.5)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr/chromium_rich_phase",
                MetallurgicalPhaseKind::Intermetallic,
                properties(520.0, 1000.0, 0.08, 1.20, 14.0, 0.80),
            )
            .limit(ComponentLimit::new("Cr", 0.30, 0.80))
            .free_energy_model(composition_phase_energy("Cr", 0.55, 0.30, 2.0)),
        )
        .property_calibration(nickel_property_calibration())
}

fn nickel_chromium_aluminum_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:ni_cr_al", ["Ni", "Cr", "Al"])
        .phase_boundary(multi_boundary(
            [("Ni", 0.80), ("Cr", 0.14), ("Al", 0.06)],
            1600.0,
            1720.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ni", 0.72), ("Cr", 0.16), ("Al", 0.12)],
            1580.0,
            1710.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ni", 0.62), ("Cr", 0.18), ("Al", 0.20)],
            1540.0,
            1690.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_al/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(68.0, 145.0, 0.42, 1.08, 18.0, 0.86),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1650.0, 3600.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_al/gamma_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(260.0, 720.0, 0.22, 1.18, 16.0, 0.94),
            )
            .limit(ComponentLimit::new("Cr", 0.05, 0.30))
            .limit(ComponentLimit::new("Al", 0.0, 0.18))
            .free_energy_model(
                composition_phase_energy("Cr", 0.16, 0.24, 4.4)
                    .composition_term(CompositionEnergyTerm::new("Al", 0.08, 0.14, 10_000.0))
                    .temperature_window(300.0, 1700.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_al/gamma_prime_ni3al",
                MetallurgicalPhaseKind::Intermetallic,
                properties(620.0, 1250.0, 0.07, 1.35, 12.0, 0.90),
            )
            .limit(ComponentLimit::new("Al", 0.06, 0.28))
            .free_energy_model(
                composition_phase_energy("Al", 0.18, 0.16, 2.2).temperature_window(700.0, 1450.0),
            )
            .kinetic_model(precipitation_kinetics(1.1e-4)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_al/chromium_rich",
                MetallurgicalPhaseKind::Intermetallic,
                properties(560.0, 1080.0, 0.06, 1.24, 14.0, 0.88),
            )
            .limit(ComponentLimit::new("Cr", 0.22, 0.55))
            .free_energy_model(composition_phase_energy("Cr", 0.35, 0.28, 1.8)),
        )
        .thermal_treatment_profile(nickel_superalloy_treatment())
        .property_calibration(nickel_property_calibration())
}

fn nickel_chromium_cobalt_molybdenum_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:ni_cr_co_mo", ["Ni", "Cr", "Co", "Mo"])
        .phase_boundary(multi_boundary(
            [("Ni", 0.62), ("Cr", 0.20), ("Co", 0.12), ("Mo", 0.06)],
            1580.0,
            1710.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ni", 0.55), ("Cr", 0.22), ("Co", 0.15), ("Mo", 0.08)],
            1560.0,
            1700.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ni", 0.48), ("Cr", 0.26), ("Co", 0.16), ("Mo", 0.10)],
            1540.0,
            1690.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_co_mo/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(75.0, 160.0, 0.42, 1.10, 18.0, 0.92),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1620.0, 3700.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_co_mo/gamma_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(300.0, 820.0, 0.20, 1.22, 15.0, 0.94),
            )
            .limit(ComponentLimit::new("Cr", 0.10, 0.35))
            .limit(ComponentLimit::new("Co", 0.0, 0.28))
            .limit(ComponentLimit::new("Mo", 0.0, 0.16))
            .free_energy_model(
                composition_phase_energy("Cr", 0.22, 0.25, 4.4)
                    .composition_term(CompositionEnergyTerm::new("Co", 0.12, 0.18, 8_000.0))
                    .composition_term(CompositionEnergyTerm::new("Mo", 0.06, 0.12, 12_000.0))
                    .temperature_window(300.0, 1720.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_co_mo/molybdenum_rich_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(720.0, 1400.0, 0.04, 1.35, 12.0, 0.88),
            )
            .limit(ComponentLimit::new("Mo", 0.04, 0.24))
            .free_energy_model(
                composition_phase_energy("Mo", 0.11, 0.14, 1.9).temperature_window(600.0, 1450.0),
            )
            .kinetic_model(precipitation_kinetics(9.0e-5)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ni_cr_co_mo/chromium_rich",
                MetallurgicalPhaseKind::Intermetallic,
                properties(580.0, 1120.0, 0.06, 1.25, 14.0, 0.90),
            )
            .limit(ComponentLimit::new("Cr", 0.24, 0.55))
            .free_energy_model(composition_phase_energy("Cr", 0.34, 0.26, 1.8)),
        )
        .thermal_treatment_profile(nickel_superalloy_treatment())
        .property_calibration(nickel_property_calibration())
}

fn copper_nickel_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:cu_ni", ["Cu", "Ni"])
        .phase_boundary(binary_boundary("Cu", "Ni", 0.0, 1358.0, 1358.0))
        .phase_boundary(binary_boundary("Cu", "Ni", 0.30, 1500.0, 1580.0))
        .phase_boundary(binary_boundary("Cu", "Ni", 0.70, 1620.0, 1700.0))
        .phase_boundary(binary_boundary("Cu", "Ni", 1.0, 1728.0, 1728.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_ni/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(48.0, 95.0, 0.48, 0.28, 48.0, 0.72),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1400.0, 3200.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_ni/continuous_solid_solution",
                MetallurgicalPhaseKind::SolidSolution,
                properties(145.0, 390.0, 0.36, 0.36, 42.0, 0.82),
            )
            .limit(ComponentLimit::new("Ni", 0.0, 1.0))
            .free_energy_model(composition_phase_energy("Ni", 0.35, 0.65, 5.0)),
        )
        .property_calibration(copper_property_calibration())
}

fn copper_beryllium_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:cu_be", ["Cu", "Be"])
        .phase_boundary(binary_boundary("Cu", "Be", 0.0, 1358.0, 1358.0))
        .phase_boundary(binary_boundary("Cu", "Be", 0.02, 1140.0, 1280.0))
        .phase_boundary(binary_boundary("Cu", "Be", 0.12, 1120.0, 1400.0))
        .phase_boundary(binary_boundary("Cu", "Be", 1.0, 1560.0, 1560.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_be/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(50.0, 100.0, 0.45, 0.08, 95.0, 0.66),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1120.0, 3000.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_be/copper_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(160.0, 520.0, 0.30, 0.045, 210.0, 0.66),
            )
            .limit(ComponentLimit::new("Be", 0.0, 0.045))
            .free_energy_model(composition_phase_energy("Be", 0.018, 0.040, 4.2)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_be/cube_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(560.0, 1200.0, 0.03, 0.18, 70.0, 0.70),
            )
            .limit(ComponentLimit::new("Be", 0.010, 0.25))
            .free_energy_model(
                composition_phase_energy("Be", 0.085, 0.16, 1.9).temperature_window(480.0, 1150.0),
            )
            .kinetic_model(precipitation_kinetics(2.0e-4)),
        )
        .thermal_treatment_profile(copper_beryllium_treatment())
        .property_calibration(copper_property_calibration())
}

fn copper_aluminum_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:cu_al", ["Cu", "Al"])
        .phase_boundary(binary_boundary("Cu", "Al", 0.0, 1358.0, 1358.0))
        .phase_boundary(binary_boundary("Cu", "Al", 0.10, 1220.0, 1320.0))
        .phase_boundary(binary_boundary("Cu", "Al", 0.25, 1180.0, 1380.0))
        .phase_boundary(binary_boundary("Cu", "Al", 1.0, 933.0, 933.0))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_al/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(48.0, 95.0, 0.48, 0.08, 120.0, 0.60),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1180.0, 3000.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_al/alpha_aluminum_bronze",
                MetallurgicalPhaseKind::SolidSolution,
                properties(180.0, 520.0, 0.30, 0.08, 95.0, 0.72),
            )
            .limit(ComponentLimit::new("Al", 0.0, 0.13))
            .free_energy_model(composition_phase_energy("Al", 0.07, 0.13, 4.2)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_al/beta_aluminum_bronze",
                MetallurgicalPhaseKind::SolidSolution,
                properties(280.0, 820.0, 0.16, 0.12, 58.0, 0.66),
            )
            .limit(ComponentLimit::new("Al", 0.10, 0.22))
            .free_energy_model(composition_phase_energy("Al", 0.16, 0.13, 2.4)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:cu_al/intermetallic",
                MetallurgicalPhaseKind::Intermetallic,
                properties(520.0, 1100.0, 0.04, 0.24, 32.0, 0.58),
            )
            .limit(ComponentLimit::new("Al", 0.18, 0.45))
            .free_energy_model(composition_phase_energy("Al", 0.28, 0.22, 1.8)),
        )
        .property_calibration(copper_property_calibration())
}

fn magnesium_aluminum_zinc_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:mg_al_zn", ["Mg", "Al", "Zn"])
        .phase_boundary(multi_boundary(
            [("Mg", 0.92), ("Al", 0.06), ("Zn", 0.02)],
            720.0,
            900.0,
        ))
        .phase_boundary(multi_boundary(
            [("Mg", 0.84), ("Al", 0.12), ("Zn", 0.04)],
            690.0,
            860.0,
        ))
        .phase_boundary(multi_boundary(
            [("Mg", 0.75), ("Al", 0.18), ("Zn", 0.07)],
            650.0,
            820.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:mg_al_zn/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(28.0, 58.0, 0.58, 0.09, 92.0, 0.42),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(820.0, 2500.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:mg_al_zn/magnesium_matrix",
                MetallurgicalPhaseKind::SolidSolution,
                properties(70.0, 230.0, 0.30, 0.12, 74.0, 0.45),
            )
            .limit(ComponentLimit::new("Al", 0.0, 0.20))
            .limit(ComponentLimit::new("Zn", 0.0, 0.10))
            .free_energy_model(
                composition_phase_energy("Al", 0.08, 0.16, 4.2)
                    .composition_term(CompositionEnergyTerm::new("Zn", 0.03, 0.08, 9_000.0)),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:mg_al_zn/beta_mg17al12",
                MetallurgicalPhaseKind::Intermetallic,
                properties(220.0, 520.0, 0.08, 0.18, 45.0, 0.30),
            )
            .limit(ComponentLimit::new("Al", 0.08, 0.35))
            .free_energy_model(
                composition_phase_energy("Al", 0.22, 0.20, 2.1).temperature_window(300.0, 740.0),
            )
            .kinetic_model(precipitation_kinetics(1.4e-4)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:mg_al_zn/mgzn_precipitate",
                MetallurgicalPhaseKind::Intermetallic,
                properties(260.0, 620.0, 0.055, 0.20, 38.0, 0.28),
            )
            .limit(ComponentLimit::new("Zn", 0.025, 0.18))
            .free_energy_model(
                composition_phase_energy("Zn", 0.08, 0.12, 1.9).temperature_window(300.0, 700.0),
            )
            .kinetic_model(precipitation_kinetics(1.8e-4)),
        )
        .thermal_treatment_profile(magnesium_precipitation_treatment())
        .property_calibration(magnesium_property_calibration())
}

fn titanium_aluminum_vanadium_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:ti_al_v", ["Ti", "Al", "V"])
        .phase_boundary(multi_boundary(
            [("Ti", 0.90), ("Al", 0.06), ("V", 0.04)],
            1850.0,
            1940.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ti", 0.82), ("Al", 0.12), ("V", 0.06)],
            1780.0,
            1910.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ti", 0.74), ("Al", 0.16), ("V", 0.10)],
            1720.0,
            1880.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_v/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(55.0, 110.0, 0.48, 0.95, 28.0, 0.75),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1900.0, 3600.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_v/alpha_titanium",
                MetallurgicalPhaseKind::SolidSolution,
                properties(320.0, 900.0, 0.15, 1.10, 18.0, 0.82),
            )
            .limit(ComponentLimit::new("Al", 0.0, 0.16))
            .limit(ComponentLimit::new("V", 0.0, 0.12))
            .free_energy_model(composition_phase_energy("Al", 0.08, 0.14, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_v/beta_titanium",
                MetallurgicalPhaseKind::SolidSolution,
                properties(280.0, 820.0, 0.22, 1.00, 22.0, 0.76),
            )
            .limit(ComponentLimit::new("V", 0.04, 0.30))
            .free_energy_model(
                composition_phase_energy("V", 0.12, 0.20, 3.5).temperature_window(900.0, 1900.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_v/titanium_aluminide",
                MetallurgicalPhaseKind::Intermetallic,
                properties(460.0, 1050.0, 0.06, 1.25, 15.0, 0.78),
            )
            .limit(ComponentLimit::new("Al", 0.12, 0.45))
            .free_energy_model(composition_phase_energy("Al", 0.28, 0.25, 1.8)),
        )
        .thermal_treatment_profile(titanium_alpha_beta_treatment())
        .property_calibration(titanium_property_calibration())
}

fn titanium_aluminum_molybdenum_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new("metallurgy:ti_al_mo", ["Ti", "Al", "Mo"])
        .phase_boundary(multi_boundary(
            [("Ti", 0.90), ("Al", 0.06), ("Mo", 0.04)],
            1840.0,
            1940.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ti", 0.82), ("Al", 0.10), ("Mo", 0.08)],
            1780.0,
            1910.0,
        ))
        .phase_boundary(multi_boundary(
            [("Ti", 0.74), ("Al", 0.14), ("Mo", 0.12)],
            1720.0,
            1880.0,
        ))
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_mo/liquid",
                MetallurgicalPhaseKind::Liquid,
                properties(58.0, 120.0, 0.46, 1.05, 24.0, 0.78),
            )
            .free_energy_model(phase_energy(5.0).temperature_window(1900.0, 3700.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_mo/alpha_titanium",
                MetallurgicalPhaseKind::SolidSolution,
                properties(330.0, 940.0, 0.14, 1.15, 17.0, 0.84),
            )
            .limit(ComponentLimit::new("Al", 0.0, 0.16))
            .limit(ComponentLimit::new("Mo", 0.0, 0.08))
            .free_energy_model(composition_phase_energy("Al", 0.08, 0.14, 4.0)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_mo/beta_titanium",
                MetallurgicalPhaseKind::SolidSolution,
                properties(300.0, 900.0, 0.20, 1.08, 21.0, 0.78),
            )
            .limit(ComponentLimit::new("Mo", 0.04, 0.24))
            .free_energy_model(
                composition_phase_energy("Mo", 0.10, 0.18, 3.2).temperature_window(850.0, 1900.0),
            ),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_mo/titanium_aluminide",
                MetallurgicalPhaseKind::Intermetallic,
                properties(470.0, 1080.0, 0.055, 1.25, 15.0, 0.78),
            )
            .limit(ComponentLimit::new("Al", 0.12, 0.42))
            .free_energy_model(composition_phase_energy("Al", 0.28, 0.25, 1.8)),
        )
        .phase_model(
            MetallurgicalPhaseModel::new(
                "metallurgy:ti_al_mo/molybdenum_rich_beta",
                MetallurgicalPhaseKind::SolidSolution,
                properties(360.0, 980.0, 0.16, 1.20, 19.0, 0.76),
            )
            .limit(ComponentLimit::new("Mo", 0.12, 0.35))
            .free_energy_model(composition_phase_energy("Mo", 0.22, 0.20, 2.0)),
        )
        .thermal_treatment_profile(titanium_alpha_beta_treatment())
        .property_calibration(titanium_property_calibration())
}

fn iron_carbon_chromium_molybdenum_vanadium_system() -> MetallurgicalSystem {
    MetallurgicalSystem::new(
        "metallurgy:fe_c_cr_mo_v",
        ["Fe", "destroy:carbon", "Cr", "Mo", "V"],
    )
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.880),
            ("destroy:carbon", 0.012),
            ("Cr", 0.050),
            ("Mo", 0.035),
            ("V", 0.023),
        ],
        1680.0,
        1800.0,
    ))
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.800),
            ("destroy:carbon", 0.020),
            ("Cr", 0.090),
            ("Mo", 0.060),
            ("V", 0.030),
        ],
        1620.0,
        1770.0,
    ))
    .phase_boundary(multi_boundary(
        [
            ("Fe", 0.720),
            ("destroy:carbon", 0.030),
            ("Cr", 0.120),
            ("Mo", 0.085),
            ("V", 0.045),
        ],
        1580.0,
        1740.0,
    ))
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_mo_v/liquid",
            MetallurgicalPhaseKind::Liquid,
            properties(68.0, 145.0, 0.48, 1.00, 28.0, 0.74),
        )
        .free_energy_model(phase_energy(6.0).temperature_window(1650.0, 3600.0)),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_mo_v/tempered_martensite",
            MetallurgicalPhaseKind::TemperedMartensite,
            properties(620.0, 1650.0, 0.10, 0.52, 24.0, 0.72),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.006, 0.06))
        .limit(ComponentLimit::new("Cr", 0.02, 0.16))
        .limit(ComponentLimit::new("Mo", 0.01, 0.12))
        .limit(ComponentLimit::new("V", 0.0, 0.08))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.025, 0.05, 3.8)
                .composition_term(CompositionEnergyTerm::new("Cr", 0.08, 0.14, 8_000.0))
                .composition_term(CompositionEnergyTerm::new("Mo", 0.04, 0.10, 12_000.0))
                .temperature_window(650.0, 1000.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_mo_v/martensite",
            MetallurgicalPhaseKind::Martensite,
            properties(780.0, 2100.0, 0.035, 0.56, 22.0, 0.66),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.010, 0.075))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.035, 0.065, 0.0)
                .temperature_window(0.0, 700.0)
                .cooling_rate_stabilization(60.0, 45_000.0),
        ),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_mo_v/chromium_carbides",
            MetallurgicalPhaseKind::Intermetallic,
            properties(1050.0, 1700.0, 0.01, 0.72, 8.0, 0.50),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.008, 0.08))
        .limit(ComponentLimit::new("Cr", 0.04, 0.25))
        .free_energy_model(
            composition_phase_energy("destroy:carbon", 0.030, 0.060, 1.6)
                .composition_term(CompositionEnergyTerm::new("Cr", 0.10, 0.18, 10_000.0))
                .temperature_window(350.0, 1350.0),
        )
        .kinetic_model(precipitation_kinetics(1.0e-4)),
    )
    .phase_model(
        MetallurgicalPhaseModel::new(
            "metallurgy:fe_c_cr_mo_v/molybdenum_vanadium_carbides",
            MetallurgicalPhaseKind::Intermetallic,
            properties(1180.0, 1850.0, 0.008, 0.60, 16.0, 0.62),
        )
        .limit(ComponentLimit::new("destroy:carbon", 0.006, 0.07))
        .limit(ComponentLimit::new("Mo", 0.015, 0.16))
        .limit(ComponentLimit::new("V", 0.005, 0.12))
        .free_energy_model(
            composition_phase_energy("Mo", 0.055, 0.12, 1.5)
                .composition_term(CompositionEnergyTerm::new("V", 0.030, 0.08, 12_000.0))
                .temperature_window(350.0, 1450.0),
        )
        .kinetic_model(precipitation_kinetics(8.0e-5)),
    )
    .thermal_treatment_profile(tool_steel_treatment())
    .property_calibration(tool_steel_property_calibration())
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

fn element(
    component: &'static str,
    melting_point_kelvin: f64,
    atomic_radius_pm: f64,
    crystal_structure: CrystalStructure,
    base_property_model: MetallurgicalPhasePropertyModel,
) -> MetallurgicalElementData {
    MetallurgicalElementData::new(
        component,
        melting_point_kelvin,
        atomic_radius_pm,
        crystal_structure,
        base_property_model,
    )
}

fn pair(
    first: &'static str,
    second: &'static str,
    solid_miscibility: SolidMiscibility,
    liquid_miscibility: LiquidMiscibility,
) -> MetallurgicalPairInteractionData {
    MetallurgicalPairInteractionData::new(first, second, solid_miscibility, liquid_miscibility)
}

fn compound<const N: usize>(
    id: &'static str,
    components: [(&'static str, f64); N],
    kind: MetallurgicalPhaseKind,
    property_model: MetallurgicalPhasePropertyModel,
    formation_energy_j_per_mol: f64,
) -> MetallurgicalCompoundPhaseData {
    MetallurgicalCompoundPhaseData::new(
        id,
        components,
        kind,
        property_model,
        formation_energy_j_per_mol,
    )
}

fn phase_energy(stability_depth: f64) -> PhaseFreeEnergyModel {
    PhaseFreeEnergyModel::new(-stability_depth * 5_000.0, 0.0)
}

fn composition_phase_energy(
    component: impl Into<crate::chemistry::metallurgy::MetallurgicalComponentId>,
    center_fraction: f64,
    width_fraction: f64,
    stability_depth: f64,
) -> PhaseFreeEnergyModel {
    phase_energy(stability_depth).composition_term(CompositionEnergyTerm::new(
        component,
        center_fraction,
        width_fraction,
        20_000.0,
    ))
}

fn precipitation_kinetics(precipitation_rate_per_second: f64) -> PhaseKineticModel {
    PhaseKineticModel::new(
        2.0e-13,
        175_000.0,
        0.01,
        0.0015,
        5.0e-5,
        precipitation_rate_per_second,
    )
}

fn steel_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(105.0, 22.0, 210.0, 340.0)
        .hardness_per_strength(0.31)
        .ductility_penalties(1250.0, 0.38, 0.62)
        .transport_penalties(0.035, 0.020, 0.18, 0.16)
}

fn stainless_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(92.0, 18.0, 160.0, 280.0)
        .hardness_per_strength(0.29)
        .ductility_penalties(1100.0, 0.30, 0.50)
        .transport_penalties(0.080, 0.040, 0.12, 0.10)
}

fn tool_steel_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(125.0, 26.0, 320.0, 420.0)
        .hardness_per_strength(0.34)
        .ductility_penalties(1500.0, 0.55, 0.72)
        .transport_penalties(0.060, 0.030, 0.20, 0.18)
}

fn aluminum_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(62.0, 12.0, 260.0, 220.0)
        .hardness_per_strength(0.27)
        .ductility_penalties(850.0, 0.34, 0.48)
        .transport_penalties(0.018, 0.012, 0.45, 0.35)
}

fn copper_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(50.0, 10.0, 180.0, 300.0)
        .hardness_per_strength(0.26)
        .ductility_penalties(750.0, 0.28, 0.45)
        .transport_penalties(0.025, 0.018, 0.55, 0.45)
}

fn nickel_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(95.0, 18.0, 300.0, 270.0)
        .hardness_per_strength(0.30)
        .ductility_penalties(1050.0, 0.36, 0.48)
        .transport_penalties(0.090, 0.045, 0.16, 0.12)
}

fn titanium_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(88.0, 17.0, 140.0, 260.0)
        .hardness_per_strength(0.29)
        .ductility_penalties(1000.0, 0.26, 0.52)
        .transport_penalties(0.060, 0.030, 0.14, 0.12)
}

fn magnesium_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(45.0, 9.0, 190.0, 180.0)
        .hardness_per_strength(0.25)
        .ductility_penalties(800.0, 0.32, 0.42)
        .transport_penalties(0.025, 0.016, 0.32, 0.25)
}

fn solder_property_calibration() -> MetallurgicalPropertyCalibration {
    MetallurgicalPropertyCalibration::neutral()
        .strength_response(22.0, 5.0, 55.0, 70.0)
        .hardness_per_strength(0.22)
        .ductility_penalties(650.0, 0.22, 0.35)
        .transport_penalties(0.020, 0.012, 0.20, 0.16)
}

fn plain_carbon_steel_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/plain_carbon_steel")
        .steel(
            1020.0,
            650.0,
            70.0,
            (520.0, 850.0),
            (3.0, 35.0),
            (650.0, 950.0),
        )
        .recovery_multiplier(1.15)
        .grain_growth_multiplier(1.20)
        .quench_vacancy_multiplier(1.15)
}

fn stainless_steel_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/stainless_steel")
        .steel(
            1120.0,
            620.0,
            95.0,
            (520.0, 840.0),
            (5.0, 45.0),
            (650.0, 980.0),
        )
        .precipitation_aging(1320.0, (700.0, 1150.0), 0.8)
        .recovery_multiplier(0.85)
        .grain_growth_multiplier(0.90)
        .quench_vacancy_multiplier(0.95)
}

fn manganese_steel_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/manganese_silicon_steel")
        .steel(
            1010.0,
            640.0,
            55.0,
            (520.0, 880.0),
            (2.0, 30.0),
            (650.0, 980.0),
        )
        .recovery_multiplier(1.05)
        .grain_growth_multiplier(1.00)
        .quench_vacancy_multiplier(1.10)
}

fn tool_steel_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/tool_steel")
        .steel(
            1080.0,
            700.0,
            50.0,
            (550.0, 880.0),
            (1.5, 25.0),
            (650.0, 1000.0),
        )
        .precipitation_aging(1120.0, (720.0, 920.0), 1.2)
        .recovery_multiplier(0.75)
        .grain_growth_multiplier(0.70)
        .quench_vacancy_multiplier(1.25)
}

fn aluminum_precipitation_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/aluminum_precipitation")
        .precipitation_aging(760.0, (360.0, 520.0), 1.6)
        .recovery_multiplier(1.35)
        .grain_growth_multiplier(1.35)
        .quench_vacancy_multiplier(0.75)
}

fn magnesium_precipitation_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/magnesium_precipitation")
        .precipitation_aging(650.0, (360.0, 500.0), 1.25)
        .recovery_multiplier(1.45)
        .grain_growth_multiplier(1.50)
        .quench_vacancy_multiplier(0.70)
}

fn copper_beryllium_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/copper_beryllium")
        .precipitation_aging(1030.0, (520.0, 700.0), 1.8)
        .recovery_multiplier(1.10)
        .grain_growth_multiplier(1.00)
        .quench_vacancy_multiplier(0.80)
}

fn nickel_superalloy_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/nickel_superalloy")
        .precipitation_aging(1450.0, (950.0, 1250.0), 1.35)
        .recovery_multiplier(0.65)
        .grain_growth_multiplier(0.55)
        .quench_vacancy_multiplier(0.65)
}

fn titanium_alpha_beta_treatment() -> ThermalTreatmentProfile {
    ThermalTreatmentProfile::new("metallurgy:thermal/titanium_alpha_beta")
        .precipitation_aging(1180.0, (720.0, 980.0), 0.65)
        .recovery_multiplier(0.95)
        .grain_growth_multiplier(0.75)
        .quench_vacancy_multiplier(0.85)
}

fn binary_boundary(
    first_component: &'static str,
    second_component: &'static str,
    second_fraction: f64,
    solidus_kelvin: f64,
    liquidus_kelvin: f64,
) -> PhaseBoundaryPoint {
    PhaseBoundaryPoint::new(
        [
            (first_component, 1.0 - second_fraction),
            (second_component, second_fraction),
        ],
        solidus_kelvin,
        liquidus_kelvin,
    )
}

fn multi_boundary<const N: usize>(
    composition: [(&'static str, f64); N],
    solidus_kelvin: f64,
    liquidus_kelvin: f64,
) -> PhaseBoundaryPoint {
    PhaseBoundaryPoint::new(composition, solidus_kelvin, liquidus_kelvin)
}
