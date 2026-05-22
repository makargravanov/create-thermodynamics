pub mod catalog;
pub mod error;
pub mod functional_group;
pub mod mixture;
pub mod molecule;
pub mod reaction;
pub mod reactions;
pub mod registry;
pub mod simulation;
pub mod substance;

pub use error::{ChemistryError, ChemistryResult};
pub use reactions::{DESTROY_EXPLICIT_REACTION_COUNT, DESTROY_REGISTERED_REACTION_COUNT};
pub use registry::{ChemistryRegistry, ChemistryRegistryBuilder};

pub fn destroy_registry_builder() -> ChemistryResult<ChemistryRegistryBuilder> {
    let builder = catalog::destroy_substances_registry_builder()?;
    reactions::destroy_reactions_registry_builder(builder)
}

#[cfg(test)]
mod tests {
    use super::destroy_registry_builder;
    use super::error::ChemistryError;
    use super::mixture::Mixture;
    use super::reaction::Reaction;
    use super::registry::ChemistryRegistryBuilder;
    use super::simulation::{react_for_tick, react_until_equilibrium};
    use super::substance::{Substance, SubstanceId};
    use super::{DESTROY_EXPLICIT_REACTION_COUNT, DESTROY_REGISTERED_REACTION_COUNT};

    fn water_id() -> SubstanceId {
        "destroy:water".into()
    }

    fn test_registry() -> super::ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:water",
                0,
                18.0,
                18_000.0,
                373.0,
                75.0,
                40_650.0,
            ))
            .substance(Substance::new(
                "destroy:hydrogen",
                0,
                2.0,
                1_000.0,
                20.0,
                28.8,
                900.0,
            ))
            .substance(Substance::new(
                "destroy:oxygen",
                0,
                32.0,
                1_140.0,
                90.0,
                29.4,
                6_820.0,
            ))
            .substance(Substance::new(
                "destroy:proton",
                1,
                1.0,
                1_000.0,
                10_000.0,
                0.0,
                0.0,
            ))
            .substance(Substance::new(
                "destroy:hydroxide",
                -1,
                17.0,
                17_000.0,
                10_000.0,
                75.0,
                0.0,
            ))
            .substance(Substance::new(
                "destroy:acid",
                0,
                10.0,
                10_000.0,
                350.0,
                80.0,
                20_000.0,
            ))
            .substance(Substance::new(
                "destroy:base",
                -1,
                9.0,
                9_000.0,
                350.0,
                70.0,
                20_000.0,
            ))
            .substance(Substance::new(
                "destroy:isomer_a",
                0,
                44.0,
                44_000.0,
                500.0,
                50.0,
                10_000.0,
            ))
            .substance(Substance::new(
                "destroy:isomer_b",
                0,
                44.0,
                44_000.0,
                500.0,
                50.0,
                10_000.0,
            ))
            .reaction(
                Reaction::builder("destroy:combustion")
                    .reactant("destroy:hydrogen", 2, 1)
                    .reactant("destroy:oxygen", 1, 1)
                    .product("destroy:water", 2)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .enthalpy_change_kj_per_mol(-240.0)
                    .build(),
            )
            .reaction(
                Reaction::builder("destroy:neutralization")
                    .reactant("destroy:proton", 1, 1)
                    .reactant("destroy:hydroxide", 1, 1)
                    .product("destroy:water", 1)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .enthalpy_change_kj_per_mol(-57.0)
                    .build(),
            )
            .reaction(
                Reaction::builder("destroy:acid_dissociation")
                    .reactant("destroy:acid", 1, 1)
                    .product("destroy:proton", 1)
                    .product("destroy:base", 1)
                    .pre_exponential_factor(1.0)
                    .activation_energy_kj_per_mol(2.0)
                    .enthalpy_change_kj_per_mol(4.0)
                    .build(),
            )
            .reaction(
                Reaction::builder("destroy:isomer_a_to_b")
                    .reactant("destroy:isomer_a", 1, 1)
                    .product("destroy:isomer_b", 1)
                    .pre_exponential_factor(1.0e8)
                    .activation_energy_kj_per_mol(20.0)
                    .enthalpy_change_kj_per_mol(5.0)
                    .reverse_reaction_id("destroy:isomer_b_to_a")
                    .build(),
            )
            .reaction(
                Reaction::builder("destroy:isomer_b_to_a")
                    .reactant("destroy:isomer_b", 1, 1)
                    .product("destroy:isomer_a", 1)
                    .pre_exponential_factor(1.0e8)
                    .activation_energy_kj_per_mol(15.0)
                    .enthalpy_change_kj_per_mol(-5.0)
                    .reverse_reaction_id("destroy:isomer_a_to_b")
                    .build(),
            )
            .build()
            .expect("test registry must be valid")
    }

    #[test]
    fn heating_empty_mixture_keeps_valid_temperature() {
        let registry = test_registry();
        let mut mixture = Mixture::empty();

        mixture.heat(&registry, 100_000.0).unwrap();

        assert!(mixture.temperature_kelvin().is_finite());
        assert_eq!(mixture.temperature_kelvin(), 298.0);
    }

    #[test]
    fn heating_water_to_boiling_consumes_latent_heat() {
        let registry = test_registry();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, water_id(), 1.0).unwrap();

        mixture
            .heat(&registry, 75.0 * 75.0 + 40_650.0 / 2.0)
            .unwrap();

        assert!((mixture.temperature_kelvin() - 373.0).abs() < 1.0e-9);
        assert!((mixture.gaseous_fraction_of(&water_id()) - 0.5).abs() < 1.0e-9);
    }

    #[test]
    fn mixing_conserves_amount_and_internal_energy() {
        let registry = test_registry();
        let mut cold = Mixture::new(300.0).unwrap();
        cold.add_substance(&registry, water_id(), 1.0).unwrap();
        let mut hot = Mixture::new(340.0).unwrap();
        hot.add_substance(&registry, water_id(), 1.0).unwrap();

        let mixed = Mixture::mix(&registry, &[(cold, 1.0), (hot, 1.0)]).unwrap();

        assert!((mixed.concentration_of(&water_id()) - 1.0).abs() < 1.0e-9);
        assert!((mixed.temperature_kelvin() - 320.0).abs() < 1.0e-9);
    }

    #[test]
    fn reaction_is_limited_by_available_reactant() {
        let registry = test_registry();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:hydrogen", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:oxygen", 1.0)
            .unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(mixture.concentration_of(&"destroy:hydrogen".into()) <= 1.0e-9);
        assert!((mixture.concentration_of(&water_id()) - 0.1).abs() < 1.0e-9);
        assert!(mixture.concentration_of(&"destroy:oxygen".into()) > 0.9);
    }

    #[test]
    fn reversible_reaction_moves_in_both_directions() {
        let registry = test_registry();
        let mut forward = Mixture::new(500.0).unwrap();
        forward
            .add_substance(&registry, "destroy:isomer_a", 1.0)
            .unwrap();
        react_for_tick(&registry, &mut forward, 1).unwrap();
        assert!(forward.concentration_of(&"destroy:isomer_b".into()) > 0.0);

        let mut reverse = Mixture::new(500.0).unwrap();
        reverse
            .add_substance(&registry, "destroy:isomer_b", 1.0)
            .unwrap();
        react_for_tick(&registry, &mut reverse, 1).unwrap();
        assert!(reverse.concentration_of(&"destroy:isomer_a".into()) > 0.0);
    }

    #[test]
    fn exothermic_reaction_heats_mixture() {
        let registry = test_registry();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:proton", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:hydroxide", 0.1)
            .unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(mixture.temperature_kelvin() > 298.0);
    }

    #[test]
    fn invalid_reaction_with_unknown_substance_fails_registry_build() {
        let error = ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:water",
                0,
                18.0,
                18_000.0,
                373.0,
                75.0,
                40_650.0,
            ))
            .reaction(
                Reaction::builder("destroy:bad")
                    .reactant("destroy:missing", 1, 1)
                    .product("destroy:water", 1)
                    .allow_mass_imbalance()
                    .build(),
            )
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::UnknownSubstance { .. }));
    }

    #[test]
    fn invalid_reaction_with_charge_imbalance_fails_registry_build() {
        let error = ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:proton",
                1,
                1.0,
                1_000.0,
                10_000.0,
                0.0,
                0.0,
            ))
            .substance(Substance::new(
                "destroy:water",
                0,
                18.0,
                18_000.0,
                373.0,
                75.0,
                40_650.0,
            ))
            .reaction(
                Reaction::builder("destroy:bad_charge")
                    .reactant("destroy:proton", 1, 1)
                    .product("destroy:water", 1)
                    .allow_mass_imbalance()
                    .build(),
            )
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::ChargeNotConserved { .. }));
    }

    #[test]
    fn destroy_reaction_catalog_builds() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();

        assert_eq!(DESTROY_EXPLICIT_REACTION_COUNT, 119);
        assert_eq!(
            DESTROY_REGISTERED_REACTION_COUNT,
            registry.reactions().count()
        );
        assert_eq!(DESTROY_REGISTERED_REACTION_COUNT, 149);
    }

    #[test]
    fn equilibrium_loop_stops_when_no_reaction_can_proceed() {
        let registry = test_registry();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, water_id(), 1.0).unwrap();

        let report = react_until_equilibrium(&registry, &mut mixture, 10, 1).unwrap();

        assert!(report.reached_equilibrium);
        assert_eq!(report.ticks, 1);
    }
}
