#[path = "core/alloy.rs"]
pub mod alloy;
#[path = "molecule/canonical.rs"]
pub mod canonical;
#[path = "data/catalog.rs"]
pub mod catalog;
#[path = "core/catalysis.rs"]
pub mod catalysis;
#[path = "core/complex.rs"]
pub mod complex;
#[path = "core/condition.rs"]
pub mod condition;
#[path = "dynamic/mod.rs"]
pub mod dynamic;
pub mod error;
#[path = "molecule/frowns.rs"]
pub mod frowns;
#[path = "molecule/functional_group.rs"]
pub mod functional_group;
#[path = "core/kinetics.rs"]
pub mod kinetics;
#[path = "core/metallurgy/mod.rs"]
pub mod metallurgy;
#[path = "data/metallurgy.rs"]
pub mod metallurgy_data;
#[path = "core/mixture.rs"]
pub mod mixture;
#[path = "molecule/mod.rs"]
pub mod molecule;
#[path = "organic/mod.rs"]
pub mod organic;
#[path = "core/reaction.rs"]
pub mod reaction;
#[path = "data/reactions.rs"]
pub mod reactions;
#[path = "molecule/reactive_site.rs"]
pub mod reactive_site;
#[path = "core/redox.rs"]
pub mod redox;
#[path = "core/registry.rs"]
pub mod registry;
#[path = "core/simulation.rs"]
pub mod simulation;
#[path = "core/solution.rs"]
pub mod solution;
#[path = "core/substance.rs"]
pub mod substance;
pub mod synthesis;
#[path = "core/thermodynamics.rs"]
pub mod thermodynamics;

#[path = "kinetics/mod.rs"]
pub mod selectivity;

pub use error::{ChemistryError, ChemistryResult};
pub use reactions::{
    DESTROY_EXPLICIT_REACTION_COUNT, DESTROY_METALLURGY_REACTION_COUNT,
    DESTROY_REGISTERED_REACTION_COUNT,
};
pub use registry::{ChemistryRegistry, ChemistryRegistryBuilder};

pub fn destroy_registry_builder() -> ChemistryResult<ChemistryRegistryBuilder> {
    let builder = catalog::destroy_substances_registry_builder()?;
    let builder = reactions::destroy_reactions_registry_builder(builder)?;
    let builder = reactions::destroy_metallurgy_reactions_registry_builder(builder)?;
    Ok(builder
        .metallurgical_systems(metallurgy_data::default_metallurgical_systems())
        .metallurgical_elements(metallurgy_data::default_metallurgical_element_data())
        .metallurgical_pair_interactions(metallurgy_data::default_metallurgical_pair_interactions())
        .metallurgical_compound_phases(metallurgy_data::default_metallurgical_compound_phases()))
}

pub fn destroy_registry_with_generated_reactions_builder(
) -> ChemistryResult<ChemistryRegistryBuilder> {
    organic::destroy_registry_with_generated_reactions_builder()
}

#[cfg(test)]
mod tests {
    use super::catalysis::{CatalystSurfaceId, CatalystSurfaceSpec};
    use super::complex::{ComplexGeometry, ComplexLigand, ComplexSpec, LigandExchangeLability};
    use super::condition::{AcidityCondition, ReactionCondition};
    use super::destroy_registry_builder;
    use super::error::ChemistryError;
    use super::kinetics::{
        ChannelConditionEffect, EnergyModel, IsomerEnergy, LightBand, ReactionChannel,
    };
    use super::mixture::{Mixture, MixturePhase};
    use super::reaction::{Reaction, StoichiometricTerm};
    use super::redox::{apply_electrolysis_cell, ElectrodeProcess, ElectrolysisCell, RedoxRole};
    use super::registry::{ChemistryRegistryBuilder, ReactionCandidateScratch};
    use super::selectivity::{
        NucleophileStrength, ReactionType, SelectivityContext, SelectivityEngine,
        SelectivityProfile, SiteDescriptorBuilder,
    };
    use super::simulation::{
        react_for_tick, react_for_tick_with_context, react_until_equilibrium, ReactionContext,
    };
    use super::solution::PrecipitationSpec;
    use super::substance::{Substance, SubstanceId, SubstancePhaseProperties};
    use super::{
        DESTROY_EXPLICIT_REACTION_COUNT, DESTROY_METALLURGY_REACTION_COUNT,
        DESTROY_REGISTERED_REACTION_COUNT,
    };

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
                    .reactant_phase_access("destroy:hydrogen", [MixturePhase::Gas])
                    .reactant_phase_access("destroy:oxygen", [MixturePhase::Gas])
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
            .reaction(
                Reaction::builder("destroy:display_reversible")
                    .reactant("destroy:acid", 1, 1)
                    .product("destroy:isomer_a", 1)
                    .display_as_reversible()
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .allow_mass_imbalance()
                    .build(),
            )
            .build()
            .expect("test registry must be valid")
    }

    fn test_substance(id: &'static str) -> Substance {
        Substance::new(id, 0, 10.0, 1_000.0, 373.0, 100.0, 20_000.0)
    }

    fn neutral_substance(id: &'static str, molar_mass: f64) -> Substance {
        Substance::new(id, 0, molar_mass, 1_000.0, 373.0, 100.0, 20_000.0)
    }

    fn surface_test_registry(reaction: Reaction) -> super::ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .reaction(reaction)
            .catalyst_surface_spec(CatalystSurfaceSpec::chemical("surface:nickel", 58.69, 0))
            .build()
            .unwrap()
    }

    fn distributed_product_registry() -> super::ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .substance(test_substance("c"))
            .reaction(
                Reaction::builder("split_stereo_products")
                    .reactant("a", 1, 1)
                    .product_distribution_variant(0.5, [("b", 1)])
                    .product_distribution_variant(0.5, [("c", 1)])
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .build()
            .unwrap()
    }

    #[test]
    fn external_products_are_counted_once_in_registry_validation() {
        ChemistryRegistryBuilder::new()
            .substance(neutral_substance("a", 10.0))
            .substance(neutral_substance("b", 5.0))
            .reaction(
                Reaction::builder("external_product_balance")
                    .reactant("a", 1, 1)
                    .product("b", 1)
                    .chemical_external_product("external:fragment", 1.0, 5.0, 0)
                    .build(),
            )
            .build()
            .expect("external product must be counted exactly once");
    }

    #[test]
    fn empty_reaction_and_substance_ids_are_rejected() {
        let reaction_error = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .reaction(
                Reaction::builder("")
                    .reactant("a", 1, 1)
                    .product("b", 1)
                    .build(),
            )
            .build()
            .unwrap_err();
        assert!(matches!(
            reaction_error,
            ChemistryError::InvalidReaction { .. }
        ));

        let substance_error = ChemistryRegistryBuilder::new()
            .substance(test_substance(""))
            .build()
            .unwrap_err();
        assert!(matches!(
            substance_error,
            ChemistryError::InvalidSubstance { .. }
        ));
    }

    #[test]
    fn external_catalyst_adds_surface_sites_without_overwriting() {
        let mut context = ReactionContext::default();
        context
            .add_external_catalyst("surface:nickel", 1.0)
            .unwrap();
        context
            .add_external_catalyst("surface:nickel", 1.0)
            .unwrap();

        assert_eq!(context.external_catalysts["surface:nickel"], 2.0);
        assert_eq!(
            context
                .surfaces
                .get(&CatalystSurfaceId::from("surface:nickel"))
                .unwrap()
                .total_sites_mol_per_bucket,
            2.0
        );
    }

    #[test]
    fn external_products_are_recorded_outside_mixture() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .reaction(
                Reaction::builder("vented_external_product")
                    .reactant("a", 1, 1)
                    .product("b", 1)
                    .chemical_external_product("external:vented", 1.0, 0.0, 0)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "a", 1.0).unwrap();
        let mut context = ReactionContext::default();

        assert!(react_for_tick_with_context(&registry, &mut mixture, &mut context, 1).unwrap());
        assert!(context.external_products["external:vented"] > 0.0);
        assert!(mixture.concentration_of(&"external:vented".into()) == 0.0);
    }

    fn channel_product_registry(
        first_activation: f64,
        second_activation: f64,
    ) -> super::ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .substance(test_substance("c"))
            .reaction(
                Reaction::builder("channel_products")
                    .reactant("a", 1, 1)
                    .channel(
                        ReactionChannel::new(
                            "channel:b",
                            [StoichiometricTerm::new("b", 1)],
                            first_activation,
                        )
                        .with_pre_exponential_factor(1.0e12),
                    )
                    .channel(
                        ReactionChannel::new(
                            "channel:c",
                            [StoichiometricTerm::new("c", 1)],
                            second_activation,
                        )
                        .with_pre_exponential_factor(1.0e12),
                    )
                    .build(),
            )
            .build()
            .unwrap()
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
    fn distributed_products_are_applied_as_concrete_substances() {
        let registry = distributed_product_registry();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "a", 1.0).unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert_eq!(mixture.concentration_of(&SubstanceId::from("a")), 0.0);
        assert!((mixture.concentration_of(&SubstanceId::from("b")) - 0.5).abs() < 1.0e-9);
        assert!((mixture.concentration_of(&SubstanceId::from("c")) - 0.5).abs() < 1.0e-9);
    }

    #[test]
    fn distributed_product_fractions_must_sum_to_one() {
        let result = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .substance(test_substance("c"))
            .reaction(
                Reaction::builder("bad_distribution")
                    .reactant("a", 1, 1)
                    .product_distribution_variant(0.5, [("b", 1)])
                    .product_distribution_variant(0.25, [("c", 1)])
                    .build(),
            )
            .build();

        assert!(matches!(
            result,
            Err(ChemistryError::InvalidReaction { .. })
        ));
    }

    #[test]
    fn distributed_products_are_checked_for_mass() {
        let result = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(Substance::new("b", 0, 8.0, 1_000.0, 373.0, 100.0, 20_000.0))
            .substance(Substance::new("c", 0, 8.0, 1_000.0, 373.0, 100.0, 20_000.0))
            .reaction(
                Reaction::builder("bad_mass_distribution")
                    .reactant("a", 1, 1)
                    .product_distribution_variant(0.5, [("b", 1)])
                    .product_distribution_variant(0.5, [("c", 1)])
                    .build(),
            )
            .build();

        assert!(matches!(
            result,
            Err(ChemistryError::MassNotConserved { .. })
        ));
    }

    #[test]
    fn equal_channel_barriers_give_equal_products() {
        let registry = channel_product_registry(0.0, 0.0);
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "a", 1.0).unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!((mixture.concentration_of(&SubstanceId::from("b")) - 0.5).abs() < 1.0e-9);
        assert!((mixture.concentration_of(&SubstanceId::from("c")) - 0.5).abs() < 1.0e-9);
    }

    #[test]
    fn lower_transition_state_energy_gives_larger_kinetic_product() {
        let registry = channel_product_registry(0.0, 5.0);
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "a", 1.0).unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(mixture.concentration_of(&SubstanceId::from("b")) > 0.85);
        assert!(mixture.concentration_of(&SubstanceId::from("c")) < 0.15);
    }

    #[test]
    fn higher_temperature_reduces_channel_selectivity() {
        let registry = channel_product_registry(0.0, 5.0);
        let mut cold = Mixture::new(250.0).unwrap();
        cold.add_substance(&registry, "a", 1.0).unwrap();
        let mut hot = Mixture::new(800.0).unwrap();
        hot.add_substance(&registry, "a", 1.0).unwrap();

        react_for_tick(&registry, &mut cold, 1).unwrap();
        react_for_tick(&registry, &mut hot, 1).unwrap();

        assert!(
            cold.concentration_of(&SubstanceId::from("b"))
                > hot.concentration_of(&SubstanceId::from("b"))
        );
        assert!(hot.concentration_of(&SubstanceId::from("b")) < 0.75);
    }

    #[test]
    fn photochemical_channel_requires_light() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .reaction(
                Reaction::builder("photo_channel")
                    .reactant("a", 1, 1)
                    .channel(
                        ReactionChannel::new("photo:b", [StoichiometricTerm::new("b", 1)], 0.0)
                            .with_pre_exponential_factor(1.0e12)
                            .with_condition_effect(ChannelConditionEffect::Light {
                                band: LightBand::Ultraviolet,
                                minimum_power: 0.1,
                                multiplier: 1.0,
                            }),
                    )
                    .build(),
            )
            .build()
            .unwrap();
        let mut dark = Mixture::new(298.0).unwrap();
        dark.add_substance(&registry, "a", 1.0).unwrap();
        let mut light = dark.clone();
        let mut context = ReactionContext::default()
            .with_light_power(LightBand::Ultraviolet, 1.0)
            .unwrap();

        assert!(!react_for_tick(&registry, &mut dark, 1).unwrap());
        assert!(react_for_tick_with_context(&registry, &mut light, &mut context, 1).unwrap());

        assert_eq!(dark.concentration_of(&SubstanceId::from("b")), 0.0);
        assert!(light.concentration_of(&SubstanceId::from("b")) > 0.0);
    }

    #[test]
    fn catalyst_surface_changes_channel_distribution() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .substance(test_substance("c"))
            .catalyst_surface_spec(CatalystSurfaceSpec::chemical("surface:nickel", 58.69, 0))
            .reaction(
                Reaction::builder("surface_channels")
                    .reactant("a", 1, 1)
                    .channel(
                        ReactionChannel::new(
                            "surface:none",
                            [StoichiometricTerm::new("b", 1)],
                            0.0,
                        )
                        .with_pre_exponential_factor(1.0e12),
                    )
                    .channel(
                        ReactionChannel::new(
                            "surface:nickel",
                            [StoichiometricTerm::new("c", 1)],
                            0.0,
                        )
                        .with_pre_exponential_factor(1.0e12)
                        .with_condition_effect(
                            ChannelConditionEffect::Surface {
                                surface_id: CatalystSurfaceId::from("surface:nickel"),
                                multiplier: 9.0,
                            },
                        ),
                    )
                    .build(),
            )
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "a", 1.0).unwrap();
        let mut context = ReactionContext::default();
        context.add_surface("surface:nickel", 1.0).unwrap();

        react_for_tick_with_context(&registry, &mut mixture, &mut context, 1).unwrap();

        assert!((mixture.concentration_of(&SubstanceId::from("c")) - 0.9).abs() < 1.0e-9);
    }

    #[test]
    fn energy_model_calculates_thermodynamic_isomer_distribution() {
        let model = EnergyModel::new()
            .with_isomer_energy(IsomerEnergy {
                substance_id: "e".into(),
                relative_gibbs_kj_per_mol: 0.0,
                phase: None,
                surface_id: None,
            })
            .unwrap()
            .with_isomer_energy(IsomerEnergy {
                substance_id: "z".into(),
                relative_gibbs_kj_per_mol: 5.0,
                phase: None,
                surface_id: None,
            })
            .unwrap();

        let distribution = model
            .equilibrium_distribution(["e".into(), "z".into()], 298.0, None)
            .unwrap();

        assert!(distribution[&SubstanceId::from("e")] > 0.85);
        assert!(distribution[&SubstanceId::from("z")] < 0.15);
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
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
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
    fn display_as_reversible_does_not_disable_reaction() {
        let registry = test_registry();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:acid", 0.1)
            .unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(mixture.concentration_of(&"destroy:isomer_a".into()) > 0.0);
    }

    #[test]
    fn registry_returns_only_reaction_candidates_for_present_substances() {
        let registry = test_registry();
        let hydrogen: SubstanceId = "destroy:hydrogen".into();
        let proton: SubstanceId = "destroy:proton".into();

        let hydrogen_candidates = registry
            .reaction_candidates_for_substances([&hydrogen])
            .into_iter()
            .map(|reaction| reaction.id.to_string())
            .collect::<Vec<_>>();
        assert_eq!(hydrogen_candidates, vec!["destroy:combustion"]);

        let proton_candidates = registry
            .reaction_candidates_for_substances([&proton])
            .into_iter()
            .map(|reaction| reaction.id.to_string())
            .collect::<Vec<_>>();
        assert!(proton_candidates.contains(&"destroy:neutralization".to_string()));
        assert!(!proton_candidates.contains(&"destroy:combustion".to_string()));
    }

    #[test]
    fn numeric_registry_indices_match_public_lookup() {
        let registry = test_registry();
        let water = water_id();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let combustion = "destroy:combustion".into();

        let water_index = registry.substance_index(&water).unwrap();
        let oxygen_index = registry.substance_index(&oxygen).unwrap();
        assert_ne!(water_index, oxygen_index);
        assert_eq!(registry.substance_by_index(water_index).unwrap().id, water);

        let reaction_index = registry.reaction_index(&combustion).unwrap();
        assert_eq!(
            registry.reaction_by_index(reaction_index).unwrap().id,
            combustion
        );
    }

    #[test]
    fn indexed_reaction_candidates_match_public_candidates() {
        let registry = test_registry();
        let hydrogen: SubstanceId = "destroy:hydrogen".into();
        let hydrogen_index = registry.substance_index(&hydrogen).unwrap();

        let public_candidates = registry
            .reaction_candidates_for_substances([&hydrogen])
            .into_iter()
            .map(|reaction| reaction.id.to_string())
            .collect::<Vec<_>>();
        let indexed_candidates = registry
            .reaction_candidate_indices_for_substance_indices([hydrogen_index])
            .into_iter()
            .map(|index| registry.reaction_by_index(index).unwrap().id.to_string())
            .collect::<Vec<_>>();

        assert_eq!(indexed_candidates, public_candidates);
    }

    #[test]
    fn reaction_candidate_scratch_deduplicates_and_keeps_unindexed_reactions() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:a",
                0,
                10.0,
                10_000.0,
                500.0,
                100.0,
                20_000.0,
            ))
            .substance(Substance::new(
                "destroy:b",
                0,
                10.0,
                10_000.0,
                500.0,
                100.0,
                20_000.0,
            ))
            .reaction(
                Reaction::builder("destroy:indexed_once")
                    .reactant("destroy:a", 1, 1)
                    .reactant("destroy:b", 1, 1)
                    .product("destroy:a", 1)
                    .product("destroy:b", 1)
                    .pre_exponential_factor(1.0)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .reaction(
                Reaction::builder("destroy:unindexed_uv")
                    .requires_uv()
                    .pre_exponential_factor(1.0)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .build()
            .unwrap();
        let a = registry.substance_index(&"destroy:a".into()).unwrap();
        let b = registry.substance_index(&"destroy:b".into()).unwrap();
        let mut scratch = ReactionCandidateScratch::new();

        registry.collect_reaction_candidate_indices_for_substance_indices([a, b], &mut scratch);
        let candidates = scratch
            .candidates()
            .iter()
            .map(|index| registry.reaction_by_index(*index).unwrap().id.to_string())
            .collect::<Vec<_>>();

        assert_eq!(
            candidates,
            vec!["destroy:unindexed_uv", "destroy:indexed_once"]
        );
    }

    #[test]
    fn uv_context_controls_reaction_rate() {
        let registry = ChemistryRegistryBuilder::new()
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
                "destroy:water",
                0,
                18.0,
                18_000.0,
                373.0,
                75.0,
                40_650.0,
            ))
            .reaction(
                Reaction::builder("destroy:uv_water")
                    .reactant("destroy:hydrogen", 2, 1)
                    .reactant("destroy:oxygen", 1, 1)
                    .product("destroy:water", 2)
                    .reactant_phase_access("destroy:hydrogen", [MixturePhase::Gas])
                    .reactant_phase_access("destroy:oxygen", [MixturePhase::Gas])
                    .requires_uv()
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .build()
            .unwrap();
        let mut dark = Mixture::new(298.0).unwrap();
        dark.add_substance(&registry, "destroy:hydrogen", 1.0)
            .unwrap();
        dark.add_substance(&registry, "destroy:oxygen", 1.0)
            .unwrap();
        react_for_tick(&registry, &mut dark, 1).unwrap();
        assert_eq!(dark.concentration_of(&"destroy:water".into()), 0.0);

        let mut lit = dark.clone();
        let mut context = ReactionContext::default().with_uv_power(1.0).unwrap();
        react_for_tick_with_context(&registry, &mut lit, &mut context, 1).unwrap();
        assert!(lit.concentration_of(&"destroy:water".into()) > 0.0);
    }

    #[test]
    fn external_reactant_is_consumed_and_catalyst_is_not() {
        let registry = ChemistryRegistryBuilder::new()
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
                "destroy:water",
                0,
                18.0,
                18_000.0,
                373.0,
                75.0,
                40_650.0,
            ))
            .reaction(
                Reaction::builder("destroy:external_water")
                    .reactant("destroy:hydrogen", 1, 1)
                    .product("destroy:water", 1)
                    .reactant_phase_access("destroy:hydrogen", [MixturePhase::Gas])
                    .chemical_external_reactant("external:oxygen_atom", 1.0, 16.0, 0)
                    .chemical_external_catalyst("external:nickel", 1.0, 58.69, 0)
                    .reaction_result("external:water_result", 1.0)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .catalyst_surface_spec(CatalystSurfaceSpec::chemical("external:nickel", 58.69, 0))
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:hydrogen", 1.0)
            .unwrap();
        let mut context = ReactionContext::default();
        context
            .add_external_reactant("external:oxygen_atom", 0.25)
            .unwrap();
        context
            .add_external_catalyst("external:nickel", 1.0)
            .unwrap();

        react_for_tick_with_context(&registry, &mut mixture, &mut context, 1).unwrap();

        assert_eq!(
            context
                .external_reactants
                .get("external:oxygen_atom")
                .copied()
                .unwrap_or(0.0),
            0.0
        );
        assert_eq!(context.external_catalysts["external:nickel"], 1.0);
        assert!(context.reaction_results["external:water_result"] > 0.0);
    }

    #[test]
    fn complex_equilibrium_binds_free_metal_and_ligands() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:copper_ii", 0.01)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:ammonia", 0.10)
            .unwrap();
        mixture
            .move_between_phases(
                &registry,
                "destroy:ammonia",
                MixturePhase::Gas,
                MixturePhase::Aqueous,
                0.10,
            )
            .unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(
            mixture.concentration_of(&"destroy:copper_ii_tetraammine".into()) > 0.0,
            "complex formation must create a real complex substance"
        );
        assert!(
            mixture.concentration_of(&"destroy:copper_ii".into()) < 0.01,
            "complex formation must reduce the free metal concentration"
        );
    }

    #[test]
    fn invalid_complex_charge_fails_registry_build() {
        let error = ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "metal",
                2,
                10.0,
                1_000.0,
                f64::MAX,
                100.0,
                20_000.0,
            ))
            .substance(Substance::new(
                "ligand", 0, 5.0, 1_000.0, 373.0, 100.0, 20_000.0,
            ))
            .complex_spec(ComplexSpec::new(
                "metal_ligand",
                "metal",
                [ComplexLigand::new("ligand", 2)],
                1,
                1.0e3,
            ))
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::ChargeNotConserved { .. }));
    }

    #[test]
    fn invalid_complex_coordination_fails_registry_build() {
        let error = ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "metal",
                2,
                10.0,
                1_000.0,
                f64::MAX,
                100.0,
                20_000.0,
            ))
            .substance(
                Substance::new("ligand", 0, 5.0, 1_000.0, 373.0, 100.0, 20_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_unlimited()),
            )
            .complex_spec(
                ComplexSpec::new(
                    "metal_ligand",
                    "metal",
                    [ComplexLigand::new("ligand", 2)],
                    2,
                    1.0e3,
                )
                .with_coordination_number(4)
                .with_geometry(ComplexGeometry::SquarePlanar),
            )
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidReaction { .. }));
    }

    #[test]
    fn complexation_competes_with_precipitation() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:water",
                0,
                18.0,
                1_000.0,
                373.0,
                75.0,
                40_000.0,
            ))
            .substance(Substance::new(
                "destroy:silver_ion",
                1,
                107.8682,
                1_000.0,
                f64::MAX,
                100.0,
                20_000.0,
            ))
            .substance(Substance::new(
                "destroy:chloride",
                -1,
                35.45,
                1_000.0,
                f64::MAX,
                100.0,
                20_000.0,
            ))
            .substance(
                Substance::new("destroy:ammonia", 0, 17.031, 1_000.0, 400.0, 80.0, 23_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_unlimited()),
            )
            .substance(
                Substance::new(
                    "destroy:silver_chloride",
                    0,
                    143.3182,
                    1_000.0,
                    728.0,
                    100.0,
                    20_000.0,
                )
                .with_phase_properties(SubstancePhaseProperties {
                    can_precipitate: true,
                    can_form_liquid_phase: false,
                    aqueous_solubility_mol_per_bucket: Some(0.0),
                    organic_solubility_mol_per_bucket: Some(0.0),
                    ..SubstancePhaseProperties::aqueous_unlimited()
                }),
            )
            .precipitation(PrecipitationSpec::new(
                "destroy:silver_chloride",
                "destroy:silver_chloride",
                [
                    (SubstanceId::from("destroy:silver_ion"), 1),
                    (SubstanceId::from("destroy:chloride"), 1),
                ],
                1.0e-4,
            ))
            .complex_spec(
                ComplexSpec::new(
                    "destroy:silver_diammine",
                    "destroy:silver_ion",
                    [ComplexLigand::new("destroy:ammonia", 2)],
                    1,
                    1.0e4,
                )
                .with_coordination_number(2)
                .with_geometry(ComplexGeometry::Linear)
                .with_ligand_exchange_lability(LigandExchangeLability::Labile),
            )
            .build()
            .unwrap();

        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:silver_ion", 0.02)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:chloride", 0.02)
            .unwrap();
        mixture.equilibrate_solution(&registry).unwrap();
        let precipitated_before =
            mixture.concentration_in_phase(&"destroy:silver_chloride".into(), MixturePhase::Solid);
        assert!(
            precipitated_before > 0.005,
            "precipitated before ammonia was {precipitated_before}, silver aq {}, chloride aq {}, silver total {}, chloride total {}",
            mixture.concentration_in_phase(&"destroy:silver_ion".into(), MixturePhase::Aqueous),
            mixture.concentration_in_phase(&"destroy:chloride".into(), MixturePhase::Aqueous),
            mixture.concentration_of(&"destroy:silver_ion".into()),
            mixture.concentration_of(&"destroy:chloride".into())
        );

        mixture
            .add_substance(&registry, "destroy:ammonia", 0.20)
            .unwrap();
        mixture.equilibrate_solution(&registry).unwrap();

        let complex = mixture.concentration_of(&"destroy:silver_diammine".into());
        let precipitated_after =
            mixture.concentration_in_phase(&"destroy:silver_chloride".into(), MixturePhase::Solid);
        assert!(
            complex > 0.001,
            "silver diammine concentration was {complex}"
        );
        assert!(
            precipitated_after < precipitated_before,
            "precipitate before {precipitated_before}, after {precipitated_after}"
        );
    }

    #[test]
    fn calculated_ph_controls_reaction_conditions_during_simulation() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1_000.0, 373.0, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(Substance::new(
                "destroy:proton",
                1,
                1.0,
                1_000.0,
                f64::MAX,
                100.0,
                20_000.0,
            ))
            .substance(test_substance("organic:acid_labile_start"))
            .substance(test_substance("organic:acid_labile_product"))
            .reaction(
                Reaction::builder("organic:acid_labile_conversion")
                    .reactant("organic:acid_labile_start", 1, 1)
                    .product("organic:acid_labile_product", 1)
                    .condition(
                        ReactionCondition::new("acid-labile conversion requires acidic solution")
                            .acidity(AcidityCondition::Acidic),
                    )
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .build()
            .unwrap();

        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        mixture
            .add_substance(&registry, "organic:acid_labile_start", 0.1)
            .unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();
        assert_eq!(
            mixture.concentration_of(&"organic:acid_labile_product".into()),
            0.0
        );

        mixture
            .add_substance(&registry, "destroy:proton", 0.01)
            .unwrap();
        react_for_tick(&registry, &mut mixture, 1).unwrap();
        assert!(
            mixture.concentration_of(&"organic:acid_labile_product".into()) > 0.0,
            "acidic pH calculated from the mixture must unblock acid-sensitive reactions"
        );
    }

    #[test]
    fn gas_pressure_from_mixture_controls_pressure_sensitive_reactions() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1_000.0, 373.0, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(test_substance("organic:pressure_start"))
            .substance(test_substance("organic:pressure_product"))
            .substance(
                Substance::new("destroy:hydrogen", 0, 2.0, 70.0, 20.0, 28.0, 900.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: super::substance::LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: Some(0.0),
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: false,
                        can_form_liquid_phase: false,
                        solvent_role: super::substance::SolventRole::NotSolvent,
                    }),
            )
            .reaction(
                Reaction::builder("organic:pressure_sensitive_conversion")
                    .reactant("organic:pressure_start", 1, 1)
                    .product("organic:pressure_product", 1)
                    .condition(
                        ReactionCondition::new("hydrogen pressure must be present")
                            .gas_pressure_atm(0.5),
                    )
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .build()
            .unwrap();

        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "organic:pressure_start", 0.1)
            .unwrap();
        react_for_tick(&registry, &mut mixture, 1).unwrap();
        assert_eq!(
            mixture.concentration_of(&"organic:pressure_product".into()),
            0.0
        );

        mixture
            .exchange_gases_with_atmosphere(
                &registry,
                &[(SubstanceId::from("destroy:hydrogen"), 1.0)],
                super::mixture::STANDARD_PRESSURE_PASCAL,
                10.0,
                1.0,
            )
            .unwrap();
        react_for_tick(&registry, &mut mixture, 1).unwrap();
        assert!(
            mixture.concentration_of(&"organic:pressure_product".into()) > 0.0,
            "gas pressure calculated from the gas phase must unblock pressure-sensitive reactions"
        );
    }

    #[test]
    fn redox_environment_from_inorganic_species_changes_organic_selectivity() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1_000.0, 373.0, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                neutral_substance("inorganic:oxidant", 50.0)
                    .with_redox_roles(vec![RedoxRole::Oxidant]),
            )
            .build()
            .unwrap();
        let profile = SelectivityProfile::new(
            ReactionType::CarbonylReduction,
            SiteDescriptorBuilder::aldehyde(),
        )
        .with_nucleophile_strength(NucleophileStrength::Strong)
        .never_suppress();

        let mut neutral = Mixture::new(298.0).unwrap();
        neutral
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        let neutral_context =
            SelectivityContext::from_mixture(&registry, &neutral, &ReactionContext::default())
                .unwrap();
        let neutral_effect = SelectivityEngine::evaluate_profile(&profile, &neutral_context);

        let mut oxidizing = neutral.clone();
        oxidizing
            .add_substance(&registry, "inorganic:oxidant", 0.5)
            .unwrap();
        let oxidizing_context =
            SelectivityContext::from_mixture(&registry, &oxidizing, &ReactionContext::default())
                .unwrap();
        let oxidizing_effect = SelectivityEngine::evaluate_profile(&profile, &oxidizing_context);

        assert!(oxidizing_context.is_oxidizing());
        assert!(
            oxidizing_effect.rate_multiplier < neutral_effect.rate_multiplier,
            "oxidizing inorganic species must slow reduction-sensitive organic profiles"
        );
        assert!(
            oxidizing_effect.activation_delta_kj_per_mol
                > neutral_effect.activation_delta_kj_per_mol
        );
    }

    #[test]
    fn surface_catalyst_is_required_and_limited_by_free_sites() {
        let registry = surface_test_registry(
            Reaction::builder("surface:isomerization")
                .reactant("a", 1, 1)
                .product("b", 1)
                .surface_requirement("surface:nickel", 1.0)
                .surface_adsorption("surface:nickel", "adsorbed", 1.0)
                .pre_exponential_factor(1.0e12)
                .activation_energy_kj_per_mol(0.0)
                .build(),
        );
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "a", 1.0).unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();
        assert_eq!(mixture.concentration_of(&"b".into()), 0.0);

        let mut context = ReactionContext::default();
        context.add_surface("surface:nickel", 0.25).unwrap();
        react_for_tick_with_context(&registry, &mut mixture, &mut context, 1).unwrap();

        assert!(mixture.concentration_of(&"b".into()) > 0.0);
        assert!(mixture.concentration_of(&"b".into()) <= 0.25 + 1.0e-9);
        assert!(context.occupied_sites(&CatalystSurfaceId::from("surface:nickel")) > 0.0);
    }

    #[test]
    fn poisoning_and_recovery_change_surface_state_explicitly() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .reaction(
                Reaction::builder("surface:poison")
                    .reactant("a", 1, 1)
                    .product("b", 1)
                    .surface_requirement("surface:nickel", 1.0)
                    .surface_poisoning("surface:nickel", "adsorbed", 1.0)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .reaction(
                Reaction::builder("surface:recover")
                    .reactant("b", 1, 1)
                    .product("a", 1)
                    .surface_recovery("surface:nickel", 1.0)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .build(),
            )
            .catalyst_surface_spec(CatalystSurfaceSpec::chemical("surface:nickel", 58.69, 0))
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "a", 0.5).unwrap();
        let mut context = ReactionContext::default();
        context.add_surface("surface:nickel", 0.25).unwrap();

        react_for_tick_with_context(&registry, &mut mixture, &mut context, 1).unwrap();
        let surface_id = CatalystSurfaceId::from("surface:nickel");
        assert!(context.poisoned_sites(&surface_id) > 0.0);
        let poisoned = context.poisoned_sites(&surface_id);

        react_for_tick_with_context(&registry, &mut mixture, &mut context, 1).unwrap();
        assert!(context.poisoned_sites(&surface_id) < poisoned);
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

        assert_eq!(DESTROY_EXPLICIT_REACTION_COUNT, 118);
        assert_eq!(
            DESTROY_REGISTERED_REACTION_COUNT + DESTROY_METALLURGY_REACTION_COUNT,
            registry.reactions().count()
        );
        assert_eq!(DESTROY_REGISTERED_REACTION_COUNT, 155);
        assert_eq!(DESTROY_METALLURGY_REACTION_COUNT, 46);
    }

    #[test]
    fn carbonate_calcination_requires_high_temperature() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let carbonate = SubstanceId::from("destroy:calcium_carbonate");
        let oxide = SubstanceId::from("destroy:calcium_oxide");

        let mut cold = Mixture::new(900.0).unwrap();
        cold.add_substance(&registry, carbonate.clone(), 0.1)
            .unwrap();
        assert!(!react_for_tick(&registry, &mut cold, 1).unwrap());
        assert_eq!(cold.concentration_of(&oxide), 0.0);

        let mut hot = Mixture::new(1_250.0).unwrap();
        hot.add_substance(&registry, carbonate, 0.1).unwrap();
        assert!(react_for_tick(&registry, &mut hot, 1).unwrap());
        assert!(hot.concentration_of(&oxide) > 0.0);
        assert!(hot.concentration_of(&SubstanceId::from("destroy:carbon_dioxide")) > 0.0);
    }

    #[test]
    fn carbon_monoxide_reduction_produces_metal_and_gas() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let oxide = SubstanceId::from("destroy:copper_ii_oxide");
        let metal = SubstanceId::from("destroy:copper_metal");
        let carbon_monoxide = SubstanceId::from("destroy:carbon_monoxide");
        let carbon_dioxide = SubstanceId::from("destroy:carbon_dioxide");
        let reaction = registry
            .reaction(&"destroy:copper_ii_oxide.carbon_monoxide_reduction".into())
            .unwrap();
        let redox = reaction
            .redox
            .as_ref()
            .expect("metallurgical oxide reduction must be redox-annotated");
        assert_eq!(redox.transferred_electrons, 2);
        assert_eq!(
            redox.oxidation_half_id.as_deref(),
            Some("destroy:carbon_monoxide_to_carbon_dioxide_in_molten_oxide")
        );
        assert_eq!(
            redox.reduction_half_id.as_deref(),
            Some("destroy:copper_ii_oxide.molten_oxide_reduction_half")
        );
        assert!(registry
            .redox_half_reaction("destroy:copper_ii_oxide.molten_oxide_reduction_half")
            .is_some());

        let mut mixture = Mixture::new(1_450.0).unwrap();
        mixture
            .add_substance(&registry, oxide.clone(), 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, carbon_monoxide.clone(), 0.1)
            .unwrap();

        assert!(react_for_tick(&registry, &mut mixture, 1).unwrap());
        assert!(mixture.concentration_of(&oxide) < 0.1);
        assert!(mixture.concentration_of(&carbon_monoxide) < 0.1);
        assert!(mixture.concentration_in_phase(&metal, MixturePhase::MoltenMetal) > 0.0);
        assert!(mixture.concentration_in_phase(&carbon_dioxide, MixturePhase::Gas) > 0.0);
    }

    #[test]
    fn molten_oxide_electrolysis_produces_metal_and_oxygen_without_free_oxide_seed() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let oxide = SubstanceId::from("destroy:copper_ii_oxide");
        let metal = SubstanceId::from("destroy:copper_metal");
        let oxygen = SubstanceId::from("destroy:oxygen");
        let oxide_ion = SubstanceId::from("destroy:oxide");
        let mut mixture = Mixture::new(1_800.0).unwrap();
        mixture
            .add_substance(&registry, oxide.clone(), 0.01)
            .unwrap();

        let cell = ElectrolysisCell::new(
            ElectrodeProcess::anode("destroy:oxide_to_oxygen_in_molten_slag")
                .in_phase(MixturePhase::MoltenSlag),
            ElectrodeProcess::cathode(
                "destroy:copper_ii_oxide.molten_oxide_electrolysis_reduction_half",
            )
            .in_phase(MixturePhase::MoltenSlag),
            10.0,
        );
        let report = apply_electrolysis_cell(
            &registry,
            &mut mixture,
            &cell,
            super::redox::FARADAY_CONSTANT_COULOMBS_PER_MOL * 0.002,
            1.0,
        )
        .unwrap();

        assert!((report.transferred_electrons_mol_per_bucket - 0.002).abs() < 1.0e-12);
        assert!(mixture.concentration_in_phase(&metal, MixturePhase::MoltenMetal) > 0.0);
        assert!(mixture.concentration_in_phase(&oxygen, MixturePhase::Gas) > 0.0);
        assert!(mixture.concentration_of(&oxide) < 0.01);
        assert!(
            mixture.concentration_of(&oxide_ion)
                <= super::mixture::TRACE_CONCENTRATION_MOL_PER_BUCKET
        );
    }

    #[test]
    fn molten_metal_oxidation_moves_metal_into_slag() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let copper = SubstanceId::from("destroy:copper_metal");
        let oxygen = SubstanceId::from("destroy:oxygen");
        let copper_i_oxide = SubstanceId::from("destroy:copper_i_oxide");
        let copper_oxide = SubstanceId::from("destroy:copper_ii_oxide");
        let mut mixture = Mixture::new(1_500.0).unwrap();
        mixture
            .add_substance(&registry, copper.clone(), 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, oxygen.clone(), 0.05)
            .unwrap();

        assert!(react_for_tick(&registry, &mut mixture, 1).unwrap());
        assert!(mixture.concentration_in_phase(&copper, MixturePhase::MoltenMetal) < 1.0);
        assert!(mixture.concentration_in_phase(&oxygen, MixturePhase::Gas) < 0.05);
        assert!(
            mixture.concentration_of(&copper_i_oxide) + mixture.concentration_of(&copper_oxide)
                > 0.0
        );
    }

    #[test]
    fn carbon_monoxide_carburization_adds_dissolved_carbon_to_metal_phase() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let iron = SubstanceId::from("destroy:iron_metal");
        let carbon_monoxide = SubstanceId::from("destroy:carbon_monoxide");
        let dissolved_carbon = SubstanceId::from("destroy:dissolved_carbon");
        let carbon_dioxide = SubstanceId::from("destroy:carbon_dioxide");
        let mut mixture = Mixture::new(1_900.0).unwrap();
        mixture.add_substance(&registry, iron, 1.0).unwrap();
        mixture
            .add_substance(&registry, carbon_monoxide.clone(), 0.4)
            .unwrap();

        assert!(react_for_tick(&registry, &mut mixture, 1).unwrap());
        assert!(mixture.concentration_in_phase(&dissolved_carbon, MixturePhase::MoltenMetal) > 0.0);
        assert!(mixture.concentration_in_phase(&carbon_dioxide, MixturePhase::Gas) > 0.0);
        assert!(mixture.concentration_in_phase(&carbon_monoxide, MixturePhase::Gas) < 0.4);
    }

    #[test]
    fn aluminum_deoxidation_moves_dissolved_oxygen_to_slag() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let iron = SubstanceId::from("destroy:iron_metal");
        let aluminum = SubstanceId::from("destroy:aluminum_metal");
        let oxygen = SubstanceId::from("destroy:dissolved_oxygen");
        let alumina = SubstanceId::from("destroy:aluminum_oxide");
        let mut mixture = Mixture::new(1_900.0).unwrap();
        mixture.add_substance(&registry, iron, 1.0).unwrap();
        mixture
            .add_substance(&registry, aluminum.clone(), 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, oxygen.clone(), 0.2)
            .unwrap();

        assert!(react_for_tick(&registry, &mut mixture, 1).unwrap());
        assert!(mixture.concentration_in_phase(&oxygen, MixturePhase::MoltenMetal) < 0.2);
        assert!(mixture.concentration_in_phase(&aluminum, MixturePhase::MoltenMetal) < 0.1);
        assert!(mixture.concentration_in_phase(&alumina, MixturePhase::MoltenSlag) > 0.0);
    }

    #[test]
    fn basic_slag_removes_sulfur_and_phosphorus_from_molten_metal() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        let iron = SubstanceId::from("destroy:iron_metal");
        let lime = SubstanceId::from("destroy:calcium_oxide");
        let sulfur = SubstanceId::from("destroy:dissolved_sulfur");
        let phosphorus = SubstanceId::from("destroy:dissolved_phosphorus");
        let oxygen = SubstanceId::from("destroy:dissolved_oxygen");
        let calcium_sulfide = SubstanceId::from("destroy:calcium_sulfide");
        let calcium_phosphate = SubstanceId::from("destroy:calcium_phosphate");
        let mut mixture = Mixture::new(1_900.0).unwrap();
        mixture.add_substance(&registry, iron, 1.0).unwrap();
        mixture.add_substance(&registry, lime.clone(), 1.0).unwrap();
        mixture
            .add_substance(&registry, sulfur.clone(), 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, phosphorus.clone(), 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, oxygen.clone(), 0.5)
            .unwrap();

        assert!(react_for_tick(&registry, &mut mixture, 1).unwrap());
        assert!(mixture.concentration_in_phase(&sulfur, MixturePhase::MoltenMetal) < 0.1);
        assert!(mixture.concentration_in_phase(&phosphorus, MixturePhase::MoltenMetal) < 0.1);
        assert!(mixture.concentration_in_phase(&lime, MixturePhase::MoltenSlag) < 1.0);
        assert!(mixture.concentration_of(&calcium_sulfide) > 0.0);
        assert!(mixture.concentration_of(&calcium_phosphate) > 0.0);
    }

    #[test]
    fn destroy_reverse_reactions_are_registered_as_real_pairs() {
        let registry = destroy_registry_builder().unwrap().build().unwrap();
        for id in [
            "destroy:chlorine_solvation",
            "destroy:hydroxide_neutralization",
            "destroy:iodine_dissolution",
            "destroy:oleum_formation",
            "destroy:sodium_amalgamization",
            "destroy:sulfur_trioxide_hydration",
            "destroy:tetraborate_equilibrium",
        ] {
            let forward = registry.reaction(&id.into()).unwrap();
            let reverse_id = forward
                .reverse_reaction_id
                .as_ref()
                .expect("forward reaction must point at reverse reaction");
            let reverse = registry.reaction(reverse_id).unwrap();
            assert_eq!(reverse.reverse_reaction_id.as_ref(), Some(&forward.id));
            assert!(!reverse.show_in_jei);
        }

        assert!(registry
            .reaction(&"destroy:iron_iii_reduction".into())
            .is_err());
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

    #[test]
    fn reaction_thermodynamics_derives_equilibrium_constant_from_gibbs_energy() {
        let thermo =
            super::thermodynamics::ReactionThermodynamics::from_equilibrium_constant_at_kelvin(
                10.0, 298.15,
            )
            .unwrap();

        assert!(
            (thermo.gibbs_free_energy_change_kj_per_mol
                - super::thermodynamics::delta_g_from_equilibrium_constant(10.0, 298.15).unwrap())
            .abs()
                < 1.0e-12
        );
        assert!(
            (thermo.equilibrium_constant_at_kelvin(0.0, 298.15).unwrap() - 10.0).abs() < 1.0e-9
        );
    }

    #[test]
    fn reaction_equilibrium_constant_tracks_temperature_from_enthalpy_and_entropy() {
        let thermo =
            super::thermodynamics::ReactionThermodynamics::from_equilibrium_constant_at_kelvin(
                1.0, 298.15,
            )
            .unwrap();

        let cold = thermo.equilibrium_constant_at_kelvin(20.0, 280.0).unwrap();
        let hot = thermo.equilibrium_constant_at_kelvin(20.0, 350.0).unwrap();

        assert!(
            hot > cold,
            "endothermic equilibrium constant should increase with temperature: cold={cold}, hot={hot}"
        );
    }

    #[test]
    fn thermodynamic_rate_factor_stops_forward_reaction_at_equilibrium() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .reaction(
                Reaction::builder("a_to_b")
                    .reactant("a", 1, 1)
                    .product("b", 1)
                    .product_phase("b", MixturePhase::Organic)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .gibbs_free_energy_change_kj_per_mol(0.0)
                    .build(),
            )
            .build()
            .unwrap();
        let reaction = registry.reaction(&"a_to_b".into()).unwrap();

        let mut far_from_equilibrium = Mixture::new(298.15).unwrap();
        far_from_equilibrium
            .add_substance(&registry, "a", 1.0)
            .unwrap();
        let forward_rate = super::simulation::reaction_rate_mol_per_bucket_per_tick(
            &registry,
            &far_from_equilibrium,
            reaction,
        )
        .unwrap();

        let mut at_equilibrium = Mixture::new(298.15).unwrap();
        at_equilibrium.add_substance(&registry, "a", 1.0).unwrap();
        at_equilibrium.add_substance(&registry, "b", 4.0).unwrap();
        let equilibrium_rate = super::simulation::reaction_rate_mol_per_bucket_per_tick(
            &registry,
            &at_equilibrium,
            reaction,
        )
        .unwrap();

        assert!(forward_rate > 0.0);
        assert_eq!(equilibrium_rate, 0.0);
    }

    #[test]
    fn reverse_reactions_must_have_mirrored_gibbs_energy() {
        let error = ChemistryRegistryBuilder::new()
            .substance(test_substance("a"))
            .substance(test_substance("b"))
            .reaction(
                Reaction::builder("a_to_b")
                    .reactant("a", 1, 1)
                    .product("b", 1)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .gibbs_free_energy_change_kj_per_mol(-5.0)
                    .reverse_reaction_id("b_to_a")
                    .build(),
            )
            .reaction(
                Reaction::builder("b_to_a")
                    .reactant("b", 1, 1)
                    .product("a", 1)
                    .pre_exponential_factor(1.0e12)
                    .activation_energy_kj_per_mol(0.0)
                    .gibbs_free_energy_change_kj_per_mol(0.0)
                    .reverse_reaction_id("a_to_b")
                    .build(),
            )
            .build()
            .unwrap_err();

        assert!(matches!(
            error,
            ChemistryError::ReversibleThermodynamicsMismatch { .. }
        ));
    }
}
