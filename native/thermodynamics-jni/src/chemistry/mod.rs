#[path = "molecule/canonical.rs"]
pub mod canonical;
#[path = "data/catalog.rs"]
pub mod catalog;
#[path = "core/catalysis.rs"]
pub mod catalysis;
#[path = "core/complex.rs"]
pub mod complex;
#[path = "dynamic/mod.rs"]
pub mod dynamic;
pub mod error;
#[path = "molecule/frowns.rs"]
pub mod frowns;
#[path = "molecule/functional_group.rs"]
pub mod functional_group;
#[path = "core/mixture.rs"]
pub mod mixture;
#[path = "molecule/graph.rs"]
pub mod molecule;
#[path = "organic/mod.rs"]
pub mod organic;
#[path = "core/reaction.rs"]
pub mod reaction;
#[path = "data/reactions.rs"]
pub mod reactions;
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

pub use error::{ChemistryError, ChemistryResult};
pub use reactions::{DESTROY_EXPLICIT_REACTION_COUNT, DESTROY_REGISTERED_REACTION_COUNT};
pub use registry::{ChemistryRegistry, ChemistryRegistryBuilder};

pub fn destroy_registry_builder() -> ChemistryResult<ChemistryRegistryBuilder> {
    let builder = catalog::destroy_substances_registry_builder()?;
    reactions::destroy_reactions_registry_builder(builder)
}

pub fn destroy_registry_with_generated_reactions_builder(
) -> ChemistryResult<ChemistryRegistryBuilder> {
    organic::destroy_registry_with_generated_reactions_builder()
}

#[cfg(test)]
mod tests {
    use super::catalysis::{CatalystSurfaceId, CatalystSurfaceSpec};
    use super::complex::{ComplexLigand, ComplexSpec};
    use super::destroy_registry_builder;
    use super::error::ChemistryError;
    use super::mixture::{Mixture, MixturePhase};
    use super::reaction::Reaction;
    use super::registry::{ChemistryRegistryBuilder, ReactionCandidateScratch};
    use super::simulation::{
        react_for_tick, react_for_tick_with_context, react_until_equilibrium, ReactionContext,
    };
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
            DESTROY_REGISTERED_REACTION_COUNT,
            registry.reactions().count()
        );
        assert_eq!(DESTROY_REGISTERED_REACTION_COUNT, 155);
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
}
