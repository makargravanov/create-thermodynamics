use super::*;
use crate::chemistry::condition::{AcidityCondition, AtmosphereCondition};
use crate::chemistry::mixture::Mixture;
use crate::chemistry::molecule::{StereoDescriptor, Stereochemistry};
use crate::chemistry::organic::generators::expand_stereo_product_distribution;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::reactive_site::ReactiveSiteKind;
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::simulation::reaction_rate_mol_per_bucket_per_tick;
use crate::chemistry::substance::SubstanceId;
use crate::chemistry::DESTROY_REGISTERED_REACTION_COUNT;
use std::collections::BTreeSet;
use std::sync::OnceLock;

const ACTIVE_DESTROY_GENERIC_REACTIONS: &[&str] = &[
    "acyl_chloride_esterification",
    "acyl_chloride_formation",
    "acyl_chloride_hydrolysis",
    "alcohol_dehydration",
    "alcohol_oxidation",
    "aldehyde_oxidation",
    "alkene_chlorination",
    "alkene_chlorohydrination",
    "alkene_hydrolysis",
    "alkene_borane_hydroboration",
    "alkene_hydrochlorination",
    "alkene_hydrogenation",
    "alkene_hydroiodination",
    "alkene_iodination",
    "alkoxide_protonation",
    "alkyne_chlorination",
    "alkyne_chlorohydrination",
    "alkyne_hydrolysis",
    "alkyne_borane_hydroboration",
    "alkyne_hydrochlorination",
    "alkyne_hydrogenation",
    "alkyne_hydroiodination",
    "alkyne_iodination",
    "amide_hydrolysis",
    "amine_phosgenation",
    "borane_oxidation",
    "borate_ester_hydrolysis",
    "cyanamide_addition",
    "carboxylic_acid_esterification",
    "cyanide_nucleophilic_addition",
    "halide_amine_substitution",
    "halide_ammonia_substitution",
    "halide_cyanide_substitution",
    "halide_hydroxide_substitution",
    "isocyanate_hydrolysis",
    "nitrile_hydrogenation",
    "nitrile_hydrolysis",
    "nitro_hydrogenation",
    "thionyl_chloride_substitution",
    "wolff_kishner_reduction",
];

const EXCLUDED_DESTROY_GENERIC_REACTIONS: &[&str] = &[
    "electrophilic_hydroboration",
    "borate_esterification",
    "borohydride_carbonyl_reduction",
    "carboxylic_acid_reduction",
];

const ACTIVE_GENERATORS_WITHOUT_CATALOG_SUBSTRATE: &[&str] = &["aldehyde_oxidation"];
const ACTIVE_GENERATORS_WITH_UNKNOWN_STEREO_DISTRIBUTION: &[&str] = &[];

fn generated_registry() -> ChemistryRegistry {
    static REGISTRY: OnceLock<ChemistryRegistry> = OnceLock::new();
    REGISTRY
        .get_or_init(|| {
            destroy_registry_with_generated_reactions_builder()
                .unwrap()
                .build()
                .unwrap()
        })
        .clone()
}

fn reaction_with_prefix<'a>(registry: &'a ChemistryRegistry, prefix: &str) -> &'a Reaction {
    registry
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with(prefix))
        .unwrap_or_else(|| panic!("missing generated reaction with prefix {prefix}"))
}

#[test]
fn generation_space_indexes_only_substances_inside_scope() {
    let registry = super::super::destroy_registry_builder()
        .unwrap()
        .build()
        .unwrap();
    let substances = registry.substances().collect::<Vec<_>>();
    let scope = GenerationScope::from_substances(&BTreeSet::from([SubstanceId::from(
        "destroy:acetic_acid",
    )]));
    let space = OrganicGenerationSpace::new(substances.iter().copied(), &scope).unwrap();

    let acids = space
        .sites_of(&ReactiveSiteKind::CarboxylicAcid)
        .collect::<Vec<_>>();
    assert_eq!(acids.len(), 1);
    assert_eq!(acids[0].substance.id.as_str(), "destroy:acetic_acid");
    assert_eq!(space.sites_of(&ReactiveSiteKind::Alcohol).count(), 0);
}

#[test]
fn organic_generation_has_no_functional_group_transition_layer() {
    let source = include_str!("mod.rs");
    assert!(!source.contains(concat!("legacy", "_group", "_from", "_site")));
    assert!(!source.contains(concat!("sites", "_of", "_legacy", "_group")));
    assert!(!source.contains(concat!("Functional", "Group")));
}

#[test]
fn acetal_and_imine_generators_create_concrete_products_with_conditions() {
    let registry = generated_registry();
    let acetal = reaction_with_prefix(&registry, "acetal_formation/");
    assert!(acetal
        .conditions
        .iter()
        .any(|condition| condition.acidity == Some(AcidityCondition::Acidic)));
    assert!(acetal
        .products
        .iter()
        .chain(
            acetal
                .channels
                .iter()
                .flat_map(|channel| channel.products.iter())
        )
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    let imine = reaction_with_prefix(&registry, "imine_formation/");
    assert!(imine
        .conditions
        .iter()
        .any(|condition| condition.max_water_activity.is_some()));
    assert!(imine.products.len() >= 2);
}

#[test]
fn reactive_site_generators_add_aromatic_nitration_and_epoxide_hydrolysis() {
    let registry = generated_registry();
    let nitration = reaction_with_prefix(&registry, "aromatic_nitration/destroy_benzene/");
    assert!(nitration
        .conditions
        .iter()
        .any(|condition| condition.acidity == Some(AcidityCondition::Acidic)));
    assert!(!nitration.channels.is_empty() || !nitration.products.is_empty());

    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let epoxide = dynamic
        .resolve_frowns(
            "destroy:graph:atoms=C.C.O.H.H.H.H;bonds=0-s-1,0-s-2,1-s-2,0-s-3,0-s-4,1-s-5,1-s-6",
        )
        .unwrap();
    let report = dynamic.generate_reactions_for(&epoxide, 1).unwrap();
    assert!(report.added_reactions > 0);
    assert!(dynamic
        .reactions()
        .any(|reaction| reaction.id.as_str().starts_with("epoxide_hydrolysis/")));
}

#[test]
fn organometallic_and_aldol_generators_create_carbon_carbon_bonds() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let methyl_magnesium_chloride = dynamic.resolve_frowns("CMgCl").unwrap();
    dynamic
        .generate_reactions_for_substances(
            [
                SubstanceId::from("destroy:acetone"),
                methyl_magnesium_chloride,
            ],
            1,
        )
        .unwrap();
    let organometallic = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("organometallic_carbonyl_addition/")
        })
        .unwrap();
    assert!(organometallic
        .conditions
        .iter()
        .any(|condition| condition.atmosphere == Some(AtmosphereCondition::Inert)));
    assert!(!organometallic.external_products.is_empty());

    let acetaldehyde = dynamic.resolve_frowns("CC=O").unwrap();
    dynamic
        .generate_reactions_for_substances([SubstanceId::from("destroy:acetone"), acetaldehyde], 1)
        .unwrap();
    assert!(dynamic
        .reactions()
        .any(|reaction| reaction.id.as_str().starts_with("aldol_addition/")));
}

#[test]
fn scoped_generation_matches_full_static_generation() {
    let registry = super::super::destroy_registry_builder()
        .unwrap()
        .build()
        .unwrap();
    let full = generate_organic_reactions(&registry).unwrap();
    let substances = registry.substances().collect::<Vec<_>>();
    let all_ids = substances
        .iter()
        .map(|substance| substance.id.clone())
        .collect::<BTreeSet<_>>();
    let scoped =
        generate_organic_reactions_for_substances(&substances, &all_ids, &all_ids).unwrap();

    let full_substance_ids = full
        .substances
        .iter()
        .map(|substance| substance.id.as_str())
        .collect::<BTreeSet<_>>();
    let scoped_substance_ids = scoped
        .substances
        .iter()
        .map(|substance| substance.id.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(full_substance_ids, scoped_substance_ids);

    let full_reaction_ids = full
        .reactions
        .iter()
        .map(|reaction| reaction.id.as_str())
        .collect::<BTreeSet<_>>();
    let scoped_reaction_ids = scoped
        .reactions
        .iter()
        .map(|reaction| reaction.id.as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(full_reaction_ids, scoped_reaction_ids);
}

#[test]
fn generated_registry_builds_without_duplicate_derived_substances() {
    let registry = generated_registry();
    let mut canonical_codes = BTreeSet::new();
    for substance in registry.substances() {
        if !substance.id.as_str().starts_with("destroy:derived_") {
            continue;
        }
        if let Some(structure) = &substance.molecular_structure {
            assert!(canonical_codes.insert(structure.canonical_code().unwrap()));
        }
    }
    assert!(registry.reactions().count() > DESTROY_REGISTERED_REACTION_COUNT);
}

#[test]
fn active_destroy_generic_reactions_are_accounted_for() {
    assert_eq!(ACTIVE_DESTROY_GENERIC_REACTIONS.len(), 40);
    assert_eq!(EXCLUDED_DESTROY_GENERIC_REACTIONS.len(), 4);

    let registry = generated_registry();
    for prefix in ACTIVE_DESTROY_GENERIC_REACTIONS {
        if ACTIVE_GENERATORS_WITHOUT_CATALOG_SUBSTRATE.contains(prefix) {
            continue;
        }
        if ACTIVE_GENERATORS_WITH_UNKNOWN_STEREO_DISTRIBUTION.contains(prefix) {
            continue;
        }
        assert!(
            registry
                .reactions()
                .any(|reaction| reaction.id.as_str().starts_with(prefix)),
            "missing generated reaction for active Destroy generator {prefix}",
        );
    }
    for prefix in EXCLUDED_DESTROY_GENERIC_REACTIONS {
        assert!(
            !registry
                .reactions()
                .any(|reaction| reaction.id.as_str().starts_with(prefix)),
            "excluded Destroy generator {prefix} should not be registered",
        );
    }
}

#[test]
fn halide_hydroxide_substitution_generates_ethanol_from_chloroethane() {
    let registry = generated_registry();
    let reaction = registry
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("halide_hydroxide_substitution/destroy_chloroethane/")
        })
        .unwrap();
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:ethanol"));
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:chloride"));
}

#[test]
fn alcohol_oxidation_generates_acetone_from_isopropanol() {
    let registry = generated_registry();
    let reaction = registry
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("alcohol_oxidation/destroy_isopropanol/")
        })
        .unwrap();
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:acetone"));
}

#[test]
fn esterification_generates_product_from_acetic_acid_and_ethanol() {
    let registry = generated_registry();
    let reaction = registry
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("carboxylic_acid_esterification/destroy_acetic_acid/destroy_ethanol/")
        })
        .unwrap();
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() != "destroy:water"));
}

#[test]
fn nitrile_hydrolysis_and_nitro_hydrogenation_are_registered() {
    let registry = generated_registry();
    assert!(registry.reactions().any(|reaction| {
        reaction
            .id
            .as_str()
            .starts_with("nitrile_hydrolysis/destroy_generic_nitrile/")
    }));
    let nitro = registry
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("nitro_hydrogenation/destroy_dinitrotoluene/")
        })
        .unwrap();
    assert!(!nitro.external_catalysts.is_empty());
}

#[test]
fn acyl_chloride_generators_are_registered() {
    let registry = generated_registry();
    let hydrolysis = reaction_with_prefix(
        &registry,
        "acyl_chloride_hydrolysis/destroy_generic_acyl_chloride/",
    );
    assert!(hydrolysis
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrochloric_acid"));

    let formation = reaction_with_prefix(&registry, "acyl_chloride_formation/destroy_acetic_acid/");
    assert!(formation
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:carbon_dioxide"));

    let esterification = reaction_with_prefix(
        &registry,
        "acyl_chloride_esterification/destroy_generic_acyl_chloride/destroy_ethanol/",
    );
    assert!(esterification
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrochloric_acid"));
}

#[test]
fn halide_substitution_generators_are_registered() {
    let registry = generated_registry();
    let ammonia = reaction_with_prefix(
        &registry,
        "halide_ammonia_substitution/destroy_chloroethane/",
    );
    assert!(ammonia
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:ammonium"));

    let cyanide = reaction_with_prefix(
        &registry,
        "halide_cyanide_substitution/destroy_chloroethane/",
    );
    assert!(cyanide
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:cyanide"));

    let amine = reaction_with_prefix(
        &registry,
        "halide_amine_substitution/destroy_chloroethane/destroy_methylamine/",
    );
    assert!(amine
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:proton"));
}

#[test]
fn electrophilic_addition_generators_are_registered() {
    let registry = generated_registry();
    for prefix in [
        "alkene_chlorination/destroy_ethene/",
        "alkene_chlorohydrination/destroy_ethene/",
        "alkene_hydrolysis/destroy_ethene/",
        "alkene_borane_hydroboration/destroy_ethene/",
        "alkene_hydrochlorination/destroy_ethene/",
        "alkene_hydrogenation/destroy_ethene/",
        "alkene_hydroiodination/destroy_ethene/",
        "alkene_iodination/destroy_ethene/",
        "alkyne_hydrogenation/destroy_acetylene/",
    ] {
        reaction_with_prefix(&registry, prefix);
    }
    let hydrogenation = reaction_with_prefix(&registry, "alkene_hydrogenation/destroy_ethene/");
    assert!(!hydrogenation.external_catalysts.is_empty());
}

#[test]
fn geometric_stereo_products_use_kinetic_channels() {
    let structure = crate::chemistry::frowns::parse_frowns(
        "destroy:graph:atoms=C.C.H.Cl.H.I;bonds=0-2-1,0-s-2,0-s-3,1-s-4,1-s-5;stereo=mix:db:0.1.3.5",
    )
    .unwrap();
    let variants = expand_stereo_product_distribution(structure).unwrap();
    let e_variant = variants
        .iter()
        .find(|variant| variant.channel_suffix.contains("db_e"))
        .unwrap();
    let z_variant = variants
        .iter()
        .find(|variant| variant.channel_suffix.contains("db_z"))
        .unwrap();

    assert!(z_variant.activation_delta_kj_per_mol > e_variant.activation_delta_kj_per_mol);
    assert!(
        z_variant.pre_exponential_factor_multiplier < e_variant.pre_exponential_factor_multiplier
    );
    assert!(e_variant
        .structure
        .stereochemistry
        .iter()
        .any(|stereo| matches!(stereo, Stereochemistry::DoubleBond(double_bond) if double_bond.descriptor == StereoDescriptor::E)));
    assert!(z_variant
        .structure
        .stereochemistry
        .iter()
        .any(|stereo| matches!(stereo, Stereochemistry::DoubleBond(double_bond) if double_bond.descriptor == StereoDescriptor::Z)));
}

#[test]
fn heteroatom_generators_are_registered() {
    let registry = generated_registry();
    reaction_with_prefix(&registry, "amide_hydrolysis/destroy_acetamide/");
    reaction_with_prefix(&registry, "amine_phosgenation/destroy_methylamine/");
    reaction_with_prefix(&registry, "cyanamide_addition/destroy_methylamine/");
    reaction_with_prefix(
        &registry,
        "isocyanate_hydrolysis/destroy_generic_isocyanate/",
    );
    reaction_with_prefix(&registry, "borane_oxidation/destroy_generic_borane/");
    reaction_with_prefix(
        &registry,
        "borate_ester_hydrolysis/destroy_generic_borate_ester/",
    );
    reaction_with_prefix(&registry, "nitrile_hydrogenation/destroy_generic_nitrile/");
    reaction_with_prefix(&registry, "thionyl_chloride_substitution/destroy_ethanol/");
    reaction_with_prefix(&registry, "wolff_kishner_reduction/destroy_acetone/");
}

#[test]
fn selectivity_engine_integration_keeps_reactions_but_suppresses_runtime_rate() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();

    let acetic_acid = SubstanceId::from("destroy:acetic_acid");
    let tert_butanol = dynamic.resolve_frowns("CC(C)(C)O").unwrap(); // 3° alcohol
    let ethanol = SubstanceId::from("destroy:ethanol"); // 1° alcohol
    let tert_butyl_chloride = dynamic.resolve_frowns("CC(C)(C)Cl").unwrap(); // 3° halide
    let chloroethane = SubstanceId::from("destroy:chloroethane"); // 1° halide

    dynamic
        .generate_reactions_for_substances(
            [
                acetic_acid.clone(),
                tert_butanol.clone(),
                ethanol.clone(),
                tert_butyl_chloride.clone(),
                chloroethane.clone(),
            ],
            1,
        )
        .unwrap();

    let registry = dynamic.to_registry().unwrap();
    let tert_esterification = reaction_for_reactants(
        &registry,
        "carboxylic_acid_esterification",
        &[acetic_acid.clone(), tert_butanol.clone()],
    );
    let eth_esterification = reaction_for_reactants(
        &registry,
        "carboxylic_acid_esterification",
        &[acetic_acid, ethanol],
    );
    assert!(
        tert_esterification.is_some(),
        "tertiary esterification remains a structural candidate"
    );
    assert!(
        eth_esterification.is_some(),
        "primary esterification remains a structural candidate"
    );

    let tert_substitution = reaction_for_reactants(
        &registry,
        "halide_hydroxide_substitution",
        &[tert_butyl_chloride.clone(), "destroy:hydroxide".into()],
    )
    .expect("tertiary halide substitution remains a structural candidate");
    let eth_substitution = reaction_for_reactants(
        &registry,
        "halide_hydroxide_substitution",
        &[chloroethane.clone(), "destroy:hydroxide".into()],
    )
    .expect("primary halide substitution remains a structural candidate");

    let mut tert_mixture = Mixture::new(298.15).unwrap();
    tert_mixture
        .add_substance(&registry, "destroy:water", 1.0)
        .unwrap();
    tert_mixture
        .add_substance(&registry, tert_butyl_chloride.clone(), 0.01)
        .unwrap();
    tert_mixture
        .set_gaseous_fraction(&registry, tert_butyl_chloride, 0.0)
        .unwrap();
    tert_mixture
        .add_substance(&registry, "destroy:hydroxide", 0.1)
        .unwrap();
    let tert_rate =
        reaction_rate_mol_per_bucket_per_tick(&registry, &tert_mixture, tert_substitution).unwrap();

    let mut eth_mixture = Mixture::new(298.15).unwrap();
    eth_mixture
        .add_substance(&registry, "destroy:water", 1.0)
        .unwrap();
    eth_mixture
        .add_substance(&registry, chloroethane.clone(), 0.01)
        .unwrap();
    eth_mixture
        .set_gaseous_fraction(&registry, chloroethane, 0.0)
        .unwrap();
    eth_mixture
        .add_substance(&registry, "destroy:hydroxide", 0.1)
        .unwrap();
    let eth_rate =
        reaction_rate_mol_per_bucket_per_tick(&registry, &eth_mixture, eth_substitution).unwrap();

    assert_eq!(tert_rate, 0.0, "tertiary SN2 must be runtime-suppressed");
    assert!(
        eth_rate > tert_rate,
        "primary SN2 should be faster than tertiary SN2"
    );
}

fn reaction_for_reactants<'a>(
    registry: &'a ChemistryRegistry,
    prefix: &str,
    reactants: &[SubstanceId],
) -> Option<&'a Reaction> {
    registry.reactions().find(|reaction| {
        reaction.id.as_str().starts_with(prefix)
            && reactants.iter().all(|id| {
                reaction
                    .reactants
                    .iter()
                    .any(|term| &term.substance_id == id)
            })
    })
}

#[test]
fn test_aromatic_eas_directing_groups_and_deactivation() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();

    // 1. Unsubstituted Benzene (already in catalog as destroy:benzene)
    let benzene = SubstanceId::from("destroy:benzene");
    dynamic
        .generate_reactions_for_substances([benzene.clone()], 1)
        .unwrap();

    // Verify all four single-site EAS reactions exist for benzene
    let nitration_exists = dynamic
        .reactions()
        .any(|r| r.id.as_str().contains("aromatic_nitration") && r.id.as_str().contains("benzene"));
    let chlorination_exists = dynamic.reactions().any(|r| {
        r.id.as_str().contains("aromatic_chlorination") && r.id.as_str().contains("benzene")
    });
    let bromination_exists = dynamic.reactions().any(|r| {
        r.id.as_str().contains("aromatic_bromination") && r.id.as_str().contains("benzene")
    });
    let sulfonation_exists = dynamic.reactions().any(|r| {
        r.id.as_str().contains("aromatic_sulfonation") && r.id.as_str().contains("benzene")
    });

    assert!(nitration_exists, "Benzene must undergo nitration");
    assert!(chlorination_exists, "Benzene must undergo chlorination");
    assert!(bromination_exists, "Benzene must undergo bromination");
    assert!(sulfonation_exists, "Benzene must undergo sulfonation");

    // 2. Friedel-Crafts on Benzene
    let chloroethane = SubstanceId::from("destroy:chloroethane");
    let acetyl_chloride = SubstanceId::from("destroy:acetyl_chloride");
    dynamic
        .generate_reactions_for_substances(
            [
                benzene.clone(),
                chloroethane.clone(),
                acetyl_chloride.clone(),
            ],
            1,
        )
        .unwrap();

    let fc_alkylation_exists = dynamic
        .reactions()
        .any(|r| r.id.as_str().contains("fc_alkylation") && r.id.as_str().contains("benzene"));
    let fc_acylation_exists = dynamic
        .reactions()
        .any(|r| r.id.as_str().contains("fc_acylation") && r.id.as_str().contains("benzene"));
    assert!(
        fc_alkylation_exists,
        "Benzene must undergo Friedel-Crafts alkylation"
    );
    assert!(
        fc_acylation_exists,
        "Benzene must undergo Friedel-Crafts acylation"
    );

    // 3. Toluene (Weakly Activating Ortho/Para-director)
    // Graph FROWNS: 7 carbons (ring 0-5 + methyl 6), 8 hydrogens (methyl 7-9 + ring 10-14)
    let toluene = dynamic.resolve_frowns(
        "destroy:graph:atoms=C.C.C.C.C.C.C.H.H.H.H.H.H.H.H;bonds=0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0,0-s-6,6-s-7,6-s-8,6-s-9,1-s-10,2-s-11,3-s-12,4-s-13,5-s-14"
    ).unwrap();
    dynamic
        .generate_reactions_for_substances([toluene.clone()], 1)
        .unwrap();

    // Get the toluene nitration reaction and inspect the channel activation energies
    let toluene_nitration = dynamic
        .reactions()
        .find(|r| r.id.as_str().contains("aromatic_nitration") && r.id.as_str().contains("toluene"))
        .expect("Toluene must undergo nitration");

    // There should be multiple channels (ortho, meta, para)
    assert!(
        toluene_nitration.channels.len() > 1,
        "Toluene nitration must have multiple regioselective channels"
    );

    let channels = &toluene_nitration.channels;
    let mut energies: Vec<f64> = channels
        .iter()
        .map(|c| c.activation_gibbs_kj_per_mol)
        .collect::<Vec<_>>();
    energies.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Lowest Ea should be 30.0 - 5.0 = 25.0 kcal/mol (para)
    // Ortho positions should be 30.0 - 4.0 = 26.0 kcal/mol
    // Meta positions should be 30.0 - 1.0 = 29.0 kcal/mol
    assert!(
        (energies[0] - 25.0).abs() < 1e-3,
        "Para-nitration Ea should be 25.0 kcal/mol, got {}",
        energies[0]
    );
    assert!(
        (energies[1] - 26.0).abs() < 1e-3,
        "Ortho-nitration Ea should be 26.0 kcal/mol, got {}",
        energies[1]
    );
    assert!(
        (energies[energies.len() - 2] - 29.0).abs() < 1e-3,
        "Meta-nitration Ea should be 29.0 kcal/mol, got {}",
        energies[energies.len() - 2]
    );

    // 4. Nitrobenzene (Strongly Deactivating Meta-director)
    // Zwitterionic representation: C-N(=O)O- with N+ and O-
    let nitrobenzene = dynamic.resolve_frowns(
        "destroy:graph:atoms=C.C.C.C.C.C.N^1.O.O^-1.H.H.H.H.H;bonds=0-s-6,6-d-7,6-s-8,1-s-9,2-s-10,3-s-11,4-s-12,5-s-13,0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0"
    ).unwrap();
    dynamic
        .generate_reactions_for_substances([nitrobenzene.clone()], 1)
        .unwrap();

    // Find nitration reaction for nitrobenzene by its substance ID
    // Find nitrobenzene nitration by matching reactant substance ID instead of
    // prefix-matching the reaction ID (which has sanitized slashes from graph FROWNS).
    let nitro_nitration = dynamic
        .reactions()
        .find(|r| {
            r.id.as_str().starts_with("aromatic_nitration/")
                && r.reactants
                    .first()
                    .is_some_and(|term| term.substance_id == nitrobenzene)
        })
        .expect("Nitrobenzene must undergo nitration");

    // Class: StronglyDeactivating (-NO2). Ortho (+15.0), Meta (+9.0), Para (+16.0)
    // Therefore, Meta-nitration has the lowest activation energy (30.0 + 9.0 = 39.0 kcal/mol)
    let mut nitro_energies: Vec<f64> = nitro_nitration
        .channels
        .iter()
        .map(|c| c.activation_gibbs_kj_per_mol)
        .collect::<Vec<_>>();
    nitro_energies.sort_by(|a, b| a.partial_cmp(b).unwrap());

    assert!(
        (nitro_energies[0] - 39.0).abs() < 1e-3,
        "Meta-nitration of nitrobenzene must have Ea of 39.0 kcal/mol, got {}",
        nitro_energies[0]
    );
    assert!(
        nitro_energies[2] > 44.0,
        "Ortho/Para-nitration must be highly deactivated, got {}",
        nitro_energies[2]
    );

    // 5. Friedel-Crafts Deactivation Blocking Rule
    dynamic
        .generate_reactions_for_substances([nitrobenzene.clone(), chloroethane], 1)
        .unwrap();

    let fc_prefix = format!("fc_alkylation/{}", nitrobenzene.as_str());
    let nitro_fc_exists = dynamic
        .reactions()
        .any(|r| r.id.as_str().starts_with(&fc_prefix));
    assert!(
        !nitro_fc_exists,
        "Friedel-Crafts alkylation must be completely blocked on deactivated nitrobenzene"
    );
}

#[test]
fn toluene_can_be_nitrated_dynamically() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();

    let toluene = dynamic.resolve_frowns(
        "destroy:graph:atoms=C.C.C.C.C.C.C.H.H.H.H.H.H.H.H;bonds=0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0,0-s-6,6-s-7,6-s-8,6-s-9,1-s-10,2-s-11,3-s-12,4-s-13,5-s-14"
    ).unwrap();

    dynamic
        .generate_reactions_for_substances([toluene.clone()], 1)
        .unwrap();

    let nitration = dynamic
        .reactions()
        .find(|r| {
            r.id.as_str().starts_with("aromatic_nitration/")
                && r.reactants
                    .first()
                    .is_some_and(|t| t.substance_id == toluene)
        })
        .expect("Toluene must undergo nitration");

    assert!(
        nitration.channels.len() >= 2,
        "Toluene should have multiple regioselective nitration channels"
    );
    assert!(nitration
        .reactants
        .iter()
        .any(|t| t.substance_id == toluene));
    // Products are stored in channels when there are multiple regioselective positions
    assert!(
        nitration.channels.iter().any(|c| !c.products.is_empty()),
        "Nitration channels should contain products"
    );
}

#[test]
fn snar_on_nitrochlorobenzene() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();

    // o-Nitrochlorobenzene: NO₂ ortho to Cl (strongest SNAr activation)
    // Ring: C0(NO₂)-C1(Cl)-C2-C3-C4-C5
    // NO₂ represented as zwitterion N⁺—O⁻ and N⁺=O
    let o_nitrochlorobenzene = dynamic.resolve_frowns(
        "destroy:graph:atoms=C.C.C.C.C.C.N^1.O.O^-1.Cl.H.H.H.H;bonds=0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0,0-s-6,6-d-7,6-s-8,1-s-9,2-s-10,3-s-11,4-s-12,5-s-13"
    ).unwrap();

    dynamic
        .generate_reactions_for_substances([o_nitrochlorobenzene.clone()], 1)
        .unwrap();

    // SNAr with hydroxide should exist (NO₂ activates the ring)
    let snar_oh = dynamic
        .reactions()
        .find(|r| {
            r.id.as_str()
                .starts_with("aryl_halide_hydroxide_substitution/")
                && r.reactants
                    .first()
                    .is_some_and(|t| t.substance_id == o_nitrochlorobenzene)
        })
        .expect("Nitrochlorobenzene must undergo SNAr with hydroxide");

    assert!(
        snar_oh
            .reactants
            .iter()
            .any(|t| t.substance_id.as_str() == "destroy:hydroxide"),
        "Hydroxide must be a reactant in SNAr"
    );
    assert!(
        snar_oh.products.len() > 0,
        "SNAr must produce a substituted product"
    );

    // SNAr with ammonia should also exist
    let snar_nh3 = dynamic
        .reactions()
        .find(|r| {
            r.id.as_str()
                .starts_with("aryl_halide_ammonia_substitution/")
                && r.reactants
                    .first()
                    .is_some_and(|t| t.substance_id == o_nitrochlorobenzene)
        })
        .expect("Nitrochlorobenzene must undergo SNAr with ammonia");

    assert!(
        snar_nh3
            .reactants
            .iter()
            .any(|t| t.substance_id.as_str() == "destroy:ammonia"),
        "Ammonia must be a reactant in SNAr"
    );
}
