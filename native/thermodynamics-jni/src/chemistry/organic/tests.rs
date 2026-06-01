use super::*;
use crate::chemistry::condition::{AcidityCondition, AtmosphereCondition};
use crate::chemistry::mixture::Mixture;
use crate::chemistry::molecule::{StereoDescriptor, Stereochemistry};
use crate::chemistry::organic::generators::expand_stereo_product_distribution;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::reactive_site::{try_find_reactive_sites, ReactiveSiteKind};
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
    "borohydride_carbonyl_reduction",
    "cyanamide_addition",
    "carboxylic_acid_esterification",
    "cyanide_nucleophilic_addition",
    "halide_amine_substitution",
    "halide_ammonia_substitution",
    "halide_cyanide_substitution",
    "halide_hydroxide_substitution",
    "isocyanate_hydrolysis",
    "lah_ester_reduction",
    "nitrile_hydrogenation",
    "nitrile_hydrolysis",
    "nitro_hydrogenation",
    "thionyl_chloride_substitution",
    "wolff_kishner_reduction",
];

const EXCLUDED_DESTROY_GENERIC_REACTIONS: &[&str] = &[
    "electrophilic_hydroboration",
    "borate_esterification",
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
fn alpha_carbon_generators_create_halogenation_dehydration_enamine_and_alkylation() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();

    let acetone = SubstanceId::from("destroy:acetone");
    let iodomethane = SubstanceId::from("destroy:iodomethane");
    let methyl_acetate = SubstanceId::from("destroy:methyl_acetate");
    let dimethylamine = dynamic.resolve_frowns("CNC").unwrap();
    let methyl_vinyl_ketone = dynamic.resolve_frowns("C=CC(=O)C").unwrap();

    dynamic
        .generate_reactions_for_substances(
            [
                acetone.clone(),
                iodomethane.clone(),
                methyl_acetate.clone(),
                dimethylamine.clone(),
                methyl_vinyl_ketone.clone(),
            ],
            1,
        )
        .unwrap();

    let alpha_chlorination = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("alpha_chlorination/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == acetone)
        })
        .expect("acetone must have alpha chlorination");
    assert!(alpha_chlorination
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrochloric_acid"));

    let enolate_alkylation = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("enolate_alkylation/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == iodomethane)
        })
        .expect("acetone enolate must alkylate iodomethane");
    assert!(enolate_alkylation
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:iodide"));
    assert!(enolate_alkylation
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    let enamine = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("enamine_formation/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == dimethylamine)
        })
        .expect("acetone and secondary amine must form an enamine");
    assert!(enamine
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("michael_addition/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == methyl_vinyl_ketone)
        })
        .expect("acetone enolate must add to a conjugated enone");

    dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("claisen_condensation/")
                && reaction
                    .reactants
                    .iter()
                    .all(|term| term.substance_id == methyl_acetate)
        })
        .expect("methyl acetate must undergo self-Claisen condensation");

    let acetaldehyde = dynamic.resolve_frowns("CC=O").unwrap();
    dynamic
        .generate_reactions_for_substances([acetaldehyde.clone()], 1)
        .unwrap();
    let aldol_product = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("aldol_addition/")
                && reaction
                    .reactants
                    .iter()
                    .all(|term| term.substance_id == acetaldehyde)
        })
        .and_then(|reaction| reaction.products.first())
        .map(|term| term.substance_id.clone())
        .expect("acetaldehyde must produce an aldol addition product");

    dynamic
        .generate_reactions_for_substances([aldol_product.clone()], 1)
        .unwrap();
    assert!(dynamic.reactions().any(|reaction| {
        reaction.id.as_str().starts_with("aldol_dehydration/")
            && reaction
                .reactants
                .iter()
                .any(|term| term.substance_id == aldol_product)
    }));
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
    assert_eq!(ACTIVE_DESTROY_GENERIC_REACTIONS.len(), 42);
    assert_eq!(EXCLUDED_DESTROY_GENERIC_REACTIONS.len(), 3);

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
fn borohydride_reduces_carbonyls_to_alcohols_with_closed_boron_stoichiometry() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    dynamic
        .generate_reactions_for_substances([SubstanceId::from("destroy:acetone")], 1)
        .unwrap();
    let registry = dynamic.to_registry().unwrap();
    let reaction = reaction_for_reactants(
        &registry,
        "borohydride_carbonyl_reduction",
        &[SubstanceId::from("destroy:acetone")],
    )
    .expect("acetone must generate borohydride reduction");
    assert!(reaction
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:borohydride"));
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:isopropanol"));
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:tetrahydroxyborate"));
}

#[test]
fn lah_ester_reduction_splits_ester_into_two_alcohols() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let acid = SubstanceId::from("destroy:acetic_acid");
    let ethanol = SubstanceId::from("destroy:ethanol");
    dynamic
        .generate_reactions_for_substances([acid.clone(), ethanol.clone()], 1)
        .unwrap();
    let registry = dynamic.to_registry().unwrap();
    let esterification = reaction_for_reactants(
        &registry,
        "carboxylic_acid_esterification",
        &[acid, ethanol],
    )
    .expect("acetic acid and ethanol must generate ethyl acetate");
    let ester = esterification
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("esterification must have an ester product")
        .substance_id
        .clone();

    dynamic
        .generate_reactions_for_substances([ester.clone()], 1)
        .unwrap();
    let registry = dynamic.to_registry().unwrap();
    let reduction = reaction_for_reactants(&registry, "lah_ester_reduction", &[ester])
        .expect("ester must generate LAH reduction");
    assert!(reduction.external_reactants.iter().any(|external| {
        external.description == "lithium aluminium hydride hydride/proton equivalents"
    }));
    let ethanol_count = reduction
        .products
        .iter()
        .filter(|term| term.substance_id.as_str() == "destroy:ethanol")
        .map(|term| term.coefficient)
        .sum::<u32>();
    assert_eq!(ethanol_count, 2);
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
fn phosphorus_generators_are_registered() {
    let mut dynamic = super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog()
        .expect("destroy catalog must load");
    let phosphine = SubstanceId::from("destroy:trimethylphosphine");
    let chloroethane = SubstanceId::from("destroy:chloroethane");
    let ethoxide = SubstanceId::from("destroy:ethoxide");
    let acetone = SubstanceId::from("destroy:acetone");

    dynamic
        .generate_reactions_for_substances([phosphine.clone(), chloroethane.clone()], 1)
        .expect("phosphorus path must generate");

    let salt = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("phosphonium_salt_formation/")
        })
        .expect("phosphonium salt formation must be registered");
    assert!(salt
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:chloride"));
    let salt_product = salt
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:chloride")
        .expect("phosphonium salt must produce a phosphonium cation")
        .substance_id
        .clone();

    dynamic
        .generate_reactions_for_substances([salt_product.clone(), ethoxide.clone()], 1)
        .expect("ylide path must generate");

    let registry = dynamic
        .to_registry()
        .expect("dynamic registry must convert");
    let ylide = reaction_for_reactants(
        &registry,
        "phosphonium_ylide_formation",
        &[salt_product.clone(), ethoxide.clone()],
    )
    .expect("phosphonium ylide formation must be registered");
    assert!(ylide
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:ethanol"));
    let ylide_product = ylide
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:ethanol")
        .expect("phosphonium ylide must produce the ylide substance")
        .substance_id
        .clone();

    dynamic
        .generate_reactions_for_substances([ylide_product.clone(), acetone.clone()], 1)
        .expect("wittig path must generate");
    let registry = dynamic
        .to_registry()
        .expect("dynamic registry must convert");
    let wittig = reaction_for_reactants(
        &registry,
        "wittig_olefination",
        &[ylide_product, acetone.clone()],
    )
    .expect("wittig olefination must be registered");
    assert!(wittig.products.len() >= 2);
    assert!(wittig
        .products
        .iter()
        .any(|term| term.substance_id.as_str() != "destroy:ethanol"));
}

#[test]
fn hwe_and_julia_olefinations_use_concrete_anionic_reagents() {
    let mut dynamic = super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog()
        .expect("destroy catalog must load");
    let acetone = SubstanceId::from("destroy:acetone");
    let phosphonate = dynamic
        .resolve_frowns("O=P(OC)(OC)C^-1")
        .expect("phosphonate carbanion must parse");
    let sulfone = dynamic
        .resolve_frowns("CS(=O)(=O)C^-1")
        .expect("sulfone carbanion must parse");

    dynamic
        .generate_reactions_for_substances(
            [phosphonate.clone(), sulfone.clone(), acetone.clone()],
            1,
        )
        .expect("olefination generators must run");
    let registry = dynamic
        .to_registry()
        .expect("dynamic registry must convert");

    let hwe = reaction_for_reactants(
        &registry,
        "horner_wadsworth_emmons_olefination",
        &[phosphonate, acetone.clone()],
    )
    .expect("Horner-Wadsworth-Emmons olefination must be registered");
    assert!(hwe.channels.len() >= 2 || hwe.products.len() >= 2);

    let julia = reaction_for_reactants(&registry, "julia_olefination", &[sulfone, acetone])
        .expect("Julia olefination must be registered");
    assert!(julia.channels.len() >= 2 || julia.products.len() >= 2);
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
fn alcohol_acyl_protection_uses_regular_ester_chemistry() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let ethanol = SubstanceId::from("destroy:ethanol");
    let acetyl_chloride = SubstanceId::from("destroy:acetyl_chloride");

    dynamic
        .generate_reactions_for_substances([ethanol.clone(), acetyl_chloride.clone()], 1)
        .expect("acetyl chloride and ethanol must generate esterification");
    let registry = dynamic
        .to_registry()
        .expect("dynamic registry must convert");
    let esterification = reaction_for_reactants(
        &registry,
        "acyl_chloride_esterification",
        &[acetyl_chloride.clone(), ethanol.clone()],
    )
    .expect("acetyl chloride must acyl-protect ethanol as a regular ester");
    let ester = esterification
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:hydrochloric_acid")
        .expect("esterification must produce an ester")
        .substance_id
        .clone();

    let ester_structure = dynamic
        .substance(&ester)
        .expect("ester product must be registered")
        .molecular_structure
        .as_ref()
        .expect("ester product must have a molecular graph");
    let ester_sites = try_find_reactive_sites(ester_structure).expect("ester sites must be valid");
    assert!(
        ester_sites
            .iter()
            .any(|site| site.kind == ReactiveSiteKind::Ester),
        "acyl-protected alcohol must be represented as an ordinary ester"
    );
    assert!(
        !ester_sites
            .iter()
            .any(|site| site.kind == ReactiveSiteKind::Alcohol),
        "acyl-protected alcohol must not expose a free alcohol center"
    );

    dynamic
        .generate_reactions_for_substances([ester.clone()], 1)
        .expect("ester product must generate hydrolysis");
    let registry = dynamic
        .to_registry()
        .expect("dynamic registry must convert");
    let hydrolysis = reaction_for_reactants(&registry, "ester_hydrolysis", &[ester])
        .expect("acyl-protected alcohol must deprotect through ester hydrolysis");
    assert!(
        hydrolysis
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:acetic_acid"),
        "ester hydrolysis must restore the acyl fragment as acetic acid"
    );
    assert!(
        hydrolysis
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:ethanol"),
        "ester hydrolysis must restore the protected alcohol"
    );
}

#[test]
fn carboxylic_acid_protection_uses_concrete_ester_products() {
    let alcohols = [
        SubstanceId::from("destroy:methanol"),
        SubstanceId::from("destroy:ethanol"),
        SubstanceId::from("destroy:tert_butanol"),
    ];

    for alcohol in alcohols {
        let mut dynamic =
            super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acid = SubstanceId::from("destroy:acetic_acid");
        dynamic
            .generate_reactions_for_substances([acid.clone(), alcohol.clone()], 1)
            .expect("carboxylic acid and alcohol must generate ester protection");
        let registry = dynamic
            .to_registry()
            .expect("dynamic registry must convert");
        let esterification = reaction_for_reactants(
            &registry,
            "carboxylic_acid_esterification",
            &[acid.clone(), alcohol.clone()],
        )
        .expect("acid protection must be represented by ordinary esterification");

        let ester = esterification
            .products
            .iter()
            .find(|term| term.substance_id.as_str() != "destroy:water")
            .expect("esterification must produce an ester")
            .substance_id
            .clone();
        let ester_structure = dynamic
            .substance(&ester)
            .expect("ester product must be registered")
            .molecular_structure
            .as_ref()
            .expect("ester product must have a molecular graph");
        let ester_sites =
            try_find_reactive_sites(ester_structure).expect("ester sites must be valid");
        assert!(
            ester_sites
                .iter()
                .any(|site| site.kind == ReactiveSiteKind::Ester),
            "protected acid must be represented as an ester"
        );
        assert!(
            !ester_sites
                .iter()
                .any(|site| site.kind == ReactiveSiteKind::CarboxylicAcid),
            "protected acid must not expose a free carboxylic acid center"
        );

        dynamic
            .generate_reactions_for_substances([ester.clone()], 1)
            .expect("acid-protecting ester must generate hydrolysis");
        let registry = dynamic
            .to_registry()
            .expect("dynamic registry must convert");
        let hydrolysis = reaction_for_reactants(&registry, "ester_hydrolysis", &[ester])
            .expect("protected acid must deprotect through ester hydrolysis");
        assert!(
            hydrolysis
                .products
                .iter()
                .any(|term| term.substance_id == acid),
            "ester hydrolysis must restore the carboxylic acid"
        );
        assert!(
            hydrolysis
                .products
                .iter()
                .any(|term| term.substance_id == alcohol),
            "ester hydrolysis must restore the protecting alcohol"
        );
    }
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

#[test]
fn tms_protection_creates_protected_ether_without_free_alcohol_site() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let ethanol = SubstanceId::from("destroy:ethanol");
    dynamic
        .generate_reactions_for_substances([ethanol.clone()], 1)
        .unwrap();

    let protection = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("alcohol_silyl_protection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == ethanol)
        })
        .expect("ethanol must generate TMS protection");
    assert!(protection
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:trimethylsilyl_chloride"));
    assert!(protection
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrochloric_acid"));

    let protected_id = protection
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:hydrochloric_acid")
        .unwrap()
        .substance_id
        .clone();
    let protected = dynamic.substance(&protected_id).unwrap();
    let protected_structure = protected.molecular_structure.as_ref().unwrap();
    let site_kinds = try_find_reactive_sites(protected_structure)
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(
        site_kinds.contains(&ReactiveSiteKind::SilylEther),
        "TMS product must expose a silyl ether site"
    );
    assert!(
        !site_kinds.contains(&ReactiveSiteKind::Alcohol),
        "TMS-protected ethanol must not expose a free alcohol site"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn tms_deprotection_restores_original_alcohol() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let ethanol = SubstanceId::from("destroy:ethanol");
    dynamic
        .generate_reactions_for_substances([ethanol.clone()], 1)
        .unwrap();

    let protected_id = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("alcohol_silyl_protection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == ethanol)
        })
        .unwrap()
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:hydrochloric_acid")
        .unwrap()
        .substance_id
        .clone();

    dynamic
        .generate_reactions_for_substances([protected_id.clone()], 1)
        .unwrap();
    let deprotection = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("silyl_ether_deprotection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == protected_id)
        })
        .expect("TMS ether must generate fluoride deprotection");
    assert!(deprotection
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:fluoride"));
    assert!(deprotection
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:proton"));
    assert!(deprotection
        .products
        .iter()
        .any(|term| term.substance_id == ethanol));
    assert!(deprotection
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:trimethylsilyl_fluoride"));

    dynamic.to_registry().unwrap();
}

#[test]
fn acetal_hydrolysis_restores_carbonyl_and_concrete_alcohols() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let acetaldehyde = dynamic.resolve_frowns("CC=O").unwrap();
    let ethanol = SubstanceId::from("destroy:ethanol");
    dynamic
        .generate_reactions_for_substances([acetaldehyde.clone(), ethanol.clone()], 1)
        .unwrap();

    let acetal_reaction = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("acetal_formation/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == acetaldehyde)
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == ethanol)
        })
        .expect("acetaldehyde and ethanol must generate an acetal");
    let acetal_id = acetal_reaction
        .products
        .iter()
        .chain(
            acetal_reaction
                .channels
                .iter()
                .flat_map(|channel| channel.products.iter()),
        )
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("acetal formation must have a concrete acetal product")
        .substance_id
        .clone();

    let acetal = dynamic.substance(&acetal_id).unwrap();
    let site_kinds = try_find_reactive_sites(acetal.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(site_kinds.contains(&ReactiveSiteKind::Acetal));
    assert!(!site_kinds.contains(&ReactiveSiteKind::Carbonyl));
    assert!(!site_kinds.contains(&ReactiveSiteKind::Aldehyde));

    dynamic
        .generate_reactions_for_substances([acetal_id.clone()], 1)
        .unwrap();
    let hydrolysis = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("acetal_deprotection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == acetal_id)
        })
        .expect("acetal must generate hydrolysis");
    assert!(hydrolysis
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));
    assert!(hydrolysis
        .orders
        .keys()
        .any(|substance| substance.as_str() == "destroy:proton"));
    assert!(hydrolysis
        .products
        .iter()
        .any(|term| term.substance_id == acetaldehyde));
    assert_eq!(
        hydrolysis
            .products
            .iter()
            .filter(|term| term.substance_id == ethanol)
            .map(|term| term.coefficient)
            .sum::<u32>(),
        2
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn boc_protection_blocks_free_amine_and_deprotection_restores_it() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let methylamine = dynamic.resolve_frowns("CN").unwrap();
    dynamic
        .generate_reactions_for_substances([methylamine.clone()], 1)
        .unwrap();

    let protection = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("amine_boc_protection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == methylamine)
        })
        .expect("methylamine must generate Boc protection");
    assert!(protection
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:di_tert_butyl_dicarbonate"));
    assert!(protection
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:tert_butanol"));
    assert!(protection
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:carbon_dioxide"));

    let boc_id = protection
        .products
        .iter()
        .find(|term| {
            !matches!(
                term.substance_id.as_str(),
                "destroy:tert_butanol" | "destroy:carbon_dioxide"
            )
        })
        .expect("Boc protection must create a protected amine product")
        .substance_id
        .clone();
    let boc = dynamic.substance(&boc_id).unwrap();
    let site_kinds = try_find_reactive_sites(boc.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(site_kinds.contains(&ReactiveSiteKind::BocCarbamate));
    assert!(!site_kinds.contains(&ReactiveSiteKind::PrimaryAmine));
    assert!(!site_kinds.contains(&ReactiveSiteKind::NonTertiaryAmine));

    dynamic
        .generate_reactions_for_substances([boc_id.clone()], 1)
        .unwrap();
    let deprotection = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("boc_deprotection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == boc_id)
        })
        .expect("Boc carbamate must generate acidic hydrolysis");
    assert!(deprotection
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));
    assert!(deprotection
        .orders
        .keys()
        .any(|substance| substance.as_str() == "destroy:proton"));
    assert!(deprotection
        .products
        .iter()
        .any(|term| term.substance_id == methylamine));

    dynamic.to_registry().unwrap();
}

#[test]
fn cbz_protection_blocks_free_amine_and_hydrogenolysis_restores_it() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let methylamine = dynamic.resolve_frowns("CN").unwrap();
    dynamic
        .generate_reactions_for_substances([methylamine.clone()], 1)
        .unwrap();

    let protection = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("amine_cbz_protection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == methylamine)
        })
        .expect("methylamine must generate Cbz protection");
    assert!(protection
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:benzyl_chloroformate"));
    assert!(protection
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrochloric_acid"));

    let cbz_id = protection
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:hydrochloric_acid")
        .expect("Cbz protection must create a protected amine product")
        .substance_id
        .clone();
    let cbz = dynamic.substance(&cbz_id).unwrap();
    let site_kinds = try_find_reactive_sites(cbz.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(site_kinds.contains(&ReactiveSiteKind::CbzCarbamate));
    assert!(!site_kinds.contains(&ReactiveSiteKind::PrimaryAmine));
    assert!(!site_kinds.contains(&ReactiveSiteKind::NonTertiaryAmine));

    dynamic
        .generate_reactions_for_substances([cbz_id.clone()], 1)
        .unwrap();
    let deprotection = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("cbz_deprotection/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == cbz_id)
        })
        .expect("Cbz carbamate must generate hydrogenolysis");
    assert!(deprotection
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrogen"));
    assert!(deprotection
        .external_catalysts
        .iter()
        .any(|catalyst| catalyst.description == "forge:dusts/palladium"));
    assert!(deprotection
        .products
        .iter()
        .any(|term| term.substance_id == methylamine));
    assert!(deprotection
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:toluene"));
    assert!(deprotection
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:carbon_dioxide"));

    dynamic.to_registry().unwrap();
}
