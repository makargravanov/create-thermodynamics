use super::*;
use crate::chemistry::condition::{AcidityCondition, AtmosphereCondition};
use crate::chemistry::mixture::Mixture;
use crate::chemistry::molecule::{StereoDescriptor, Stereochemistry};
use crate::chemistry::organic::generators::expand_stereo_product_distribution;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::reactive_site::{try_find_reactive_sites, ReactiveSiteKind};
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::simulation::reaction_rate_mol_per_bucket_per_tick;
use crate::chemistry::substance::{SubstanceId, SubstanceRepresentation};
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

const ORGANIC_MODEL_SOURCES: &[(&str, &str)] = &[
    ("organic/engine.rs", include_str!("engine.rs")),
    ("dynamic/mod.rs", include_str!("../dynamic/mod.rs")),
    ("synthesis.rs", include_str!("../synthesis.rs")),
    ("generators/acid_derivatives.rs", include_str!("generators/acid_derivatives.rs")),
    ("generators/addition.rs", include_str!("generators/addition.rs")),
    ("generators/alcohol.rs", include_str!("generators/alcohol.rs")),
    ("generators/aromatic.rs", include_str!("generators/aromatic.rs")),
    ("generators/boron.rs", include_str!("generators/boron.rs")),
    ("generators/c1_nitrogen.rs", include_str!("generators/c1_nitrogen.rs")),
    ("generators/carbonyl.rs", include_str!("generators/carbonyl.rs")),
    ("generators/combustion.rs", include_str!("generators/combustion.rs")),
    ("generators/cracking.rs", include_str!("generators/cracking.rs")),
    ("generators/cyclization.rs", include_str!("generators/cyclization.rs")),
    ("generators/enolate.rs", include_str!("generators/enolate.rs")),
    ("generators/heteroatom.rs", include_str!("generators/heteroatom.rs")),
    ("generators/heterocycle.rs", include_str!("generators/heterocycle.rs")),
    ("generators/organic_redox.rs", include_str!("generators/organic_redox.rs")),
    ("generators/organometallic.rs", include_str!("generators/organometallic.rs")),
    ("generators/phosphorus.rs", include_str!("generators/phosphorus.rs")),
    ("generators/polycondensation.rs", include_str!("generators/polycondensation.rs")),
    ("generators/protecting_groups.rs", include_str!("generators/protecting_groups.rs")),
    ("generators/radical.rs", include_str!("generators/radical.rs")),
    ("generators/rearrangement.rs", include_str!("generators/rearrangement.rs")),
    ("generators/ring_closure.rs", include_str!("generators/ring_closure.rs")),
    ("generators/substitution.rs", include_str!("generators/substitution.rs")),
];

const FORBIDDEN_TARGETED_ORGANIC_FRAGMENTS: &[&str] = &[
    "generate_caffeine",
    "generate_indole",
    "generate_uracil",
    "generate_xanthine",
    "\"destroy:caffeine\"",
    "\"destroy:indole\"",
    "\"destroy:uracil\"",
    "\"destroy:xanthine\"",
];

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

#[test]
fn organic_model_does_not_use_target_molecule_shortcuts() {
    for (path, source) in ORGANIC_MODEL_SOURCES {
        for forbidden in FORBIDDEN_TARGETED_ORGANIC_FRAGMENTS {
            assert!(
                !source.contains(forbidden),
                "{path} contains targeted organic shortcut {forbidden}"
            );
        }
    }
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
fn isocyanate_amine_addition_creates_urea_like_product() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let methyl_isocyanate = dynamic.resolve_frowns("CN=C=O").unwrap();
    let methylamine = dynamic.resolve_frowns("CN").unwrap();

    dynamic
        .generate_reactions_for_substances([methyl_isocyanate.clone(), methylamine.clone()], 1)
        .unwrap();

    let reaction = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("isocyanate_amine_addition/")
        })
        .expect("isocyanate and amine must add to a urea-like product");
    let product_id = reaction
        .products
        .iter()
        .find(|term| term.substance_id != methyl_isocyanate && term.substance_id != methylamine)
        .expect("addition must create an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap();
    let product_sites = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(product_sites.contains(&ReactiveSiteKind::UreaLike));
}

#[test]
fn isocyanate_ammonolysis_creates_urea_like_product() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let methyl_isocyanate = dynamic.resolve_frowns("CN=C=O").unwrap();

    dynamic
        .generate_reactions_for_substances([methyl_isocyanate.clone()], 1)
        .unwrap();

    let reaction = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("isocyanate_ammonolysis/"))
        .expect("isocyanate must react with ammonia into a urea-like product");
    assert!(reaction
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:ammonia"));
    let product_id = reaction.products[0].substance_id.clone();
    let product = dynamic.substance(&product_id).unwrap();
    let product_sites = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(product_sites.contains(&ReactiveSiteKind::UreaLike));
}

#[test]
fn amine_formylation_transfers_formyl_group_from_real_donor() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let methylamine = dynamic.resolve_frowns("CN").unwrap();
    let formic_acid = dynamic.resolve_frowns("O=CO").unwrap();

    dynamic
        .generate_reactions_for_substances([methylamine.clone(), formic_acid.clone()], 1)
        .unwrap();

    let reaction = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("amine_formylation/"))
        .expect("amine and formyl donor must form a formamide");
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));
    let product_id = reaction
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("formylation must create an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap();
    let product_sites = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(product_sites.contains(&ReactiveSiteKind::Amide));
}

#[test]
fn organic_redox_generates_graph_based_oxidation_paths() {
    let registry = generated_registry();
    let alcohol = reaction_with_prefix(&registry, "alcohol_oxidation/destroy_ethanol/");
    assert!(alcohol
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:dichromate"));
    assert!(alcohol.selectivity_profile.is_some());

    let peroxide =
        reaction_with_prefix(&registry, "alcohol_peroxide_overoxidation/destroy_ethanol/");
    assert!(peroxide
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrogen_peroxide"));
    assert!(peroxide
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:acetic_acid"));
}

#[test]
fn alkene_epoxidation_creates_epoxide_center_from_double_bond() {
    let registry = generated_registry();
    let epoxidation = reaction_with_prefix(&registry, "alkene_epoxidation/destroy_ethene/");
    assert!(epoxidation
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrogen_peroxide"));
    let product_id = epoxidation
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("epoxidation must have an organic product")
        .substance_id
        .clone();
    let product = registry.substance(&product_id).unwrap();
    let site_kinds = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(site_kinds.contains(&ReactiveSiteKind::Epoxide));
}

#[test]
fn rearrangement_generators_use_graph_migration_rules() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let acetone = SubstanceId::from("destroy:acetone");
    let acetone_oxime = dynamic.resolve_frowns("CC(=NO)C").unwrap();
    dynamic
        .generate_reactions_for_substances([acetone.clone(), acetone_oxime.clone()], 1)
        .unwrap();

    let baeyer_villiger = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("baeyer_villiger_rearrangement/")
        })
        .expect("ketone should generate Baeyer-Villiger oxygen insertion");
    assert!(baeyer_villiger
        .reactants
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:hydrogen_peroxide"));
    let ester_product = baeyer_villiger
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("Baeyer-Villiger must create an organic ester")
        .substance_id
        .clone();
    let ester = dynamic.substance(&ester_product).unwrap();
    let ester_sites = try_find_reactive_sites(ester.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(ester_sites.contains(&ReactiveSiteKind::Ester));

    let beckmann = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("beckmann_rearrangement/"))
        .expect("oxime should generate Beckmann migration");
    let amide_product = beckmann.products[0].substance_id.clone();
    let amide = dynamic.substance(&amide_product).unwrap();
    let amide_sites = try_find_reactive_sites(amide.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(amide_sites.contains(&ReactiveSiteKind::Amide));
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
fn organometallic_formation_creates_dynamic_reagent_from_halide() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    dynamic
        .generate_reactions_for_substances([SubstanceId::from("destroy:iodomethane")], 1)
        .unwrap();

    let organomagnesium = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("organomagnesium_formation/")
        })
        .expect("alkyl iodide should form an organomagnesium reagent");
    let product_id = organomagnesium.products[0].substance_id.clone();
    let product = dynamic.substance(&product_id).unwrap();
    let site_kinds = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(site_kinds.contains(&ReactiveSiteKind::Organomagnesium));

    let organolithium = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("organolithium_formation/"))
        .expect("alkyl iodide should form an organolithium reagent");
    assert!(organolithium.external_products.iter().any(|external| {
        external.description.as_str() == "external:lithium_I_salt"
            && external.molar_mass_grams.is_some()
    }));
}

#[test]
fn organometallic_reagent_adds_to_nitrile_and_opens_epoxide() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let methyl_magnesium_chloride = dynamic.resolve_frowns("CMgCl").unwrap();
    let acetonitrile = dynamic.resolve_frowns("CC#N").unwrap();
    let ethylene_oxide = dynamic
        .resolve_frowns(
            "destroy:graph:atoms=C.C.O.H.H.H.H;bonds=0-s-1,0-s-2,1-s-2,0-s-3,0-s-4,1-s-5,1-s-6",
        )
        .unwrap();

    dynamic
        .generate_reactions_for_substances(
            [methyl_magnesium_chloride, acetonitrile, ethylene_oxide],
            1,
        )
        .unwrap();

    assert!(dynamic.reactions().any(|reaction| reaction
        .id
        .as_str()
        .starts_with("organometallic_nitrile_addition/")));
    assert!(dynamic.reactions().any(|reaction| reaction
        .id
        .as_str()
        .starts_with("organometallic_epoxide_opening/")));
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
fn activated_methylene_condenses_with_carbonyl_to_alkene() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let malonaldehyde = dynamic
        .resolve_frowns(
            "destroy:graph:atoms=C.O.C.C.O.H.H.H.H;\
             bonds=0-d-1,0-s-2,2-s-3,3-d-4,0-s-5,2-s-6,2-s-7,3-s-8",
        )
        .unwrap();
    let acetaldehyde = dynamic.resolve_frowns("CC=O").unwrap();

    dynamic
        .generate_reactions_for_substances([malonaldehyde.clone(), acetaldehyde.clone()], 1)
        .unwrap();

    let reaction = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("knoevenagel_condensation/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == malonaldehyde)
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == acetaldehyde)
        })
        .expect("activated methylene and carbonyl must condense");
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));
    let product_id = reaction
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("condensation must create an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap();
    let product_sites = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(product_sites.contains(&ReactiveSiteKind::Alkene));
}

#[test]
fn bis_nucleophile_and_dicarbonyl_condense_to_n_heterocycle() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let urea_like = dynamic
        .resolve_frowns(
            "destroy:graph:atoms=C.O.N.N.H.H.H.H;\
             bonds=0-d-1,0-s-2,0-s-3,2-s-4,2-s-5,3-s-6,3-s-7",
        )
        .unwrap();
    let dicarbonyl = dynamic
        .resolve_frowns(
            "destroy:graph:atoms=C.O.C.C.O.H.H.H.H;\
             bonds=0-d-1,0-s-2,2-s-3,3-d-4,0-s-5,2-s-6,2-s-7,3-s-8",
        )
        .unwrap();
    let urea_like_sites = try_find_reactive_sites(
        dynamic
            .substance(&urea_like)
            .unwrap()
            .molecular_structure
            .as_ref()
            .unwrap(),
    )
    .unwrap()
    .into_iter()
    .map(|site| site.kind)
    .collect::<Vec<_>>();
    let dicarbonyl_sites = try_find_reactive_sites(
        dynamic
            .substance(&dicarbonyl)
            .unwrap()
            .molecular_structure
            .as_ref()
            .unwrap(),
    )
    .unwrap()
    .into_iter()
    .map(|site| site.kind)
    .collect::<Vec<_>>();
    assert!(urea_like_sites.contains(&ReactiveSiteKind::UreaLike));
    assert!(dicarbonyl_sites.contains(&ReactiveSiteKind::DicarbonylElectrophile));
    let urea_substance = dynamic.substance(&urea_like).unwrap();
    let urea_structure = urea_substance.molecular_structure.as_ref().unwrap();
    let urea_site = try_find_reactive_sites(urea_structure)
        .unwrap()
        .into_iter()
        .find(|site| site.kind == ReactiveSiteKind::UreaLike)
        .unwrap();
    let urea_center = SiteParticipant {
        substance: urea_substance,
        structure: urea_structure,
        site: urea_site,
    }
    .bis_nucleophile_center()
    .unwrap();
    assert_eq!(urea_center.class, BisNucleophileClass::UreaLike);

    let report = dynamic
        .generate_reactions_for_substances([urea_like.clone(), dicarbonyl.clone()], 1)
        .unwrap();
    assert!(
        report.generator_errors.is_empty(),
        "generation errors: {:?}",
        report.generator_errors
    );

    let reaction = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("bis_nucleophile_dicarbonyl_condensation/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == urea_like)
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == dicarbonyl)
        })
        .expect("bis-nucleophile and activated 1,3-dicarbonyl must condense");
    assert_eq!(
        reaction
            .products
            .iter()
            .filter(|term| term.substance_id.as_str() == "destroy:water")
            .map(|term| term.coefficient)
            .sum::<u32>(),
        2
    );
    let product_id = reaction
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("condensation must create an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap();
    let product_sites = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(product_sites.contains(&ReactiveSiteKind::UreaLike));
    assert!(product_sites.contains(&ReactiveSiteKind::Alkene));
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
fn chain_growth_polymerization_registers_a_polymer_material_from_an_alkene() {
    let registry = generated_registry();
    let reaction = reaction_with_prefix(&registry, "chain_growth_polymerization/destroy_ethene/");
    assert_eq!(reaction.reactants.len(), 1);
    assert_eq!(reaction.products.len(), 1);
    let monomer = registry.substance(&"destroy:ethene".into()).unwrap();
    let product = registry
        .substance(&reaction.products[0].substance_id)
        .unwrap();
    assert!(matches!(
        product.representation,
        SubstanceRepresentation::Polymer { .. }
    ));
    assert!(product.molecular_structure.is_none());
    assert!(
        (product.molar_mass_grams
            - monomer.molar_mass_grams * f64::from(reaction.reactants[0].coefficient))
        .abs()
            < 1.0e-6,
        "the polymer material mass must match the consumed monomer count"
    );
    assert_ne!(
        reaction.products[0].substance_id, monomer.id,
        "the polymer material is a distinct substance from the monomer"
    );
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

    assert!(
        (tert_rate - 0.0).abs() < 0.01,
        "tertiary SN2 must be runtime-suppressed"
    );
    assert!(
        eth_rate >= tert_rate,
        "primary SN2 should be faster than or equal to tertiary SN2"
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

#[test]
fn lactonization_closes_a_five_membered_lactone_from_a_hydroxy_acid() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // 4-hydroxybutanoic acid: HOOC-CH2-CH2-CH2-OH. The acyl carbon and the
    // alcohol oxygen are separated by a 4-edge path, so closure forms a
    // 5-membered ring (gamma-butyrolactone).
    let hydroxy_acid = dynamic.resolve_frowns("OC(=O)CCCO").unwrap();
    dynamic.generate_reactions_for(&hydroxy_acid, 1).unwrap();

    let lactonization = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("lactonization_5/"))
        .expect("hydroxy acid must close to a 5-membered lactone");
    assert!(lactonization
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    let product_id = lactonization
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("lactonization must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let site_kinds = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    // A lactone is a cyclic ester.
    assert!(site_kinds.contains(&ReactiveSiteKind::Ester));

    dynamic.to_registry().unwrap();
}

#[test]
fn lactonization_fuses_a_lactone_onto_an_existing_aromatic_ring() {
    // EMERGENCE CHECK: the generic ring-closure core uses `would_form_ring_of_size`,
    // a graph BFS blind to whether the closing atoms already sit in another ring.
    // So a closure that FUSES a new ring onto an existing one must fall out for free,
    // with no fused-ring-specific code. 2-(hydroxymethyl)benzoic acid (a benzene ring
    // bearing ortho -COOH and -CH2OH) closes into phthalide: a 5-membered lactone
    // sharing its C–C edge with the benzene ring.
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let benzoic = dynamic
        .resolve_frowns("destroy:benzene:C(=O)O,CO,,,,")
        .unwrap();
    dynamic.generate_reactions_for(&benzoic, 1).unwrap();

    let lactonization = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("lactonization_5/"))
        .expect("ortho hydroxymethyl benzoic acid must close to a 5-membered lactone");

    let product_id = lactonization
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("lactonization must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();

    // The new lactone is present (cyclic ester)...
    let site_kinds = try_find_reactive_sites(structure)
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(
        site_kinds.contains(&ReactiveSiteKind::Ester),
        "the fused product is a lactone (cyclic ester)"
    );

    // ...AND the benzene ring survived intact: six aromatic-order (1.5) C–C bonds
    // remain, proving the lactone FUSED onto the ring rather than consuming it.
    let aromatic_cc_bonds = structure
        .bonds
        .iter()
        .filter(|bond| {
            crate::chemistry::molecule::bond_order_matches(bond.order, 1.5)
                && structure.atoms[bond.from].element == "C"
                && structure.atoms[bond.to].element == "C"
        })
        .count();
    assert_eq!(
        aromatic_cc_bonds, 6,
        "the benzene ring stays aromatic after a lactone fuses onto it"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn lactamization_closes_a_five_membered_lactam_from_an_amino_acid() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // 4-aminobutanoic acid (GABA): HOOC-CH2-CH2-CH2-NH2 closes to a 5-membered
    // lactam (2-pyrrolidinone).
    let amino_acid = dynamic.resolve_frowns("OC(=O)CCCN").unwrap();
    dynamic.generate_reactions_for(&amino_acid, 1).unwrap();

    let lactamization = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("lactamization_5/"))
        .expect("amino acid must close to a 5-membered lactam");
    assert!(lactamization
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    let product_id = lactamization
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("lactamization must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let site_kinds = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    // The closure consumed the free acid into the ring.
    assert!(!site_kinds.contains(&ReactiveSiteKind::CarboxylicAcid));
    // The cyclic secondary amide is correctly perceived as an amide, and its
    // delocalised nitrogen is NOT surfaced as a basic/nucleophilic amine.
    assert!(site_kinds.contains(&ReactiveSiteKind::Amide));
    assert!(!site_kinds.contains(&ReactiveSiteKind::PrimaryAmine));
    assert!(!site_kinds.contains(&ReactiveSiteKind::NonTertiaryAmine));
    // The product keeps its single nitrogen, now in the lactam ring.
    let structure = product.molecular_structure.as_ref().unwrap();
    let nitrogen_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "N")
        .count();
    assert_eq!(nitrogen_count, 1, "lactam keeps its single nitrogen");

    dynamic.to_registry().unwrap();
}

#[test]
fn lactonization_is_rejected_for_a_strained_three_membered_ring() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // 2-hydroxyacetic acid (glycolic acid): HOOC-CH2-OH. Closure would form a
    // strained 3-membered alpha-lactone, which the ring-size gate rejects.
    let glycolic_acid = dynamic.resolve_frowns("OC(=O)CO").unwrap();
    dynamic.generate_reactions_for(&glycolic_acid, 1).unwrap();
    assert!(
        dynamic
            .reactions()
            .all(|reaction| !reaction.id.as_str().starts_with("lactonization_")),
        "alpha-lactone closure must be rejected as too strained"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn paal_knorr_closes_a_furan_from_a_1_4_diketone() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Hexane-2,5-dione: CH3-CO-CH2-CH2-CO-CH3. The two carbonyl carbons are
    // 1,4-related, so acid-catalysed cyclodehydration closes a furan ring
    // (2,5-dimethylfuran) with loss of water.
    let diketone = dynamic.resolve_frowns("CC(=O)CCC(=O)C").unwrap();
    dynamic.generate_reactions_for(&diketone, 1).unwrap();

    let paal_knorr = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("paal_knorr_furan/"))
        .expect("1,4-diketone must close to a furan");
    // Furan closure is intramolecular: the diketone is the sole reactant and
    // exactly one water leaves (only one of the two carbonyl oxygens departs).
    assert!(paal_knorr
        .reactants
        .iter()
        .any(|term| term.substance_id == diketone));
    let water_term = paal_knorr
        .products
        .iter()
        .find(|term| term.substance_id.as_str() == "destroy:water")
        .expect("furan closure must expel water");
    assert_eq!(
        water_term.coefficient, 1,
        "one carbonyl oxygen leaves as water"
    );

    let product_id = paal_knorr
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("Paal–Knorr must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    // The product furan is aromatic: its ring bonds were aromatised on finish().
    let aromatic_bonds = structure
        .bonds
        .iter()
        .filter(|bond| crate::chemistry::molecule::bond_order_matches(bond.order, 1.5))
        .count();
    assert!(
        aromatic_bonds >= 5,
        "furan ring must be aromatic (got {aromatic_bonds} aromatic bonds)"
    );
    // The closure consumed both carbonyls; no ketone site should remain.
    let site_kinds = try_find_reactive_sites(structure)
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(!site_kinds.contains(&ReactiveSiteKind::Ketone));
    // The single furan oxygen is retained (one carbonyl O left as water).
    let oxygen_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "O")
        .count();
    assert_eq!(oxygen_count, 1, "furan keeps a single ring oxygen");

    dynamic.to_registry().unwrap();
}

#[test]
fn paal_knorr_closes_a_pyrrole_from_a_1_4_diketone_and_amine() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Hexane-2,5-dione + methylamine condense to an N-methylpyrrole, losing two
    // waters. The closure is intermolecular: the amine is a separate donor.
    let diketone = dynamic.resolve_frowns("CC(=O)CCC(=O)C").unwrap();
    let methylamine = dynamic.resolve_frowns("CN").unwrap();
    dynamic
        .generate_reactions_for_substances([diketone.clone(), methylamine.clone()], 1)
        .unwrap();

    let pyrrole = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("paal_knorr_pyrrole/"))
        .expect("1,4-diketone and amine must close to a pyrrole");
    // Two waters leave (one per carbonyl oxygen).
    let water_term = pyrrole
        .products
        .iter()
        .find(|term| term.substance_id.as_str() == "destroy:water")
        .expect("pyrrole closure must expel water");
    assert_eq!(
        water_term.coefficient, 2,
        "two carbonyl oxygens leave as water"
    );
    // Both the diketone and the amine are consumed.
    assert!(pyrrole
        .reactants
        .iter()
        .any(|term| term.substance_id == diketone));
    assert!(pyrrole
        .reactants
        .iter()
        .any(|term| term.substance_id == methylamine));

    let product_id = pyrrole
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("Paal–Knorr must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    // The pyrrole ring is aromatic.
    let aromatic_bonds = structure
        .bonds
        .iter()
        .filter(|bond| crate::chemistry::molecule::bond_order_matches(bond.order, 1.5))
        .count();
    assert!(
        aromatic_bonds >= 5,
        "pyrrole ring must be aromatic (got {aromatic_bonds} aromatic bonds)"
    );
    // No carbonyl survived the condensation.
    let site_kinds = try_find_reactive_sites(structure)
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(!site_kinds.contains(&ReactiveSiteKind::Ketone));
    // The ring nitrogen is the amine's nitrogen, now aromatic and oxygen-free.
    let oxygen_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "O")
        .count();
    assert_eq!(oxygen_count, 0, "pyrrole has no oxygen");
    // Exactly one nitrogen — the amine fragment was spliced in once, not
    // duplicated or dropped (a 2-N or 0-N product would still be aromatic and
    // oxygen-free, so this guards a blind spot the other assertions miss).
    let nitrogen_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "N")
        .count();
    assert_eq!(nitrogen_count, 1, "pyrrole keeps a single ring nitrogen");

    dynamic.to_registry().unwrap();
}

#[test]
fn paal_knorr_closes_a_thiophene_from_a_1_4_diketone_and_hydrogen_sulfide() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Hexane-2,5-dione + hydrogen sulfide condense to 2,5-dimethylthiophene,
    // losing two waters. The sulfur donor (H2S) bridges both ring carbons.
    let diketone = dynamic.resolve_frowns("CC(=O)CCC(=O)C").unwrap();
    let hydrogen_sulfide = dynamic.resolve_frowns("S").unwrap();
    dynamic
        .generate_reactions_for_substances([diketone.clone(), hydrogen_sulfide.clone()], 1)
        .unwrap();

    let thiophene = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("paal_knorr_thiophene/"))
        .expect("1,4-diketone and H2S must close to a thiophene");
    let water_term = thiophene
        .products
        .iter()
        .find(|term| term.substance_id.as_str() == "destroy:water")
        .expect("thiophene closure must expel water");
    assert_eq!(
        water_term.coefficient, 2,
        "two carbonyl oxygens leave as water"
    );
    assert!(thiophene
        .reactants
        .iter()
        .any(|term| term.substance_id == hydrogen_sulfide));

    let product_id = thiophene
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("Paal–Knorr must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    // The thiophene ring is aromatic.
    let aromatic_bonds = structure
        .bonds
        .iter()
        .filter(|bond| crate::chemistry::molecule::bond_order_matches(bond.order, 1.5))
        .count();
    assert!(
        aromatic_bonds >= 5,
        "thiophene ring must be aromatic (got {aromatic_bonds} aromatic bonds)"
    );
    // No carbonyl survived; the product carries exactly one ring sulfur and no oxygen.
    let site_kinds = try_find_reactive_sites(structure)
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(!site_kinds.contains(&ReactiveSiteKind::Ketone));
    let sulfur_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "S")
        .count();
    assert_eq!(sulfur_count, 1, "thiophene keeps a single ring sulfur");
    let oxygen_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "O")
        .count();
    assert_eq!(oxygen_count, 0, "thiophene has no oxygen");

    dynamic.to_registry().unwrap();
}

#[test]
fn paal_knorr_thiophene_rejects_a_monothiol_donor() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Ethanethiol (CCS) has a single S-H, so it cannot bridge both ring carbons.
    // The generator's `thiol.hydrogens.len() < 2` guard must reject it: pairing a
    // monothiol with the 1,4-diketone yields no thiophene closure.
    let diketone = dynamic.resolve_frowns("CC(=O)CCC(=O)C").unwrap();
    let ethanethiol = dynamic.resolve_frowns("CCS").unwrap();
    dynamic
        .generate_reactions_for_substances([diketone.clone(), ethanethiol.clone()], 1)
        .unwrap();

    assert!(
        dynamic
            .reactions()
            .all(|reaction| !reaction.id.as_str().starts_with("paal_knorr_thiophene/")),
        "a monothiol with one S-H must not be accepted as a thiophene sulfur donor"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn diels_alder_closes_a_cyclohexene_from_butadiene_and_ethylene() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // 1,3-butadiene (C=CC=C) + ethylene (C=C) cycloadd to cyclohexene. No atoms
    // leave: the six ring carbons are the four diene carbons plus the two
    // dienophile carbons; one ring double bond remains.
    let butadiene = dynamic.resolve_frowns("C=CC=C").unwrap();
    let ethylene = dynamic.resolve_frowns("C=C").unwrap();
    dynamic
        .generate_reactions_for_substances([butadiene.clone(), ethylene.clone()], 1)
        .unwrap();

    let diels_alder = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("diels_alder/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == ethylene)
        })
        .expect("butadiene and ethylene must cycloadd to a cyclohexene");
    // Both partners are consumed; the cycloaddition is atom-economical (no byproduct).
    assert!(diels_alder
        .reactants
        .iter()
        .any(|term| term.substance_id == butadiene));
    assert!(diels_alder
        .reactants
        .iter()
        .any(|term| term.substance_id == ethylene));
    assert_eq!(
        diels_alder.products.len(),
        1,
        "Diels–Alder produces a single product with no byproduct"
    );

    let product_id = diels_alder.products[0].substance_id.clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    // Cyclohexene: six carbons, exactly one C=C double bond, all atoms retained.
    let carbon_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "C")
        .count();
    assert_eq!(carbon_count, 6, "cyclohexene has six ring carbons");
    let double_bonds = structure
        .bonds
        .iter()
        .filter(|bond| crate::chemistry::molecule::bond_order_matches(bond.order, 2.0))
        .count();
    assert_eq!(
        double_bonds, 1,
        "cyclohexene retains exactly one C=C double bond"
    );
    // The product is a ring: a six-membered carbocycle has 6 C–C ring bonds.
    let carbon_carbon_bonds = structure
        .bonds
        .iter()
        .filter(|bond| {
            structure.atoms[bond.from].element == "C" && structure.atoms[bond.to].element == "C"
        })
        .count();
    assert_eq!(
        carbon_carbon_bonds, 6,
        "cyclohexene ring has six carbon-carbon bonds"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn diels_alder_rejects_a_non_conjugated_diene() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Penta-1,4-diene (C=CCC=C) has its two double bonds separated by an sp3 CH2,
    // so they are not conjugated. conjugated_diene_carbons walks a single C2–C3
    // bond to a C3=C4 double bond and finds only a saturated carbon here, so no
    // [4+2] may fire.
    let pentadiene = dynamic.resolve_frowns("C=CCC=C").unwrap();
    let ethylene = dynamic.resolve_frowns("C=C").unwrap();
    dynamic
        .generate_reactions_for_substances([pentadiene.clone(), ethylene.clone()], 1)
        .unwrap();

    assert!(
        dynamic
            .reactions()
            .all(|reaction| !reaction.id.as_str().starts_with("diels_alder/")),
        "a non-conjugated 1,4-diene must not undergo a Diels–Alder cycloaddition"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn diels_alder_rejects_an_alkyne_dienophile() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Acetylene (C#C) is an alkyne, not an alkene; the generator's is_alkyne guard
    // must keep it out of the dienophile role for this [4+2] generator.
    let butadiene = dynamic.resolve_frowns("C=CC=C").unwrap();
    let acetylene = dynamic.resolve_frowns("C#C").unwrap();
    dynamic
        .generate_reactions_for_substances([butadiene.clone(), acetylene.clone()], 1)
        .unwrap();

    // Butadiene homodimerizes, so a diels_alder reaction does exist — but none may
    // involve acetylene, which the is_alkyne guard keeps out of the dienophile role.
    assert!(
        dynamic.reactions().all(|reaction| {
            !reaction.id.as_str().starts_with("diels_alder/")
                || reaction
                    .reactants
                    .iter()
                    .all(|term| term.substance_id != acetylene)
        }),
        "an alkyne must not be accepted as a dienophile by this generator"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn diels_alder_homodimerizes_butadiene() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Two molecules of one species react: butadiene acts as both diene and
    // dienophile (2 C4H6 → 4-vinylcyclohexene, C8H12). This intermolecular
    // homodimerization is a real Diels–Alder and must be generated even though
    // both partners share a substance id.
    let butadiene = dynamic.resolve_frowns("C=CC=C").unwrap();
    dynamic.generate_reactions_for(&butadiene, 1).unwrap();

    let dimer = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("diels_alder/"))
        .expect("butadiene must homodimerize via Diels–Alder");
    // The one species is consumed twice: a single coefficient-2 reactant term
    // (second-order), not two duplicate terms.
    assert_eq!(
        dimer.reactants.len(),
        1,
        "homodimerization lists the species once, not twice"
    );
    assert_eq!(dimer.reactants[0].substance_id, butadiene);
    assert_eq!(
        dimer.reactants[0].coefficient, 2,
        "two molecules of butadiene are consumed"
    );

    let product_id = dimer.products[0].substance_id.clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    // 4-vinylcyclohexene: eight carbons (six-membered ring + vinyl), and two C=C
    // double bonds (one in the ring, one in the pendant vinyl group).
    let carbon_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "C")
        .count();
    assert_eq!(carbon_count, 8, "4-vinylcyclohexene has eight carbons");
    let double_bonds = structure
        .bonds
        .iter()
        .filter(|bond| crate::chemistry::molecule::bond_order_matches(bond.order, 2.0))
        .count();
    assert_eq!(
        double_bonds, 2,
        "ring double bond plus the pendant vinyl double bond"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn retro_diels_alder_splits_a_cyclohexene_adduct_back_to_diene_and_dienophile() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let butadiene = dynamic.resolve_frowns("C=CC=C").unwrap();
    let ethylene = dynamic.resolve_frowns("C=C").unwrap();
    dynamic
        .generate_reactions_for_substances([butadiene.clone(), ethylene.clone()], 1)
        .unwrap();
    let adduct = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("diels_alder/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == ethylene)
        })
        .unwrap()
        .products[0]
        .substance_id
        .clone();

    dynamic.generate_reactions_for(&adduct, 1).unwrap();

    let retro = dynamic
        .reactions()
        .find(|reaction| {
            reaction.id.as_str().starts_with("retro_diels_alder/")
                && reaction
                    .reactants
                    .iter()
                    .any(|term| term.substance_id == adduct)
        })
        .expect("a Diels-Alder adduct must have a thermal cycloreversion path");
    assert!(retro
        .products
        .iter()
        .any(|term| term.substance_id == butadiene));
    assert!(retro
        .products
        .iter()
        .any(|term| term.substance_id == ethylene));

    dynamic.to_registry().unwrap();
}

#[test]
fn alkene_photoisomerization_creates_separate_e_and_z_channels() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let butene = dynamic.resolve_frowns("CC=CC").unwrap();
    dynamic.generate_reactions_for(&butene, 1).unwrap();

    let photo = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("alkene_photoisomerization/")
        })
        .expect("an unsymmetrical alkene must have light-driven E/Z isomerization channels");
    assert_eq!(photo.channels.len(), 2);

    let mut has_e = false;
    let mut has_z = false;
    for channel in &photo.channels {
        let product_id = &channel.products[0].substance_id;
        let product = dynamic.substance(product_id).unwrap();
        let structure = product.molecular_structure.as_ref().unwrap();
        has_e |= structure.stereochemistry.iter().any(|stereo| {
            matches!(stereo, Stereochemistry::DoubleBond(double) if double.descriptor == StereoDescriptor::E)
        });
        has_z |= structure.stereochemistry.iter().any(|stereo| {
            matches!(stereo, Stereochemistry::DoubleBond(double) if double.descriptor == StereoDescriptor::Z)
        });
    }
    assert!(has_e, "one channel must produce the E isomer");
    assert!(has_z, "one channel must produce the Z isomer");

    dynamic.to_registry().unwrap();
}

#[test]
fn intramolecular_n_alkylation_closes_a_pyrrolidine() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // 4-chlorobutan-1-amine (N-C-C-C-C-Cl): the amine nitrogen displaces the
    // terminal chloride in an internal SN2, closing a five-membered pyrrolidine
    // ring and expelling HCl (chloride + proton).
    let amino_halide = dynamic.resolve_frowns("NCCCCCl").unwrap();
    dynamic.generate_reactions_for(&amino_halide, 1).unwrap();

    let closure = dynamic
        .reactions()
        .find(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("intramolecular_n_alkylation_5/")
        })
        .expect("4-chlorobutan-1-amine must close to a pyrrolidine");
    // Intramolecular: the single substance is the only reactant.
    assert_eq!(closure.reactants.len(), 1);
    assert_eq!(closure.reactants[0].substance_id, amino_halide);
    // HX leaves as a halide ion plus a proton.
    assert!(closure
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:chloride"));
    assert!(closure
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:proton"));

    let product_id = closure
        .products
        .iter()
        .find(|term| {
            term.substance_id.as_str() != "destroy:chloride"
                && term.substance_id.as_str() != "destroy:proton"
        })
        .expect("must have an organic ring product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    // Pyrrolidine: one nitrogen, four carbons, no halogen remaining.
    let nitrogen_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "N")
        .count();
    assert_eq!(nitrogen_count, 1, "pyrrolidine has one ring nitrogen");
    let chlorine_count = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "Cl")
        .count();
    assert_eq!(chlorine_count, 0, "the chloride left the molecule");
    // The nitrogen is now bonded to two ring carbons (it closed the ring).
    let nitrogen_index = structure
        .atoms
        .iter()
        .position(|atom| atom.element == "N")
        .unwrap();
    let nitrogen_carbon_bonds = structure
        .neighbors(nitrogen_index)
        .into_iter()
        .filter(|(n, _)| structure.atoms[*n].element == "C")
        .count();
    assert_eq!(
        nitrogen_carbon_bonds, 2,
        "ring nitrogen bonds to two carbons after closure"
    );
    // Prove a RING actually closed (not an open-chain secondary amine, which would
    // also have 2 N–C bonds): pyrrolidine has 4 carbons all retained, and the ring
    // means C-count + N-count == bond count among those ring atoms forms a cycle.
    let total_carbons = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "C")
        .count();
    assert_eq!(total_carbons, 4, "pyrrolidine ring has four carbons");
    // A 5-membered ring of 4 C + 1 N has exactly 5 ring bonds; an open chain of the
    // same atoms would have only 4. Count bonds among the {C,N} heavy atoms.
    let heavy_bonds = structure
        .bonds
        .iter()
        .filter(|bond| {
            let a = &structure.atoms[bond.from].element;
            let b = &structure.atoms[bond.to].element;
            (a == "C" || a == "N") && (b == "C" || b == "N")
        })
        .count();
    assert_eq!(
        heavy_bonds, 5,
        "a closed pyrrolidine ring has five heavy-atom ring bonds (open chain would have four)"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn intramolecular_n_alkylation_rejects_a_strained_three_membered_ring() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // 2-chloroethan-1-amine (N-C-C-Cl) could only close a 3-membered aziridine,
    // below MIN_CLOSABLE_RING (4); the ring-size guard must reject it.
    let amino_halide = dynamic.resolve_frowns("NCCCl").unwrap();
    dynamic.generate_reactions_for(&amino_halide, 1).unwrap();

    assert!(
        dynamic.reactions().all(|reaction| !reaction
            .id
            .as_str()
            .starts_with("intramolecular_n_alkylation_")),
        "a strained three-membered aziridine closure must be rejected"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn intramolecular_n_alkylation_does_not_cross_substance_boundaries() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Two DIFFERENT amino-halide molecules in one registry. The intramolecular
    // closure must only fire within a single molecule (same substance id), never
    // pairing molecule A's amine with molecule B's halide.
    let a = dynamic.resolve_frowns("NCCCCCl").unwrap();
    let b = dynamic.resolve_frowns("NCCCCCCl").unwrap();
    dynamic
        .generate_reactions_for_substances([a.clone(), b.clone()], 1)
        .unwrap();

    for reaction in dynamic.reactions() {
        if reaction
            .id
            .as_str()
            .starts_with("intramolecular_n_alkylation_")
        {
            // Each closure has exactly one reactant — the single molecule it closes.
            assert_eq!(
                reaction.reactants.len(),
                1,
                "an intramolecular closure consumes one substance, not a cross-substance pair"
            );
        }
    }

    dynamic.to_registry().unwrap();
}

#[test]
fn intramolecular_n_alkylation_rejects_a_vinyl_halide() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // 4-amino-1-chlorobut-1-ene (N-C-C-C=C-Cl): the chloride sits on an sp2 alkene
    // carbon. An internal SN2 needs back-side attack on a tetrahedral sp3 carbon,
    // impossible at a planar vinyl carbon, so the sp2 guard must reject it.
    let vinyl_halide = dynamic.resolve_frowns("NCCC=CCl").unwrap();
    dynamic.generate_reactions_for(&vinyl_halide, 1).unwrap();

    assert!(
        dynamic.reactions().all(|reaction| !reaction
            .id
            .as_str()
            .starts_with("intramolecular_n_alkylation_")),
        "an internal SN2 on an sp2 vinyl-halide carbon must be rejected"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn amidation_condenses_an_acid_and_an_amine_into_an_amide() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    // Acetic acid + methylamine → N-methylacetamide + water (intermolecular).
    let acid = dynamic.resolve_frowns("CC(=O)O").unwrap();
    let amine = dynamic.resolve_frowns("CN").unwrap();
    dynamic
        .generate_reactions_for_substances([acid.clone(), amine.clone()], 1)
        .unwrap();

    let amidation = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("amidation/"))
        .expect("acid and amine on separate molecules must condense to an amide");
    assert!(amidation
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));
    assert!(amidation.reactants.iter().any(|t| t.substance_id == acid));
    assert!(amidation.reactants.iter().any(|t| t.substance_id == amine));

    let product_id = amidation
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("amidation must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    let site_kinds = try_find_reactive_sites(structure)
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(
        site_kinds.contains(&ReactiveSiteKind::Amide),
        "the condensation product is an amide"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn amidine_cyclization_closes_a_benzimidazole_fused_onto_benzene() {
    // EMERGENCE CHECK (generic imidazole closure): o-formamidoaniline is a benzene
    // ring bearing an ortho -NH2 and an ortho -NHCHO (formamide). The free amine
    // nitrogen attacks the formamide carbon, forming a C=N while the carbonyl
    // oxygen leaves as water — closing a 5-membered imidazole that shares its C–C
    // edge with the benzene ring (benzimidazole). The ring-closure core is blind to
    // the pre-existing benzene, so the fused bicycle falls out with no benzimidazole
    // /purine-specific code. This is the same closure that builds xanthine's imidazole.
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let aniline = dynamic
        .resolve_frowns("destroy:benzene:N,NC=O,,,,")
        .unwrap();
    dynamic.generate_reactions_for(&aniline, 1).unwrap();

    let closure = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("amidine_cyclization_5/"))
        .expect("o-formamidoaniline must close to a 5-membered imidazole (benzimidazole)");
    assert!(closure
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    let product_id = closure
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("benzimidazole closure must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();

    // The benzene ring survives (6 aromatic C–C bonds) AND the new imidazole is
    // aromatic too: a fused benzimidazole has more aromatic bonds than benzene alone.
    let aromatic_bonds = structure
        .bonds
        .iter()
        .filter(|bond| crate::chemistry::molecule::bond_order_matches(bond.order, 1.5))
        .count();
    assert!(
        aromatic_bonds > 6,
        "the fused benzimidazole is aromatic across both rings (got {aromatic_bonds})"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn amidine_cyclization_closes_xanthine_imidazole_onto_the_pyrimidinedione() {
    // CAFFEINE PATH, ring-closure half: 6-amino-5-formamidouracil is the Traube
    // intermediate — a pyrimidine-2,4-dione (uracil) bearing an ortho amino group
    // (C6-NH2, the nucleophile) and an ortho formamido group (C5-NH-CHO, the amide).
    // The amino nitrogen attacks the formyl carbon, closing a 5-membered imidazole
    // that shares its C5–C6 edge with the pyrimidine ring: that fused bicycle IS
    // xanthine's purine skeleton. The generic closure builds it with no purine code.
    //
    // Atom map: 0 N1(H) 1 C2(=O2) 3 N3(H) 4 C4(=O5) 6 C5 7 C6 8 N-amino(2H)
    //           9 N-amido(H) 10 C-formyl(H) 11 O-formyl. H: 12,13,14,15,16,17.
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let precursor = dynamic
        .resolve_frowns(
            "destroy:graph:atoms=N.C.O.N.C.O.C.C.N.N.C.O.H.H.H.H.H.H;\
             bonds=0-s-1,1-d-2,1-s-3,3-s-4,4-d-5,4-s-6,6-d-7,7-s-0,7-s-8,6-s-9,\
             9-s-10,10-d-11,0-s-12,3-s-13,8-s-14,8-s-15,9-s-16,10-s-17",
        )
        .unwrap();
    dynamic.generate_reactions_for(&precursor, 1).unwrap();

    let closure = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("amidine_cyclization_5/"))
        .expect("the formamido intermediate must close xanthine's 5-membered imidazole");
    assert!(closure
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    let product_id = closure
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("xanthine closure must have an organic product")
        .substance_id
        .clone();
    assert_ne!(
        product_id, precursor,
        "the product is a new (cyclized) substance"
    );
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();

    // The pyrimidinedione's two carbonyls survive the closure: the imidazole FUSED
    // onto the ring, it did not consume it. (Exocyclic lactam C=O stay double bonds.)
    let carbonyls = structure
        .bonds
        .iter()
        .filter(|bond| {
            crate::chemistry::molecule::bond_order_matches(bond.order, 2.0)
                && ((structure.atoms[bond.from].element == "C"
                    && structure.atoms[bond.to].element == "O")
                    || (structure.atoms[bond.from].element == "O"
                        && structure.atoms[bond.to].element == "C"))
        })
        .count();
    assert_eq!(
        carbonyls, 2,
        "xanthine keeps both pyrimidinedione carbonyls"
    );

    dynamic.to_registry().unwrap();
}

#[test]
fn carbonyl_and_hydrazine_like_bis_nucleophile_form_hydrazone() {
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let acetone = dynamic.resolve_frowns("CC(=O)C").unwrap();
    let hydrazine = SubstanceId::from("destroy:hydrazine");
    dynamic
        .generate_reactions_for_substances([acetone.clone(), hydrazine], 1)
        .unwrap();

    let reaction = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("hydrazone_formation/"))
        .expect("carbonyl plus hydrazine-like bis-nucleophile must form a hydrazone");
    assert!(reaction
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:water"));

    let product_id = reaction
        .products
        .iter()
        .find(|term| term.substance_id.as_str() != "destroy:water")
        .expect("hydrazone formation must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap();
    let site_kinds = try_find_reactive_sites(product.molecular_structure.as_ref().unwrap())
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();
    assert!(site_kinds.contains(&ReactiveSiteKind::Hydrazone));
}

/// Uracil: pyrimidine-2,4-dione. Ring N1 sits between C2=O and C6... (N1 is bonded
/// to C2=O and to C6; N3 sits between C2=O and C4=O). Both ring nitrogens are imide
/// N-H. Graph: 0 N1(H) 1 C2 2 O 3 N3(H) 4 C4 5 O 6 C5(H) 7 C6(H).
const URACIL_FROWNS: &str = "destroy:graph:atoms=N.C.O.N.C.O.C.C.H.H.H.H;\
     bonds=0-s-1,1-d-2,1-s-3,3-s-4,4-d-5,4-s-6,6-d-7,7-s-0,0-s-8,3-s-9,6-s-10,7-s-11";

#[test]
fn ring_amide_nitrogen_is_not_perceived_as_a_basic_amine() {
    // CONTAINMENT CHECK: a ring imide N-H must surface as the dedicated
    // AmideNitrogen site (alkylation-only), NOT as a basic PrimaryAmine /
    // NonTertiaryAmine — otherwise it would wrongly feed esterification, imine
    // formation, Paal–Knorr, etc.
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let uracil = dynamic.resolve_frowns(URACIL_FROWNS).unwrap();
    let structure = dynamic
        .substance(&uracil)
        .unwrap()
        .molecular_structure
        .as_ref()
        .unwrap()
        .clone();
    let site_kinds = try_find_reactive_sites(&structure)
        .unwrap()
        .into_iter()
        .map(|site| site.kind)
        .collect::<Vec<_>>();

    assert!(
        site_kinds.contains(&ReactiveSiteKind::AmideNitrogen),
        "uracil's ring N-H must surface as AmideNitrogen"
    );
    assert!(
        !site_kinds.contains(&ReactiveSiteKind::PrimaryAmine)
            && !site_kinds.contains(&ReactiveSiteKind::NonTertiaryAmine),
        "an amide/imide nitrogen must never be perceived as a basic amine"
    );
}

#[test]
fn amide_n_alkylation_methylates_a_ring_amide_n_h() {
    // Uracil + iodomethane → N-methyluracil + iodide + proton. The ring imide N-H
    // is alkylated by the methyl halide (generic over any amide N-H).
    let mut dynamic =
        super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
    let uracil = dynamic.resolve_frowns(URACIL_FROWNS).unwrap();
    let iodomethane = SubstanceId::from("destroy:iodomethane");
    dynamic
        .generate_reactions_for_substances([uracil.clone(), iodomethane.clone()], 1)
        .unwrap();

    let alkylation = dynamic
        .reactions()
        .find(|reaction| reaction.id.as_str().starts_with("amide_n_alkylation/"))
        .expect("iodomethane must N-methylate uracil's ring amide N-H");
    assert!(alkylation
        .reactants
        .iter()
        .any(|t| t.substance_id == uracil));
    assert!(alkylation
        .reactants
        .iter()
        .any(|t| t.substance_id == iodomethane));
    assert!(alkylation
        .products
        .iter()
        .any(|term| term.substance_id.as_str() == "destroy:iodide"));

    let product_id = alkylation
        .products
        .iter()
        .find(|term| {
            !matches!(
                term.substance_id.as_str(),
                "destroy:iodide" | "destroy:proton"
            )
        })
        .expect("alkylation must have an organic product")
        .substance_id
        .clone();
    let product = dynamic.substance(&product_id).unwrap().clone();
    let structure = product.molecular_structure.as_ref().unwrap();
    // One more carbon than uracil (the added methyl).
    let carbons = structure
        .atoms
        .iter()
        .filter(|atom| atom.element == "C")
        .count();
    assert_eq!(
        carbons, 5,
        "N-methyluracil has uracil's 4 carbons plus a methyl"
    );

    dynamic.to_registry().unwrap();
}
