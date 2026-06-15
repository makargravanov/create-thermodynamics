use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};

pub(crate) fn generate_borate_esterification(
    boric_acid: &BoricAcidSite<'_>,
    alcohol: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let boron_substance = boric_acid.participant.substance;
    let alcohol_substance = alcohol.participant.substance;

    let mut boron_editor = MolecularEditor::new(boric_acid.participant.structure);
    let boron_mapping =
        boron_editor.remove_atoms(&[boric_acid.hydroxyl_oxygen, boric_acid.hydroxyl_hydrogen])?;
    let boron = mapped_atom(&boron_mapping, boric_acid.boron, "boric acid boron")?;
    let boron_fragment = boron_editor.finish()?;

    let mut alcohol_editor = MolecularEditor::new(alcohol.participant.structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol.hydrogen])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol.oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &boron_fragment,
        boron,
        &alcohol_fragment,
        alcohol_oxygen,
        1.0,
    )?)?;

    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "borate_esterification",
        &boric_acid.participant,
        &alcohol.participant,
    ))
    .reactant(boron_substance.id.clone(), 1, 1)
    .reactant(alcohol_substance.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .condition(
        ReactionCondition::new("borate esterification is favored by acidic, water-poor conditions")
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.35),
    )
    .activation_energy_kj_per_mol(32.0)
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::FischerEsterification,
            SiteDescriptorBuilder::build(
                crate::chemistry::reactive_site::ReactiveSiteKind::BoricAcid,
                crate::chemistry::selectivity::types::SubstitutionDegree::Primary,
                0,
                0,
                0,
                false,
                false,
                false,
            ),
        )
        .with_secondary_site(SiteDescriptorBuilder::from_alcohol_site(alcohol))
        .never_suppress(),
    )
    .build())
}

pub(crate) fn generate_borane_oxidation(
    site: &BoraneSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let boron = site.boron;
    let mut editor = MolecularEditor::new(structure);
    editor.insert_bridging_atom(carbon, boron, "O", 0.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "borane_oxidation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen_peroxide", 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .build())
}

pub(crate) fn generate_borate_ester_hydrolysis(
    site: &BorateEsterSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let oxygen = site.oxygen;
    let boron = site.boron;
    let (first, first_mapping, second, second_mapping) =
        MolecularEditor::split_at_bond(structure, oxygen, boron)?;
    let (boron_fragment, boron_mapping, alcohol_fragment, oxygen_mapping) =
        if first_mapping[boron].is_some() {
            (first, first_mapping, second, second_mapping)
        } else {
            (second, second_mapping, first, first_mapping)
        };

    let mut boron_editor = MolecularEditor::new(&boron_fragment);
    let boron = mapped_atom(&boron_mapping, boron, "borate boron")?;
    add_hydroxyl(&mut boron_editor, boron)?;
    let boron_product = resolver.resolve(boron_editor.finish()?)?;

    let mut alcohol_editor = MolecularEditor::new(&alcohol_fragment);
    let oxygen = mapped_atom(&oxygen_mapping, oxygen, "borate ester oxygen")?;
    alcohol_editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let alcohol_product = resolver.resolve(alcohol_editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "borate_ester_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product(boron_product, 1)
    .product(alcohol_product, 1)
    .build())
}
