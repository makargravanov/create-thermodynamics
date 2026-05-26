use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;
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
