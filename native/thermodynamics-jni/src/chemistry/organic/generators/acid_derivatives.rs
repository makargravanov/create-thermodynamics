use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::{SelectivityEngine, SiteDescriptorBuilder},
    types::{SelectivityContext, SelectivityRecommendation},
};

pub(crate) fn generate_carboxylic_acid_esterification(
    acid_site: &CarboxylicAcidSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    let acid_desc = SiteDescriptorBuilder::from_carboxylic_acid_site(acid_site);
    let alcohol_desc = SiteDescriptorBuilder::from_alcohol_site(alcohol_site);
    
    let evaluation = SelectivityEngine::fischer_esterification(&acid_desc, &alcohol_desc, context);
    
    if matches!(
        evaluation.recommendation,
        SelectivityRecommendation::Suppressed | SelectivityRecommendation::None
    ) {
        return Ok(None);
    }

    let base_ea = 25.0;
    let adjusted_ea = base_ea + evaluation.primary.activation_delta;

    let acid = acid_site.participant.substance;
    let acid_structure = acid_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let alcohol_structure = alcohol_site.participant.structure;
    let acid_carbon = acid_site.carbon;
    let acid_hydroxyl_oxygen = acid_site.hydroxyl_oxygen;
    let acid_proton = acid_site.hydroxyl_hydrogen;
    let alcohol_oxygen = alcohol_site.oxygen;
    let alcohol_proton = alcohol_site.hydrogen;

    let mut acid_editor = MolecularEditor::new(acid_structure);
    let acid_mapping = acid_editor.remove_atoms(&[acid_proton, acid_hydroxyl_oxygen])?;
    let acid_carbon = mapped_atom(&acid_mapping, acid_carbon, "acid carbon")?;
    let acid_fragment = acid_editor.finish()?;

    let mut alcohol_editor = MolecularEditor::new(alcohol_structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol_proton])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol_oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let product_structure = MolecularEditor::join_structures(
        &acid_fragment,
        acid_carbon,
        &alcohol_fragment,
        alcohol_oxygen,
        1.0,
    )?;
    let product = resolver.resolve(product_structure)?;
    Ok(Some(Reaction::builder(generated_pair_site_reaction_id(
        "carboxylic_acid_esterification",
        &acid_site.participant,
        &alcohol_site.participant,
    ))
    .reactant(acid.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 1, 0)
    .catalyst_order("destroy:sulfuric_acid", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .activation_energy_kj_per_mol(adjusted_ea)
    .build()))
}

pub(crate) fn generate_acyl_chloride_formation(
    site: &CarboxylicAcidSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let hydroxyl_oxygen = site.hydroxyl_oxygen;
    let proton = site.hydroxyl_hydrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydroxyl_oxygen, proton])?;
    let carbon = mapped_atom(&mapping, carbon, "carboxylic acid carbon")?;
    editor.add_atom(carbon, "Cl", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "acyl_chloride_formation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:phosgene", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .product("destroy:carbon_dioxide", 1)
    .build())
}

pub(crate) fn generate_acyl_chloride_hydrolysis(
    site: &AcylChlorideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let chlorine = site.chlorine;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[chlorine])?;
    let carbon = mapped_atom(&mapping, carbon, "acyl chloride carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "acyl_chloride_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
}

pub(crate) fn generate_acyl_chloride_esterification(
    acyl_chloride_site: &AcylChlorideSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let acyl_chloride = acyl_chloride_site.participant.substance;
    let acyl_chloride_structure = acyl_chloride_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let alcohol_structure = alcohol_site.participant.structure;
    let acyl_carbon = acyl_chloride_site.carbon;
    let chlorine = acyl_chloride_site.chlorine;
    let alcohol_oxygen = alcohol_site.oxygen;
    let alcohol_proton = alcohol_site.hydrogen;
    let mut acyl_editor = MolecularEditor::new(acyl_chloride_structure);
    let acyl_mapping = acyl_editor.remove_atoms(&[chlorine])?;
    let acyl_carbon = mapped_atom(&acyl_mapping, acyl_carbon, "acyl chloride carbon")?;
    let acyl_fragment = acyl_editor.finish()?;

    let mut alcohol_editor = MolecularEditor::new(alcohol_structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol_proton])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol_oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &acyl_fragment,
        acyl_carbon,
        &alcohol_fragment,
        alcohol_oxygen,
        1.0,
    )?)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "acyl_chloride_esterification",
        &acyl_chloride_site.participant,
        &alcohol_site.participant,
    ))
    .reactant(acyl_chloride.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
}

pub(crate) fn generate_amide_hydrolysis(
    site: &AmideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let nitrogen = site.nitrogen;
    let hydrogens = &site.nitrogen_hydrogens;
    if hydrogens.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("amide_hydrolysis", &site.participant),
            reason: "unsubstituted amide must have exactly two explicit nitrogen hydrogens"
                .to_string(),
        });
    }
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen, hydrogens[0], hydrogens[1]])?;
    let carbon = mapped_atom(&mapping, carbon, "amide carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "amide_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(product, 1)
    .product("destroy:ammonia", 1)
    .build())
}
