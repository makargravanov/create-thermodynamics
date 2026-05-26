use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
pub(crate) fn generate_epoxide_hydrolysis(
    epoxide: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let oxygen = epoxide
        .site
        .atoms
        .iter()
        .copied()
        .find(|atom| epoxide.structure.atoms[*atom].element == "O")
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("epoxide_hydrolysis", &epoxide),
            reason: "epoxide site has no oxygen atom".to_string(),
        })?;
    let carbons = epoxide
        .site
        .atoms
        .iter()
        .copied()
        .filter(|atom| epoxide.structure.atoms[*atom].element == "C")
        .collect::<Vec<_>>();
    if carbons.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("epoxide_hydrolysis", &epoxide),
            reason: "epoxide site must contain exactly two carbon atoms".to_string(),
        });
    }
    let mut editor = MolecularEditor::new(epoxide.structure);
    editor.remove_bond(oxygen, carbons[0])?;
    add_hydroxyl(&mut editor, carbons[0])?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(
        Reaction::builder(generated_site_reaction_id("epoxide_hydrolysis", &epoxide))
            .reactant(epoxide.substance.id.clone(), 1, 1)
            .reactant("destroy:water", 1, 1)
            .catalyst_order("destroy:proton", 1)
            .product(product, 1)
            .condition(
                ReactionCondition::new("epoxide hydrolysis requires aqueous acid")
                    .acidity(AcidityCondition::Acidic)
                    .min_water_activity(0.1),
            )
            .build(),
    )
}

pub(crate) fn deprotonated_alcohol_fragment(
    site: &AlcoholSite<'_>,
    _role: &str,
) -> ChemistryResult<(MolecularStructure, usize)> {
    let structure = site.participant.structure;
    let oxygen = site.oxygen;
    let hydrogen = site.hydrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let oxygen = mapped_atom(&mapping, oxygen, "alcohol oxygen")?;
    Ok((editor.finish()?, oxygen))
}

pub(crate) fn generate_nitrile_hydrolysis(
    site: &NitrileSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let nitrogen = site.nitrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "nitrile carbon")?;
    editor.add_atom(carbon, "O", 0.0, 2.0)?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "nitrile_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(product, 1)
    .build())
}

pub(crate) fn generate_nitro_hydrogenation(
    site: &NitroSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let [first_oxygen, second_oxygen] = site.oxygens;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[first_oxygen, second_oxygen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "nitro nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "nitro_hydrogenation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen", 3, 1)
    .product("destroy:water", 2)
    .product(product, 1)
    .external_catalyst("forge:dusts/palladium", 1.0)
    .build())
}

pub(crate) fn generate_amine_phosgenation(
    site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let hydrogens = &site.hydrogens;
    if hydrogens.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("amine_phosgenation", &site.participant),
            reason: "primary amine must have exactly two explicit hydrogens".to_string(),
        });
    }
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogens[0], hydrogens[1]])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "amine nitrogen")?;
    let carbon = editor.add_atom(nitrogen, "C", 0.0, 2.0)?;
    editor.add_atom(carbon, "O", 0.0, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "amine_phosgenation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:phosgene", 1, 1)
    .product("destroy:hydrochloric_acid", 2)
    .product(product, 1)
    .build())
}

pub(crate) fn generate_cyanamide_addition(
    site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let hydrogen = *site
        .hydrogens
        .first()
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("cyanamide_addition", &site.participant),
            reason: "amine has no explicit hydrogen".to_string(),
        })?;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "amine nitrogen")?;
    let carbon = editor.add_atom(nitrogen, "C", 0.0, 1.0)?;
    let imine_nitrogen = editor.add_atom(carbon, "N", 0.0, 2.0)?;
    editor.add_atom(imine_nitrogen, "H", 0.0, 1.0)?;
    let amine_nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(amine_nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(amine_nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "cyanamide_addition",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:cyanamide", 1, 1)
    .product(product, 1)
    .build())
}

pub(crate) fn generate_isocyanate_hydrolysis(
    site: &IsocyanateSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let functional_carbon = site.functional_carbon;
    let oxygen = site.oxygen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[functional_carbon, oxygen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "isocyanate nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "isocyanate_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product("destroy:carbon_dioxide", 1)
    .product(product, 1)
    .build())
}

pub(crate) fn generate_nitrile_hydrogenation(
    site: &NitrileSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let nitrogen = site.nitrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "nitrile carbon")?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "nitrile_hydrogenation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen", 2, 1)
    .product(product, 1)
    .external_catalyst("forge:dusts/nickel", 1.0)
    .build())
}
