use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::ReactionCondition;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor};
use crate::chemistry::reaction::Reaction;

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

pub(crate) fn generate_isocyanate_ammonolysis(
    site: &IsocyanateSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let mut editor = MolecularEditor::new(site.participant.structure);
    editor.set_bond_order(site.nitrogen, site.functional_carbon, 1.0)?;
    editor.add_atom(site.nitrogen, "H", 0.0, 1.0)?;
    let ammonia_nitrogen = editor.add_atom(site.functional_carbon, "N", 0.0, 1.0)?;
    editor.add_atom(ammonia_nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(ammonia_nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "isocyanate_ammonolysis",
        &site.participant,
    ))
    .reactant(site.participant.substance.id.clone(), 1, 1)
    .reactant("destroy:ammonia", 1, 1)
    .product(product, 1)
    .condition(
        ReactionCondition::new("isocyanate ammonolysis requires a dry medium")
            .max_water_activity(0.1),
    )
    .activation_energy_kj_per_mol(20.0)
    .build())
}

pub(crate) fn generate_isocyanate_amine_addition(
    isocyanate_site: &IsocyanateSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if isocyanate_site.participant.substance.id == amine_site.participant.substance.id {
        return Ok(None);
    }
    let Some(amine_hydrogen) = amine_site.hydrogens.first().copied() else {
        return Ok(None);
    };

    let mut isocyanate_editor = MolecularEditor::new(isocyanate_site.participant.structure);
    isocyanate_editor.set_bond_order(
        isocyanate_site.nitrogen,
        isocyanate_site.functional_carbon,
        1.0,
    )?;
    isocyanate_editor.add_atom(isocyanate_site.nitrogen, "H", 0.0, 1.0)?;
    let isocyanate_fragment = isocyanate_editor.finish()?;

    let mut amine_editor = MolecularEditor::new(amine_site.participant.structure);
    let amine_mapping = amine_editor.remove_atoms(&[amine_hydrogen])?;
    let amine_nitrogen = mapped_atom(&amine_mapping, amine_site.nitrogen, "amine nitrogen")?;
    let amine_fragment = amine_editor.finish()?;

    let mut product_editor = MolecularEditor::new(&isocyanate_fragment);
    product_editor.add_group(
        isocyanate_site.functional_carbon,
        &amine_fragment,
        amine_nitrogen,
        1.0,
    )?;
    let product = resolver.resolve(product_editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "isocyanate_amine_addition",
            &isocyanate_site.participant,
            &amine_site.participant,
        ))
        .reactant(isocyanate_site.participant.substance.id.clone(), 1, 1)
        .reactant(amine_site.participant.substance.id.clone(), 1, 1)
        .product(product, 1)
        .condition(
            ReactionCondition::new("isocyanate amine addition requires a dry medium")
                .max_water_activity(0.1),
        )
        .activation_energy_kj_per_mol(18.0)
        .build(),
    ))
}

pub(crate) fn generate_amine_formylation(
    amine_site: &AmineSite<'_>,
    donor_site: &FormylationDonorCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if amine_site.participant.substance.id == donor_site.participant.substance.id {
        return Ok(None);
    }
    let Some(amine_hydrogen) = amine_site.hydrogens.first().copied() else {
        return Ok(None);
    };
    let leaving_atom = donor_site
        .participant
        .structure
        .neighbors(donor_site.carbon)
        .into_iter()
        .find_map(|(neighbor, order)| {
            (neighbor != donor_site.oxygen
                && neighbor != donor_site.hydrogen
                && bond_order_matches(order, 1.0)
                && matches!(
                    donor_site.participant.structure.atoms[neighbor]
                        .element
                        .as_str(),
                    "O" | "N" | "Cl"
                ))
            .then_some(neighbor)
        });
    let Some(leaving_atom) = leaving_atom else {
        return Ok(None);
    };

    let mut amine_editor = MolecularEditor::new(amine_site.participant.structure);
    let amine_mapping = amine_editor.remove_atoms(&[amine_hydrogen])?;
    let nitrogen = mapped_atom(
        &amine_mapping,
        amine_site.nitrogen,
        "formylated amine nitrogen",
    )?;
    let formyl_carbon = amine_editor.add_atom(nitrogen, "C", 0.0, 1.0)?;
    amine_editor.add_atom(formyl_carbon, "O", 0.0, 2.0)?;
    amine_editor.add_atom(formyl_carbon, "H", 0.0, 1.0)?;
    let formamide = resolver.resolve(amine_editor.finish()?)?;

    let leaving_product = if donor_site.participant.structure.atoms[leaving_atom].element == "Cl" {
        "destroy:hydrochloric_acid".into()
    } else {
        let mut donor_editor = MolecularEditor::new(donor_site.participant.structure);
        let donor_mapping = donor_editor.remove_atoms(&[
            donor_site.carbon,
            donor_site.oxygen,
            donor_site.hydrogen,
        ])?;
        let leaving_atom = mapped_atom(&donor_mapping, leaving_atom, "formyl leaving atom")?;
        donor_editor.add_atom(leaving_atom, "H", 0.0, 1.0)?;
        resolver.resolve(donor_editor.finish()?)?
    };

    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "amine_formylation",
            &amine_site.participant,
            &donor_site.participant,
        ))
        .reactant(amine_site.participant.substance.id.clone(), 1, 1)
        .reactant(donor_site.participant.substance.id.clone(), 1, 1)
        .product(formamide, 1)
        .product(leaving_product, 1)
        .condition(
            ReactionCondition::new("amine formylation requires removal of excess water")
                .max_water_activity(0.2),
        )
        .activation_energy_kj_per_mol(35.0)
        .build(),
    ))
}
