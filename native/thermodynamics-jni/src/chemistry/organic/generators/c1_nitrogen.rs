use std::collections::BTreeSet;

use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::ReactionCondition;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::substance::SubstanceId;

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
    let phosgene = known_reagent_structure(resolver, "destroy:phosgene", "amine_phosgenation")?;
    let (phosgene_fragment, phosgene_mapping, phosgene_carbon) =
        phosgene_carbonyl_fragment(phosgene)?;
    let phosgene_carbon = mapped_atom(
        &phosgene_mapping,
        phosgene_carbon,
        "phosgene carbonyl carbon",
    )?;

    let mut amine_editor = MolecularEditor::new(structure);
    let amine_mapping = amine_editor.remove_atoms(&[hydrogens[0], hydrogens[1]])?;
    let nitrogen = mapped_atom(&amine_mapping, nitrogen, "amine nitrogen")?;
    let product = resolver.resolve(MolecularEditor::join_structures(
        &amine_editor.finish()?,
        nitrogen,
        &phosgene_fragment,
        phosgene_carbon,
        2.0,
    )?)?;
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
    let cyanamide = known_reagent_structure(resolver, "destroy:cyanamide", "cyanamide_addition")?;
    let (cyanamide_carbon, imine_nitrogen, _) = cyanamide_atoms(cyanamide)?;
    let amine_hydrogen = atom_fragment(structure, hydrogen)?;

    let mut cyanamide_editor = MolecularEditor::new(cyanamide);
    cyanamide_editor.set_bond_order(cyanamide_carbon, imine_nitrogen, 2.0)?;
    cyanamide_editor.add_group(imine_nitrogen, &amine_hydrogen, 0, 1.0)?;
    let cyanamide_fragment = cyanamide_editor.finish()?;

    let mut amine_editor = MolecularEditor::new(structure);
    let amine_mapping = amine_editor.remove_atoms(&[hydrogen])?;
    let nitrogen = mapped_atom(&amine_mapping, nitrogen, "amine nitrogen")?;
    let product = resolver.resolve(MolecularEditor::join_structures(
        &amine_editor.finish()?,
        nitrogen,
        &cyanamide_fragment,
        cyanamide_carbon,
        1.0,
    )?)?;
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
    let water = known_reagent_structure(resolver, "destroy:water", "isocyanate_hydrolysis")?;
    let water_hydrogens = water_hydrogens(water)?;
    let first_hydrogen = atom_fragment(water, water_hydrogens[0])?;
    let second_hydrogen = atom_fragment(water, water_hydrogens[1])?;

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[functional_carbon, oxygen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "isocyanate nitrogen")?;
    editor.add_group(nitrogen, &first_hydrogen, 0, 1.0)?;
    editor.add_group(nitrogen, &second_hydrogen, 0, 1.0)?;
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
    let ammonia = known_reagent_structure(resolver, "destroy:ammonia", "isocyanate_ammonolysis")?;
    let (ammonia_fragment, ammonia_mapping, ammonia_nitrogen, transferred_hydrogen) =
        ammonia_nh2_fragment(ammonia)?;
    let ammonia_nitrogen = mapped_atom(&ammonia_mapping, ammonia_nitrogen, "ammonia nitrogen")?;
    let transferred_hydrogen = atom_fragment(ammonia, transferred_hydrogen)?;

    let mut editor = MolecularEditor::new(site.participant.structure);
    editor.set_bond_order(site.nitrogen, site.functional_carbon, 1.0)?;
    editor.add_group(site.nitrogen, &transferred_hydrogen, 0, 1.0)?;
    editor.add_group(
        site.functional_carbon,
        &ammonia_fragment,
        ammonia_nitrogen,
        1.0,
    )?;
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

    let amine_hydrogen_fragment = atom_fragment(amine_site.participant.structure, amine_hydrogen)?;
    let leaving_branch = leaving_branch(
        donor_site.participant.structure,
        leaving_atom,
        donor_site.carbon,
        "amine_formylation",
    )?;
    let (formyl_fragment, formyl_mapping) =
        fragment_without_atoms(donor_site.participant.structure, &leaving_branch)?;
    let formyl_carbon = mapped_atom(&formyl_mapping, donor_site.carbon, "formyl carbon")?;

    let mut amine_editor = MolecularEditor::new(amine_site.participant.structure);
    let amine_mapping = amine_editor.remove_atoms(&[amine_hydrogen])?;
    let nitrogen = mapped_atom(
        &amine_mapping,
        amine_site.nitrogen,
        "formylated amine nitrogen",
    )?;
    let formamide = resolver.resolve(MolecularEditor::join_structures(
        &amine_editor.finish()?,
        nitrogen,
        &formyl_fragment,
        formyl_carbon,
        1.0,
    )?)?;

    let leaving_product = if donor_site.participant.structure.atoms[leaving_atom].element == "Cl" {
        "destroy:hydrochloric_acid".into()
    } else {
        let (leaving_fragment, leaving_mapping) =
            fragment_from_atoms(donor_site.participant.structure, &leaving_branch)?;
        let leaving_atom = mapped_atom(&leaving_mapping, leaving_atom, "formyl leaving atom")?;
        let mut leaving_editor = MolecularEditor::new(&leaving_fragment);
        leaving_editor.add_group(leaving_atom, &amine_hydrogen_fragment, 0, 1.0)?;
        resolver.resolve(leaving_editor.finish()?)?
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

fn known_reagent_structure<'a>(
    resolver: &'a DerivedSubstanceResolver,
    id: &str,
    reaction_id: &str,
) -> ChemistryResult<&'a MolecularStructure> {
    let substance_id = SubstanceId::from(id);
    resolver
        .known_structure(&substance_id)
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: format!("required reagent '{id}' has no known molecular structure"),
        })
}

fn phosgene_carbonyl_fragment(
    structure: &MolecularStructure,
) -> ChemistryResult<(MolecularStructure, Vec<Option<usize>>, usize)> {
    let (carbon, _) = carbonyl_in_structure(structure, "phosgene")?;
    let chlorines = structure
        .neighbors(carbon)
        .into_iter()
        .filter_map(|(neighbor, order)| {
            (structure.atoms[neighbor].element == "Cl" && bond_order_matches(order, 1.0))
                .then_some(neighbor)
        })
        .collect::<Vec<_>>();
    if chlorines.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "amine_phosgenation".to_string(),
            reason: "phosgene structure must have two chloride leaving atoms".to_string(),
        });
    }
    let (fragment, mapping) = fragment_without_atoms(structure, &chlorines)?;
    Ok((fragment, mapping, carbon))
}

fn cyanamide_atoms(structure: &MolecularStructure) -> ChemistryResult<(usize, usize, usize)> {
    for (carbon, atom) in structure.atoms.iter().enumerate() {
        if atom.element != "C" {
            continue;
        }
        let mut imine_nitrogen = None;
        let mut amino_nitrogen = None;
        for (neighbor, order) in structure.neighbors(carbon) {
            if structure.atoms[neighbor].element != "N" {
                continue;
            }
            if bond_order_matches(order, 3.0) {
                imine_nitrogen = Some(neighbor);
            } else if bond_order_matches(order, 1.0) {
                amino_nitrogen = Some(neighbor);
            }
        }
        if let (Some(imine), Some(amino)) = (imine_nitrogen, amino_nitrogen) {
            return Ok((carbon, imine, amino));
        }
    }
    Err(ChemistryError::InvalidReaction {
        reaction_id: "cyanamide_addition".to_string(),
        reason: "cyanamide structure must contain N#C-N connectivity".to_string(),
    })
}

fn water_hydrogens(structure: &MolecularStructure) -> ChemistryResult<[usize; 2]> {
    let oxygen = structure
        .atoms
        .iter()
        .position(|atom| atom.element == "O")
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: "isocyanate_hydrolysis".to_string(),
            reason: "water structure has no oxygen".to_string(),
        })?;
    let hydrogens = bonded_hydrogens(structure, oxygen);
    hydrogens
        .try_into()
        .map_err(|_| ChemistryError::InvalidReaction {
            reaction_id: "isocyanate_hydrolysis".to_string(),
            reason: "water structure must have exactly two explicit hydrogens".to_string(),
        })
}

fn ammonia_nh2_fragment(
    structure: &MolecularStructure,
) -> ChemistryResult<(MolecularStructure, Vec<Option<usize>>, usize, usize)> {
    let nitrogen = structure
        .atoms
        .iter()
        .position(|atom| atom.element == "N")
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: "isocyanate_ammonolysis".to_string(),
            reason: "ammonia structure has no nitrogen".to_string(),
        })?;
    let hydrogens = bonded_hydrogens(structure, nitrogen);
    if hydrogens.len() != 3 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "isocyanate_ammonolysis".to_string(),
            reason: "ammonia structure must have exactly three explicit hydrogens".to_string(),
        });
    }
    let transferred_hydrogen = hydrogens[0];
    let (fragment, mapping) = fragment_without_atoms(structure, &[transferred_hydrogen])?;
    Ok((fragment, mapping, nitrogen, transferred_hydrogen))
}

fn carbonyl_in_structure(
    structure: &MolecularStructure,
    reaction_id: &str,
) -> ChemistryResult<(usize, usize)> {
    for (carbon, atom) in structure.atoms.iter().enumerate() {
        if atom.element != "C" {
            continue;
        }
        if let Some((oxygen, _)) =
            structure
                .neighbors(carbon)
                .into_iter()
                .find(|(neighbor, order)| {
                    structure.atoms[*neighbor].element == "O" && bond_order_matches(*order, 2.0)
                })
        {
            return Ok((carbon, oxygen));
        }
    }
    Err(ChemistryError::InvalidReaction {
        reaction_id: reaction_id.to_string(),
        reason: "structure has no carbonyl group".to_string(),
    })
}

fn atom_fragment(
    structure: &MolecularStructure,
    atom: usize,
) -> ChemistryResult<MolecularStructure> {
    fragment_from_atoms(structure, &[atom]).map(|(fragment, _)| fragment)
}

fn fragment_without_atoms(
    structure: &MolecularStructure,
    removed_atoms: &[usize],
) -> ChemistryResult<(MolecularStructure, Vec<Option<usize>>)> {
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(removed_atoms)?;
    Ok((editor.finish()?, mapping))
}

fn fragment_from_atoms(
    structure: &MolecularStructure,
    kept_atoms: &[usize],
) -> ChemistryResult<(MolecularStructure, Vec<Option<usize>>)> {
    let kept = kept_atoms.iter().copied().collect::<BTreeSet<_>>();
    let removed = (0..structure.atoms.len())
        .filter(|atom| !kept.contains(atom))
        .collect::<Vec<_>>();
    fragment_without_atoms(structure, &removed)
}

fn leaving_branch(
    structure: &MolecularStructure,
    leaving_atom: usize,
    blocked_atom: usize,
    reaction_id: &str,
) -> ChemistryResult<Vec<usize>> {
    if leaving_atom >= structure.atoms.len() || blocked_atom >= structure.atoms.len() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: "leaving branch references an atom outside the structure".to_string(),
        });
    }
    let mut stack = vec![leaving_atom];
    let mut visited = vec![false; structure.atoms.len()];
    visited[blocked_atom] = true;
    while let Some(atom) = stack.pop() {
        if visited[atom] {
            continue;
        }
        visited[atom] = true;
        for (neighbor, _) in structure.neighbors(atom) {
            if !visited[neighbor] {
                stack.push(neighbor);
            }
        }
    }
    let branch = visited
        .into_iter()
        .enumerate()
        .filter_map(|(atom, seen)| (seen && atom != blocked_atom).then_some(atom))
        .collect::<Vec<_>>();
    if !branch.contains(&leaving_atom) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: "leaving branch does not contain the leaving atom".to_string(),
        });
    }
    Ok(branch)
}
