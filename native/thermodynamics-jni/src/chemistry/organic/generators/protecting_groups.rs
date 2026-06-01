//! Generators for protecting group reactions
//!
//! This module implements chemical reactions for protection and deprotection
//! of functional groups used in organic synthesis.
//!
//! NOTE: This is a minimal implementation focused on silyl ether (TMS) chemistry.
//! Other protecting groups (Boc, Cbz, acetals) are not yet fully implemented.

use super::super::centers::{AlcoholSite, SilylEtherCenter};
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::ReactionCondition;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor};
use crate::chemistry::reaction::Reaction;

/// Add TMS (trimethylsilyl) protecting group to an alcohol
///
/// Reaction: R-OH + Me3SiCl -> R-OSiMe3 + HCl
/// Note: In real synthesis, a base is used to scavenge HCl, but for simplicity
/// we model this as the direct reaction.
pub(crate) fn generate_alcohol_silyl_protection(
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = alcohol_site.participant.substance;
    let structure = alcohol_site.participant.structure;
    let oxygen = alcohol_site.oxygen;
    let hydrogen = alcohol_site.hydrogen;

    let mut editor = MolecularEditor::new(structure);
    // Remove the hydroxyl hydrogen
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let oxygen = mapped_atom(&mapping, oxygen, "alcohol oxygen")?;

    // Add silicon atom attached to oxygen
    let silicon = editor.add_atom(oxygen, "Si", 0.0, 1.0)?;

    // Add three methyl groups to silicon
    for _ in 0..3 {
        let methyl_carbon = editor.add_atom(silicon, "C", 0.0, 1.0)?;
        editor.add_atom(methyl_carbon, "H", 0.0, 1.0)?;
        editor.add_atom(methyl_carbon, "H", 0.0, 1.0)?;
        editor.add_atom(methyl_carbon, "H", 0.0, 1.0)?;
    }

    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "alcohol_silyl_protection",
        &alcohol_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:trimethylsilyl_chloride", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .condition(
        ReactionCondition::new("silyl protection requires dry conditions")
            .max_water_activity(0.1),
    )
    .activation_energy_kj_per_mol(15.0)
    .build())
}

/// Remove TMS protecting group from a silyl ether
///
/// Reaction: R-OSiMe3 + F- -> R-OH + Me3SiF
/// Fluoride cleavage is the standard method for removing TMS groups.
pub(crate) fn generate_silyl_ether_deprotection(
    silyl_site: &SilylEtherCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = silyl_site.participant.substance;
    let structure = silyl_site.participant.structure;
    let oxygen = silyl_site.oxygen;
    let silicon = silyl_site.silicon;

    // Find and remove the silicon with its three methyl groups
    let mut editor = MolecularEditor::new(structure);

    // Collect atoms to remove: silicon and its methyl groups
    let mut atoms_to_remove = vec![silicon];
    for (neighbor, order) in structure.neighbors(silicon) {
        if structure.atoms[neighbor].element == "C" && bond_order_matches(order, 1.0) {
            atoms_to_remove.push(neighbor);
            // Also remove hydrogens on these methyl carbons
            for (h_neighbor, h_order) in structure.neighbors(neighbor) {
                if structure.atoms[h_neighbor].element == "H" && bond_order_matches(h_order, 1.0) {
                    atoms_to_remove.push(h_neighbor);
                }
            }
        }
    }

    let mapping = editor.remove_atoms(&atoms_to_remove)?;
    let oxygen = mapped_atom(&mapping, oxygen, "silyl ether oxygen")?;

    // Add hydrogen to the oxygen
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;

    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "silyl_ether_deprotection",
        &silyl_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:fluoride", 1, 1) // F- from TBAF or similar
    .product(product, 1)
    .product("destroy:trimethylsilyl_fluoride", 1)
    .condition(
        ReactionCondition::new("fluoride deprotection requires fluoride source"),
    )
    .activation_energy_kj_per_mol(20.0)
    .build())
}

// TODO: Implement other protecting groups:
// - Acetal/ketal formation and hydrolysis
// - Boc protection and deprotection
// - Cbz protection and hydrogenolysis
// - Ester protection and hydrolysis
