use super::super::space::SiteParticipant;
use super::*;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{MolecularEditor, MolecularStructure};
use crate::chemistry::reactive_site::{ReactiveSite, ReactiveSiteKind};
use crate::chemistry::substance::Substance;
pub(crate) fn add_hydroxyl(editor: &mut MolecularEditor, parent: usize) -> ChemistryResult<usize> {
    let oxygen = editor.add_atom(parent, "O", 0.0, 1.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    Ok(oxygen)
}

pub(crate) fn add_addition_group(
    editor: &mut MolecularEditor,
    parent: usize,
    group: AdditionGroup,
) -> ChemistryResult<()> {
    match group {
        AdditionGroup::Atom(element) => {
            editor.add_atom(parent, element, 0.0, 1.0)?;
        }
        AdditionGroup::Hydroxyl => {
            add_hydroxyl(editor, parent)?;
        }
        AdditionGroup::Borane => {
            let boron = editor.add_atom(parent, "B", 0.0, 1.0)?;
            editor.add_atom(boron, "H", 0.0, 1.0)?;
            editor.add_atom(boron, "H", 0.0, 1.0)?;
        }
    }
    Ok(())
}

pub(crate) fn bonded_hydrogens(structure: &MolecularStructure, parent: usize) -> Vec<usize> {
    structure
        .neighbors(parent)
        .into_iter()
        .map(|(neighbor, _)| neighbor)
        .filter(|neighbor| structure.atoms[*neighbor].element == "H")
        .collect()
}

pub(crate) fn halide_ion(
    structure: &MolecularStructure,
    halogen: usize,
    prefix: &str,
    participant: &SiteParticipant<'_>,
) -> ChemistryResult<&'static str> {
    match structure.atoms[halogen].element.as_str() {
        "Cl" => Ok("destroy:chloride"),
        "F" => Ok("destroy:fluoride"),
        "I" => Ok("destroy:iodide"),
        _ => Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id(prefix, participant),
            reason: "halide group does not contain a supported halogen".to_string(),
        }),
    }
}

pub(crate) fn carbonyl_atoms_from_site(
    structure: &MolecularStructure,
    site: &ReactiveSite,
    role: &str,
) -> ChemistryResult<(usize, usize)> {
    for carbon in site
        .atoms
        .iter()
        .copied()
        .filter(|atom| structure.atoms[*atom].element == "C")
    {
        if let Some((oxygen, _)) =
            structure
                .neighbors(carbon)
                .into_iter()
                .find(|(neighbor, order)| {
                    structure.atoms[*neighbor].element == "O"
                        && crate::chemistry::molecule::bond_order_matches(*order, 2.0)
                })
        {
            return Ok((carbon, oxygen));
        }
    }
    Err(ChemistryError::InvalidReaction {
        reaction_id: "<generated-organic>".to_string(),
        reason: format!("{role} site does not contain a carbonyl bond"),
    })
}

pub(crate) fn enol_atoms(
    structure: &MolecularStructure,
    site: &ReactiveSite,
) -> ChemistryResult<(usize, usize)> {
    let (carbonyl, _) = carbonyl_atoms_from_site(structure, site, "enol")?;
    let alpha = site
        .atoms
        .iter()
        .copied()
        .find(|atom| {
            *atom != carbonyl
                && structure.atoms[*atom].element == "C"
                && structure
                    .neighbors(carbonyl)
                    .iter()
                    .any(|(neighbor, order)| {
                        *neighbor == *atom
                            && crate::chemistry::molecule::bond_order_matches(*order, 1.0)
                    })
        })
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: "<generated-organic>".to_string(),
            reason: "enol site does not contain an alpha carbon".to_string(),
        })?;
    Ok((carbonyl, alpha))
}

pub(crate) fn organometallic_atoms(
    structure: &MolecularStructure,
    site: &ReactiveSite,
) -> ChemistryResult<(usize, usize, Vec<usize>)> {
    let mut organo_carbon = None;
    let mut metal = None;
    for atom in &site.atoms {
        match structure.atoms[*atom].element.as_str() {
            "C" => organo_carbon = Some(*atom),
            "Mg" | "Li" | "Cu" => metal = Some(*atom),
            _ => {}
        }
    }
    let organo_carbon = organo_carbon.ok_or_else(|| ChemistryError::InvalidReaction {
        reaction_id: "<generated-organic>".to_string(),
        reason: "organometallic site has no carbon atom".to_string(),
    })?;
    let metal = metal.ok_or_else(|| ChemistryError::InvalidReaction {
        reaction_id: "<generated-organic>".to_string(),
        reason: "organometallic site has no metal atom".to_string(),
    })?;
    let mut residue_atoms = vec![metal];
    for (neighbor, order) in structure.neighbors(metal) {
        if neighbor != organo_carbon && crate::chemistry::molecule::bond_order_matches(order, 1.0) {
            residue_atoms.push(neighbor);
        }
    }
    residue_atoms.sort_unstable();
    residue_atoms.dedup();
    Ok((organo_carbon, metal, residue_atoms))
}

pub(crate) fn atom_mass_sum(
    structure: &MolecularStructure,
    atoms: &[usize],
) -> ChemistryResult<f64> {
    atoms.iter().try_fold(0.0, |sum, atom| {
        Ok(sum + crate::chemistry::molecule::element_mass(&structure.atoms[*atom].element)?)
    })
}

pub(crate) fn atom_charge_sum(
    structure: &MolecularStructure,
    atoms: &[usize],
) -> ChemistryResult<i32> {
    let charge = atoms
        .iter()
        .map(|atom| structure.atoms[*atom].charge)
        .sum::<f64>();
    if (charge - charge.round()).abs() > 1.0e-9 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<generated-organic>".to_string(),
            reason: "external residue has non-integral charge".to_string(),
        });
    }
    Ok(charge.round() as i32)
}

pub(crate) fn first_bonded_hydrogen(structure: &MolecularStructure, atom: usize) -> Option<usize> {
    structure
        .neighbors(atom)
        .into_iter()
        .map(|(neighbor, _)| neighbor)
        .find(|neighbor| structure.atoms[*neighbor].element == "H")
}

pub(crate) fn mapped_atom(
    mapping: &[Option<usize>],
    old_index: usize,
    role: &str,
) -> ChemistryResult<usize> {
    mapping
        .get(old_index)
        .and_then(|value| *value)
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: "<generated-organic>".to_string(),
            reason: format!("{role} was removed during graph transformation"),
        })
}

pub(crate) fn generated_pair_reaction_id(
    prefix: &str,
    first: &Substance,
    second: &Substance,
) -> String {
    format!(
        "{prefix}/{}/{}",
        sanitize_id(first.id.as_str()),
        sanitize_id(second.id.as_str())
    )
}

pub(crate) fn generated_site_reaction_id(
    prefix: &str,
    participant: &SiteParticipant<'_>,
) -> String {
    format!(
        "{}/{}/{}",
        prefix,
        sanitize_id(participant.substance.id.as_str()),
        participant
            .site
            .atoms
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("_")
    )
}

pub(crate) fn generated_pair_site_reaction_id(
    prefix: &str,
    first: &SiteParticipant<'_>,
    second: &SiteParticipant<'_>,
) -> String {
    format!(
        "{}/{}/{}/{}",
        generated_pair_reaction_id(prefix, first.substance, second.substance),
        first
            .site
            .atoms
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("_"),
        second
            .site
            .atoms
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("_"),
        site_kind_suffix(&first.site.kind)
    )
}

pub(crate) fn site_kind_suffix(kind: &ReactiveSiteKind) -> &'static str {
    match kind {
        ReactiveSiteKind::AcidAnhydride => "acid_anhydride",
        ReactiveSiteKind::AcylChloride => "acyl_chloride",
        ReactiveSiteKind::Alcohol => "alcohol",
        ReactiveSiteKind::Alkene => "alkene",
        ReactiveSiteKind::Alkoxide => "alkoxide",
        ReactiveSiteKind::Alkyne => "alkyne",
        ReactiveSiteKind::Aldehyde => "aldehyde",
        ReactiveSiteKind::Amide => "amide",
        ReactiveSiteKind::AromaticCarbon => "aromatic_carbon",
        ReactiveSiteKind::AromaticRing => "aromatic_ring",
        ReactiveSiteKind::ArylHalide => "aryl_halide",
        ReactiveSiteKind::Azide => "azide",
        ReactiveSiteKind::Borane => "borane",
        ReactiveSiteKind::BoricAcid => "boric_acid",
        ReactiveSiteKind::BorateEster => "borate_ester",
        ReactiveSiteKind::Carbonyl => "carbonyl",
        ReactiveSiteKind::CarboxylicAcid => "carboxylic_acid",
        ReactiveSiteKind::Diazonium => "diazonium",
        ReactiveSiteKind::Enol => "enol",
        ReactiveSiteKind::Enolate => "enolate",
        ReactiveSiteKind::Epoxide => "epoxide",
        ReactiveSiteKind::Ester => "ester",
        ReactiveSiteKind::Ether => "ether",
        ReactiveSiteKind::Halide => "halide",
        ReactiveSiteKind::Imine => "imine",
        ReactiveSiteKind::Isocyanate => "isocyanate",
        ReactiveSiteKind::Ketone => "ketone",
        ReactiveSiteKind::Nitrile => "nitrile",
        ReactiveSiteKind::Nitro => "nitro",
        ReactiveSiteKind::NonTertiaryAmine => "non_tertiary_amine",
        ReactiveSiteKind::Organocopper => "organocopper",
        ReactiveSiteKind::Organolithium => "organolithium",
        ReactiveSiteKind::Organomagnesium => "organomagnesium",
        ReactiveSiteKind::Phenol => "phenol",
        ReactiveSiteKind::PrimaryAmine => "primary_amine",
        ReactiveSiteKind::Sulfide => "sulfide",
        ReactiveSiteKind::SulfonylChloride => "sulfonyl_chloride",
        ReactiveSiteKind::Thiol => "thiol",
    }
}

pub(crate) fn sanitize_id(value: &str) -> String {
    value
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect()
}
