use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{bond_order_matches, MolecularStructure};
use crate::chemistry::reactive_site::ReactiveSiteKind;

use super::space::SiteParticipant;

#[derive(Clone)]
pub(crate) enum TypedReactionCenter<'a> {
    AlkylHalide(AlkylHalideCenter<'a>),
    Alcohol(AlcoholCenter<'a>),
    Alkoxide(AlkoxideCenter<'a>),
    Carbonyl(CarbonylCenter<'a>),
    CarboxylicAcid(CarboxylicAcidCenter<'a>),
    AcylChloride(AcylChlorideCenter<'a>),
    Amide(AmideCenter<'a>),
    Amine(AmineCenter<'a>),
    Nitrile(NitrileCenter<'a>),
    Nitro(NitroCenter<'a>),
    UnsaturatedBond(UnsaturatedBondCenter<'a>),
    Borane(BoraneCenter<'a>),
    BorateEster(BorateEsterCenter<'a>),
    Isocyanate(IsocyanateCenter<'a>),
    Generic(SiteParticipant<'a>),
}

#[derive(Clone)]
pub(crate) struct AlkylHalideCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) halogen: usize,
    pub(crate) degree: usize,
}

#[derive(Clone)]
pub(crate) struct AlcoholCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) oxygen: usize,
    pub(crate) hydrogen: usize,
    pub(crate) degree: usize,
}

#[derive(Clone)]
pub(crate) struct AlkoxideCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) oxygen: usize,
}

#[derive(Clone)]
pub(crate) struct CarbonylCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) oxygen: usize,
    pub(crate) is_ketone: bool,
}

#[derive(Clone)]
pub(crate) struct CarboxylicAcidCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) hydroxyl_oxygen: usize,
    pub(crate) hydroxyl_hydrogen: usize,
}

#[derive(Clone)]
pub(crate) struct AcylChlorideCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) chlorine: usize,
}

#[derive(Clone)]
pub(crate) struct AmideCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) nitrogen: usize,
    pub(crate) nitrogen_hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct AmineCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) hydrogens: Vec<usize>,
}

#[derive(Clone)]
pub(crate) struct NitrileCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) nitrogen: usize,
}

#[derive(Clone)]
pub(crate) struct NitroCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) oxygens: [usize; 2],
}

#[derive(Clone)]
pub(crate) struct UnsaturatedBondCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) high_degree_carbon: usize,
    pub(crate) low_degree_carbon: usize,
    pub(crate) is_alkyne: bool,
}

#[derive(Clone)]
pub(crate) struct BoraneCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) carbon: usize,
    pub(crate) boron: usize,
}

#[derive(Clone)]
pub(crate) struct BorateEsterCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) oxygen: usize,
    pub(crate) boron: usize,
}

#[derive(Clone)]
pub(crate) struct IsocyanateCenter<'a> {
    pub(crate) participant: SiteParticipant<'a>,
    pub(crate) nitrogen: usize,
    pub(crate) functional_carbon: usize,
    pub(crate) oxygen: usize,
}

pub(crate) fn typed_center_from_site<'a>(
    participant: SiteParticipant<'a>,
) -> ChemistryResult<TypedReactionCenter<'a>> {
    match participant.site.kind {
        ReactiveSiteKind::Halide => Ok(TypedReactionCenter::AlkylHalide(alkyl_halide_center(
            participant,
        )?)),
        ReactiveSiteKind::Alcohol => Ok(TypedReactionCenter::Alcohol(alcohol_center(participant)?)),
        ReactiveSiteKind::Alkoxide => {
            Ok(TypedReactionCenter::Alkoxide(alkoxide_center(participant)?))
        }
        ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => {
            Ok(TypedReactionCenter::Carbonyl(carbonyl_center(participant)?))
        }
        ReactiveSiteKind::CarboxylicAcid => Ok(TypedReactionCenter::CarboxylicAcid(
            carboxylic_acid_center(participant)?,
        )),
        ReactiveSiteKind::AcylChloride => Ok(TypedReactionCenter::AcylChloride(
            acyl_chloride_center(participant)?,
        )),
        ReactiveSiteKind::Amide => Ok(TypedReactionCenter::Amide(amide_center(participant)?)),
        ReactiveSiteKind::PrimaryAmine | ReactiveSiteKind::NonTertiaryAmine => {
            Ok(TypedReactionCenter::Amine(amine_center(participant)?))
        }
        ReactiveSiteKind::Nitrile => Ok(TypedReactionCenter::Nitrile(nitrile_center(participant)?)),
        ReactiveSiteKind::Nitro => Ok(TypedReactionCenter::Nitro(nitro_center(participant)?)),
        ReactiveSiteKind::Alkene | ReactiveSiteKind::Alkyne => Ok(
            TypedReactionCenter::UnsaturatedBond(unsaturated_bond_center(participant)?),
        ),
        ReactiveSiteKind::Borane => Ok(TypedReactionCenter::Borane(borane_center(participant)?)),
        ReactiveSiteKind::BorateEster => Ok(TypedReactionCenter::BorateEster(borate_ester_center(
            participant,
        )?)),
        ReactiveSiteKind::Isocyanate => Ok(TypedReactionCenter::Isocyanate(isocyanate_center(
            participant,
        )?)),
        _ => Ok(TypedReactionCenter::Generic(participant)),
    }
}

fn alkyl_halide_center<'a>(
    participant: SiteParticipant<'a>,
) -> ChemistryResult<AlkylHalideCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Halide)?;
    let carbon = site_atom_by_element(&participant, "C", "halide carbon")?;
    let halogen = participant
        .site
        .leaving_atom
        .or_else(|| {
            participant.site.atoms.iter().copied().find(|atom| {
                matches!(
                    participant.structure.atoms[*atom].element.as_str(),
                    "F" | "Cl" | "Br" | "I"
                )
            })
        })
        .ok_or_else(|| site_error(&participant, "halide site has no supported halogen"))?;
    Ok(AlkylHalideCenter {
        degree: participant
            .site
            .substitution_degree
            .unwrap_or_else(|| participant.structure.carbon_degree(carbon)),
        participant,
        carbon,
        halogen,
    })
}

fn alcohol_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<AlcoholCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Alcohol)?;
    let oxygen = site_atom_by_element(&participant, "O", "alcohol oxygen")?;
    let carbon = bonded_site_atom(&participant, oxygen, "C", 1.0, "alcohol carbon")?;
    let hydrogen = first_bonded_hydrogen(participant.structure, oxygen)
        .ok_or_else(|| site_error(&participant, "alcohol oxygen has no explicit hydrogen"))?;
    Ok(AlcoholCenter {
        degree: participant
            .site
            .substitution_degree
            .unwrap_or_else(|| participant.structure.carbon_degree(carbon)),
        participant,
        carbon,
        oxygen,
        hydrogen,
    })
}

fn alkoxide_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<AlkoxideCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Alkoxide)?;
    let oxygen = site_atom_by_element(&participant, "O", "alkoxide oxygen")?;
    bonded_site_atom(&participant, oxygen, "C", 1.0, "alkoxide carbon")?;
    Ok(AlkoxideCenter {
        participant,
        oxygen,
    })
}

fn carbonyl_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<CarbonylCenter<'a>> {
    let (carbon, oxygen) = carbonyl_atoms_from_site(participant.structure, &participant)?;
    let carbon_neighbors = participant
        .structure
        .neighbors(carbon)
        .into_iter()
        .filter(|(neighbor, order)| {
            *neighbor != oxygen
                && participant.structure.atoms[*neighbor].element == "C"
                && bond_order_matches(*order, 1.0)
        })
        .count();
    Ok(CarbonylCenter {
        is_ketone: participant.site.kind == ReactiveSiteKind::Ketone
            || (participant.site.kind == ReactiveSiteKind::Carbonyl && carbon_neighbors >= 2),
        participant,
        carbon,
        oxygen,
    })
}

fn carboxylic_acid_center<'a>(
    participant: SiteParticipant<'a>,
) -> ChemistryResult<CarboxylicAcidCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::CarboxylicAcid)?;
    let (carbon, carbonyl_oxygen) = carbonyl_atoms_from_site(participant.structure, &participant)?;
    let hydroxyl_oxygen = participant
        .structure
        .neighbors(carbon)
        .into_iter()
        .find_map(|(neighbor, order)| {
            (neighbor != carbonyl_oxygen
                && participant.structure.atoms[neighbor].element == "O"
                && bond_order_matches(order, 1.0))
            .then_some(neighbor)
        })
        .ok_or_else(|| site_error(&participant, "carboxylic acid has no hydroxyl oxygen"))?;
    let hydroxyl_hydrogen = first_bonded_hydrogen(participant.structure, hydroxyl_oxygen)
        .ok_or_else(|| {
            site_error(
                &participant,
                "carboxylic acid has no explicit hydroxyl hydrogen",
            )
        })?;
    Ok(CarboxylicAcidCenter {
        participant,
        carbon,
        hydroxyl_oxygen,
        hydroxyl_hydrogen,
    })
}

fn acyl_chloride_center<'a>(
    participant: SiteParticipant<'a>,
) -> ChemistryResult<AcylChlorideCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::AcylChloride)?;
    let (carbon, _) = carbonyl_atoms_from_site(participant.structure, &participant)?;
    let chlorine = bonded_site_atom(&participant, carbon, "Cl", 1.0, "acyl chloride chlorine")?;
    Ok(AcylChlorideCenter {
        participant,
        carbon,
        chlorine,
    })
}

fn amide_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<AmideCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Amide)?;
    let (carbon, oxygen) = carbonyl_atoms_from_site(participant.structure, &participant)?;
    let nitrogen = participant
        .structure
        .neighbors(carbon)
        .into_iter()
        .find_map(|(neighbor, order)| {
            (neighbor != oxygen
                && participant.structure.atoms[neighbor].element == "N"
                && bond_order_matches(order, 1.0))
            .then_some(neighbor)
        })
        .ok_or_else(|| site_error(&participant, "amide has no nitrogen atom"))?;
    let nitrogen_hydrogens = bonded_hydrogens(participant.structure, nitrogen);
    Ok(AmideCenter {
        participant,
        carbon,
        nitrogen,
        nitrogen_hydrogens,
    })
}

fn amine_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<AmineCenter<'a>> {
    let nitrogen = site_atom_by_element(&participant, "N", "amine nitrogen")?;
    let hydrogens = bonded_hydrogens(participant.structure, nitrogen);
    if hydrogens.is_empty() {
        return Err(site_error(
            &participant,
            "amine has no explicit nitrogen hydrogen",
        ));
    }
    Ok(AmineCenter {
        participant,
        nitrogen,
        hydrogens,
    })
}

fn nitrile_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<NitrileCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Nitrile)?;
    let carbon = site_atom_by_element(&participant, "C", "nitrile carbon")?;
    let nitrogen = bonded_site_atom(&participant, carbon, "N", 3.0, "nitrile nitrogen")?;
    Ok(NitrileCenter {
        participant,
        carbon,
        nitrogen,
    })
}

fn nitro_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<NitroCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Nitro)?;
    let nitrogen = site_atom_by_element(&participant, "N", "nitro nitrogen")?;
    let oxygens = participant
        .site
        .atoms
        .iter()
        .copied()
        .filter(|atom| participant.structure.atoms[*atom].element == "O")
        .collect::<Vec<_>>();
    let oxygens: [usize; 2] = oxygens.try_into().map_err(|_| {
        site_error(
            &participant,
            "nitro center must contain exactly two oxygens",
        )
    })?;
    Ok(NitroCenter {
        participant,
        nitrogen,
        oxygens,
    })
}

fn unsaturated_bond_center<'a>(
    participant: SiteParticipant<'a>,
) -> ChemistryResult<UnsaturatedBondCenter<'a>> {
    let is_alkyne = match participant.site.kind {
        ReactiveSiteKind::Alkene => false,
        ReactiveSiteKind::Alkyne => true,
        _ => return Err(site_error(&participant, "site is not an unsaturated bond")),
    };
    let carbons = participant
        .site
        .atoms
        .iter()
        .copied()
        .filter(|atom| participant.structure.atoms[*atom].element == "C")
        .collect::<Vec<_>>();
    if carbons.len() != 2 {
        return Err(site_error(
            &participant,
            "unsaturated bond must contain exactly two carbons",
        ));
    }
    let first_degree = participant
        .structure
        .carbon_degree(carbons[0])
        .saturating_sub(1);
    let second_degree = participant
        .structure
        .carbon_degree(carbons[1])
        .saturating_sub(1);
    let (high_degree_carbon, low_degree_carbon) = if second_degree > first_degree {
        (carbons[1], carbons[0])
    } else {
        (carbons[0], carbons[1])
    };
    Ok(UnsaturatedBondCenter {
        participant,
        high_degree_carbon,
        low_degree_carbon,
        is_alkyne,
    })
}

fn borane_center<'a>(participant: SiteParticipant<'a>) -> ChemistryResult<BoraneCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Borane)?;
    let carbon = site_atom_by_element(&participant, "C", "borane carbon")?;
    let boron = bonded_site_atom(&participant, carbon, "B", 1.0, "borane boron")?;
    Ok(BoraneCenter {
        participant,
        carbon,
        boron,
    })
}

fn borate_ester_center<'a>(
    participant: SiteParticipant<'a>,
) -> ChemistryResult<BorateEsterCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::BorateEster)?;
    let oxygen = site_atom_by_element(&participant, "O", "borate ester oxygen")?;
    let boron = bonded_site_atom(&participant, oxygen, "B", 1.0, "borate ester boron")?;
    Ok(BorateEsterCenter {
        participant,
        oxygen,
        boron,
    })
}

fn isocyanate_center<'a>(
    participant: SiteParticipant<'a>,
) -> ChemistryResult<IsocyanateCenter<'a>> {
    require_kind(&participant, ReactiveSiteKind::Isocyanate)?;
    let nitrogen = site_atom_by_element(&participant, "N", "isocyanate nitrogen")?;
    let functional_carbon =
        bonded_site_atom(&participant, nitrogen, "C", 2.0, "isocyanate carbon")?;
    let oxygen = bonded_site_atom(
        &participant,
        functional_carbon,
        "O",
        2.0,
        "isocyanate oxygen",
    )?;
    Ok(IsocyanateCenter {
        participant,
        nitrogen,
        functional_carbon,
        oxygen,
    })
}

fn require_kind(
    participant: &SiteParticipant<'_>,
    expected: ReactiveSiteKind,
) -> ChemistryResult<()> {
    if participant.site.kind == expected {
        Ok(())
    } else {
        Err(site_error(
            participant,
            &format!(
                "expected {:?} reactive site, got {:?}",
                expected, participant.site.kind
            ),
        ))
    }
}

fn site_atom_by_element(
    participant: &SiteParticipant<'_>,
    element: &str,
    label: &str,
) -> ChemistryResult<usize> {
    participant
        .site
        .atoms
        .iter()
        .copied()
        .find(|atom| participant.structure.atoms[*atom].element == element)
        .ok_or_else(|| site_error(participant, &format!("reactive site is missing {label}")))
}

fn bonded_site_atom(
    participant: &SiteParticipant<'_>,
    parent: usize,
    element: &str,
    order: f64,
    label: &str,
) -> ChemistryResult<usize> {
    participant
        .structure
        .neighbors(parent)
        .into_iter()
        .find_map(|(neighbor, bond_order)| {
            (participant.site.atoms.contains(&neighbor)
                && participant.structure.atoms[neighbor].element == element
                && bond_order_matches(bond_order, order))
            .then_some(neighbor)
        })
        .ok_or_else(|| site_error(participant, &format!("reactive site is missing {label}")))
}

fn carbonyl_atoms_from_site(
    structure: &MolecularStructure,
    participant: &SiteParticipant<'_>,
) -> ChemistryResult<(usize, usize)> {
    for atom in &participant.site.atoms {
        if structure.atoms[*atom].element != "C" {
            continue;
        }
        for (neighbor, order) in structure.neighbors(*atom) {
            if participant.site.atoms.contains(&neighbor)
                && structure.atoms[neighbor].element == "O"
                && bond_order_matches(order, 2.0)
            {
                return Ok((*atom, neighbor));
            }
        }
    }
    Err(site_error(participant, "site has no carbonyl C=O bond"))
}

fn bonded_hydrogens(structure: &MolecularStructure, parent: usize) -> Vec<usize> {
    structure
        .neighbors(parent)
        .into_iter()
        .filter_map(|(neighbor, order)| {
            (structure.atoms[neighbor].element == "H" && bond_order_matches(order, 1.0))
                .then_some(neighbor)
        })
        .collect()
}

fn first_bonded_hydrogen(structure: &MolecularStructure, atom: usize) -> Option<usize> {
    bonded_hydrogens(structure, atom).into_iter().next()
}

fn site_error(participant: &SiteParticipant<'_>, reason: &str) -> ChemistryError {
    ChemistryError::InvalidReaction {
        reaction_id: format!(
            "typed_center/{}/{}",
            participant.substance.id.as_str().replace([':', '/'], "_"),
            participant
                .site
                .atoms
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join("_")
        ),
        reason: reason.to_string(),
    }
}
