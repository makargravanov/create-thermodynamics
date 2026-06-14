use std::collections::BTreeSet;

use super::error::{ChemistryError, ChemistryResult};
use super::functional_group::{find_functional_groups, FunctionalGroupType};
use super::molecule::{bond_order_matches, MolecularStructure};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReactiveSiteKind {
    AcidAnhydride,
    AcylChloride,
    Alcohol,
    Alkene,
    Alkoxide,
    AlkylHydrogen,
    Alkyne,
    Aldehyde,
    Amide,
    AmideNitrogen,
    AromaticCarbon,
    AromaticRing,
    ArylHalide,
    Azide,
    Borane,
    BoricAcid,
    BorateEster,
    Carbonyl,
    CarboxylicAcid,
    Diazonium,
    Enol,
    Enolate,
    Epoxide,
    Ester,
    Ether,
    Halide,
    Imine,
    Isocyanate,
    Ketone,
    Nitrile,
    Nitro,
    NonTertiaryAmine,
    NucleophilicPhosphorus,
    Organocopper,
    Organolithium,
    Organomagnesium,
    Oxime,
    Phenol,
    PrimaryAmine,
    Phosphine,
    PhosphonateCarbanion,
    PhosphoniumSalt,
    PhosphorusYlide,
    SulfoneCarbanion,
    Sulfide,
    Sulfoxide,
    SulfonylChloride,
    Thiol,
    SilylEther,
    Acetal,
    Ketal,
    BocCarbamate,
    CbzCarbamate,
    Hydrazone,
    BisNucleophile,
    DicarbonylElectrophile,
    UreaLike,
    FormylationDonor,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReactiveRole {
    AcidicProton,
    AlphaCarbon,
    AromaticDirector,
    Electrophile,
    LeavingGroup,
    Nucleophile,
    Oxidizable,
    Reducible,
    UnsaturatedBond,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReactiveSiteKey {
    pub kind: ReactiveSiteKind,
    pub atoms: Vec<usize>,
    pub roles: Vec<ReactiveRole>,
    pub primary_atom: Option<usize>,
    pub anchor_atoms: Vec<usize>,
    pub leaving_atom: Option<usize>,
    pub bond_order: Option<i32>,
    pub substitution_degree: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReactiveSite {
    pub kind: ReactiveSiteKind,
    pub atoms: Vec<usize>,
    pub roles: Vec<ReactiveRole>,
    pub primary_atom: Option<usize>,
    pub anchor_atoms: Vec<usize>,
    pub leaving_atom: Option<usize>,
    pub bond_order: Option<i32>,
    pub substitution_degree: Option<usize>,
}

impl ReactiveSite {
    fn new(
        kind: ReactiveSiteKind,
        atoms: impl Into<Vec<usize>>,
        roles: impl Into<Vec<ReactiveRole>>,
    ) -> Self {
        let mut atoms = atoms.into();
        atoms.sort_unstable();
        atoms.dedup();
        let mut roles = roles.into();
        roles.sort();
        roles.dedup();
        Self {
            kind,
            atoms,
            roles,
            primary_atom: None,
            anchor_atoms: Vec::new(),
            leaving_atom: None,
            bond_order: None,
            substitution_degree: None,
        }
    }

    pub fn key(&self) -> ReactiveSiteKey {
        ReactiveSiteKey {
            kind: self.kind.clone(),
            atoms: self.atoms.clone(),
            roles: self.roles.clone(),
            primary_atom: self.primary_atom,
            anchor_atoms: self.anchor_atoms.clone(),
            leaving_atom: self.leaving_atom,
            bond_order: self.bond_order,
            substitution_degree: self.substitution_degree,
        }
    }

    pub fn validate_against(&self, structure: &MolecularStructure) -> ChemistryResult<()> {
        for atom in &self.atoms {
            if *atom >= structure.atoms.len() {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: "<reactive-site>".to_string(),
                    reason: format!(
                        "reactive site {:?} references missing atom {atom}",
                        self.kind
                    ),
                });
            }
        }
        if let Some(atom) = self.primary_atom {
            if !self.atoms.contains(&atom) {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: "<reactive-site>".to_string(),
                    reason: format!(
                        "reactive site {:?} primary atom is outside site atoms",
                        self.kind
                    ),
                });
            }
        }
        for atom in &self.anchor_atoms {
            if !self.atoms.contains(atom) {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: "<reactive-site>".to_string(),
                    reason: format!(
                        "reactive site {:?} anchor atom is outside site atoms",
                        self.kind
                    ),
                });
            }
        }
        if let Some(atom) = self.leaving_atom {
            if !self.atoms.contains(&atom) {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: "<reactive-site>".to_string(),
                    reason: format!(
                        "reactive site {:?} leaving atom is outside site atoms",
                        self.kind
                    ),
                });
            }
            if !self.roles.contains(&ReactiveRole::LeavingGroup) {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: "<reactive-site>".to_string(),
                    reason: format!(
                        "reactive site {:?} has leaving atom without leaving role",
                        self.kind
                    ),
                });
            }
        }
        if self.bond_order.is_some() && self.anchor_atoms.len() < 2 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: "<reactive-site>".to_string(),
                reason: format!(
                    "reactive site {:?} has bond order without two anchors",
                    self.kind
                ),
            });
        }
        Ok(())
    }
}

pub fn find_reactive_sites(structure: &MolecularStructure) -> Vec<ReactiveSite> {
    try_find_reactive_sites(structure).expect("reactive site search produced an invalid site")
}

pub fn try_find_reactive_sites(
    structure: &MolecularStructure,
) -> ChemistryResult<Vec<ReactiveSite>> {
    let mut sites = Vec::new();
    for group in find_functional_groups(structure) {
        let (kind, roles) = match group.group_type {
            FunctionalGroupType::AcidAnhydride => (
                ReactiveSiteKind::AcidAnhydride,
                vec![ReactiveRole::Electrophile],
            ),
            FunctionalGroupType::AcylChloride => (
                ReactiveSiteKind::AcylChloride,
                vec![ReactiveRole::Electrophile, ReactiveRole::LeavingGroup],
            ),
            FunctionalGroupType::Alcohol => (
                ReactiveSiteKind::Alcohol,
                vec![ReactiveRole::Nucleophile, ReactiveRole::AcidicProton],
            ),
            FunctionalGroupType::Alkene => (
                ReactiveSiteKind::Alkene,
                vec![ReactiveRole::UnsaturatedBond, ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::Alkoxide => {
                (ReactiveSiteKind::Alkoxide, vec![ReactiveRole::Nucleophile])
            }
            FunctionalGroupType::Alkyne => (
                ReactiveSiteKind::Alkyne,
                vec![ReactiveRole::UnsaturatedBond, ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::Borane => {
                (ReactiveSiteKind::Borane, vec![ReactiveRole::Electrophile])
            }
            FunctionalGroupType::BorateEster => (
                ReactiveSiteKind::BorateEster,
                vec![ReactiveRole::Electrophile],
            ),
            FunctionalGroupType::BoricAcid => (
                ReactiveSiteKind::BoricAcid,
                vec![ReactiveRole::AcidicProton],
            ),
            FunctionalGroupType::Carbonyl => {
                if group.is_ketone == Some(true) {
                    (
                        ReactiveSiteKind::Ketone,
                        vec![ReactiveRole::Electrophile, ReactiveRole::Reducible],
                    )
                } else {
                    (
                        ReactiveSiteKind::Aldehyde,
                        vec![
                            ReactiveRole::Electrophile,
                            ReactiveRole::Oxidizable,
                            ReactiveRole::Reducible,
                        ],
                    )
                }
            }
            FunctionalGroupType::CarboxylicAcid => (
                ReactiveSiteKind::CarboxylicAcid,
                vec![ReactiveRole::AcidicProton, ReactiveRole::Electrophile],
            ),
            FunctionalGroupType::Ester => {
                (ReactiveSiteKind::Ester, vec![ReactiveRole::Electrophile])
            }
            FunctionalGroupType::Halide => (
                ReactiveSiteKind::Halide,
                vec![ReactiveRole::Electrophile, ReactiveRole::LeavingGroup],
            ),
            FunctionalGroupType::Isocyanate => (
                ReactiveSiteKind::Isocyanate,
                vec![ReactiveRole::Electrophile],
            ),
            FunctionalGroupType::Nitrile => {
                (ReactiveSiteKind::Nitrile, vec![ReactiveRole::Electrophile])
            }
            FunctionalGroupType::Nitro => (ReactiveSiteKind::Nitro, vec![ReactiveRole::Reducible]),
            FunctionalGroupType::NonTertiaryAmine => (
                ReactiveSiteKind::NonTertiaryAmine,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::NonTertiaryBorane => {
                (ReactiveSiteKind::Borane, vec![ReactiveRole::Electrophile])
            }
            FunctionalGroupType::PrimaryAmine => (
                ReactiveSiteKind::PrimaryAmine,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::Phosphine => {
                (ReactiveSiteKind::Phosphine, vec![ReactiveRole::Nucleophile])
            }
            FunctionalGroupType::NucleophilicPhosphorus => (
                ReactiveSiteKind::NucleophilicPhosphorus,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::PhosphonateCarbanion => (
                ReactiveSiteKind::PhosphonateCarbanion,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::PhosphoniumSalt => (
                ReactiveSiteKind::PhosphoniumSalt,
                vec![ReactiveRole::AcidicProton],
            ),
            FunctionalGroupType::PhosphorusYlide => (
                ReactiveSiteKind::PhosphorusYlide,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::SulfoneCarbanion => (
                ReactiveSiteKind::SulfoneCarbanion,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::UnsubstitutedAmide => {
                (ReactiveSiteKind::Amide, vec![ReactiveRole::Electrophile])
            }
            FunctionalGroupType::SubstitutedAmide => {
                (ReactiveSiteKind::Amide, vec![ReactiveRole::Electrophile])
            }
            FunctionalGroupType::AmideNitrogen => (
                ReactiveSiteKind::AmideNitrogen,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::SilylEther => (ReactiveSiteKind::SilylEther, vec![]),
            FunctionalGroupType::Acetal => (ReactiveSiteKind::Acetal, vec![]),
            FunctionalGroupType::Ketal => (ReactiveSiteKind::Ketal, vec![]),
            FunctionalGroupType::BocCarbamate => (ReactiveSiteKind::BocCarbamate, vec![]),
            FunctionalGroupType::CbzCarbamate => (ReactiveSiteKind::CbzCarbamate, vec![]),
            FunctionalGroupType::Oxime => (
                ReactiveSiteKind::Oxime,
                vec![ReactiveRole::Electrophile, ReactiveRole::LeavingGroup],
            ),
            FunctionalGroupType::Hydrazone => (
                ReactiveSiteKind::Hydrazone,
                vec![ReactiveRole::Electrophile, ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::BisNucleophile => (
                ReactiveSiteKind::BisNucleophile,
                vec![ReactiveRole::Nucleophile],
            ),
            FunctionalGroupType::DicarbonylElectrophile => (
                ReactiveSiteKind::DicarbonylElectrophile,
                vec![ReactiveRole::Electrophile],
            ),
            FunctionalGroupType::UreaLike => (
                ReactiveSiteKind::UreaLike,
                vec![ReactiveRole::Nucleophile, ReactiveRole::Electrophile],
            ),
            FunctionalGroupType::FormylationDonor => (
                ReactiveSiteKind::FormylationDonor,
                vec![ReactiveRole::Electrophile],
            ),
        };
        sites.push(ReactiveSite::new(kind, group.atoms, roles));
    }

    // Remove conflicting sites: protected functional groups should not have their
    // original reactive sites detected
    remove_conflicting_protected_sites(structure, &mut sites);

    add_aromatic_sites(structure, &mut sites);
    add_oxygen_sites(structure, &mut sites);
    add_sulfur_sites(structure, &mut sites);
    add_nitrogen_sites(structure, &mut sites);
    add_organometallic_sites(structure, &mut sites);
    add_alpha_sites(structure, &mut sites);
    enrich_sites(structure, &mut sites);
    deduplicate_sites(&mut sites);
    for site in &sites {
        site.validate_against(structure)?;
    }
    Ok(sites)
}

fn add_aromatic_sites(structure: &MolecularStructure, sites: &mut Vec<ReactiveSite>) {
    let aromatic_carbons = (0..structure.atoms.len())
        .filter(|atom| {
            structure.atoms[*atom].element == "C"
                && structure
                    .neighbors(*atom)
                    .iter()
                    .filter(|(_, order)| bond_order_matches(*order, 1.5))
                    .count()
                    >= 2
        })
        .collect::<Vec<_>>();
    if aromatic_carbons.len() >= 5 {
        sites.push(ReactiveSite::new(
            ReactiveSiteKind::AromaticRing,
            aromatic_carbons.clone(),
            [ReactiveRole::AromaticDirector, ReactiveRole::Nucleophile],
        ));
        for carbon in aromatic_carbons {
            let kind = if structure.neighbors(carbon).iter().any(|(neighbor, order)| {
                bond_order_matches(*order, 1.0)
                    && matches!(
                        structure.atoms[*neighbor].element.as_str(),
                        "Cl" | "Br" | "I" | "F"
                    )
            }) {
                ReactiveSiteKind::ArylHalide
            } else {
                ReactiveSiteKind::AromaticCarbon
            };
            sites.push(ReactiveSite::new(
                kind,
                [carbon],
                [ReactiveRole::AromaticDirector],
            ));
        }
    }
}

fn add_oxygen_sites(structure: &MolecularStructure, sites: &mut Vec<ReactiveSite>) {
    for oxygen in 0..structure.atoms.len() {
        if structure.atoms[oxygen].element != "O" {
            continue;
        }
        let carbon_neighbors = structure
            .neighbors(oxygen)
            .into_iter()
            .filter(|(neighbor, order)| {
                bond_order_matches(*order, 1.0) && structure.atoms[*neighbor].element == "C"
            })
            .map(|(neighbor, _)| neighbor)
            .collect::<Vec<_>>();
        if carbon_neighbors.len() == 2 {
            if carbon_neighbors.iter().all(|carbon| {
                structure.neighbors(*carbon).iter().any(|(other, order)| {
                    carbon_neighbors.contains(other) && bond_order_matches(*order, 1.0)
                })
            }) {
                let mut atoms = carbon_neighbors.clone();
                atoms.push(oxygen);
                sites.push(ReactiveSite::new(
                    ReactiveSiteKind::Epoxide,
                    atoms,
                    [ReactiveRole::Electrophile],
                ));
            } else {
                let mut atoms = carbon_neighbors.clone();
                atoms.push(oxygen);
                sites.push(ReactiveSite::new(
                    ReactiveSiteKind::Ether,
                    atoms,
                    [ReactiveRole::Nucleophile],
                ));
            }
        }
        if carbon_neighbors.len() == 1
            && structure.hydrogen_count(oxygen) == 1
            && is_aromatic_carbon(structure, carbon_neighbors[0])
        {
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Phenol,
                [carbon_neighbors[0], oxygen],
                [ReactiveRole::AcidicProton, ReactiveRole::AromaticDirector],
            ));
        }
    }
}

fn add_sulfur_sites(structure: &MolecularStructure, sites: &mut Vec<ReactiveSite>) {
    for sulfur in 0..structure.atoms.len() {
        if structure.atoms[sulfur].element != "S" {
            continue;
        }
        let neighbors = structure.neighbors(sulfur);
        let chlorines = neighbors
            .iter()
            .filter(|(neighbor, order)| {
                bond_order_matches(*order, 1.0) && structure.atoms[*neighbor].element == "Cl"
            })
            .count();
        let oxygens = neighbors
            .iter()
            .filter(|(neighbor, order)| structure.atoms[*neighbor].element == "O" && *order >= 1.5)
            .count();
        if chlorines > 0 && oxygens >= 2 {
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::SulfonylChloride,
                [sulfur],
                [ReactiveRole::Electrophile, ReactiveRole::LeavingGroup],
            ));
        }
        for (carbon, order) in &neighbors {
            if bond_order_matches(*order, 1.0)
                && structure.atoms[*carbon].element == "C"
                && structure.atoms[*carbon].charge < -0.1
                && oxygens >= 2
            {
                let mut atoms = vec![sulfur, *carbon];
                atoms.extend(
                    neighbors
                        .iter()
                        .filter(|(neighbor, order)| {
                            structure.atoms[*neighbor].element == "O" && *order >= 1.5
                        })
                        .map(|(neighbor, _)| *neighbor),
                );
                sites.push(ReactiveSite::new(
                    ReactiveSiteKind::SulfoneCarbanion,
                    atoms,
                    [ReactiveRole::Nucleophile],
                ));
            }
        }
        if neighbors
            .iter()
            .any(|(neighbor, _)| structure.atoms[*neighbor].element == "H")
        {
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Thiol,
                [sulfur],
                [ReactiveRole::AcidicProton, ReactiveRole::Nucleophile],
            ));
        }
        let carbon_substituents: Vec<usize> = neighbors
            .iter()
            .filter(|(neighbor, order)| {
                bond_order_matches(*order, 1.0) && structure.atoms[*neighbor].element == "C"
            })
            .map(|(neighbor, _)| *neighbor)
            .collect();
        if oxygens == 1 && carbon_substituents.len() == 2 {
            let mut atoms = vec![sulfur];
            atoms.extend(carbon_substituents.iter().copied());
            let oxygen_atom = neighbors
                .iter()
                .find(|(neighbor, order)| {
                    structure.atoms[*neighbor].element == "O" && *order >= 1.5
                })
                .map(|(neighbor, _)| *neighbor);
            if let Some(o) = oxygen_atom {
                atoms.push(o);
            }
            atoms.sort_unstable();
            atoms.dedup();
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Sulfoxide,
                atoms,
                [ReactiveRole::Oxidizable],
            ));
        }
        if carbon_substituents.len() >= 2 && oxygens == 0 {
            let mut atoms = vec![sulfur];
            atoms.extend(carbon_substituents.iter().copied());
            atoms.sort_unstable();
            atoms.dedup();
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Sulfide,
                atoms,
                [ReactiveRole::Nucleophile, ReactiveRole::Oxidizable],
            ));
        }
    }
}

fn add_nitrogen_sites(structure: &MolecularStructure, sites: &mut Vec<ReactiveSite>) {
    for nitrogen in 0..structure.atoms.len() {
        if structure.atoms[nitrogen].element != "N" {
            continue;
        }
        let neighbors = structure.neighbors(nitrogen);
        if structure.atoms[nitrogen].charge > 0.5
            && neighbors
                .iter()
                .any(|(neighbor, _)| structure.atoms[*neighbor].element == "N")
        {
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Diazonium,
                [nitrogen],
                [ReactiveRole::Electrophile, ReactiveRole::LeavingGroup],
            ));
        }
        if neighbors
            .iter()
            .filter(|(neighbor, _)| structure.atoms[*neighbor].element == "N")
            .count()
            >= 2
        {
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Azide,
                [nitrogen],
                [ReactiveRole::Nucleophile],
            ));
        }
        if neighbors.iter().any(|(neighbor, order)| {
            bond_order_matches(*order, 2.0) && structure.atoms[*neighbor].element == "C"
        }) {
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Imine,
                [nitrogen],
                [ReactiveRole::Electrophile],
            ));
        }
    }
}

fn add_organometallic_sites(structure: &MolecularStructure, sites: &mut Vec<ReactiveSite>) {
    for bond in &structure.bonds {
        if !bond_order_matches(bond.order, 1.0) {
            continue;
        }
        let left = structure.atoms[bond.from].element.as_str();
        let right = structure.atoms[bond.to].element.as_str();
        let metal = match (left, right) {
            ("C", "Mg") | ("Mg", "C") => Some(ReactiveSiteKind::Organomagnesium),
            ("C", "Li") | ("Li", "C") => Some(ReactiveSiteKind::Organolithium),
            ("C", "Cu") | ("Cu", "C") => Some(ReactiveSiteKind::Organocopper),
            _ => None,
        };
        if let Some(kind) = metal {
            sites.push(ReactiveSite::new(
                kind,
                [bond.from, bond.to],
                [ReactiveRole::Nucleophile],
            ));
        }
    }
}

fn add_alpha_sites(structure: &MolecularStructure, sites: &mut Vec<ReactiveSite>) {
    let carbonyl_carbons = sites
        .iter()
        .filter(|site| {
            matches!(
                site.kind,
                ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Ester
            )
        })
        .filter_map(|site| carbonyl_carbon_in_site(structure, site))
        .collect::<BTreeSet<_>>();
    for carbonyl in carbonyl_carbons {
        for (neighbor, order) in structure.neighbors(carbonyl) {
            if bond_order_matches(order, 1.0)
                && structure.atoms[neighbor].element == "C"
                && structure.hydrogen_count(neighbor) > 0
            {
                sites.push(ReactiveSite::new(
                    ReactiveSiteKind::Enol,
                    [carbonyl, neighbor],
                    [ReactiveRole::AlphaCarbon, ReactiveRole::Nucleophile],
                ));
            }
        }
    }
}

fn is_aromatic_carbon(structure: &MolecularStructure, atom: usize) -> bool {
    structure.atoms[atom].element == "C"
        && structure
            .neighbors(atom)
            .iter()
            .filter(|(_, order)| bond_order_matches(*order, 1.5))
            .count()
            >= 2
}

/// Remove conflicting sites when a protected functional group is present.
/// A protected group should not have its original reactive site detected.
fn remove_conflicting_protected_sites(
    structure: &MolecularStructure,
    sites: &mut Vec<ReactiveSite>,
) {
    let mut protected_alcohol_oxygens: Vec<usize> = Vec::new();
    let mut protected_amine_nitrogens: Vec<usize> = Vec::new();
    let mut protected_carbamate_carbons: Vec<usize> = Vec::new();

    for site in sites.iter() {
        if site.kind == ReactiveSiteKind::SilylEther {
            if let Some(&oxygen) = site
                .atoms
                .iter()
                .find(|&&atom| structure.atoms[atom].element == "O")
            {
                protected_alcohol_oxygens.push(oxygen);
            }
        }
        if matches!(
            site.kind,
            ReactiveSiteKind::BocCarbamate | ReactiveSiteKind::CbzCarbamate
        ) {
            if let Some(&nitrogen) = site
                .atoms
                .iter()
                .find(|&&atom| structure.atoms[atom].element == "N")
            {
                protected_amine_nitrogens.push(nitrogen);
            }
            if let Some(carbon) = carbonyl_carbon_in_site(structure, site) {
                protected_carbamate_carbons.push(carbon);
            }
        }
    }

    sites.retain(|site| {
        if site.kind == ReactiveSiteKind::Alcohol
            && site.atoms.iter().any(|&atom| {
                structure.atoms[atom].element == "O" && protected_alcohol_oxygens.contains(&atom)
            })
        {
            return false;
        }
        if matches!(
            site.kind,
            ReactiveSiteKind::PrimaryAmine | ReactiveSiteKind::NonTertiaryAmine
        ) && site.atoms.iter().any(|&atom| {
            structure.atoms[atom].element == "N" && protected_amine_nitrogens.contains(&atom)
        }) {
            return false;
        }
        if site.kind == ReactiveSiteKind::Ester
            && site.atoms.iter().any(|&atom| {
                structure.atoms[atom].element == "C" && protected_carbamate_carbons.contains(&atom)
            })
        {
            return false;
        }
        true
    });
}

fn carbonyl_carbon_in_site(structure: &MolecularStructure, site: &ReactiveSite) -> Option<usize> {
    site.atoms.iter().copied().find(|atom| {
        structure.atoms[*atom].element == "C"
            && structure.neighbors(*atom).iter().any(|(neighbor, order)| {
                site.atoms.contains(neighbor)
                    && structure.atoms[*neighbor].element == "O"
                    && bond_order_matches(*order, 2.0)
            })
    })
}

fn enrich_sites(structure: &MolecularStructure, sites: &mut [ReactiveSite]) {
    for site in sites {
        site.primary_atom = primary_atom_for_site(site);
        site.anchor_atoms = anchor_atoms_for_site(site);
        site.leaving_atom = leaving_atom_for_site(structure, site);
        site.bond_order = bond_order_for_site(structure, site);
        site.substitution_degree = site.primary_atom.and_then(|atom| {
            (structure.atoms[atom].element == "C").then(|| structure.carbon_degree(atom))
        });
    }
}

fn primary_atom_for_site(site: &ReactiveSite) -> Option<usize> {
    match site.kind {
        ReactiveSiteKind::Alcohol
        | ReactiveSiteKind::Alkoxide
        | ReactiveSiteKind::Aldehyde
        | ReactiveSiteKind::Ketone
        | ReactiveSiteKind::Carbonyl
        | ReactiveSiteKind::CarboxylicAcid
        | ReactiveSiteKind::AcylChloride
        | ReactiveSiteKind::AcidAnhydride
        | ReactiveSiteKind::Ester
        | ReactiveSiteKind::Halide
        | ReactiveSiteKind::Nitrile
        | ReactiveSiteKind::Nitro
        | ReactiveSiteKind::Amide
        | ReactiveSiteKind::Borane
        | ReactiveSiteKind::BoricAcid
        | ReactiveSiteKind::BorateEster
        | ReactiveSiteKind::Alkene
        | ReactiveSiteKind::Alkyne
        | ReactiveSiteKind::Enol
        | ReactiveSiteKind::AromaticCarbon
        | ReactiveSiteKind::ArylHalide => site.atoms.first().copied(),
        ReactiveSiteKind::Hydrazone
        | ReactiveSiteKind::BisNucleophile
        | ReactiveSiteKind::DicarbonylElectrophile
        | ReactiveSiteKind::UreaLike
        | ReactiveSiteKind::FormylationDonor => site.atoms.first().copied(),
        ReactiveSiteKind::AromaticRing => None,
        ReactiveSiteKind::SilylEther => site.atoms.first().copied(),
        ReactiveSiteKind::Acetal | ReactiveSiteKind::Ketal => site.atoms.first().copied(),
        _ => site.atoms.first().copied(),
    }
}

fn anchor_atoms_for_site(site: &ReactiveSite) -> Vec<usize> {
    match site.kind {
        ReactiveSiteKind::Alkene
        | ReactiveSiteKind::Alkyne
        | ReactiveSiteKind::Carbonyl
        | ReactiveSiteKind::Aldehyde
        | ReactiveSiteKind::Ketone
        | ReactiveSiteKind::Enol
        | ReactiveSiteKind::Organomagnesium
        | ReactiveSiteKind::Organolithium
        | ReactiveSiteKind::Organocopper
        | ReactiveSiteKind::Oxime
        | ReactiveSiteKind::Hydrazone
        | ReactiveSiteKind::DicarbonylElectrophile => site.atoms.iter().copied().take(2).collect(),
        ReactiveSiteKind::Epoxide => site.atoms.clone(),
        _ => site
            .primary_atom
            .map(|atom| vec![atom])
            .unwrap_or_else(|| site.atoms.clone()),
    }
}

fn leaving_atom_for_site(structure: &MolecularStructure, site: &ReactiveSite) -> Option<usize> {
    if !site.roles.contains(&ReactiveRole::LeavingGroup) {
        return None;
    }
    site.atoms.iter().copied().find(|atom| {
        matches!(
            structure.atoms[*atom].element.as_str(),
            "F" | "Cl" | "Br" | "I"
        )
    })
}

fn bond_order_for_site(structure: &MolecularStructure, site: &ReactiveSite) -> Option<i32> {
    let first = *site.anchor_atoms.first()?;
    let second = *site.anchor_atoms.get(1)?;
    structure
        .neighbors(first)
        .into_iter()
        .find_map(|(neighbor, order)| {
            (neighbor == second).then_some((order * 1000.0).round() as i32)
        })
}

fn deduplicate_sites(sites: &mut Vec<ReactiveSite>) {
    let mut seen = BTreeSet::new();
    sites.retain(|site| seen.insert(site.key()));
    sites.sort_by(|left, right| left.key().cmp(&right.key()));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;

    fn site_kinds(code: &str) -> Vec<ReactiveSiteKind> {
        find_reactive_sites(&parse_frowns(code).unwrap())
            .into_iter()
            .map(|site| site.kind)
            .collect()
    }

    #[test]
    fn carbonyl_has_multiple_reactive_roles_without_duplicate_sites() {
        let sites = find_reactive_sites(&parse_frowns("CC(=O)C").unwrap());
        assert!(sites
            .iter()
            .any(|site| site.kind == ReactiveSiteKind::Ketone));
        assert!(sites.iter().any(|site| site.kind == ReactiveSiteKind::Enol));
        let unique = sites.iter().collect::<BTreeSet<_>>();
        assert_eq!(unique.len(), sites.len());
    }

    #[test]
    fn detects_aromatic_phenol_ether_and_epoxide_sites() {
        let phenol = site_kinds("destroy:benzene:O,,,,,");
        assert!(phenol.contains(&ReactiveSiteKind::AromaticRing));
        assert!(phenol.contains(&ReactiveSiteKind::Phenol));

        let ether = site_kinds("COC");
        assert!(ether.contains(&ReactiveSiteKind::Ether));

        let epoxide = site_kinds(
            "destroy:graph:atoms=C.C.O.H.H.H.H;bonds=0-s-1,0-s-2,1-s-2,0-s-3,0-s-4,1-s-5,1-s-6",
        );
        assert!(epoxide.contains(&ReactiveSiteKind::Epoxide));
    }

    #[test]
    fn detects_generalized_cn_heterocycle_building_sites() {
        let hydrazone = site_kinds(
            "destroy:graph:atoms=C.N.N.H.H.H.H;\
             bonds=0-d-1,1-s-2,0-s-3,0-s-4,2-s-5,2-s-6",
        );
        assert!(hydrazone.contains(&ReactiveSiteKind::Hydrazone));

        let urea = site_kinds(
            "destroy:graph:atoms=C.O.N.N.H.H.H.H;\
             bonds=0-d-1,0-s-2,0-s-3,2-s-4,2-s-5,3-s-6,3-s-7",
        );
        assert!(urea.contains(&ReactiveSiteKind::UreaLike));
        assert!(urea.contains(&ReactiveSiteKind::BisNucleophile));

        let dicarbonyl = site_kinds(
            "destroy:graph:atoms=C.O.C.C.O.H.H.H.H;\
             bonds=0-d-1,0-s-2,2-s-3,3-d-4,0-s-5,2-s-6,2-s-7,3-s-8",
        );
        assert!(dicarbonyl.contains(&ReactiveSiteKind::DicarbonylElectrophile));

        let formamide = site_kinds(
            "destroy:graph:atoms=C.O.N.H.H.H;\
             bonds=0-d-1,0-s-2,0-s-3,2-s-4,2-s-5",
        );
        assert!(formamide.contains(&ReactiveSiteKind::FormylationDonor));
    }
}
