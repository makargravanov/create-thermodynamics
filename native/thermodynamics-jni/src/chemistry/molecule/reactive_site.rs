use std::collections::BTreeSet;

use super::functional_group::{find_functional_groups, FunctionalGroupType};
use super::molecule::{bond_order_matches, MolecularStructure};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReactiveSiteKind {
    AcidAnhydride,
    AcylChloride,
    Alcohol,
    Alkene,
    Alkoxide,
    Alkyne,
    Aldehyde,
    Amide,
    AromaticCarbon,
    AromaticRing,
    ArylHalide,
    Azide,
    Borane,
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
    Ketone,
    Nitrile,
    Nitro,
    NonTertiaryAmine,
    Organocopper,
    Organolithium,
    Organomagnesium,
    Phenol,
    PrimaryAmine,
    Sulfide,
    SulfonylChloride,
    Thiol,
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
pub struct ReactiveSite {
    pub kind: ReactiveSiteKind,
    pub atoms: Vec<usize>,
    pub roles: Vec<ReactiveRole>,
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
        Self { kind, atoms, roles }
    }
}

pub fn find_reactive_sites(structure: &MolecularStructure) -> Vec<ReactiveSite> {
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
                ReactiveSiteKind::CarboxylicAcid,
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
            FunctionalGroupType::Isocyanate => {
                (ReactiveSiteKind::Imine, vec![ReactiveRole::Electrophile])
            }
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
            FunctionalGroupType::UnsubstitutedAmide => {
                (ReactiveSiteKind::Amide, vec![ReactiveRole::Electrophile])
            }
        };
        sites.push(ReactiveSite::new(kind, group.atoms, roles));
    }

    add_aromatic_sites(structure, &mut sites);
    add_oxygen_sites(structure, &mut sites);
    add_sulfur_sites(structure, &mut sites);
    add_nitrogen_sites(structure, &mut sites);
    add_organometallic_sites(structure, &mut sites);
    add_alpha_sites(structure, &mut sites);
    deduplicate_sites(&mut sites);
    sites
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
        if neighbors
            .iter()
            .filter(|(neighbor, order)| {
                bond_order_matches(*order, 1.0) && structure.atoms[*neighbor].element == "C"
            })
            .count()
            >= 2
        {
            sites.push(ReactiveSite::new(
                ReactiveSiteKind::Sulfide,
                [sulfur],
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
        .flat_map(|site| site.atoms.iter().copied())
        .filter(|atom| structure.atoms[*atom].element == "C")
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

fn deduplicate_sites(sites: &mut Vec<ReactiveSite>) {
    let mut seen = BTreeSet::new();
    sites.retain(|site| seen.insert((site.kind.clone(), site.atoms.clone(), site.roles.clone())));
    sites.sort_by(|left, right| {
        left.kind
            .cmp(&right.kind)
            .then_with(|| left.atoms.cmp(&right.atoms))
            .then_with(|| left.roles.cmp(&right.roles))
    });
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
}
