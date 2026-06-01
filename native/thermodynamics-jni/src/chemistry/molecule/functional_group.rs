use std::collections::BTreeSet;

use super::molecule::{bond_order_matches, MolecularStructure};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FunctionalGroupType {
    AcidAnhydride,
    AcylChloride,
    Alcohol,
    Alkene,
    Alkoxide,
    Alkyne,
    Borane,
    BorateEster,
    BoricAcid,
    Carbonyl,
    CarboxylicAcid,
    Ester,
    Halide,
    Isocyanate,
    Nitrile,
    Nitro,
    NonTertiaryAmine,
    NonTertiaryBorane,
    PrimaryAmine,
    Phosphine,
    PhosphonateCarbanion,
    PhosphoniumSalt,
    PhosphorusYlide,
    SulfoneCarbanion,
    UnsubstitutedAmide,
    // Protecting groups
    SilylEther,
    Acetal,
    Ketal,
    BocCarbamate,
    CbzCarbamate,
    FmocCarbamate,
    AcylProtectedAmine,
    EsterProtectedAcid,
    ProtectedThiol,
    Thioacetal,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FunctionalGroup {
    pub group_type: FunctionalGroupType,
    pub atoms: Vec<usize>,
    pub degree: Option<usize>,
    pub is_ketone: Option<bool>,
}

impl FunctionalGroup {
    fn new(group_type: FunctionalGroupType, atoms: impl Into<Vec<usize>>) -> Self {
        Self {
            group_type,
            atoms: atoms.into(),
            degree: None,
            is_ketone: None,
        }
    }

    fn with_degree(mut self, degree: usize) -> Self {
        self.degree = Some(degree);
        self
    }

    fn with_ketone(mut self, is_ketone: bool) -> Self {
        self.is_ketone = Some(is_ketone);
        self
    }
}

pub fn find_functional_groups(structure: &MolecularStructure) -> Vec<FunctionalGroup> {
    let mut groups = Vec::new();
    let mut carbonyl_carbons_to_ignore = BTreeSet::new();
    let mut alkenes_to_ignore = BTreeSet::new();
    let mut alkynes_to_ignore = BTreeSet::new();

    for carbon in 0..structure.atoms.len() {
        if structure.atoms[carbon].element != "C" || carbonyl_carbons_to_ignore.contains(&carbon) {
            continue;
        }

        let carbonyl_oxygens = bonded(structure, carbon, "O", Some(2.0));
        let single_oxygens = bonded(structure, carbon, "O", Some(1.0));
        let chlorines = bonded(structure, carbon, "Cl", Some(1.0));
        let fluorines = bonded(structure, carbon, "F", Some(1.0));
        let mut halogens = chlorines.clone();
        halogens.extend(bonded(structure, carbon, "I", Some(1.0)));
        let nitrogens = bonded(structure, carbon, "N", Some(1.0));
        let borons = bonded(structure, carbon, "B", Some(1.0));
        let carbon_hydrogens = structure.hydrogen_count(carbon);
        let carbons = bonded(structure, carbon, "C", Some(1.0));
        let alkene_carbons = bonded(structure, carbon, "C", Some(2.0));
        let alkyne_carbons = bonded(structure, carbon, "C", Some(3.0));
        let double_nitrogens = bonded(structure, carbon, "N", Some(2.0));
        let nitrile_nitrogens = bonded(structure, carbon, "N", Some(3.0));
        let r_groups = bonded(structure, carbon, "R", None);

        if carbonyl_oxygens.len() == 1 {
            let carbonyl_oxygen = carbonyl_oxygens[0];
            if double_nitrogens.len() == 1 {
                continue;
            }

            if single_oxygens.len() == 1 {
                let alcohol_oxygen = single_oxygens[0];
                let oxygen_carbons = bonded(structure, alcohol_oxygen, "C", Some(1.0));
                if oxygen_carbons.len() == 2 {
                    if let Some(other_carbon) = oxygen_carbons
                        .into_iter()
                        .find(|candidate| *candidate != carbon)
                    {
                        let other_carbonyl_oxygens =
                            bonded(structure, other_carbon, "O", Some(2.0));
                        if other_carbonyl_oxygens.len() == 1 {
                            groups.push(FunctionalGroup::new(
                                FunctionalGroupType::AcidAnhydride,
                                vec![
                                    carbon,
                                    carbonyl_oxygen,
                                    other_carbon,
                                    other_carbonyl_oxygens[0],
                                    alcohol_oxygen,
                                ],
                            ));
                        } else {
                            groups.push(FunctionalGroup::new(
                                FunctionalGroupType::Ester,
                                vec![carbon, other_carbon, carbonyl_oxygen, alcohol_oxygen],
                            ));
                        }
                        carbonyl_carbons_to_ignore.insert(other_carbon);
                        continue;
                    }
                } else if structure.hydrogen_count(alcohol_oxygen) == 1 {
                    let hydrogen = bonded(structure, alcohol_oxygen, "H", Some(1.0))[0];
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::CarboxylicAcid,
                        vec![carbon, carbonyl_oxygen, alcohol_oxygen, hydrogen],
                    ));
                    continue;
                }
            } else if nitrogens.len() == 1 {
                let nitrogen = nitrogens[0];
                let hydrogens = bonded(structure, nitrogen, "H", Some(1.0));
                if hydrogens.len() == 2 {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::UnsubstitutedAmide,
                        vec![
                            carbon,
                            carbonyl_oxygen,
                            nitrogen,
                            hydrogens[0],
                            hydrogens[1],
                        ],
                    ));
                    continue;
                }
            } else if chlorines.len() == 1 {
                groups.push(FunctionalGroup::new(
                    FunctionalGroupType::AcylChloride,
                    vec![carbon, carbonyl_oxygen, chlorines[0]],
                ));
                continue;
            } else if carbons.len() == 2 {
                groups.push(
                    FunctionalGroup::new(
                        FunctionalGroupType::Carbonyl,
                        vec![carbon, carbonyl_oxygen],
                    )
                    .with_ketone(true),
                );
            } else if carbons.len() + carbon_hydrogens + r_groups.len() == 2 {
                groups.push(
                    FunctionalGroup::new(
                        FunctionalGroupType::Carbonyl,
                        vec![carbon, carbonyl_oxygen],
                    )
                    .with_ketone(false),
                );
            }
        } else {
            for halogen in halogens {
                if chlorines.len() < 3 && fluorines.is_empty() {
                    groups.push(
                        FunctionalGroup::new(FunctionalGroupType::Halide, vec![carbon, halogen])
                            .with_degree(carbons.len()),
                    );
                }
            }

            for oxygen in single_oxygens {
                let oxygen_hydrogens = bonded(structure, oxygen, "H", Some(1.0));
                if oxygen_hydrogens.len() == 1 {
                    groups.push(
                        FunctionalGroup::new(
                            FunctionalGroupType::Alcohol,
                            vec![carbon, oxygen, oxygen_hydrogens[0]],
                        )
                        .with_degree(carbons.len()),
                    );
                } else if structure.atoms[oxygen].charge == -1.0 {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::Alkoxide,
                        vec![carbon, oxygen],
                    ));
                }
                let borate_borons = bonded(structure, oxygen, "B", None);
                if borate_borons.len() == 1 {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::BorateEster,
                        vec![carbon, oxygen, borate_borons[0]],
                    ));
                }
            }

            for nitrogen in nitrogens {
                let double_carbons = bonded(structure, nitrogen, "C", Some(2.0));
                let aromatic_oxygens = bonded(structure, nitrogen, "O", Some(1.5));
                if double_carbons.len() == 1 {
                    let isocyanate_carbon = double_carbons[0];
                    let isocyanate_oxygens = bonded(structure, isocyanate_carbon, "O", Some(2.0));
                    if isocyanate_oxygens.len() == 1 {
                        groups.push(FunctionalGroup::new(
                            FunctionalGroupType::Isocyanate,
                            vec![carbon, nitrogen, isocyanate_carbon, isocyanate_oxygens[0]],
                        ));
                    }
                } else if aromatic_oxygens.len() == 2 {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::Nitro,
                        vec![carbon, nitrogen, aromatic_oxygens[0], aromatic_oxygens[1]],
                    ));
                } else if nitrile_nitrogens.is_empty() && double_nitrogens.is_empty() {
                    let hydrogens = bonded(structure, nitrogen, "H", Some(1.0));
                    for hydrogen in &hydrogens {
                        groups.push(FunctionalGroup::new(
                            FunctionalGroupType::NonTertiaryAmine,
                            vec![carbon, nitrogen, *hydrogen],
                        ));
                    }
                    if hydrogens.len() == 2 {
                        groups.push(FunctionalGroup::new(
                            FunctionalGroupType::PrimaryAmine,
                            vec![carbon, nitrogen, hydrogens[0], hydrogens[1]],
                        ));
                    }
                }
            }

            if nitrile_nitrogens.len() == 1 && carbons.len() == 1 {
                let second_carbon = carbons[0];
                let second_carbon_nitrogens = bonded(structure, second_carbon, "N", Some(1.0));
                if second_carbon_nitrogens.len() == 1 {
                    let first_nitrogen = second_carbon_nitrogens[0];
                    if bonded(structure, first_nitrogen, "N", Some(2.0)).is_empty() {
                        groups.push(FunctionalGroup::new(
                            FunctionalGroupType::Nitrile,
                            vec![carbon, nitrile_nitrogens[0]],
                        ));
                    }
                } else if second_carbon_nitrogens.is_empty() {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::Nitrile,
                        vec![carbon, nitrile_nitrogens[0]],
                    ));
                }
            }

            for boron in borons {
                for hydrogen in bonded(structure, boron, "H", Some(1.0)) {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::NonTertiaryBorane,
                        vec![carbon, boron, hydrogen],
                    ));
                }
                groups.push(FunctionalGroup::new(
                    FunctionalGroupType::Borane,
                    vec![carbon, boron],
                ));
            }
        }

        for alkene_carbon in alkene_carbons {
            if alkenes_to_ignore.contains(&alkene_carbon) {
                continue;
            }
            let first_degree = structure.carbon_degree(carbon).saturating_sub(1);
            let second_degree = structure.carbon_degree(alkene_carbon).saturating_sub(1);
            if first_degree >= second_degree {
                groups.push(FunctionalGroup::new(
                    FunctionalGroupType::Alkene,
                    vec![carbon, alkene_carbon],
                ));
            }
            if second_degree >= first_degree {
                groups.push(FunctionalGroup::new(
                    FunctionalGroupType::Alkene,
                    vec![alkene_carbon, carbon],
                ));
            }
            alkenes_to_ignore.insert(carbon);
        }

        if alkyne_carbons.len() == 1 {
            let alkyne_carbon = alkyne_carbons[0];
            if !alkynes_to_ignore.contains(&alkyne_carbon) {
                let first_degree = structure.carbon_degree(carbon).saturating_sub(1);
                let second_degree = structure.carbon_degree(alkyne_carbon).saturating_sub(1);
                if first_degree >= second_degree {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::Alkyne,
                        vec![carbon, alkyne_carbon],
                    ));
                }
                if second_degree >= first_degree {
                    groups.push(FunctionalGroup::new(
                        FunctionalGroupType::Alkyne,
                        vec![alkyne_carbon, carbon],
                    ));
                }
                alkynes_to_ignore.insert(carbon);
                alkynes_to_ignore.insert(alkyne_carbon);
            }
        }
    }

    for boron in 0..structure.atoms.len() {
        if structure.atoms[boron].element != "B" {
            continue;
        }
        for oxygen in bonded(structure, boron, "O", None) {
            let hydrogens = bonded(structure, oxygen, "H", Some(1.0));
            if hydrogens.len() == 1 {
                groups.push(FunctionalGroup::new(
                    FunctionalGroupType::BoricAcid,
                    vec![boron, oxygen, hydrogens[0]],
                ));
            }
        }
    }

    for phosphorus in 0..structure.atoms.len() {
        if structure.atoms[phosphorus].element != "P" {
            continue;
        }
        let hydrogens = bonded(structure, phosphorus, "H", Some(1.0));
        let carbons = bonded(structure, phosphorus, "C", Some(1.0));
        let double_oxygens = bonded(structure, phosphorus, "O", Some(2.0));
        let positive = structure.atoms[phosphorus].charge > 0.1;
        let negative_carbon = carbons
            .iter()
            .copied()
            .find(|carbon| structure.atoms[*carbon].charge < -0.1);

        if !positive && !double_oxygens.is_empty() {
            if let Some(alpha_carbon) = negative_carbon {
                let mut atoms = vec![phosphorus, alpha_carbon];
                atoms.extend(double_oxygens.iter().copied());
                atoms.extend(carbons.iter().copied().filter(|atom| *atom != alpha_carbon));
                atoms.sort_unstable();
                atoms.dedup();
                groups.push(FunctionalGroup::new(
                    FunctionalGroupType::PhosphonateCarbanion,
                    atoms,
                ));
                continue;
            }
        }

        if positive {
            if let Some(alpha_carbon) = negative_carbon {
                let mut atoms = vec![phosphorus, alpha_carbon];
                atoms.extend(bonded(structure, alpha_carbon, "H", Some(1.0)));
                atoms.extend(carbons.iter().copied().filter(|atom| *atom != alpha_carbon));
                atoms.extend(hydrogens.iter().copied());
                atoms.sort_unstable();
                atoms.dedup();
                groups.push(FunctionalGroup::new(FunctionalGroupType::PhosphorusYlide, atoms));
                continue;
            }
            if let Some(alpha_carbon) = carbons
                .iter()
                .copied()
                .find(|carbon| !bonded(structure, *carbon, "H", Some(1.0)).is_empty())
            {
                let mut atoms = vec![phosphorus, alpha_carbon];
                atoms.extend(bonded(structure, alpha_carbon, "H", Some(1.0)));
                atoms.extend(carbons.iter().copied().filter(|atom| *atom != alpha_carbon));
                atoms.extend(hydrogens.iter().copied());
                atoms.sort_unstable();
                atoms.dedup();
                groups.push(FunctionalGroup::new(FunctionalGroupType::PhosphoniumSalt, atoms));
                continue;
            }
        }

        if !positive
            && structure.atoms[phosphorus].charge.abs() < 0.1
            && hydrogens.is_empty()
            && double_oxygens.is_empty()
            && carbons.len() == 3
        {
            let mut atoms = vec![phosphorus];
            atoms.extend(carbons.iter().copied());
            atoms.sort_unstable();
            atoms.dedup();
            groups.push(FunctionalGroup::new(FunctionalGroupType::Phosphine, atoms));
        }
    }

    // Add protecting group detection
    add_protecting_groups(structure, &mut groups);

    groups
}

fn add_protecting_groups(structure: &MolecularStructure, groups: &mut Vec<FunctionalGroup>) {
    // Detect silyl ethers: C-O-Si
    for oxygen in 0..structure.atoms.len() {
        if structure.atoms[oxygen].element != "O" {
            continue;
        }
        let neighbors = structure.neighbors(oxygen);
        let carbon_neighbors: Vec<usize> = neighbors
            .iter()
            .filter(|(atom, order)| {
                structure.atoms[*atom].element == "C"
                    && bond_order_matches(*order, 1.0)
            })
            .map(|(atom, _)| *atom)
            .collect();
        let silicon_neighbors: Vec<usize> = neighbors
            .iter()
            .filter(|(atom, order)| {
                structure.atoms[*atom].element == "Si"
                    && bond_order_matches(*order, 1.0)
            })
            .map(|(atom, _)| *atom)
            .collect();
        
        if !carbon_neighbors.is_empty() && !silicon_neighbors.is_empty() {
            // This is a silyl ether - the oxygen is bonded to both carbon and silicon
            let mut atoms = vec![oxygen];
            atoms.extend(carbon_neighbors.iter().copied());
            atoms.extend(silicon_neighbors.iter().copied());
            groups.push(FunctionalGroup::new(FunctionalGroupType::SilylEther, atoms));
        }
    }

    // Detect acetals/ketals: C with two ether oxygens (not part of carbonyl)
    for carbon in 0..structure.atoms.len() {
        if structure.atoms[carbon].element != "C" {
            continue;
        }
        let neighbors = structure.neighbors(carbon);
        
        // Check if this carbon has two single-bonded oxygens that are not carbonyl oxygens
        let ether_oxygens: Vec<usize> = neighbors
            .iter()
            .filter(|(atom, order)| {
                if structure.atoms[*atom].element != "O" || !bond_order_matches(*order, 1.0) {
                    return false;
                }
                // Verify oxygen is not part of a carbonyl (not double-bonded to any carbon)
                !structure.neighbors(*atom).iter().any(|(n, o)| {
                    structure.atoms[*n].element == "C" && bond_order_matches(*o, 2.0)
                })
            })
            .map(|(atom, _)| *atom)
            .collect();
        
        if ether_oxygens.len() >= 2 {
            // This is an acetal/ketal
            let mut atoms = vec![carbon];
            atoms.extend(ether_oxygens.iter().copied().take(2));
            
            // Determine if it's an acetal (one hydrogen on central carbon) or ketal (two carbons)
            let carbon_neighbors_count = neighbors
                .iter()
                .filter(|(atom, order)| {
                    structure.atoms[*atom].element == "C" && bond_order_matches(*order, 1.0)
                })
                .count();
            
            // Count hydrogens on the central carbon
            let hydrogen_count = structure.hydrogen_count(carbon);
            
            if hydrogen_count >= 1 && carbon_neighbors_count == 1 {
                groups.push(FunctionalGroup::new(FunctionalGroupType::Acetal, atoms));
            } else if carbon_neighbors_count == 2 {
                groups.push(FunctionalGroup::new(FunctionalGroupType::Ketal, atoms));
            }
        }
    }

    // Detect carbamate-protected amines: N-C(=O)-O-R pattern
    for nitrogen in 0..structure.atoms.len() {
        if structure.atoms[nitrogen].element != "N" {
            continue;
        }
        // Find carbonyl carbon attached to nitrogen
        let carbamate_carbonyl = structure.neighbors(nitrogen).iter().find_map(|(atom, order)| {
            if structure.atoms[*atom].element == "C" && bond_order_matches(*order, 1.0) {
                // Check if this carbon has a double-bonded oxygen (carbonyl)
                let has_carbonyl_oxygen = structure.neighbors(*atom).iter().any(|(n, o)| {
                    structure.atoms[*n].element == "O" && bond_order_matches(*o, 2.0)
                });
                // And a single-bonded oxygen
                let has_ether_oxygen = structure.neighbors(*atom).iter().any(|(n, o)| {
                    structure.atoms[*n].element == "O" && bond_order_matches(*o, 1.0)
                });
                if has_carbonyl_oxygen && has_ether_oxygen {
                    return Some(*atom);
                }
            }
            None
        });
        
        if let Some(carbonyl_carbon) = carbamate_carbonyl {
            // Determine protecting group kind based on the alkyl group on oxygen
            let oxygen = structure.neighbors(carbonyl_carbon).iter().find_map(|(atom, order)| {
                if structure.atoms[*atom].element == "O" && bond_order_matches(*order, 1.0) {
                    Some(*atom)
                } else {
                    None
                }
            });
            
            let protecting_kind = if let Some(oxy) = oxygen {
                // Check what's attached to the oxygen
                let oxy_neighbors = structure.neighbors(oxy);
                // Look for tert-butyl pattern (central carbon with 3 methyl groups)
                let has_tert_butyl = oxy_neighbors.iter().any(|(atom, order)| {
                    if structure.atoms[*atom].element != "C" || !bond_order_matches(*order, 1.0) {
                        return false;
                    }
                    let carbon_neighbors = structure.neighbors(*atom)
                        .iter()
                        .filter(|(n, o)| {
                            structure.atoms[*n].element == "C" && bond_order_matches(*o, 1.0)
                        })
                        .count();
                    carbon_neighbors >= 3
                });
                
                // Look for benzyl pattern (carbon attached to aromatic ring)
                let has_benzyl = oxy_neighbors.iter().any(|(atom, order)| {
                    if structure.atoms[*atom].element != "C" || !bond_order_matches(*order, 1.0) {
                        return false;
                    }
                    structure.neighbors(*atom).iter().any(|(n, o)| {
                        bond_order_matches(*o, 1.5) && structure.atoms[*n].element == "C"
                    })
                });
                
                // Look for fluorenyl pattern (complex polycyclic)
                let has_fluorenyl = oxy_neighbors.iter().any(|(atom, order)| {
                    if structure.atoms[*atom].element != "C" || !bond_order_matches(*order, 1.0) {
                        return false;
                    }
                    // Fluorenyl has multiple aromatic carbons attached
                    let aromatic_neighbors = structure.neighbors(*atom)
                        .iter()
                        .filter(|(n, o)| {
                            bond_order_matches(*o, 1.5) && structure.atoms[*n].element == "C"
                        })
                        .count();
                    aromatic_neighbors >= 2
                });
                
                if has_tert_butyl {
                    FunctionalGroupType::BocCarbamate
                } else if has_fluorenyl {
                    FunctionalGroupType::FmocCarbamate
                } else if has_benzyl {
                    FunctionalGroupType::CbzCarbamate
                } else {
                    FunctionalGroupType::AcylProtectedAmine
                }
            } else {
                FunctionalGroupType::AcylProtectedAmine
            };
            
            groups.push(FunctionalGroup::new(protecting_kind, vec![nitrogen, carbonyl_carbon]));
        }
    }

    // Detect ester-protected acids (already detected as Ester, but mark specifically as protected)
    // These are already found as regular esters, we just need to make sure they don't 
    // also get detected as carboxylic acids
    // The existing logic handles this - esters are separate from acids

    // Detect thioacetals: C with two sulfur atoms
    for carbon in 0..structure.atoms.len() {
        if structure.atoms[carbon].element != "C" {
            continue;
        }
        let sulfur_neighbors: Vec<usize> = structure.neighbors(carbon)
            .iter()
            .filter(|(atom, order)| {
                structure.atoms[*atom].element == "S" && bond_order_matches(*order, 1.0)
            })
            .map(|(atom, _)| *atom)
            .collect();
        
        if sulfur_neighbors.len() >= 2 {
            let mut atoms = vec![carbon];
            atoms.extend(sulfur_neighbors.iter().copied().take(2));
            groups.push(FunctionalGroup::new(FunctionalGroupType::Thioacetal, atoms));
        }
    }

    // Detect protected thiols: C-S-C (thioether) - marked as protected thiol
    for sulfur in 0..structure.atoms.len() {
        if structure.atoms[sulfur].element != "S" {
            continue;
        }
        let carbon_neighbors: Vec<usize> = structure.neighbors(sulfur)
            .iter()
            .filter(|(atom, order)| {
                structure.atoms[*atom].element == "C" && bond_order_matches(*order, 1.0)
            })
            .map(|(atom, _)| *atom)
            .collect();
        
        if carbon_neighbors.len() >= 2 && !structure.neighbors(sulfur).iter().any(|(atom, _)| structure.atoms[*atom].element == "H") {
            // This is a protected thiol (thioether) - no hydrogen on sulfur
            let mut atoms = vec![sulfur];
            atoms.extend(carbon_neighbors.iter().copied().take(2));
            groups.push(FunctionalGroup::new(FunctionalGroupType::ProtectedThiol, atoms));
        }
    }
}

fn bonded(
    structure: &MolecularStructure,
    atom_index: usize,
    element: &str,
    order: Option<f64>,
) -> Vec<usize> {
    structure
        .neighbors(atom_index)
        .into_iter()
        .filter_map(|(neighbor, bond_order)| {
            let atom = &structure.atoms[neighbor];
            let order_matches = order
                .map(|expected| bond_order_matches(bond_order, expected))
                .unwrap_or(true);
            (atom.element == element && order_matches).then_some(neighbor)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::molecule::parse_legacy_structure;

    fn group_types(code: &str) -> Vec<FunctionalGroupType> {
        find_functional_groups(&parse_legacy_structure(code).unwrap())
            .into_iter()
            .map(|group| group.group_type)
            .collect()
    }

    #[test]
    fn detects_carboxylic_acid_and_ester() {
        assert_eq!(
            group_types("destroy:linear:CC(=O)O"),
            vec![FunctionalGroupType::CarboxylicAcid]
        );
        assert_eq!(
            group_types("destroy:linear:CC(=O)OC"),
            vec![FunctionalGroupType::Ester]
        );
    }

    #[test]
    fn detects_alcohol_halide_and_amine() {
        let ethanol = group_types("destroy:linear:CCO");
        assert!(ethanol.contains(&FunctionalGroupType::Alcohol));

        let chloroethane = group_types("destroy:linear:CCCl");
        assert!(chloroethane.contains(&FunctionalGroupType::Halide));

        let ethylamine = group_types("destroy:linear:CCN");
        assert!(ethylamine.contains(&FunctionalGroupType::PrimaryAmine));
        assert_eq!(
            ethylamine
                .iter()
                .filter(|group| **group == FunctionalGroupType::NonTertiaryAmine)
                .count(),
            2
        );
    }

    #[test]
    fn detects_unsaturated_and_nitro_groups() {
        assert!(group_types("destroy:linear:C=C").contains(&FunctionalGroupType::Alkene));
        assert!(group_types("destroy:linear:C#C").contains(&FunctionalGroupType::Alkyne));
        assert!(group_types("destroy:linear:CN(~O)(~O)").contains(&FunctionalGroupType::Nitro));
    }
}
