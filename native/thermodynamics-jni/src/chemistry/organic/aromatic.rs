use crate::chemistry::molecule::{bond_order_matches, MolecularStructure};
use std::collections::{BTreeSet, VecDeque};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SubstituentClass {
    StronglyActivating,
    ModeratelyActivating,
    WeaklyActivating,
    WeaklyDeactivating,
    ModeratelyDeactivating,
    StronglyDeactivating,
}

impl SubstituentClass {
    pub(crate) fn effect_for_distance(self, distance: usize) -> f64 {
        match self {
            Self::StronglyActivating => match distance {
                1 => -15.0,
                2 => -3.0,
                _ => -16.0,
            },
            Self::ModeratelyActivating => match distance {
                1 => -10.0,
                2 => -2.0,
                _ => -11.0,
            },
            Self::WeaklyActivating => match distance {
                1 => -4.0,
                2 => -1.0,
                _ => -5.0,
            },
            Self::WeaklyDeactivating => match distance {
                1 => 2.0,
                2 => 6.0,
                _ => 1.5,
            },
            Self::ModeratelyDeactivating => match distance {
                1 => 8.0,
                2 => 4.0,
                _ => 8.5,
            },
            Self::StronglyDeactivating => match distance {
                1 => 15.0,
                2 => 9.0,
                _ => 16.0,
            },
        }
    }

    pub(crate) fn is_deactivating(self) -> bool {
        matches!(self, Self::ModeratelyDeactivating | Self::StronglyDeactivating)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct AromaticSubstituent {
    pub(crate) ring_carbon: usize,
    pub(crate) substituent_atom: usize,
    pub(crate) class: SubstituentClass,
}

#[derive(Debug, Clone)]
pub(crate) struct AromaticRingDescriptor<'a> {
    pub(crate) structure: &'a MolecularStructure,
    pub(crate) ring_atoms: Vec<usize>,
}

impl<'a> AromaticRingDescriptor<'a> {
    pub(crate) fn from_start_carbon(structure: &'a MolecularStructure, start: usize) -> Self {
        let mut visited = BTreeSet::new();
        let mut queue = VecDeque::new();
        queue.push_back(start);
        visited.insert(start);
        while let Some(curr) = queue.pop_front() {
            for (neighbor, order) in structure.neighbors(curr) {
                if bond_order_matches(order, 1.5) && visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
        Self {
            structure,
            ring_atoms: visited.into_iter().collect(),
        }
    }

    pub(crate) fn from_ring_atoms(structure: &'a MolecularStructure, ring_atoms: &[usize]) -> Self {
        Self {
            structure,
            ring_atoms: ring_atoms.to_vec(),
        }
    }

    pub(crate) fn ring_distance(&self, start: usize, end: usize) -> Option<usize> {
        if start == end {
            return Some(0);
        }
        let mut queue = VecDeque::new();
        let mut visited = BTreeSet::new();
        queue.push_back((start, 0));
        visited.insert(start);

        while let Some((curr, dist)) = queue.pop_front() {
            if curr == end {
                return Some(dist);
            }
            for (neighbor, order) in self.structure.neighbors(curr) {
                if self.ring_atoms.contains(&neighbor) && bond_order_matches(order, 1.5) {
                    if visited.insert(neighbor) {
                        queue.push_back((neighbor, dist + 1));
                    }
                }
            }
        }
        None
    }

    pub(crate) fn substituents(&self) -> Vec<AromaticSubstituent> {
        let mut substituents = Vec::new();
        for &ring_car in &self.ring_atoms {
            for (neighbor, order) in self.structure.neighbors(ring_car) {
                if !self.ring_atoms.contains(&neighbor)
                    && self.structure.atoms[neighbor].element != "H"
                    && !bond_order_matches(order, 1.5)
                {
                    let class = classify_substituent(self.structure, ring_car, neighbor);
                    substituents.push(AromaticSubstituent {
                        ring_carbon: ring_car,
                        substituent_atom: neighbor,
                        class,
                    });
                }
            }
        }
        substituents
    }

    pub(crate) fn compute_eas_activation_delta(&self, target_carbon: usize) -> f64 {
        let mut total_delta = 0.0;
        for sub in self.substituents() {
            if let Some(dist) = self.ring_distance(sub.ring_carbon, target_carbon) {
                let mut effect = sub.class.effect_for_distance(dist);
                if dist == 1 && is_substituent_bulky(self.structure, sub.substituent_atom) {
                    effect += 3.0;
                }
                total_delta += effect;
            }
        }
        total_delta
    }

    pub(crate) fn compute_snar_activation_delta(&self, halogen_carbon: usize) -> f64 {
        let mut delta = 0.0;
        for sub in self.substituents() {
            if matches!(
                sub.class,
                SubstituentClass::ModeratelyDeactivating | SubstituentClass::StronglyDeactivating
            ) {
                if let Some(dist) = self.ring_distance(sub.ring_carbon, halogen_carbon) {
                    if dist == 1 || dist == 3 {
                        delta -= match sub.class {
                            SubstituentClass::StronglyDeactivating => 12.0,
                            SubstituentClass::ModeratelyDeactivating => 6.0,
                            _ => 0.0,
                        };
                    }
                }
            }
        }
        delta
    }

    pub(crate) fn is_deactivated_for_fc(&self) -> bool {
        for sub in self.substituents() {
            if sub.class.is_deactivating()
                || matches!(sub.class, SubstituentClass::WeaklyDeactivating)
            {
                return true;
            }
        }
        false
    }
}

fn classify_substituent(
    structure: &MolecularStructure,
    _ring_carbon: usize,
    sub_atom: usize,
) -> SubstituentClass {
    let atom_element = structure.atoms[sub_atom].element.as_str();
    match atom_element {
        "O" => {
            let has_h = structure.neighbors(sub_atom).iter().any(|(n, _)| {
                structure.atoms[*n].element == "H"
            });
            if has_h {
                SubstituentClass::StronglyActivating
            } else {
                let is_ester = structure.neighbors(sub_atom).iter().any(|(n, _)| {
                    structure.atoms[*n].element == "C" && has_double_bonded_oxygen(structure, *n)
                });
                if is_ester {
                    SubstituentClass::ModeratelyActivating
                } else {
                    SubstituentClass::ModeratelyActivating
                }
            }
        }
        "N" => {
            if structure.atoms[sub_atom].charge > 0.1 {
                SubstituentClass::StronglyDeactivating
            } else {
                let is_amide = structure.neighbors(sub_atom).iter().any(|(n, _)| {
                    structure.atoms[*n].element == "C" && has_double_bonded_oxygen(structure, *n)
                });
                if is_amide {
                    SubstituentClass::ModeratelyActivating
                } else {
                    SubstituentClass::StronglyActivating
                }
            }
        }
        "C" => {
            if has_double_bonded_oxygen(structure, sub_atom) {
                SubstituentClass::ModeratelyDeactivating
            } else if has_triple_bonded_nitrogen(structure, sub_atom) {
                SubstituentClass::StronglyDeactivating
            } else if halogen_count(structure, sub_atom) >= 3 {
                SubstituentClass::StronglyDeactivating
            } else {
                SubstituentClass::WeaklyActivating
            }
        }
        "F" | "Cl" | "Br" | "I" => SubstituentClass::WeaklyDeactivating,
        "S" => {
            let double_bonded_oxygens = structure
                .neighbors(sub_atom)
                .iter()
                .filter(|(n, order)| {
                    structure.atoms[*n].element == "O" && bond_order_matches(*order, 2.0)
                })
                .count();
            if double_bonded_oxygens >= 2 {
                SubstituentClass::ModeratelyDeactivating
            } else {
                SubstituentClass::WeaklyActivating
            }
        }
        _ => SubstituentClass::WeaklyActivating,
    }
}

fn has_double_bonded_oxygen(structure: &MolecularStructure, carbon: usize) -> bool {
    structure
        .neighbors(carbon)
        .iter()
        .any(|(n, order)| structure.atoms[*n].element == "O" && bond_order_matches(*order, 2.0))
}

fn has_triple_bonded_nitrogen(structure: &MolecularStructure, carbon: usize) -> bool {
    structure
        .neighbors(carbon)
        .iter()
        .any(|(n, order)| structure.atoms[*n].element == "N" && bond_order_matches(*order, 3.0))
}

fn halogen_count(structure: &MolecularStructure, carbon: usize) -> usize {
    structure
        .neighbors(carbon)
        .iter()
        .filter(|(n, _)| {
            matches!(
                structure.atoms[*n].element.as_str(),
                "F" | "Cl" | "Br" | "I"
            )
        })
        .count()
}

fn is_substituent_bulky(structure: &MolecularStructure, sub_atom: usize) -> bool {
    let element = structure.atoms[sub_atom].element.as_str();
    if element == "Br" || element == "I" {
        return true;
    }
    if element == "C" {
        let non_h_neighbors = structure
            .neighbors(sub_atom)
            .iter()
            .filter(|(n, _)| structure.atoms[*n].element != "H")
            .count();
        return non_h_neighbors >= 3;
    }
    false
}
