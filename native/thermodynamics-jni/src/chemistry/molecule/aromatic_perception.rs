use std::collections::{BTreeMap, BTreeSet};

use super::MolecularStructure;

pub fn aromatize(mut structure: MolecularStructure) -> super::ChemistryResult<MolecularStructure> {
    let cycles = find_simple_cycles_up_to_len(&structure, 10);
    if cycles.is_empty() {
        return Ok(structure);
    }

    let systems = build_fused_ring_systems(&structure, cycles);

    let mut aromatic_bonds = BTreeSet::new();

    for system in systems {
        if !is_conjugated_ring_system(&structure, &system) {
            continue;
        }
        let Some(pi_electrons) = pi_electron_count(&structure, &system) else {
            continue;
        };
        if follows_huckel_rule(pi_electrons) {
            aromatic_bonds.extend(system.bonds.iter().copied());
        }
    }

    for bond in &mut structure.bonds {
        let key = ordered_pair(bond.from, bond.to);
        if aromatic_bonds.contains(&key) {
            bond.order = 1.5;
        }
    }

    Ok(structure)
}

fn ordered_pair(a: usize, b: usize) -> (usize, usize) {
    if a <= b {
        (a, b)
    } else {
        (b, a)
    }
}

fn bond_order_at(structure: &MolecularStructure, a: usize, b: usize) -> Option<f64> {
    structure
        .bonds
        .iter()
        .find(|bond| (bond.from == a && bond.to == b) || (bond.from == b && bond.to == a))
        .map(|bond| bond.order)
}

fn find_simple_cycles_up_to_len(structure: &MolecularStructure, max_len: usize) -> Vec<Vec<usize>> {
    let n = structure.atoms.len();
    let mut raw_cycles = Vec::new();

    for start in 0..n {
        let mut path = Vec::new();
        let mut visited = vec![false; n];
        dfs_find_cycles(
            structure,
            start,
            start,
            max_len,
            &mut path,
            &mut visited,
            &mut raw_cycles,
        );
    }

    dedupe_cycles(raw_cycles)
}

fn dfs_find_cycles(
    structure: &MolecularStructure,
    start: usize,
    current: usize,
    max_len: usize,
    path: &mut Vec<usize>,
    visited: &mut [bool],
    cycles: &mut Vec<Vec<usize>>,
) {
    path.push(current);
    visited[current] = true;

    for (neighbor, _) in structure.neighbors(current) {
        if neighbor < start {
            continue;
        }
        if neighbor == start {
            if path.len() >= 3 && path.len() <= max_len {
                cycles.push(path.clone());
            }
        } else if !visited[neighbor] && path.len() < max_len {
            dfs_find_cycles(structure, start, neighbor, max_len, path, visited, cycles);
        }
    }

    visited[current] = false;
    path.pop();
}

fn dedupe_cycles(cycles: Vec<Vec<usize>>) -> Vec<Vec<usize>> {
    let mut seen = BTreeSet::new();
    let mut unique = Vec::new();

    for cycle in cycles {
        let mut key = cycle.clone();
        key.sort_unstable();
        if seen.insert(key) {
            unique.push(cycle);
        }
    }

    unique
}

#[derive(Debug, Clone)]
#[allow(dead_code)]
struct RingSystem {
    atoms: BTreeSet<usize>,
    bonds: BTreeSet<(usize, usize)>,
    cycles: Vec<Vec<usize>>,
}

fn cycle_bonds(cycle: &[usize]) -> BTreeSet<(usize, usize)> {
    let mut bonds = BTreeSet::new();
    for i in 0..cycle.len() {
        let a = cycle[i];
        let b = cycle[(i + 1) % cycle.len()];
        bonds.insert(ordered_pair(a, b));
    }
    bonds
}

fn cycles_share_bond(a: &BTreeSet<(usize, usize)>, b: &BTreeSet<(usize, usize)>) -> bool {
    a.iter().any(|bond| b.contains(bond))
}

#[derive(Debug, Clone)]
struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSet {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, x: usize) -> usize {
        if self.parent[x] != x {
            let root = self.find(self.parent[x]);
            self.parent[x] = root;
        }
        self.parent[x]
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        match self.rank[ra].cmp(&self.rank[rb]) {
            std::cmp::Ordering::Less => {
                self.parent[ra] = rb;
            }
            std::cmp::Ordering::Greater => {
                self.parent[rb] = ra;
            }
            std::cmp::Ordering::Equal => {
                self.parent[rb] = ra;
                self.rank[ra] += 1;
            }
        }
    }
}

fn build_fused_ring_systems(
    structure: &MolecularStructure,
    cycles: Vec<Vec<usize>>,
) -> Vec<RingSystem> {
    let cycle_bond_sets: Vec<BTreeSet<(usize, usize)>> =
        cycles.iter().map(|cycle| cycle_bonds(cycle)).collect();

    let mut dsu = DisjointSet::new(cycles.len());

    for i in 0..cycles.len() {
        for j in (i + 1)..cycles.len() {
            if cycles_share_bond(&cycle_bond_sets[i], &cycle_bond_sets[j]) {
                dsu.union(i, j);
            }
        }
    }

    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..cycles.len() {
        groups.entry(dsu.find(i)).or_default().push(i);
    }

    let mut systems = Vec::new();

    for (_, indices) in groups {
        let mut atoms = BTreeSet::new();
        let mut bonds = BTreeSet::new();
        let mut grouped_cycles = Vec::new();

        for idx in indices {
            for &atom in &cycles[idx] {
                atoms.insert(atom);
            }
            for &bond in &cycle_bond_sets[idx] {
                if bond_order_at(structure, bond.0, bond.1).is_some() {
                    bonds.insert(bond);
                }
            }
            grouped_cycles.push(cycles[idx].clone());
        }

        systems.push(RingSystem {
            atoms,
            bonds,
            cycles: grouped_cycles,
        });
    }

    systems
}

fn approx_eq(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-6
}

fn is_conjugated_ring_system(structure: &MolecularStructure, system: &RingSystem) -> bool {
    system
        .atoms
        .iter()
        .all(|&atom_idx| has_p_orbital_in_system(structure, atom_idx, system))
}

fn has_p_orbital_in_system(
    structure: &MolecularStructure,
    atom_idx: usize,
    system: &RingSystem,
) -> bool {
    let atom = &structure.atoms[atom_idx];
    let element = atom.element.as_str();
    let charge = atom.charge;

    let neighbors = structure.neighbors(atom_idx);
    let ring_neighbors = neighbors
        .iter()
        .filter(|(n, _)| system.atoms.contains(n))
        .count();

    let has_endocyclic_pi_bond = neighbors
        .iter()
        .any(|(neighbor, order)| system.atoms.contains(neighbor) && *order >= 1.5);

    if has_endocyclic_pi_bond {
        return matches!(element, "C" | "N" | "O" | "S" | "B" | "P");
    }

    match element {
        "C" => approx_eq(charge, 1.0) || approx_eq(charge, -1.0),
        "B" => approx_eq(charge, 0.0) && ring_neighbors >= 2,
        "N" => (approx_eq(charge, 0.0) && ring_neighbors >= 2) || approx_eq(charge, -1.0),
        "O" | "S" => approx_eq(charge, 0.0) && ring_neighbors >= 2,
        _ => false,
    }
}

fn pi_electron_contribution(
    structure: &MolecularStructure,
    atom_idx: usize,
    system: &RingSystem,
) -> Option<u32> {
    let atom = &structure.atoms[atom_idx];
    let element = atom.element.as_str();
    let charge = atom.charge;

    let neighbors = structure.neighbors(atom_idx);
    let has_endocyclic_pi_bond = neighbors
        .iter()
        .any(|(neighbor, order)| system.atoms.contains(neighbor) && *order >= 1.5);

    match element {
        "C" => {
            if approx_eq(charge, 1.0) {
                Some(0)
            } else if approx_eq(charge, -1.0) && !has_endocyclic_pi_bond {
                Some(2)
            } else {
                Some(1)
            }
        }
        "N" => {
            if approx_eq(charge, 1.0) {
                Some(1)
            } else if has_endocyclic_pi_bond {
                Some(1)
            } else if approx_eq(charge, 0.0) || approx_eq(charge, -1.0) {
                Some(2)
            } else {
                None
            }
        }
        "O" | "S" => {
            if has_endocyclic_pi_bond {
                Some(1)
            } else if approx_eq(charge, 0.0) {
                Some(2)
            } else {
                None
            }
        }
        "B" => {
            if approx_eq(charge, 0.0) || approx_eq(charge, 1.0) {
                Some(0)
            } else {
                None
            }
        }
        _ => None,
    }
}

fn pi_electron_count(structure: &MolecularStructure, system: &RingSystem) -> Option<u32> {
    let mut total = 0;
    for &atom_idx in &system.atoms {
        let contribution = pi_electron_contribution(structure, atom_idx, system)?;
        total += contribution;
    }
    Some(total)
}

fn follows_huckel_rule(pi_electrons: u32) -> bool {
    pi_electrons >= 2 && (pi_electrons - 2) % 4 == 0
}

#[allow(dead_code)]
fn mark_system_aromatic(structure: &mut MolecularStructure, system: &RingSystem) {
    for bond in &mut structure.bonds {
        let key = ordered_pair(bond.from, bond.to);
        if system.bonds.contains(&key) {
            bond.order = 1.5;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::canonical::canonical_structure_code;
    use crate::chemistry::frowns::parse_frowns;
    use crate::chemistry::molecule::bond_order_matches;

    fn count_aromatic_bonds(structure: &MolecularStructure, ring_indices: &[usize]) -> usize {
        structure
            .bonds
            .iter()
            .filter(|b| {
                ring_indices.contains(&b.from)
                    && ring_indices.contains(&b.to)
                    && bond_order_matches(b.order, 1.5)
            })
            .count()
    }

    #[test]
    fn aromatizes_benzene_from_alternating_bonds() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-0,\
                   0-s-6,1-s-7,2-s-8,3-s-9,4-s-10,5-s-11",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4, 5]), 6);
    }

    #[test]
    fn does_not_aromatize_cyclobutadiene() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-0,\
                   0-s-4,1-s-5,2-s-6,3-s-7",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3]), 0);
    }

    #[test]
    fn aromatizes_pyridine() {
        let structure = parse_frowns(
            "destroy:graph:atoms=N.C.C.C.C.C.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-0,\
                   1-s-6,2-s-7,3-s-8,4-s-9,5-s-10",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4, 5]), 6);
    }

    #[test]
    fn aromatizes_pyrrole() {
        let structure = parse_frowns(
            "destroy:graph:atoms=N.C.C.C.C.H.H.H.H.H;\
             bonds=0-s-1,1-d-2,2-s-3,3-d-4,4-s-0,\
                   0-s-5,1-s-6,2-s-7,3-s-8,4-s-9",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4]), 5);
    }

    #[test]
    fn aromatizes_furan() {
        let structure = parse_frowns(
            "destroy:graph:atoms=O.C.C.C.C.H.H.H.H;\
             bonds=0-s-1,1-d-2,2-s-3,3-d-4,4-s-0,\
                   1-s-5,2-s-6,3-s-7,4-s-8",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4]), 5);
    }

    #[test]
    fn aromatizes_thiophene() {
        let structure = parse_frowns(
            "destroy:graph:atoms=S.C.C.C.C.H.H.H.H;\
             bonds=0-s-1,1-d-2,2-s-3,3-d-4,4-s-0,\
                   1-s-5,2-s-6,3-s-7,4-s-8",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4]), 5);
    }

    #[test]
    fn aromatizes_naphthalene_as_fused_system() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.C.C.C.C.H.H.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-0,\
                   5-s-6,6-d-7,7-s-8,8-d-9,9-s-4,\
                   0-s-10,1-s-11,2-s-12,3-s-13,6-s-14,7-s-15,8-s-16,9-s-17",
        )
        .unwrap();

        assert_eq!(
            count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
            11
        );
    }

    #[test]
    fn aromatizes_anthracene() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.C.C.C.C.C.C.C.C.H.H.H.H.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-0,\
                   5-s-6,6-d-7,7-s-8,8-d-9,9-s-4,\
                   9-s-10,10-d-11,11-s-12,12-d-13,13-s-8,\
                   0-s-14,1-s-15,2-s-16,3-s-17,\
                   6-s-18,7-s-19,10-s-20,11-s-21,12-s-22,13-s-23",
        )
        .unwrap();

        let rings: Vec<usize> = (0..14).collect();
        let count = count_aromatic_bonds(&structure, &rings);
        assert_eq!(count, 16);
    }

    #[test]
    fn aromatizes_phenanthrene() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.C.C.C.C.C.C.C.C.H.H.H.H.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-0,\
                   3-s-6,6-d-7,7-s-8,8-d-9,9-s-2,\
                   5-s-10,10-d-11,11-s-12,12-d-13,13-s-0,\
                   0-s-14,1-s-15,2-s-16,4-s-17,6-s-18,7-s-19,8-s-20,9-s-21,\
                   10-s-22,11-s-23",
        )
        .unwrap();

        let rings: Vec<usize> = (0..14).collect();
        let count = count_aromatic_bonds(&structure, &rings);
        assert_eq!(count, 16);
    }

    #[test]
    fn aromatizes_cyclopentadienyl_anion() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C^-1.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-s-0,\
                   0-s-5,1-s-6,2-s-7,3-s-8",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4]), 5);
    }

    #[test]
    fn aromatizes_tropylium_cation() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.C^1.H.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-6,6-s-0,\
                   0-s-7,1-s-8,2-s-9,3-s-10,4-s-11,5-s-12,6-s-13",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4, 5, 6]), 7);
    }

    #[test]
    fn does_not_aromatize_cyclooctatetraene() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.C.C.H.H.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-6,6-d-7,7-s-0,\
                   0-s-8,1-s-9,2-s-10,3-s-11,4-s-12,5-s-13,6-s-14,7-s-15",
        )
        .unwrap();

        assert_eq!(
            count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4, 5, 6, 7]),
            0
        );
    }

    #[test]
    fn does_not_aromatize_cyclopentadienyl_cation() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C^1.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-s-0,\
                   0-s-5,1-s-6,2-s-7,3-s-8",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4]), 0);
    }

    #[test]
    fn does_not_aromatize_non_conjugated_cyclohexadiene() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-s-5,5-s-0,\
                   0-s-6,1-s-7,2-s-8,3-s-9,4-s-10,5-s-11",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4, 5]), 0);
    }

    #[test]
    fn does_not_aromatize_ring_with_silicon() {
        let structure = parse_frowns(
            "destroy:graph:atoms=Si.C.C.C.C.C.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-0,\
                   1-s-6,2-s-7,3-s-8,4-s-9,5-s-10",
        )
        .unwrap();

        assert_eq!(count_aromatic_bonds(&structure, &[0, 1, 2, 3, 4, 5]), 0);
    }

    #[test]
    fn aromatize_is_idempotent() {
        let structure = parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.H.H.H.H.H.H;\
             bonds=0-d-1,1-s-2,2-d-3,3-s-4,4-d-5,5-s-0,\
                   0-s-6,1-s-7,2-s-8,3-s-9,4-s-10,5-s-11",
        )
        .unwrap();

        let once = aromatize(structure.clone()).unwrap();
        let twice = aromatize(once.clone()).unwrap();

        assert_eq!(
            canonical_structure_code(&once).unwrap(),
            canonical_structure_code(&twice).unwrap()
        );
    }
}
