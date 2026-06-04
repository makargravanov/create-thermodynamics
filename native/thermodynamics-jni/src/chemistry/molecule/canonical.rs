use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};

use super::error::{ChemistryError, ChemistryResult};
use super::molecule::{
    bond_order_matches, DoubleBondStereo, MolecularAtom, MolecularStructure, StereoDescriptor,
    StereoMixtureKind, Stereochemistry, TetrahedralStereo, ValenceSaturation,
};

const DESTROY_NAMESPACE: &str = "destroy";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CanonicalizationProgress {
    pub stage: CanonicalizationStage,
    pub processed_branches: usize,
    pub pending_branches: usize,
    pub best_code_found: bool,
    pub estimated_remaining_branches: Option<usize>,
    pub completed: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CanonicalizationStage {
    Prepare,
    Tree,
    Refine,
    Search,
    Complete,
}

#[derive(Debug, Clone)]
pub struct CanonicalizationJob {
    graph: CanonicalGraph,
    stack: Vec<Vec<Vec<usize>>>,
    best_body: Option<String>,
    processed_branches: usize,
    stage: CanonicalizationStage,
    completed: bool,
    result: Option<String>,
}

pub struct CanonicalizationHandle {
    progress: Arc<Mutex<CanonicalizationProgress>>,
    result: Arc<Mutex<Option<ChemistryResult<String>>>>,
    worker: Option<JoinHandle<()>>,
}

impl CanonicalizationHandle {
    pub fn progress(&self) -> CanonicalizationProgress {
        self.progress
            .lock()
            .expect("canonicalization progress mutex must not be poisoned")
            .clone()
    }

    pub fn result(&self) -> Option<ChemistryResult<String>> {
        self.result
            .lock()
            .expect("canonicalization result mutex must not be poisoned")
            .clone()
    }

    pub fn join(mut self) -> ChemistryResult<String> {
        if let Some(worker) = self.worker.take() {
            worker
                .join()
                .map_err(|_| ChemistryError::InvalidSubstance {
                    substance_id: "<canonicalization>".to_string(),
                    reason: "canonicalization worker panicked".to_string(),
                })?;
        }
        self.result()
            .ok_or_else(|| ChemistryError::InvalidSubstance {
                substance_id: "<canonicalization>".to_string(),
                reason: "canonicalization worker finished without a result".to_string(),
            })?
    }
}

pub fn spawn_canonicalization(structure: MolecularStructure) -> CanonicalizationHandle {
    let progress = Arc::new(Mutex::new(CanonicalizationProgress {
        stage: CanonicalizationStage::Prepare,
        processed_branches: 0,
        pending_branches: 0,
        best_code_found: false,
        estimated_remaining_branches: None,
        completed: false,
    }));
    let result = Arc::new(Mutex::new(None));
    let worker_progress = Arc::clone(&progress);
    let worker_result = Arc::clone(&result);
    let worker = thread::spawn(move || {
        let outcome = (|| {
            let mut job = CanonicalizationJob::new(&structure)?;
            loop {
                let current = job.progress();
                *worker_progress
                    .lock()
                    .expect("canonicalization progress mutex must not be poisoned") = current;
                if job.progress().completed {
                    return job.run_to_completion();
                }
                job.step()?;
            }
        })();
        *worker_result
            .lock()
            .expect("canonicalization result mutex must not be poisoned") = Some(outcome);
        let mut progress = worker_progress
            .lock()
            .expect("canonicalization progress mutex must not be poisoned");
        progress.completed = true;
        progress.stage = CanonicalizationStage::Complete;
    });

    CanonicalizationHandle {
        progress,
        result,
        worker: Some(worker),
    }
}

impl CanonicalizationJob {
    pub fn new(structure: &MolecularStructure) -> ChemistryResult<Self> {
        structure.validate()?;
        let graph = CanonicalGraph::new(structure)?;
        if !graph.has_cycle() && graph.stereochemistry.is_empty() {
            let result = format!(
                "{DESTROY_NAMESPACE}:linear:{}",
                graph.canonical_tree_body()?
            );
            return Ok(Self {
                graph,
                stack: Vec::new(),
                best_body: None,
                processed_branches: 0,
                stage: CanonicalizationStage::Tree,
                completed: true,
                result: Some(result),
            });
        }

        Ok(Self {
            stack: vec![graph.initial_partition()],
            graph,
            best_body: None,
            processed_branches: 0,
            stage: CanonicalizationStage::Prepare,
            completed: false,
            result: None,
        })
    }

    pub fn progress(&self) -> CanonicalizationProgress {
        CanonicalizationProgress {
            stage: if self.completed {
                CanonicalizationStage::Complete
            } else {
                self.stage.clone()
            },
            processed_branches: self.processed_branches,
            pending_branches: self.stack.len(),
            best_code_found: self.best_body.is_some() || self.result.is_some(),
            estimated_remaining_branches: (!self.completed).then_some(self.stack.len()),
            completed: self.completed,
        }
    }

    pub fn result(&self) -> Option<&str> {
        self.result.as_deref()
    }

    pub fn run_to_completion(&mut self) -> ChemistryResult<String> {
        while !self.completed {
            self.step()?;
        }
        self.result
            .clone()
            .ok_or_else(|| ChemistryError::InvalidSubstance {
                substance_id: "<canonicalization>".to_string(),
                reason: "canonicalization finished without a result".to_string(),
            })
    }

    pub fn step(&mut self) -> ChemistryResult<CanonicalizationProgress> {
        if self.completed {
            return Ok(self.progress());
        }
        self.stage = CanonicalizationStage::Search;
        let Some(partition) = self.stack.pop() else {
            let body = self
                .best_body
                .clone()
                .ok_or_else(|| ChemistryError::InvalidSubstance {
                    substance_id: self.graph.source_code.clone(),
                    reason: "cyclic graph has no canonical ordering".to_string(),
                })?;
            self.result = Some(format!("{DESTROY_NAMESPACE}:graph:{body}"));
            self.completed = true;
            return Ok(self.progress());
        };

        self.processed_branches += 1;
        self.stage = CanonicalizationStage::Refine;
        let refined = self.graph.refine_partition(partition)?;
        if self.graph.partition_is_discrete(&refined) {
            let order = refined.iter().map(|cell| cell[0]).collect::<Vec<_>>();
            let body = self.graph.graph_body_for_order(&order)?;
            if self
                .best_body
                .as_ref()
                .map(|current| body < *current)
                .unwrap_or(true)
            {
                self.best_body = Some(body);
            }
            return Ok(self.progress());
        }

        self.stage = CanonicalizationStage::Search;
        let cell_index = refined
            .iter()
            .enumerate()
            .filter(|(_, cell)| cell.len() > 1)
            .min_by_key(|(_, cell)| cell.len())
            .map(|(index, _)| index)
            .ok_or_else(|| ChemistryError::InvalidSubstance {
                substance_id: self.graph.source_code.clone(),
                reason: "non-discrete partition has no ambiguous cell".to_string(),
            })?;
        let mut candidates = refined[cell_index].clone();
        candidates.sort_by_key(|atom| self.graph.vertex_sort_key(*atom));
        for atom in candidates.into_iter().rev() {
            let mut next = refined.clone();
            let remainder = next[cell_index]
                .iter()
                .copied()
                .filter(|candidate| *candidate != atom)
                .collect::<Vec<_>>();
            next.splice(cell_index..=cell_index, [vec![atom], remainder]);
            self.stack.push(next);
        }
        Ok(self.progress())
    }
}

pub fn canonical_structure_code(structure: &MolecularStructure) -> ChemistryResult<String> {
    CanonicalizationJob::new(structure)?.run_to_completion()
}

#[derive(Debug, Clone)]
struct CanonicalGraph {
    source_code: String,
    vertices: Vec<CanonicalVertex>,
    edges: Vec<CanonicalEdge>,
    adjacency: Vec<Vec<(usize, String)>>,
    stereochemistry: Vec<Stereochemistry>,
}

#[derive(Debug, Clone)]
struct CanonicalVertex {
    token: String,
    neutral_hydrogens: usize,
    original_index: usize,
}

#[derive(Debug, Clone)]
struct CanonicalEdge {
    from: usize,
    to: usize,
    order: String,
}

impl CanonicalGraph {
    fn new(structure: &MolecularStructure) -> ChemistryResult<Self> {
        let included = included_atoms(structure);
        let mut old_to_new = vec![None; structure.atoms.len()];
        let mut vertices = Vec::new();
        for old_index in included.iter().copied() {
            old_to_new[old_index] = Some(vertices.len());
            vertices.push(CanonicalVertex {
                token: atom_token(&structure.atoms[old_index]),
                neutral_hydrogens: 0,
                original_index: old_index,
            });
        }

        for (atom_index, atom) in structure.atoms.iter().enumerate() {
            if atom.element != "H" || atom.charge != 0.0 || included.contains(&atom_index) {
                continue;
            }
            let included_neighbors = structure
                .neighbors(atom_index)
                .into_iter()
                .filter_map(|(neighbor, _)| old_to_new[neighbor])
                .collect::<Vec<_>>();
            if included_neighbors.len() == 1 {
                vertices[included_neighbors[0]].neutral_hydrogens += 1;
            }
        }

        let mut edges = Vec::new();
        let mut adjacency = vec![Vec::new(); vertices.len()];
        for bond in &structure.bonds {
            let Some(from) = old_to_new[bond.from] else {
                continue;
            };
            let Some(to) = old_to_new[bond.to] else {
                continue;
            };
            let order = graph_bond_token(bond.order)?.to_string();
            let (from, to) = if from <= to { (from, to) } else { (to, from) };
            edges.push(CanonicalEdge {
                from,
                to,
                order: order.clone(),
            });
            adjacency[from].push((to, order.clone()));
            adjacency[to].push((from, order));
        }
        for neighbors in &mut adjacency {
            neighbors.sort();
        }
        let stereochemistry = structure
            .stereochemistry
            .iter()
            .filter_map(|stereo| remap_stereochemistry(stereo, &old_to_new))
            .collect::<Vec<_>>();

        Ok(Self {
            source_code: structure.source_code.clone(),
            vertices,
            edges,
            adjacency,
            stereochemistry,
        })
    }

    fn has_cycle(&self) -> bool {
        self.edges.len() >= self.vertices.len()
    }

    fn initial_partition(&self) -> Vec<Vec<usize>> {
        let mut classes: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for vertex in 0..self.vertices.len() {
            classes
                .entry(self.initial_label(vertex))
                .or_default()
                .push(vertex);
        }
        let mut partition = classes.into_values().collect::<Vec<_>>();
        for cell in &mut partition {
            cell.sort_by_key(|atom| self.vertex_sort_key(*atom));
        }
        partition
    }

    fn initial_label(&self, vertex: usize) -> String {
        let mut incident_orders = self.adjacency[vertex]
            .iter()
            .map(|(_, order)| order.clone())
            .collect::<Vec<_>>();
        incident_orders.sort();
        format!(
            "{}|h{}|d{}|{}",
            self.vertices[vertex].token,
            self.vertices[vertex].neutral_hydrogens,
            self.adjacency[vertex].len(),
            incident_orders.join(".")
        )
    }

    fn refine_partition(&self, mut partition: Vec<Vec<usize>>) -> ChemistryResult<Vec<Vec<usize>>> {
        loop {
            let class_by_vertex = class_by_vertex(self.vertices.len(), &partition)?;
            let mut changed = false;
            let mut next = Vec::new();
            for cell in partition {
                let mut split: BTreeMap<String, Vec<usize>> = BTreeMap::new();
                for vertex in cell {
                    split
                        .entry(self.refinement_signature(vertex, &class_by_vertex))
                        .or_default()
                        .push(vertex);
                }
                if split.len() > 1 {
                    changed = true;
                }
                for mut split_cell in split.into_values() {
                    split_cell.sort_by_key(|atom| self.vertex_sort_key(*atom));
                    next.push(split_cell);
                }
            }
            partition = next;
            if !changed {
                return Ok(partition);
            }
        }
    }

    fn refinement_signature(&self, vertex: usize, class_by_vertex: &[usize]) -> String {
        let mut neighbors = self.adjacency[vertex]
            .iter()
            .map(|(neighbor, order)| format!("{order}:{}", class_by_vertex[*neighbor]))
            .collect::<Vec<_>>();
        neighbors.sort();
        format!("{}[{}]", self.initial_label(vertex), neighbors.join(","))
    }

    fn partition_is_discrete(&self, partition: &[Vec<usize>]) -> bool {
        partition.len() == self.vertices.len() && partition.iter().all(|cell| cell.len() == 1)
    }

    fn vertex_sort_key(&self, vertex: usize) -> String {
        format!(
            "{}|h{}|{}",
            self.vertices[vertex].token,
            self.vertices[vertex].neutral_hydrogens,
            self.vertices[vertex].original_index
        )
    }

    fn graph_body_for_order(&self, order: &[usize]) -> ChemistryResult<String> {
        let mut remap = vec![None; self.vertices.len()];
        for (new_index, old_index) in order.iter().enumerate() {
            remap[*old_index] = Some(new_index);
        }
        let atoms = order
            .iter()
            .map(|vertex| self.vertices[*vertex].token.clone())
            .collect::<Vec<_>>()
            .join(".");
        let mut bonds = Vec::new();
        for edge in &self.edges {
            let from = remap[edge.from].ok_or_else(|| ChemistryError::InvalidSubstance {
                substance_id: self.source_code.clone(),
                reason: "canonical order misses an edge endpoint".to_string(),
            })?;
            let to = remap[edge.to].ok_or_else(|| ChemistryError::InvalidSubstance {
                substance_id: self.source_code.clone(),
                reason: "canonical order misses an edge endpoint".to_string(),
            })?;
            let (from, to) = if from <= to { (from, to) } else { (to, from) };
            bonds.push(format!("{from}-{}-{to}", edge.order));
        }
        bonds.sort();
        let stereo = self.stereochemistry_body_for_order(&remap)?;
        if stereo.is_empty() {
            Ok(format!("atoms={atoms};bonds={}", bonds.join(",")))
        } else {
            Ok(format!(
                "atoms={atoms};bonds={};stereo={}",
                bonds.join(","),
                stereo.join(",")
            ))
        }
    }

    fn stereochemistry_body_for_order(
        &self,
        remap: &[Option<usize>],
    ) -> ChemistryResult<Vec<String>> {
        let mut tokens = self
            .stereochemistry
            .iter()
            .map(|stereo| canonical_stereo_token(stereo, remap, &self.source_code))
            .collect::<ChemistryResult<Vec<_>>>()?;
        tokens.sort();
        Ok(tokens)
    }

    fn canonical_tree_body(&self) -> ChemistryResult<String> {
        let centers = self.tree_centers()?;
        centers
            .into_iter()
            .map(|center| self.rooted_tree_code(center, None))
            .collect::<ChemistryResult<Vec<_>>>()?
            .into_iter()
            .min()
            .ok_or_else(|| ChemistryError::InvalidSubstance {
                substance_id: self.source_code.clone(),
                reason: "tree has no canonical center".to_string(),
            })
    }

    fn tree_centers(&self) -> ChemistryResult<Vec<usize>> {
        if self.vertices.is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: self.source_code.clone(),
                reason: "graph has no vertices".to_string(),
            });
        }
        if self.vertices.len() <= 2 {
            return Ok((0..self.vertices.len()).collect());
        }
        let mut degree = self.adjacency.iter().map(Vec::len).collect::<Vec<_>>();
        let mut leaves = degree
            .iter()
            .enumerate()
            .filter_map(|(vertex, degree)| (*degree <= 1).then_some(vertex))
            .collect::<VecDeque<_>>();
        let mut remaining = self.vertices.len();
        while remaining > 2 {
            let leaf_count = leaves.len();
            if leaf_count == 0 {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.source_code.clone(),
                    reason: "cyclic graph reached tree center calculation".to_string(),
                });
            }
            remaining -= leaf_count;
            for _ in 0..leaf_count {
                let leaf = leaves.pop_front().expect("leaf count was measured");
                degree[leaf] = 0;
                for (neighbor, _) in &self.adjacency[leaf] {
                    if degree[*neighbor] == 0 {
                        continue;
                    }
                    degree[*neighbor] -= 1;
                    if degree[*neighbor] == 1 {
                        leaves.push_back(*neighbor);
                    }
                }
            }
        }
        let mut centers = leaves.into_iter().collect::<Vec<_>>();
        if centers.is_empty() {
            centers = degree
                .iter()
                .enumerate()
                .filter_map(|(vertex, degree)| (*degree > 0).then_some(vertex))
                .collect();
        }
        centers.sort_unstable();
        centers.dedup();
        Ok(centers)
    }

    fn rooted_tree_code(&self, vertex: usize, parent: Option<usize>) -> ChemistryResult<String> {
        let mut branches = Vec::new();
        for (neighbor, order) in &self.adjacency[vertex] {
            if Some(*neighbor) == parent {
                continue;
            }
            let child = self.rooted_tree_code(*neighbor, Some(vertex))?;
            branches.push(format!("{}{}", linear_bond_token(order)?, child));
        }
        branches.sort();
        if branches.is_empty() {
            return Ok(self.vertices[vertex].token.clone());
        }
        Ok(format!(
            "{}{}",
            self.vertices[vertex].token,
            branches
                .into_iter()
                .map(|branch| format!("({branch})"))
                .collect::<Vec<_>>()
                .join("")
        ))
    }
}

fn class_by_vertex(vertex_count: usize, partition: &[Vec<usize>]) -> ChemistryResult<Vec<usize>> {
    let mut result = vec![usize::MAX; vertex_count];
    for (class, cell) in partition.iter().enumerate() {
        for vertex in cell {
            if *vertex >= vertex_count || result[*vertex] != usize::MAX {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: "<canonicalization>".to_string(),
                    reason: "invalid canonical partition".to_string(),
                });
            }
            result[*vertex] = class;
        }
    }
    if result.iter().any(|class| *class == usize::MAX) {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: "<canonicalization>".to_string(),
            reason: "canonical partition does not cover every vertex".to_string(),
        });
    }
    Ok(result)
}

fn included_atoms(structure: &MolecularStructure) -> BTreeSet<usize> {
    let mut included = structure
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| (atom.element != "H" || atom.charge != 0.0).then_some(index))
        .collect::<BTreeSet<_>>();
    for stereo in &structure.stereochemistry {
        for atom in stereochemistry_atoms(stereo) {
            included.insert(atom);
        }
    }
    if included.is_empty() {
        included.extend(0..structure.atoms.len());
    }
    included
}

fn stereochemistry_atoms(stereo: &Stereochemistry) -> Vec<usize> {
    match stereo {
        Stereochemistry::Tetrahedral(tetrahedral) => {
            let mut atoms = vec![tetrahedral.center];
            atoms.extend(tetrahedral.substituents);
            atoms
        }
        Stereochemistry::DoubleBond(double_bond) => vec![
            double_bond.first,
            double_bond.second,
            double_bond.first_substituent,
            double_bond.second_substituent,
        ],
        Stereochemistry::Mixture { atoms, .. } => atoms.clone(),
    }
}

fn remap_stereochemistry(
    stereo: &Stereochemistry,
    mapping: &[Option<usize>],
) -> Option<Stereochemistry> {
    match stereo {
        Stereochemistry::Tetrahedral(tetrahedral) => {
            let center = mapping.get(tetrahedral.center).copied().flatten()?;
            let mut substituents = [0usize; 4];
            for (slot, substituent) in substituents
                .iter_mut()
                .zip(tetrahedral.substituents.iter().copied())
            {
                *slot = mapping.get(substituent).copied().flatten()?;
            }
            Some(Stereochemistry::Tetrahedral(TetrahedralStereo {
                center,
                substituents,
                descriptor: tetrahedral.descriptor,
            }))
        }
        Stereochemistry::DoubleBond(double_bond) => {
            Some(Stereochemistry::DoubleBond(DoubleBondStereo {
                first: mapping.get(double_bond.first).copied().flatten()?,
                second: mapping.get(double_bond.second).copied().flatten()?,
                first_substituent: mapping
                    .get(double_bond.first_substituent)
                    .copied()
                    .flatten()?,
                second_substituent: mapping
                    .get(double_bond.second_substituent)
                    .copied()
                    .flatten()?,
                descriptor: double_bond.descriptor,
            }))
        }
        Stereochemistry::Mixture { atoms, kind } => {
            let atoms = atoms
                .iter()
                .map(|atom| mapping.get(*atom).copied().flatten())
                .collect::<Option<Vec<_>>>()?;
            Some(Stereochemistry::Mixture { atoms, kind: *kind })
        }
    }
}

fn canonical_stereo_token(
    stereo: &Stereochemistry,
    remap: &[Option<usize>],
    source_code: &str,
) -> ChemistryResult<String> {
    match stereo {
        Stereochemistry::Tetrahedral(tetrahedral) => {
            let center = remapped_atom(remap, tetrahedral.center, source_code)?;
            let ordered = tetrahedral
                .substituents
                .iter()
                .map(|atom| remapped_atom(remap, *atom, source_code))
                .collect::<ChemistryResult<Vec<_>>>()?;
            let mut sorted = ordered.clone();
            sorted.sort_unstable();
            let descriptor = if permutation_is_odd(&ordered, &sorted)? {
                opposite_tetrahedral_descriptor(tetrahedral.descriptor)?
            } else {
                tetrahedral.descriptor
            };
            Ok(format!(
                "t:{center}:{}:{}",
                sorted
                    .into_iter()
                    .map(|atom| atom.to_string())
                    .collect::<Vec<_>>()
                    .join("."),
                stereo_descriptor_token(descriptor)
            ))
        }
        Stereochemistry::DoubleBond(double_bond) => {
            let first = remapped_atom(remap, double_bond.first, source_code)?;
            let second = remapped_atom(remap, double_bond.second, source_code)?;
            let first_substituent =
                remapped_atom(remap, double_bond.first_substituent, source_code)?;
            let second_substituent =
                remapped_atom(remap, double_bond.second_substituent, source_code)?;
            let (first, second, first_substituent, second_substituent) = if first <= second {
                (first, second, first_substituent, second_substituent)
            } else {
                (second, first, second_substituent, first_substituent)
            };
            Ok(format!(
                "db:{first}={second}:{first_substituent}-{second_substituent}:{}",
                stereo_descriptor_token(double_bond.descriptor)
            ))
        }
        Stereochemistry::Mixture { atoms, kind } => {
            let mut atoms = atoms
                .iter()
                .map(|atom| remapped_atom(remap, *atom, source_code))
                .collect::<ChemistryResult<Vec<_>>>()?;
            atoms.sort_unstable();
            Ok(format!(
                "mix:{}:{}",
                stereo_mixture_kind_token(*kind),
                atoms
                    .into_iter()
                    .map(|atom| atom.to_string())
                    .collect::<Vec<_>>()
                    .join(".")
            ))
        }
    }
}

fn remapped_atom(
    remap: &[Option<usize>],
    atom: usize,
    source_code: &str,
) -> ChemistryResult<usize> {
    remap
        .get(atom)
        .copied()
        .flatten()
        .ok_or_else(|| ChemistryError::InvalidSubstance {
            substance_id: source_code.to_string(),
            reason: "stereochemistry references an atom outside canonical order".to_string(),
        })
}

fn permutation_is_odd(ordered: &[usize], sorted: &[usize]) -> ChemistryResult<bool> {
    if ordered.len() != sorted.len() {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: "<stereochemistry>".to_string(),
            reason: "stereo permutation has inconsistent length".to_string(),
        });
    }
    let mut positions = ordered
        .iter()
        .map(|atom| {
            sorted
                .iter()
                .position(|candidate| candidate == atom)
                .ok_or_else(|| ChemistryError::InvalidSubstance {
                    substance_id: "<stereochemistry>".to_string(),
                    reason: "stereo permutation contains an unknown atom".to_string(),
                })
        })
        .collect::<ChemistryResult<Vec<_>>>()?;
    let mut swaps = 0usize;
    for index in 0..positions.len() {
        while positions[index] != index {
            let target = positions[index];
            positions.swap(index, target);
            swaps += 1;
        }
    }
    Ok(swaps % 2 == 1)
}

fn opposite_tetrahedral_descriptor(
    descriptor: StereoDescriptor,
) -> ChemistryResult<StereoDescriptor> {
    match descriptor {
        StereoDescriptor::Clockwise => Ok(StereoDescriptor::CounterClockwise),
        StereoDescriptor::CounterClockwise => Ok(StereoDescriptor::Clockwise),
        _ => Err(ChemistryError::InvalidSubstance {
            substance_id: "<stereochemistry>".to_string(),
            reason: "tetrahedral stereo has a non-tetrahedral descriptor".to_string(),
        }),
    }
}

fn stereo_descriptor_token(descriptor: StereoDescriptor) -> &'static str {
    match descriptor {
        StereoDescriptor::Clockwise => "cw",
        StereoDescriptor::CounterClockwise => "ccw",
        StereoDescriptor::Cis => "cis",
        StereoDescriptor::Trans => "trans",
        StereoDescriptor::E => "e",
        StereoDescriptor::Z => "z",
    }
}

fn stereo_mixture_kind_token(kind: StereoMixtureKind) -> &'static str {
    match kind {
        StereoMixtureKind::Tetrahedral => "tetra",
        StereoMixtureKind::DoubleBond => "db",
        StereoMixtureKind::General => "general",
    }
}

fn atom_token(atom: &MolecularAtom) -> String {
    let mut token = atom.element.clone();
    if atom.element == "R" && atom.r_group_number != 0 {
        token.push_str(&atom.r_group_number.to_string());
    }
    if atom.valence_saturation == ValenceSaturation::UnsaturatedAllowed {
        token.push('!');
    }
    if atom.charge != 0.0 {
        token.push('^');
        if (atom.charge - atom.charge.round()).abs() <= 1.0e-9 {
            token.push_str(&(atom.charge.round() as i64).to_string());
        } else {
            token.push_str(&atom.charge.to_string());
        }
    }
    token
}

fn graph_bond_token(order: f64) -> ChemistryResult<&'static str> {
    if bond_order_matches(order, 1.0) {
        Ok("1")
    } else if bond_order_matches(order, 2.0) {
        Ok("2")
    } else if bond_order_matches(order, 3.0) {
        Ok("3")
    } else if bond_order_matches(order, 1.5) {
        Ok("1.5")
    } else {
        Err(ChemistryError::InvalidSubstance {
            substance_id: "<structure>".to_string(),
            reason: format!("unsupported bond order {order}"),
        })
    }
}

fn linear_bond_token(order: &str) -> ChemistryResult<&'static str> {
    match order {
        "1" => Ok(""),
        "2" => Ok("="),
        "3" => Ok("#"),
        "1.5" => Ok("~"),
        _ => Err(ChemistryError::InvalidSubstance {
            substance_id: "<structure>".to_string(),
            reason: format!("unsupported bond order {order}"),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;
    use crate::chemistry::molecule::{MolecularAtom, MolecularBond};

    #[test]
    fn stepwise_job_reports_progress_until_completion() {
        let structure = parse_frowns("destroy:benzene:C,,,,,").unwrap();
        let mut job = CanonicalizationJob::new(&structure).unwrap();
        let mut previous = job.progress().processed_branches;
        while !job.progress().completed {
            let progress = job.step().unwrap();
            assert!(progress.processed_branches >= previous);
            previous = progress.processed_branches;
        }
        assert_eq!(
            job.result().unwrap(),
            canonical_structure_code(&structure).unwrap()
        );
    }

    #[test]
    fn large_symmetric_cycle_canonicalizes_without_factorial_storage() {
        let atom_count = 24;
        let atoms = (0..atom_count)
            .map(|_| MolecularAtom {
                element: "C".to_string(),
                charge: 0.0,
                r_group_number: 0,
                valence_saturation: ValenceSaturation::Saturate,
            })
            .collect::<Vec<_>>();
        let bonds = (0..atom_count)
            .map(|index| MolecularBond {
                from: index,
                to: (index + 1) % atom_count,
                order: 1.0,
            })
            .collect::<Vec<_>>();
        let structure = MolecularStructure {
            source_code: "test:large-cycle".to_string(),
            atoms,
            bonds,
            stereochemistry: Vec::new(),
        };

        let code = canonical_structure_code(&structure).unwrap();
        assert!(code.starts_with("destroy:graph:"));
    }

    #[test]
    fn charged_hydrogen_and_r_group_change_graph_code() {
        let protonated = MolecularStructure {
            source_code: "test:charged-hydrogen".to_string(),
            atoms: vec![
                MolecularAtom {
                    element: "C".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 1.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
            ],
            bonds: vec![MolecularBond {
                from: 0,
                to: 1,
                order: 1.0,
            }],
            stereochemistry: Vec::new(),
        };
        let r_group = MolecularStructure {
            source_code: "test:r-group".to_string(),
            atoms: vec![
                MolecularAtom {
                    element: "C".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
                MolecularAtom {
                    element: "R".to_string(),
                    charge: 0.0,
                    r_group_number: 1,
                    valence_saturation: ValenceSaturation::Saturate,
                },
            ],
            bonds: vec![MolecularBond {
                from: 0,
                to: 1,
                order: 1.0,
            }],
            stereochemistry: Vec::new(),
        };

        assert_ne!(
            canonical_structure_code(&protonated).unwrap(),
            canonical_structure_code(&r_group).unwrap()
        );
    }

    #[test]
    fn symmetric_destroy_topologies_are_canonicalized_stably() {
        for (first, second) in [
            ("destroy:cubane:C,,,,,,", "destroy:cubane:,,C,,,,"),
            ("destroy:octasulfur:hello", "destroy:octasulfur:anything"),
            (
                "destroy:cyclohexene:C,,,,,,,,",
                "destroy:cyclohexene:,C,,,,,,,",
            ),
        ] {
            assert_eq!(
                canonical_structure_code(&parse_frowns(first).unwrap()).unwrap(),
                canonical_structure_code(&parse_frowns(second).unwrap()).unwrap()
            );
        }
    }

    #[test]
    fn background_canonicalization_returns_same_result_as_sync() {
        let structure = parse_frowns("destroy:benzene:C,,,,,").unwrap();
        let expected = canonical_structure_code(&structure).unwrap();
        let handle = spawn_canonicalization(structure);
        let result = handle.join().unwrap();

        assert_eq!(result, expected);
    }
}
