use std::collections::VecDeque;

use super::MolecularStructure;

/// If atoms `a` and `b` are connected through the molecular graph by a path that
/// does NOT use a direct `a`–`b` bond, returns the number of atoms in the ring
/// that would form when a new `a`–`b` bond is added (i.e. the smallest such
/// ring). Returns `None` when no alternate path exists, when the indices are out
/// of range, or when `a == b` — in those cases closing the bond would not create
/// a ring.
///
/// Ring size (atom count) equals the number of edges on the shortest alternate
/// path plus one for the new bond. Example: a 4-edge path `C–C–C–C–O` closes
/// into a 5-membered lactone ring.
pub(crate) fn would_form_ring_of_size(
    structure: &MolecularStructure,
    a: usize,
    b: usize,
) -> Option<usize> {
    if a == b || a >= structure.atoms.len() || b >= structure.atoms.len() {
        return None;
    }
    let shortest_edges = shortest_path_len_excluding_direct_edge(structure, a, b)?;
    Some(shortest_edges + 1)
}

/// Shortest path length in edges from `a` to `b`, forbidding the direct `a`–`b`
/// edge as the first step so the result measures the *alternate* route that a
/// new `a`–`b` bond would close into a ring. Returns `None` if `b` is
/// unreachable without the direct edge.
fn shortest_path_len_excluding_direct_edge(
    structure: &MolecularStructure,
    a: usize,
    b: usize,
) -> Option<usize> {
    let n = structure.atoms.len();
    let mut dist = vec![usize::MAX; n];
    dist[a] = 0;
    let mut queue = VecDeque::new();
    queue.push_back(a);
    while let Some(current) = queue.pop_front() {
        for (neighbor, _) in structure.neighbors(current) {
            // Forbid stepping directly across the bond we are testing for closure.
            if current == a && neighbor == b {
                continue;
            }
            if dist[neighbor] == usize::MAX {
                dist[neighbor] = dist[current] + 1;
                queue.push_back(neighbor);
            }
        }
    }
    (dist[b] != usize::MAX).then_some(dist[b])
}
