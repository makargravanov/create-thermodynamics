use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::molecule::{
    bond_order_matches, parse_legacy_structure, MolecularAtom, MolecularBond, MolecularStructure,
};

const DESTROY_NAMESPACE: &str = "destroy";

pub fn parse_frowns(input: &str) -> ChemistryResult<MolecularStructure> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err(invalid_frowns(input, "FROWNS code must not be empty"));
    }
    validate_branch_balance(trimmed)?;
    let full_code = if trimmed.split(':').count() == 3 {
        trimmed.to_string()
    } else if trimmed.contains(':') {
        return Err(invalid_frowns(
            input,
            "FROWNS code must be either <namespace>:<topology>:<body> or a linear body",
        ));
    } else {
        format!("{DESTROY_NAMESPACE}:linear:{trimmed}")
    };
    let parts = full_code.splitn(3, ':').collect::<Vec<_>>();
    if parts[1] == "graph" {
        parse_graph_structure(&full_code, parts[2])
    } else {
        parse_legacy_structure(&full_code)
    }
}

pub fn write_frowns(structure: &MolecularStructure) -> ChemistryResult<String> {
    canonical_structure_code(structure)
}

pub fn canonical_structure_code(structure: &MolecularStructure) -> ChemistryResult<String> {
    structure.validate()?;
    let included = included_atoms(structure);
    if has_cycle(structure, &included) {
        Ok(format!(
            "{DESTROY_NAMESPACE}:graph:{}",
            canonical_graph_body(structure, &included)?
        ))
    } else {
        Ok(format!(
            "{DESTROY_NAMESPACE}:linear:{}",
            canonical_linear_body(structure, &included)?
        ))
    }
}

fn canonical_linear_body(
    structure: &MolecularStructure,
    included: &BTreeSet<usize>,
) -> ChemistryResult<String> {
    let roots = included.iter().copied().collect::<Vec<_>>();
    let mut candidates = Vec::new();
    for root in roots {
        candidates.push(serialize_from(
            structure,
            root,
            None,
            &included,
            &mut BTreeSet::new(),
        )?);
    }
    candidates.into_iter().min().ok_or_else(|| {
        invalid_frowns(
            &structure.source_code,
            "structure has no serializable atoms",
        )
    })
}

fn included_atoms(structure: &MolecularStructure) -> BTreeSet<usize> {
    let mut included = structure
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| (atom.element != "H" || atom.charge != 0.0).then_some(index))
        .collect::<BTreeSet<_>>();
    if included.is_empty() {
        included.extend(0..structure.atoms.len());
    }
    included
}

fn has_cycle(structure: &MolecularStructure, included: &BTreeSet<usize>) -> bool {
    let bond_count = structure
        .bonds
        .iter()
        .filter(|bond| included.contains(&bond.from) && included.contains(&bond.to))
        .count();
    bond_count >= included.len()
}

fn canonical_graph_body(
    structure: &MolecularStructure,
    included: &BTreeSet<usize>,
) -> ChemistryResult<String> {
    let labels = refined_atom_labels(structure, included)?;
    let mut classes: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for atom in included {
        classes
            .entry(labels[*atom].clone())
            .or_default()
            .push(*atom);
    }
    for atoms in classes.values_mut() {
        atoms.sort_by_key(|atom| atom_token(&structure.atoms[*atom]));
    }

    let mut class_values = classes.into_values().collect::<Vec<_>>();
    let permutation_count = class_values
        .iter()
        .map(|class| factorial(class.len()))
        .try_fold(1usize, |acc, value| acc.checked_mul(value))
        .ok_or_else(|| {
            invalid_frowns(
                &structure.source_code,
                "graph has too many equivalent atom permutations to canonicalize safely",
            )
        })?;
    if permutation_count > 100_000 {
        let mut order = included.iter().copied().collect::<Vec<_>>();
        order.sort_by_key(|atom| {
            (
                labels[*atom].clone(),
                atom_token(&structure.atoms[*atom]),
                *atom,
            )
        });
        return graph_code_for_order(structure, &order);
    }

    let mut best: Option<String> = None;
    enumerate_class_orders(&mut class_values, 0, &mut |order| {
        let code = graph_code_for_order(structure, order)?;
        if best.as_ref().map(|current| code < *current).unwrap_or(true) {
            best = Some(code);
        }
        Ok(())
    })?;
    best.ok_or_else(|| invalid_frowns(&structure.source_code, "structure has no graph code"))
}

fn refined_atom_labels(
    structure: &MolecularStructure,
    included: &BTreeSet<usize>,
) -> ChemistryResult<Vec<String>> {
    let mut labels = structure.atoms.iter().map(atom_token).collect::<Vec<_>>();
    for _ in 0..structure.atoms.len() {
        let mut signatures = BTreeMap::new();
        let mut atom_signatures = Vec::new();
        for atom in included {
            let mut neighbors = structure
                .neighbors(*atom)
                .into_iter()
                .filter(|(neighbor, _)| included.contains(neighbor))
                .map(|(neighbor, order)| {
                    Ok(format!("{}{}", graph_bond_token(order)?, labels[neighbor]))
                })
                .collect::<ChemistryResult<Vec<_>>>()?;
            neighbors.sort();
            let signature = format!(
                "{}[{}]",
                atom_token(&structure.atoms[*atom]),
                neighbors.join(",")
            );
            signatures.insert(signature.clone(), String::new());
            atom_signatures.push((*atom, signature));
        }
        for (rank, (_, value)) in signatures.iter_mut().enumerate() {
            *value = rank.to_string();
        }
        let mut next = labels.clone();
        for (atom, signature) in atom_signatures {
            next[atom] = signatures[&signature].clone();
        }
        labels = next;
    }
    Ok(labels)
}

fn enumerate_class_orders<F>(
    classes: &mut [Vec<usize>],
    index: usize,
    callback: &mut F,
) -> ChemistryResult<()>
where
    F: FnMut(&[usize]) -> ChemistryResult<()>,
{
    if index == classes.len() {
        let order = classes.iter().flatten().copied().collect::<Vec<_>>();
        callback(&order)?;
        return Ok(());
    }
    let mut permutations = Vec::new();
    collect_permutations(&classes[index], 0, &mut permutations);
    for permutation in permutations {
        classes[index] = permutation;
        enumerate_class_orders(classes, index + 1, callback)?;
    }
    Ok(())
}

fn collect_permutations(values: &[usize], index: usize, permutations: &mut Vec<Vec<usize>>) {
    let mut current = values.to_vec();
    collect_permutations_inner(&mut current, index, permutations);
}

fn collect_permutations_inner(
    values: &mut Vec<usize>,
    index: usize,
    permutations: &mut Vec<Vec<usize>>,
) {
    if index == values.len() {
        permutations.push(values.clone());
        return;
    }
    for swap_index in index..values.len() {
        values.swap(index, swap_index);
        collect_permutations_inner(values, index + 1, permutations);
        values.swap(index, swap_index);
    }
}

fn factorial(value: usize) -> usize {
    (1..=value).product()
}

fn graph_code_for_order(
    structure: &MolecularStructure,
    order: &[usize],
) -> ChemistryResult<String> {
    let mut remap = vec![None; structure.atoms.len()];
    for (new_index, old_index) in order.iter().enumerate() {
        remap[*old_index] = Some(new_index);
    }

    let atoms = order
        .iter()
        .map(|atom| atom_token(&structure.atoms[*atom]))
        .collect::<Vec<_>>()
        .join(".");
    let mut bonds = Vec::new();
    for bond in &structure.bonds {
        let Some(from) = remap[bond.from] else {
            continue;
        };
        let Some(to) = remap[bond.to] else {
            continue;
        };
        let (from, to) = if from <= to { (from, to) } else { (to, from) };
        bonds.push(format!("{from}-{}-{to}", graph_bond_token(bond.order)?));
    }
    bonds.sort();
    Ok(format!("atoms={atoms};bonds={}", bonds.join(",")))
}

fn serialize_from(
    structure: &MolecularStructure,
    atom: usize,
    parent: Option<usize>,
    included: &BTreeSet<usize>,
    visited: &mut BTreeSet<usize>,
) -> ChemistryResult<String> {
    if !visited.insert(atom) {
        return Ok(format!("@{}", atom_token(&structure.atoms[atom])));
    }
    let mut children = structure
        .neighbors(atom)
        .into_iter()
        .filter(|(neighbor, _)| Some(*neighbor) != parent && included.contains(neighbor))
        .collect::<Vec<_>>();
    children.sort_by_key(|(neighbor, order)| {
        let token = atom_token(&structure.atoms[*neighbor]);
        format!("{}{}", bond_token(*order).unwrap_or("?"), token)
    });
    if children.is_empty() {
        return Ok(atom_token(&structure.atoms[atom]));
    }

    let mut branches = Vec::new();
    for (child, order) in children {
        let child_code = if visited.contains(&child) {
            format!("@{}", atom_token(&structure.atoms[child]))
        } else {
            serialize_from(structure, child, Some(atom), included, visited)?
        };
        branches.push(format!("({}{child_code})", bond_token(order)?));
    }
    branches.sort();
    Ok(format!(
        "{}{}",
        atom_token(&structure.atoms[atom]),
        branches.join("")
    ))
}

fn atom_token(atom: &MolecularAtom) -> String {
    let mut token = atom.element.clone();
    if atom.element == "R" && atom.r_group_number != 0 {
        token.push_str(&atom.r_group_number.to_string());
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

fn parse_graph_structure(source_code: &str, body: &str) -> ChemistryResult<MolecularStructure> {
    let (atoms_part, bonds_part) = body
        .split_once(";bonds=")
        .ok_or_else(|| invalid_frowns(source_code, "graph FROWNS must contain bonds section"))?;
    let atoms_part = atoms_part
        .strip_prefix("atoms=")
        .ok_or_else(|| invalid_frowns(source_code, "graph FROWNS must contain atoms section"))?;
    let atoms = if atoms_part.is_empty() {
        Vec::new()
    } else {
        atoms_part
            .split('.')
            .map(parse_graph_atom)
            .collect::<ChemistryResult<Vec<_>>>()?
    };
    if atoms.is_empty() {
        return Err(invalid_frowns(
            source_code,
            "graph FROWNS must contain atoms",
        ));
    }
    let bonds = if bonds_part.is_empty() {
        Vec::new()
    } else {
        bonds_part
            .split(',')
            .map(|token| parse_graph_bond(source_code, token, atoms.len()))
            .collect::<ChemistryResult<Vec<_>>>()?
    };
    let structure = MolecularStructure {
        source_code: source_code.to_string(),
        atoms,
        bonds,
    };
    structure.validate()?;
    Ok(structure)
}

fn parse_graph_atom(token: &str) -> ChemistryResult<MolecularAtom> {
    let (element_part, charge) = if let Some((element, charge)) = token.split_once('^') {
        (
            element,
            charge
                .parse::<f64>()
                .map_err(|_| invalid_frowns(token, "invalid atom charge"))?,
        )
    } else {
        (token, 0.0)
    };
    if element_part.is_empty() {
        return Err(invalid_frowns(token, "atom token must not be empty"));
    }
    let (element, r_group_number) = if let Some(rest) = element_part.strip_prefix('R') {
        let number = if rest.is_empty() {
            0
        } else {
            rest.parse::<u8>()
                .map_err(|_| invalid_frowns(token, "invalid R group number"))?
        };
        ("R".to_string(), number)
    } else {
        (element_part.to_string(), 0)
    };
    Ok(MolecularAtom {
        element,
        charge,
        r_group_number,
    })
}

fn parse_graph_bond(
    source_code: &str,
    token: &str,
    atom_count: usize,
) -> ChemistryResult<MolecularBond> {
    let parts = token.split('-').collect::<Vec<_>>();
    if parts.len() != 3 {
        return Err(invalid_frowns(source_code, "invalid graph bond token"));
    }
    let from = parts[0]
        .parse::<usize>()
        .map_err(|_| invalid_frowns(source_code, "invalid graph bond start"))?;
    let order = parse_graph_bond_order(parts[1])
        .ok_or_else(|| invalid_frowns(source_code, "invalid graph bond order"))?;
    let to = parts[2]
        .parse::<usize>()
        .map_err(|_| invalid_frowns(source_code, "invalid graph bond end"))?;
    if from >= atom_count || to >= atom_count || from == to {
        return Err(invalid_frowns(
            source_code,
            "graph bond references invalid atom",
        ));
    }
    Ok(MolecularBond { from, to, order })
}

fn bond_token(order: f64) -> ChemistryResult<&'static str> {
    if bond_order_matches(order, 1.0) {
        Ok("")
    } else if bond_order_matches(order, 2.0) {
        Ok("=")
    } else if bond_order_matches(order, 3.0) {
        Ok("#")
    } else if bond_order_matches(order, 1.5) {
        Ok("~")
    } else {
        Err(invalid_frowns(
            "<structure>",
            &format!("unsupported bond order {order}"),
        ))
    }
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
        Err(invalid_frowns(
            "<structure>",
            &format!("unsupported bond order {order}"),
        ))
    }
}

fn parse_graph_bond_order(token: &str) -> Option<f64> {
    match token {
        "s" => Some(1.0),
        "d" => Some(2.0),
        "t" => Some(3.0),
        "a" => Some(1.5),
        "1" => Some(1.0),
        "2" => Some(2.0),
        "3" => Some(3.0),
        "1.5" => Some(1.5),
        _ => None,
    }
}

fn validate_branch_balance(input: &str) -> ChemistryResult<()> {
    let mut depth = 0usize;
    for byte in input.bytes() {
        match byte {
            b'(' => depth += 1,
            b')' => {
                if depth == 0 {
                    return Err(invalid_frowns(
                        input,
                        "branch is closed before it is opened",
                    ));
                }
                depth -= 1;
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(invalid_frowns(input, "branch is not closed"));
    }
    Ok(())
}

fn invalid_frowns(source: &str, reason: &str) -> ChemistryError {
    ChemistryError::InvalidSubstance {
        substance_id: source.to_string(),
        reason: reason.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_short_linear_frowns_and_adds_explicit_hydrogens() {
        let ethanol = parse_frowns("CCO").unwrap();
        assert_eq!(
            ethanol
                .atoms
                .iter()
                .filter(|atom| atom.element == "H" && atom.charge == 0.0)
                .count(),
            6
        );
    }

    #[test]
    fn canonicalizes_linear_direction_and_branches() {
        let first = write_frowns(&parse_frowns("CCO").unwrap()).unwrap();
        let second = write_frowns(&parse_frowns("OCC").unwrap()).unwrap();
        let third = write_frowns(&parse_frowns("C(C)O").unwrap()).unwrap();
        assert_eq!(first, second);
        assert_eq!(first, third);
    }

    #[test]
    fn supports_charges_aromatic_bonds_and_numbered_r_groups() {
        let structure = parse_frowns("R2C(~O^-0.5)O^-0.5").unwrap();
        assert!(structure
            .atoms
            .iter()
            .any(|atom| { atom.element == "R" && atom.r_group_number == 2 }));
        let code = write_frowns(&structure).unwrap();
        assert!(code.contains("R2"));
        assert!(code.contains("^-0.5"));
    }

    #[test]
    fn canonicalizes_topology_reflections_from_source() {
        let first = parse_frowns("destroy:benzene:C,,,,,").unwrap();
        let second = parse_frowns("destroy:benzene:,C,,,,").unwrap();
        assert_eq!(
            write_frowns(&first).unwrap(),
            write_frowns(&second).unwrap()
        );
    }

    #[test]
    fn graph_code_is_independent_of_atom_order() {
        let first = parse_frowns("destroy:graph:atoms=C.C.O;bonds=0-s-1,1-s-2").unwrap();
        let second = parse_frowns("destroy:graph:atoms=O.C.C;bonds=2-s-1,1-s-0").unwrap();

        assert_eq!(
            write_frowns(&first).unwrap(),
            write_frowns(&second).unwrap()
        );
    }

    #[test]
    fn different_connectivity_does_not_collapse() {
        let ethanol_like = parse_frowns("destroy:graph:atoms=C.C.O;bonds=0-s-1,1-s-2").unwrap();
        let ether_like = parse_frowns("destroy:graph:atoms=C.O.C;bonds=0-s-1,1-s-2").unwrap();

        assert_ne!(
            write_frowns(&ethanol_like).unwrap(),
            write_frowns(&ether_like).unwrap()
        );
    }

    #[test]
    fn rejects_bad_frowns() {
        assert!(parse_frowns("destroy:missing:C").is_err());
        assert!(parse_frowns("C(").is_err());
        assert!(parse_frowns("C^bad").is_err());
        assert!(parse_frowns("C(C)(C)(C)(C)(C)").is_err());
    }
}
