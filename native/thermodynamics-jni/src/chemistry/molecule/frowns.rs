use super::error::{ChemistryError, ChemistryResult};
use super::molecule::{parse_legacy_structure, MolecularAtom, MolecularBond, MolecularStructure};

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
    super::canonical::canonical_structure_code(structure)
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
