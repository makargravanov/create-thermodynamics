use std::collections::VecDeque;

use super::error::{ChemistryError, ChemistryResult};

#[derive(Debug, Clone, PartialEq)]
pub struct MolecularSummary {
    pub molar_mass_grams: f64,
    pub charge: i32,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MolecularAtom {
    pub element: String,
    pub charge: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MolecularBond {
    pub from: usize,
    pub to: usize,
    pub order: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MolecularStructure {
    pub source_code: String,
    pub atoms: Vec<MolecularAtom>,
    pub bonds: Vec<MolecularBond>,
}

impl MolecularStructure {
    pub fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    pub fn summary(&self) -> ChemistryResult<MolecularSummary> {
        self.validate()?;
        let mut mass = 0.0;
        let mut charge = 0.0;
        for atom in &self.atoms {
            mass += element_mass(&atom.element)?;
            charge += atom.charge;
        }
        Ok(MolecularSummary {
            molar_mass_grams: mass,
            charge: charge.round() as i32,
        })
    }

    pub fn validate(&self) -> ChemistryResult<()> {
        if self.atoms.is_empty() {
            return Err(invalid_structure(
                &self.source_code,
                "structure must contain at least one atom",
            ));
        }
        for (index, atom) in self.atoms.iter().enumerate() {
            element_mass(&atom.element)?;
            if !atom.charge.is_finite() {
                return Err(invalid_structure(
                    &self.source_code,
                    &format!("atom {index} has non-finite charge"),
                ));
            }
        }
        for bond in &self.bonds {
            if bond.from >= self.atoms.len() || bond.to >= self.atoms.len() || bond.from == bond.to
            {
                return Err(invalid_structure(
                    &self.source_code,
                    "bond refers to invalid atom index",
                ));
            }
            if !bond.order.is_finite() || bond.order <= 0.0 {
                return Err(invalid_structure(
                    &self.source_code,
                    "bond order must be positive and finite",
                ));
            }
        }
        if !self.is_connected() {
            return Err(invalid_structure(
                &self.source_code,
                "structure is disconnected",
            ));
        }
        let bond_orders = self.bond_orders_by_atom();
        for (index, atom) in self.atoms.iter().enumerate() {
            if atom.element == "R" {
                continue;
            }
            let max_valency = max_valency(&atom.element);
            if bond_orders[index] - atom.charge.abs() > max_valency + 1.0e-6 {
                return Err(invalid_structure(
                    &self.source_code,
                    &format!("atom {index} exceeds valency for {}", atom.element),
                ));
            }
        }
        Ok(())
    }

    pub fn neighbors(&self, atom_index: usize) -> Vec<(usize, f64)> {
        self.bonds
            .iter()
            .filter_map(|bond| {
                if bond.from == atom_index {
                    Some((bond.to, bond.order))
                } else if bond.to == atom_index {
                    Some((bond.from, bond.order))
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn bonded_atoms_by_element(&self, atom_index: usize, element: &str) -> Vec<usize> {
        self.neighbors(atom_index)
            .into_iter()
            .filter_map(|(neighbor, _)| {
                (self.atoms[neighbor].element == element).then_some(neighbor)
            })
            .collect()
    }

    pub fn bonded_atoms_by_element_and_order(
        &self,
        atom_index: usize,
        element: &str,
        order: f64,
    ) -> Vec<usize> {
        self.neighbors(atom_index)
            .into_iter()
            .filter_map(|(neighbor, bond_order)| {
                (self.atoms[neighbor].element == element && bond_order_matches(bond_order, order))
                    .then_some(neighbor)
            })
            .collect()
    }

    pub fn explicit_hydrogen_count(&self, atom_index: usize) -> usize {
        self.bonded_atoms_by_element(atom_index, "H").len()
    }

    pub fn hydrogen_count(&self, atom_index: usize) -> usize {
        self.explicit_hydrogen_count(atom_index)
    }

    pub fn carbon_degree(&self, atom_index: usize) -> usize {
        self.bonded_atoms_by_element(atom_index, "C").len()
    }

    pub fn bond_orders_by_atom(&self) -> Vec<f64> {
        let mut orders = vec![0.0; self.atoms.len()];
        for bond in &self.bonds {
            orders[bond.from] += bond.order;
            orders[bond.to] += bond.order;
        }
        orders
    }

    fn is_connected(&self) -> bool {
        if self.atoms.is_empty() {
            return false;
        }
        let mut seen = vec![false; self.atoms.len()];
        let mut queue = VecDeque::from([0usize]);
        seen[0] = true;
        while let Some(index) = queue.pop_front() {
            for bond in self
                .bonds
                .iter()
                .filter(|bond| bond.from == index || bond.to == index)
            {
                let other = if bond.from == index {
                    bond.to
                } else {
                    bond.from
                };
                if !seen[other] {
                    seen[other] = true;
                    queue.push_back(other);
                }
            }
        }
        seen.into_iter().all(|value| value)
    }
}

pub fn bond_order_matches(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= 1.0e-6
}

#[derive(Debug, Default)]
struct StructureBuilder {
    source_code: String,
    atoms: Vec<MolecularAtom>,
    bonds: Vec<MolecularBond>,
}

impl StructureBuilder {
    fn new(source_code: impl Into<String>) -> Self {
        Self {
            source_code: source_code.into(),
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }

    fn add_atom(&mut self, element: &str, charge: f64) -> ChemistryResult<usize> {
        element_mass(element)?;
        let index = self.atoms.len();
        self.atoms.push(MolecularAtom {
            element: element.to_string(),
            charge,
        });
        Ok(index)
    }

    fn add_bond(&mut self, from: usize, to: usize, order: f64) {
        self.bonds.push(MolecularBond { from, to, order });
    }

    fn finish(self) -> ChemistryResult<MolecularStructure> {
        self.finish_with_normalization(true)
    }

    fn finish_without_normalization(self) -> ChemistryResult<MolecularStructure> {
        self.finish_with_normalization(false)
    }

    fn finish_with_normalization(
        mut self,
        normalize_hydrogens: bool,
    ) -> ChemistryResult<MolecularStructure> {
        if normalize_hydrogens {
            self.add_missing_hydrogens();
        }
        let structure = MolecularStructure {
            source_code: self.source_code,
            atoms: self.atoms,
            bonds: self.bonds,
        };
        structure.validate()?;
        Ok(structure)
    }

    fn add_missing_hydrogens(&mut self) {
        let bond_orders = {
            let structure = MolecularStructure {
                source_code: self.source_code.clone(),
                atoms: self.atoms.clone(),
                bonds: self.bonds.clone(),
            };
            structure.bond_orders_by_atom()
        };
        let original_atom_count = self.atoms.len();
        let mut additions = Vec::new();
        for (index, atom) in self.atoms[..original_atom_count].iter().enumerate() {
            if atom.element == "H" || atom.element == "R" {
                continue;
            }
            let hydrogens = hydrogens_to_add(&atom.element, bond_orders[index], atom.charge);
            for _ in 0..hydrogens as usize {
                additions.push(index);
            }
        }
        for parent in additions {
            let hydrogen = self.atoms.len();
            self.atoms.push(MolecularAtom {
                element: "H".to_string(),
                charge: 0.0,
            });
            self.bonds.push(MolecularBond {
                from: parent,
                to: hydrogen,
                order: 1.0,
            });
        }
    }
}

pub fn parse_legacy_structure(structure_code: &str) -> ChemistryResult<MolecularStructure> {
    let parts = structure_code.splitn(3, ':').collect::<Vec<_>>();
    if parts.len() != 3 {
        return Err(invalid_structure(
            structure_code,
            &format!("bad legacy structure '{structure_code}'"),
        ));
    }
    if parts[1] == "linear" {
        parse_linear_structure(structure_code, parts[2])
    } else {
        parse_topology_structure(structure_code, parts[1], parts[2])
    }
}

pub fn parse_java_structure(code: &str) -> ChemistryResult<MolecularStructure> {
    parse_java_structure_with_normalization(code, true)
}

fn parse_java_structure_with_normalization(
    code: &str,
    normalize_hydrogens: bool,
) -> ChemistryResult<MolecularStructure> {
    let mut builder = StructureBuilder::new(code);
    let mut offset = 0usize;
    let mut root = None;
    let groups = extract_add_group_calls(code)?;
    while let Some(relative) = code[offset..].find("LegacyElement.") {
        let absolute = offset + relative;
        if let Some((_, end, _)) = groups
            .iter()
            .find(|(start, end, _)| absolute >= *start && absolute < *end)
        {
            offset = *end;
            continue;
        }
        let start = offset + relative + "LegacyElement.".len();
        let mut end = start;
        while end < code.len()
            && (code.as_bytes()[end].is_ascii_uppercase() || code.as_bytes()[end] == b'_')
        {
            end += 1;
        }
        let legacy_name = &code[start..end];
        let symbol = legacy_element_symbol(legacy_name)?;
        let tail = &code[end..code.len().min(end + 96)];
        let charge = parse_java_charge(tail, legacy_name)?;
        let bond_order = parse_java_bond_order(tail);
        let index = builder.add_atom(symbol, charge)?;
        if let Some(parent) = root {
            builder.add_bond(parent, index, bond_order);
        } else {
            root = Some(index);
        }
        offset = end;
    }
    let root = root.ok_or_else(|| {
        invalid_structure(code, "java structure must contain at least one root atom")
    })?;
    for (_, _, group_code) in groups {
        let group = parse_java_structure_with_normalization(&group_code, false)?;
        let offset = builder.atoms.len();
        builder.atoms.extend(group.atoms);
        builder
            .bonds
            .extend(group.bonds.into_iter().map(|bond| MolecularBond {
                from: bond.from + offset,
                to: bond.to + offset,
                order: bond.order,
            }));
        builder.add_bond(root, offset, 1.0);
    }
    if normalize_hydrogens {
        builder.finish()
    } else {
        builder.finish_without_normalization()
    }
}

fn extract_add_group_calls(code: &str) -> ChemistryResult<Vec<(usize, usize, String)>> {
    let mut groups = Vec::new();
    let mut offset = 0usize;
    while let Some(relative) = code[offset..].find(".addGroup(") {
        let start = offset + relative;
        let args_start = start + ".addGroup(".len();
        let mut depth = 1usize;
        let mut cursor = args_start;
        while cursor < code.len() && depth > 0 {
            match code.as_bytes()[cursor] {
                b'(' => depth += 1,
                b')' => depth -= 1,
                _ => {}
            }
            cursor += 1;
        }
        if depth != 0 {
            return Err(invalid_structure(code, "unclosed addGroup call"));
        }
        let args = &code[args_start..cursor - 1];
        let group_code = first_top_level_arg(args).ok_or_else(|| {
            invalid_structure(code, "addGroup call must contain a molecular structure")
        })?;
        groups.push((start, cursor, group_code.trim().to_string()));
        offset = cursor;
    }
    Ok(groups)
}

fn first_top_level_arg(args: &str) -> Option<&str> {
    let mut depth = 0usize;
    for (index, byte) in args.bytes().enumerate() {
        match byte {
            b'(' => depth += 1,
            b')' => depth = depth.saturating_sub(1),
            b',' if depth == 0 => return Some(&args[..index]),
            _ => {}
        }
    }
    if args.trim().is_empty() {
        None
    } else {
        Some(args)
    }
}

fn parse_linear_structure(source_code: &str, group: &str) -> ChemistryResult<MolecularStructure> {
    let mut builder = StructureBuilder::new(source_code);
    parse_linear_group_into(&mut builder, group, None, 1.0)?;
    builder.finish()
}

fn parse_linear_group_into(
    builder: &mut StructureBuilder,
    group: &str,
    attachment: Option<usize>,
    external_bond_order: f64,
) -> ChemistryResult<Option<usize>> {
    let mut current = attachment;
    let mut first = None;
    let mut stack: Vec<(usize, f64)> = Vec::new();
    let mut pending_bond = external_bond_order.max(1.0);
    let chars = group.chars().collect::<Vec<_>>();
    let mut i = 0usize;
    while i < chars.len() {
        match chars[i] {
            '(' => {
                if let Some(index) = current {
                    stack.push((index, pending_bond));
                    pending_bond = 1.0;
                }
                i += 1;
            }
            ')' => {
                if let Some((index, restored_bond)) = stack.pop() {
                    current = Some(index);
                    pending_bond = restored_bond;
                } else {
                    current = None;
                }
                i += 1;
            }
            '=' => {
                pending_bond = 2.0;
                i += 1;
            }
            '#' => {
                pending_bond = 3.0;
                i += 1;
            }
            '~' => {
                pending_bond = 1.5;
                i += 1;
            }
            '-' => {
                pending_bond = 1.0;
                i += 1;
            }
            c if c.is_ascii_uppercase() => {
                let start = i;
                i += 1;
                while i < chars.len() && chars[i].is_ascii_lowercase() {
                    i += 1;
                }
                let symbol = &group[start..i];
                let mut charge = 0.0;
                if i < chars.len() && chars[i] == '^' {
                    i += 1;
                    let charge_start = i;
                    if i < chars.len() && (chars[i] == '-' || chars[i] == '+') {
                        i += 1;
                    }
                    while i < chars.len() && (chars[i].is_ascii_digit() || chars[i] == '.') {
                        i += 1;
                    }
                    charge = group[charge_start..i].parse::<f64>().map_err(|_| {
                        invalid_structure(group, &format!("bad charge in '{group}'"))
                    })?;
                }
                let new_index = builder.add_atom(symbol, charge)?;
                if let Some(parent) = current {
                    builder.add_bond(parent, new_index, pending_bond);
                }
                if first.is_none() {
                    first = Some(new_index);
                }
                current = Some(new_index);
                pending_bond = 1.0;
            }
            ',' | ' ' => {
                i += 1;
            }
            other => {
                return Err(invalid_structure(
                    group,
                    &format!("unexpected character '{other}' in '{group}'"),
                ));
            }
        }
    }
    Ok(first)
}

fn parse_topology_structure(
    source_code: &str,
    topology: &str,
    groups: &str,
) -> ChemistryResult<MolecularStructure> {
    let mut builder = StructureBuilder::new(source_code);
    let attachment_sites = add_topology(&mut builder, topology)?;
    for (index, group) in groups.split(',').enumerate() {
        if index >= attachment_sites.len() {
            break;
        }
        if !group.is_empty() {
            parse_linear_group_into(&mut builder, group, Some(attachment_sites[index]), 1.0)?;
        }
    }
    builder.finish()
}

fn add_topology(builder: &mut StructureBuilder, topology: &str) -> ChemistryResult<Vec<usize>> {
    match topology {
        "benzene" => add_ring(builder, "C", 6, 1.5),
        "cubane" => add_cubane(builder),
        "cyclohexene" => add_cyclohexene(builder),
        "cyclopentadienide" => {
            let sites = add_ring(builder, "C", 5, 1.5)?;
            builder.atoms[0].charge = -1.0;
            Ok(sites)
        }
        "diborane" => {
            let b1 = builder.add_atom("B", 0.0)?;
            let h1 = builder.add_atom("H", 0.0)?;
            let b2 = builder.add_atom("B", 0.0)?;
            let h2 = builder.add_atom("H", 0.0)?;
            builder.add_bond(b1, h1, 0.5);
            builder.add_bond(h1, b2, 0.5);
            builder.add_bond(b2, h2, 0.5);
            builder.add_bond(h2, b1, 0.5);
            Ok(vec![b1, h1, b2, h2])
        }
        "octasulfur" => add_ring(builder, "S", 8, 1.0).map(|_| Vec::new()),
        "tetraborate" => add_tetraborate(builder),
        "anthracene" | "anthraquinone" => add_fused_aromatic_14(builder),
        "isohydrobenzofuran" => add_isohydrobenzofuran(builder),
        other => Err(invalid_structure(
            topology,
            &format!("unknown topology '{other}'"),
        )),
    }
}

fn add_ring(
    builder: &mut StructureBuilder,
    element: &str,
    count: usize,
    order: f64,
) -> ChemistryResult<Vec<usize>> {
    let mut sites = Vec::with_capacity(count);
    for _ in 0..count {
        sites.push(builder.add_atom(element, 0.0)?);
    }
    for index in 0..count {
        builder.add_bond(sites[index], sites[(index + 1) % count], order);
    }
    Ok(sites)
}

fn add_cubane(builder: &mut StructureBuilder) -> ChemistryResult<Vec<usize>> {
    let mut sites = Vec::with_capacity(8);
    for _ in 0..8 {
        sites.push(builder.add_atom("C", 0.0)?);
    }
    for (from, to) in [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ] {
        builder.add_bond(sites[from], sites[to], 1.0);
    }
    Ok(sites)
}

fn add_cyclohexene(builder: &mut StructureBuilder) -> ChemistryResult<Vec<usize>> {
    let sites = add_ring(builder, "C", 6, 1.0)?;
    builder.bonds[0].order = 2.0;
    Ok(vec![
        sites[0], sites[1], sites[2], sites[2], sites[3], sites[3], sites[4], sites[4], sites[5],
        sites[5],
    ])
}

fn add_tetraborate(builder: &mut StructureBuilder) -> ChemistryResult<Vec<usize>> {
    let atoms = [
        ("B", -1.0),
        ("B", -1.0),
        ("O", 0.0),
        ("O", 0.0),
        ("O", 0.0),
        ("B", 0.0),
        ("O", 0.0),
        ("O", 0.0),
        ("B", 0.0),
    ];
    let mut indexes = Vec::new();
    for (element, charge) in atoms {
        indexes.push(builder.add_atom(element, charge)?);
    }
    for (from, to) in [
        (0, 2),
        (2, 1),
        (1, 3),
        (3, 5),
        (5, 4),
        (5, 6),
        (0, 7),
        (0, 8),
        (8, 1),
    ] {
        builder.add_bond(indexes[from], indexes[to], 1.0);
    }
    Ok(vec![indexes[4], indexes[6], indexes[7], indexes[8]])
}

fn add_fused_aromatic_14(builder: &mut StructureBuilder) -> ChemistryResult<Vec<usize>> {
    let sites = add_ring(builder, "C", 14, 1.5)?;
    Ok(sites.into_iter().take(10).collect())
}

fn add_isohydrobenzofuran(builder: &mut StructureBuilder) -> ChemistryResult<Vec<usize>> {
    let sites = add_ring(builder, "C", 8, 1.5)?;
    let oxygen = builder.add_atom("O", 0.0)?;
    builder.add_bond(sites[0], oxygen, 1.0);
    builder.add_bond(oxygen, sites[1], 1.0);
    Ok(sites.into_iter().take(6).collect())
}

fn parse_java_charge(tail: &str, legacy_name: &str) -> ChemistryResult<f64> {
    if let Some(stripped) = tail.strip_prefix(",") {
        let number = stripped
            .trim_start()
            .chars()
            .take_while(|c| c.is_ascii_digit() || *c == '-' || *c == '+' || *c == '.')
            .collect::<String>();
        if !number.is_empty() {
            return number.parse::<f64>().map_err(|_| {
                invalid_structure(
                    "<java-structure>",
                    &format!("bad charge near '{legacy_name}'"),
                )
            });
        }
    }
    Ok(0.0)
}

fn parse_java_bond_order(tail: &str) -> f64 {
    if tail.contains("BondType.DOUBLE") {
        2.0
    } else if tail.contains("BondType.TRIPLE") {
        3.0
    } else if tail.contains("BondType.AROMATIC") {
        1.5
    } else {
        1.0
    }
}

fn invalid_structure(source: &str, reason: &str) -> ChemistryError {
    ChemistryError::InvalidSubstance {
        substance_id: source.to_string(),
        reason: reason.to_string(),
    }
}

fn hydrogens_to_add(element: &str, bonds: f64, charge: f64) -> f64 {
    let target = next_lowest_valency(element, bonds);
    (target - charge.abs() - bonds).max(0.0).floor()
}

fn next_lowest_valency(element: &str, bonds: f64) -> f64 {
    let valencies: &[f64] = match element {
        "R" => &[1.0, 2.0, 3.0],
        "C" => &[4.0],
        "H" => &[1.0],
        "S" => &[2.0, 0.0, 4.0, 6.0],
        "N" => &[3.0, 4.0],
        "O" => &[0.0, 2.0],
        "B" => &[3.0],
        "F" | "Na" | "Cl" | "K" | "Ni" | "Zn" | "Zr" | "I" | "Pt" => &[1.0],
        "Ca" | "Hg" => &[2.0],
        "Cr" => &[2.0, 3.0, 6.0],
        "Fe" => &[0.0, 2.0, 3.0],
        "Cu" => &[1.0, 2.0],
        "Au" => &[0.0, 4.0],
        "Pb" => &[2.0, 4.0],
        "Ar" => &[0.0],
        _ => &[0.0],
    };
    valencies
        .iter()
        .copied()
        .find(|v| *v >= bonds)
        .unwrap_or(0.0)
}

fn max_valency(element: &str) -> f64 {
    match element {
        "C" => 4.0,
        "H" => 1.0,
        "S" => 6.0,
        "N" => 4.0,
        "O" => 3.0,
        "B" => 3.0,
        "F" | "Na" | "Cl" | "K" | "Ni" | "Zn" | "Zr" | "I" => 1.0,
        "Pt" => 4.0,
        "Ca" | "Hg" => 2.0,
        "Cr" => 6.0,
        "Fe" => 3.0,
        "Cu" => 2.0,
        "Au" => 4.0,
        "Pb" => 4.0,
        "Ar" => 0.0,
        _ => 0.0,
    }
}

pub fn element_mass(symbol: &str) -> ChemistryResult<f64> {
    let mass = match symbol {
        "R" => 0.0001,
        "C" => 12.01,
        "H" => 1.01,
        "S" => 32.07,
        "N" => 14.01,
        "O" => 16.00,
        "B" => 10.81,
        "F" => 19.00,
        "Na" => 23.00,
        "Cl" => 35.45,
        "K" => 39.10,
        "Ca" => 40.08,
        "Cr" => 52.00,
        "Fe" => 55.85,
        "Ni" => 58.69,
        "Cu" => 63.55,
        "Zn" => 65.38,
        "Zr" => 91.22,
        "I" => 126.90,
        "Pt" => 195.08,
        "Au" => 196.97,
        "Hg" => 200.59,
        "Pb" => 207.20,
        "Ar" => 39.95,
        _ => {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: "<structure>".to_string(),
                reason: format!("unknown element '{symbol}'"),
            })
        }
    };
    Ok(mass)
}

pub fn legacy_element_symbol(name: &str) -> ChemistryResult<&'static str> {
    let symbol = match name {
        "R_GROUP" => "R",
        "CARBON" => "C",
        "HYDROGEN" => "H",
        "SULFUR" => "S",
        "NITROGEN" => "N",
        "OXYGEN" => "O",
        "BORON" => "B",
        "FLUORINE" => "F",
        "SODIUM" => "Na",
        "CHLORINE" => "Cl",
        "POTASSIUM" => "K",
        "CALCIUM" => "Ca",
        "CHROMIUM" => "Cr",
        "IRON" => "Fe",
        "NICKEL" => "Ni",
        "COPPER" => "Cu",
        "ZINC" => "Zn",
        "ZIRCONIUM" => "Zr",
        "IODINE" => "I",
        "PLATINUM" => "Pt",
        "GOLD" => "Au",
        "MERCURY" => "Hg",
        "LEAD" => "Pb",
        "ARGON" => "Ar",
        _ => {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: "<java-structure>".to_string(),
                reason: format!("unknown legacy element '{name}'"),
            })
        }
    };
    Ok(symbol)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_structure_builds_connected_graph() {
        let structure = parse_legacy_structure("destroy:linear:CC(=O)O").unwrap();
        assert_eq!(structure.atom_count(), 8);
        assert_eq!(structure.bond_count(), 7);
        let summary = structure.summary().unwrap();
        assert_eq!(summary.charge, 0);
        assert!((summary.molar_mass_grams - 60.06).abs() < 0.001);
    }

    #[test]
    fn ethanol_has_explicit_hydrogens() {
        let structure = parse_legacy_structure("destroy:linear:CCO").unwrap();
        assert_eq!(structure.atom_count(), 9);
        assert_eq!(structure.bond_count(), 8);
        assert_eq!(
            structure
                .atoms
                .iter()
                .filter(|atom| atom.element == "H")
                .count(),
            6
        );
        let summary = structure.summary().unwrap();
        assert!((summary.molar_mass_grams - 46.08).abs() < 0.001);
    }

    #[test]
    fn benzene_topology_tracks_ring_and_substituent() {
        let structure = parse_legacy_structure("destroy:benzene:C,,,,,").unwrap();
        assert_eq!(structure.atom_count(), 15);
        assert_eq!(structure.bond_count(), 15);
        let summary = structure.summary().unwrap();
        assert!((summary.molar_mass_grams - 92.15).abs() < 0.001);
    }

    #[test]
    fn java_structure_keeps_atoms_and_bonds() {
        let structure = parse_java_structure(
            "LegacyMolecularStructure.atom(LegacyElement.NITROGEN, 1) .addAtom(LegacyElement.OXYGEN, BondType.DOUBLE) .addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1))",
        )
        .unwrap();
        assert_eq!(structure.atom_count(), 3);
        assert_eq!(structure.bond_count(), 2);
        let summary = structure.summary().unwrap();
        assert_eq!(summary.charge, 0);
    }

    #[test]
    fn diborane_topology_keeps_bridge_hydrogens_and_adds_terminal_hydrogens() {
        let structure = parse_legacy_structure("destroy:diborane:,,,").unwrap();
        assert_eq!(
            structure
                .atoms
                .iter()
                .filter(|atom| atom.element == "H")
                .count(),
            6
        );
        let summary = structure.summary().unwrap();
        assert_eq!(summary.charge, 0);
        assert!((summary.molar_mass_grams - 27.68).abs() < 0.001);
    }
}
