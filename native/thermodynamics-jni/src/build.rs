use std::env;
use std::fs;
use std::collections::HashSet;
use std::path::PathBuf;

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let destroy_reactions = manifest_dir
        .join("..")
        .join("..")
        .join("Destroy")
        .join("src")
        .join("main")
        .join("java")
        .join("com")
        .join("petrolpark")
        .join("destroy")
        .join("chemistry")
        .join("legacy")
        .join("index")
        .join("DestroyReactions.java");
    let destroy_molecules = manifest_dir
        .join("..")
        .join("..")
        .join("Destroy")
        .join("src")
        .join("main")
        .join("java")
        .join("com")
        .join("petrolpark")
        .join("destroy")
        .join("chemistry")
        .join("legacy")
        .join("index")
        .join("DestroyMolecules.java");

    println!("cargo:rerun-if-changed={}", destroy_reactions.display());
    println!("cargo:rerun-if-changed={}", destroy_molecules.display());
    println!("cargo:rerun-if-changed=build.rs");

    let source = fs::read_to_string(&destroy_reactions)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", destroy_reactions.display()));
    let molecule_source = fs::read_to_string(&destroy_molecules)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", destroy_molecules.display()));

    let explicit_blocks = split_reaction_blocks(&source);
    let parsed_id_counts = count_parsed_ids(&explicit_blocks);
    let molecule_id_map = parse_molecule_id_map(&molecule_source);
    let acid_entries = parse_acid_entries(&source, &molecule_id_map);
    let mut emitted_ids = HashSet::new();

    let mut generated = String::new();
    generated.push_str("use super::error::ChemistryResult;\n");
    generated.push_str("use super::reaction::Reaction;\n");
    generated.push_str("use super::registry::ChemistryRegistryBuilder;\n\n");
    generated.push_str(&format!(
        "pub const DESTROY_EXPLICIT_REACTION_COUNT: usize = {};\n",
        explicit_blocks.len()
    ));
    generated.push_str(&format!(
        "pub const DESTROY_REGISTERED_REACTION_COUNT: usize = {};\n\n",
        explicit_blocks.len() + acid_entries.len() * 3
    ));
    generated.push_str(
        "pub fn destroy_reactions_registry_builder(mut builder: ChemistryRegistryBuilder) -> ChemistryResult<ChemistryRegistryBuilder> {\n",
    );

    for block in explicit_blocks {
        let line = emit_reaction_builder(
            &block,
            &parsed_id_counts,
            &molecule_id_map,
            &mut emitted_ids,
        );
        generated.push_str(&line);
    }
    for acid in acid_entries {
        generated.push_str(&emit_acid_reactions(&acid));
    }

    generated.push_str("    Ok(builder)\n}\n");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    let out_file = out_dir.join("destroy_reactions.rs");
    fs::write(&out_file, generated)
        .unwrap_or_else(|error| panic!("failed to write {}: {error}", out_file.display()));
}

#[derive(Debug, Clone)]
struct ReactionBlock {
    text: String,
}

#[derive(Debug, Clone)]
struct AcidEntry {
    acid: String,
    conjugate_base: String,
    pka: f64,
}

fn split_reaction_blocks(source: &str) -> Vec<ReactionBlock> {
    let mut blocks = Vec::new();
    let mut current: Option<String> = None;

    for line in source.lines() {
        if line.contains("= builder()") {
            current = Some(String::new());
        }
        if let Some(block) = current.as_mut() {
            block.push_str(line);
            block.push('\n');
            if line.trim_start().contains(".build()") {
                blocks.push(ReactionBlock {
                    text: block.clone(),
                });
                current = None;
            }
        }
    }

    blocks
}

fn parse_acid_entries(
    source: &str,
    molecule_id_map: &std::collections::HashMap<String, String>,
) -> Vec<AcidEntry> {
    let mut entries = Vec::new();
    for line in source.lines() {
        let trimmed = line.trim();
        if !trimmed.starts_with("builder().acid(") {
            continue;
        }
        let args = trimmed
            .strip_prefix("builder().acid(")
            .and_then(|value| value.strip_suffix(");"))
            .unwrap_or(trimmed);
        let parts = split_top_level_commas(args);
        if parts.len() != 3 {
            panic!("failed to parse acid declaration: {trimmed}");
        }
        entries.push(AcidEntry {
            acid: destroy_constant_to_id(parts[0].trim(), molecule_id_map),
            conjugate_base: destroy_constant_to_id(parts[1].trim(), molecule_id_map),
            pka: parse_float(parts[2].trim()),
        });
    }
    entries
}

fn count_parsed_ids(blocks: &[ReactionBlock]) -> std::collections::HashMap<String, usize> {
    let mut counts = std::collections::HashMap::new();
    for block in blocks {
        let text = strip_call(&strip_call(&block.text, ".reverseReaction("), ".withResult(");
        if let Some(id) = extract_single_call_arg(&text, ".id(") {
            let id = unquote_java_string(&id);
            *counts.entry(id).or_insert(0) += 1;
        }
    }
    counts
}

fn emit_reaction_builder(
    block: &ReactionBlock,
    parsed_id_counts: &std::collections::HashMap<String, usize>,
    molecule_id_map: &std::collections::HashMap<String, String>,
    emitted_ids: &mut HashSet<String>,
) -> String {
    let text = &block.text;
    let without_reverse = strip_call(text, ".reverseReaction(");
    let cleaned = strip_call(&without_reverse, ".withResult(");
    let id = extract_single_call_arg(&cleaned, ".id(")
        .unwrap_or_else(|| panic!("missing reaction id in block:\n{cleaned}"));
    let id = unquote_java_string(&id);
    let constant_name = reaction_constant_name(text)
        .unwrap_or_else(|| panic!("missing reaction constant in block:\n{cleaned}"));
    let parsed_reaction_id = format!("destroy:{id}");
    let constant_reaction_id = format!("destroy:{}", constant_name.to_ascii_lowercase());
    let reaction_id = if parsed_id_counts.get(&id).copied().unwrap_or(0) > 1
        && constant_reaction_id != parsed_reaction_id
    {
        constant_reaction_id
    } else {
        parsed_reaction_id
    };
    if emitted_ids.contains(&reaction_id) {
        panic!("duplicate reaction id '{reaction_id}'");
    }
    emitted_ids.insert(reaction_id.clone());

    let mut out = String::new();
    out.push_str("    builder = builder.reaction(\n");
    out.push_str(&format!(
        "        Reaction::builder({:?})\n",
        reaction_id
    ));

    for args in extract_call_args(&cleaned, ".addReactant(") {
        let parts = split_top_level_commas(&args);
        if parts.is_empty() {
            continue;
        }
        let substance = destroy_molecule_to_id(parts[0].trim(), molecule_id_map);
        match parts.len() {
            1 => {
                out.push_str(&format!(
                    "            .reactant({:?}, 1, 1)\n",
                    substance
                ));
            }
            2 => {
                let coefficient = parse_u32(parts[1].trim());
                if coefficient == 0 {
                    out.push_str(&format!(
                        "            .catalyst_order({:?}, 0)\n",
                        substance
                    ));
                } else {
                    out.push_str(&format!(
                        "            .reactant({:?}, {}, {})\n",
                        substance, coefficient, coefficient
                    ));
                }
            }
            _ => {
                let coefficient = parse_u32(parts[1].trim());
                let order = parse_u32(parts[2].trim());
                if coefficient == 0 {
                    out.push_str(&format!(
                        "            .catalyst_order({:?}, {})\n",
                        substance, order
                    ));
                } else {
                    out.push_str(&format!(
                        "            .reactant({:?}, {}, {})\n",
                        substance, coefficient, order
                    ));
                }
            }
        }
    }

    for args in extract_call_args(&cleaned, ".addCatalyst(") {
        let parts = split_top_level_commas(&args);
        if parts.len() != 2 {
            panic!("invalid catalyst declaration in block:\n{cleaned}");
        }
        let substance = destroy_molecule_to_id(parts[0].trim(), molecule_id_map);
        let order = parse_u32(parts[1].trim());
        out.push_str(&format!(
            "            .catalyst_order({:?}, {})\n",
            substance, order
        ));
    }

    for args in extract_call_args(&cleaned, ".addProduct(") {
        let parts = split_top_level_commas(&args);
        if parts.is_empty() {
            continue;
        }
        let substance = destroy_molecule_to_id(parts[0].trim(), molecule_id_map);
        let coefficient = if parts.len() == 1 {
            1
        } else {
            parse_u32(parts[1].trim())
        };
        out.push_str(&format!(
            "            .product({:?}, {})\n",
            substance, coefficient
        ));
    }

    for args in extract_call_args(&cleaned, ".addSimpleItemReactant(") {
        let parts = split_top_level_commas(&args);
        let moles = parse_float(parts.last().unwrap().trim());
        out.push_str(&format!(
            "            .external_reactant({:?}, {:?})\n",
            format!("addSimpleItemReactant({args})"),
            moles
        ));
    }

    for args in extract_call_args(&cleaned, ".addSimpleItemTagReactant(") {
        let parts = split_top_level_commas(&args);
        let moles = parse_float(parts.last().unwrap().trim());
        out.push_str(&format!(
            "            .external_reactant({:?}, {:?})\n",
            format!("addSimpleItemTagReactant({args})"),
            moles
        ));
    }

    for args in extract_call_args(&cleaned, ".addSimpleItemCatalyst(") {
        let parts = split_top_level_commas(&args);
        let moles = parse_float(parts.last().unwrap().trim());
        out.push_str(&format!(
            "            .external_catalyst({:?}, {:?})\n",
            format!("addSimpleItemCatalyst({args})"),
            moles
        ));
    }

    for args in extract_call_args(&cleaned, ".addSimpleItemTagCatalyst(") {
        let parts = split_top_level_commas(&args);
        let moles = parse_float(parts.last().unwrap().trim());
        out.push_str(&format!(
            "            .external_catalyst({:?}, {:?})\n",
            format!("addSimpleItemTagCatalyst({args})"),
            moles
        ));
    }

    for args in extract_call_args(&without_reverse, ".withResult(") {
        let parts = split_top_level_commas(&args);
        let moles = parse_float(parts.first().unwrap().trim());
        out.push_str(&format!(
            "            .reaction_result({:?}, {:?})\n",
            format!("withResult({args})"),
            moles
        ));
    }

    if cleaned.contains(".requireUV()") {
        out.push_str("            .requires_uv()\n");
    }
    if cleaned.contains(".reversible()") || text.contains(".reverseReaction(") {
        out.push_str("            .display_as_reversible()\n");
    }
    if cleaned.contains(".dontIncludeInJei()") {
        out.push_str("            .show_in_jei(false)\n");
    }
    for args in extract_call_args(text, ".includeInJeiIf(") {
        out.push_str(&format!(
            "            .show_in_jei_condition({:?})\n",
            format!("includeInJeiIf({args})")
        ));
    }

    if let Some(value) = extract_single_call_arg(&cleaned, ".preexponentialFactor(") {
        out.push_str(&format!(
            "            .pre_exponential_factor({:?})\n",
            parse_float(&value)
        ));
    }
    if let Some(value) = extract_single_call_arg(&cleaned, ".activationEnergy(") {
        out.push_str(&format!(
            "            .activation_energy_kj_per_mol({:?})\n",
            parse_float(&value)
        ));
    }
    if let Some(value) = extract_single_call_arg(&cleaned, ".enthalpyChange(") {
        out.push_str(&format!(
            "            .enthalpy_change_kj_per_mol({:?})\n",
            parse_float(&value)
        ));
    }

    out.push_str("            .allow_mass_imbalance()\n");

    out.push_str("            .build(),\n");
    out.push_str("    );\n");
    out
}

fn emit_acid_reactions(acid: &AcidEntry) -> String {
    let mut out = String::new();
    let acid_id = acid.acid.strip_prefix("destroy:").unwrap_or(&acid.acid);
    let dissociation_id = format!("destroy:{acid_id}.dissociation");
    let neutralization_id = format!("destroy:{acid_id}.neutralization");
    let association_id = format!("destroy:{acid_id}.association");
    let rate = 0.5 * 10f64.powf(-acid.pka);
    let room_temperature_energy = 8.314_462_618_153_24 * 0.298;

    out.push_str("    builder = builder.reaction(\n");
    out.push_str(&format!(
        "        Reaction::builder({:?})\n",
        dissociation_id
    ));
    out.push_str(&format!("            .reactant({:?}, 1, 1)\n", acid.acid));
    out.push_str("            .catalyst_order(\"destroy:water\", 1)\n");
    out.push_str("            .product(\"destroy:proton\", 1)\n");
    out.push_str(&format!("            .product({:?}, 1)\n", acid.conjugate_base));
    out.push_str(&format!("            .activation_energy_kj_per_mol({room_temperature_energy:?})\n"));
    out.push_str(&format!("            .pre_exponential_factor({rate:?})\n"));
    out.push_str("            .show_in_jei(false)\n");
    out.push_str("            .allow_mass_imbalance()\n");
    out.push_str("            .build(),\n");
    out.push_str("    );\n");

    out.push_str("    builder = builder.reaction(\n");
    out.push_str(&format!(
        "        Reaction::builder({:?})\n",
        neutralization_id
    ));
    out.push_str(&format!("            .reactant({:?}, 1, 1)\n", acid.acid));
    out.push_str("            .reactant(\"destroy:hydroxide\", 1, 1)\n");
    out.push_str(&format!("            .product({:?}, 1)\n", acid.conjugate_base));
    out.push_str("            .product(\"destroy:water\", 1)\n");
    out.push_str(&format!("            .activation_energy_kj_per_mol({room_temperature_energy:?})\n"));
    out.push_str(&format!("            .pre_exponential_factor({rate:?})\n"));
    out.push_str("            .show_in_jei(false)\n");
    out.push_str("            .allow_mass_imbalance()\n");
    out.push_str("            .build(),\n");
    out.push_str("    );\n");

    out.push_str("    builder = builder.reaction(\n");
    out.push_str(&format!(
        "        Reaction::builder({:?})\n",
        association_id
    ));
    out.push_str(&format!("            .reactant({:?}, 1, 1)\n", acid.conjugate_base));
    out.push_str("            .reactant(\"destroy:proton\", 1, 1)\n");
    out.push_str(&format!("            .product({:?}, 1)\n", acid.acid));
    out.push_str(&format!("            .activation_energy_kj_per_mol({room_temperature_energy:?})\n"));
    out.push_str("            .pre_exponential_factor(1.0)\n");
    out.push_str("            .show_in_jei(false)\n");
    out.push_str("            .allow_mass_imbalance()\n");
    out.push_str("            .build(),\n");
    out.push_str("    );\n");

    out
}

fn strip_call(text: &str, needle: &str) -> String {
    let mut output = String::with_capacity(text.len());
    let mut cursor = 0usize;
    while let Some(relative) = text[cursor..].find(needle) {
        let start = cursor + relative;
        output.push_str(&text[cursor..start]);
        let open = start + needle.len() - 1;
        let close = matching_paren(text, open)
            .unwrap_or_else(|| panic!("unbalanced call for {needle} in:\n{text}"));
        cursor = close + 1;
    }
    output.push_str(&text[cursor..]);
    output
}

fn extract_call_args(text: &str, needle: &str) -> Vec<String> {
    let mut values = Vec::new();
    let mut cursor = 0usize;
    while let Some(relative) = text[cursor..].find(needle) {
        let start = cursor + relative;
        let open = start + needle.len() - 1;
        let close = matching_paren(text, open)
            .unwrap_or_else(|| panic!("unbalanced call for {needle} in:\n{text}"));
        values.push(text[open + 1..close].trim().to_string());
        cursor = close + 1;
    }
    values
}

fn extract_single_call_arg(text: &str, needle: &str) -> Option<String> {
    extract_call_args(text, needle).into_iter().next()
}

fn matching_paren(text: &str, open_index: usize) -> Option<usize> {
    let bytes = text.as_bytes();
    if bytes.get(open_index) != Some(&b'(') {
        return None;
    }
    let mut depth = 0i32;
    for index in open_index..bytes.len() {
        match bytes[index] {
            b'(' => depth += 1,
            b')' => {
                depth -= 1;
                if depth == 0 {
                    return Some(index);
                }
            }
            _ => {}
        }
    }
    None
}

fn split_top_level_commas(value: &str) -> Vec<String> {
    let mut parts = Vec::new();
    let mut depth = 0i32;
    let mut start = 0usize;
    for (index, ch) in value.char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => depth -= 1,
            ',' if depth == 0 => {
                parts.push(value[start..index].trim().to_string());
                start = index + ch.len_utf8();
            }
            _ => {}
        }
    }
    parts.push(value[start..].trim().to_string());
    parts
}

fn parse_float(value: &str) -> f64 {
    let trimmed = value.trim().trim_end_matches(['f', 'F']);
    trimmed
        .parse::<f64>()
        .unwrap_or_else(|error| panic!("failed to parse float '{value}': {error}"))
}

fn parse_u32(value: &str) -> u32 {
    value
        .trim()
        .parse::<u32>()
        .unwrap_or_else(|error| panic!("failed to parse integer '{value}': {error}"))
}

fn destroy_molecule_to_id(
    value: &str,
    molecule_id_map: &std::collections::HashMap<String, String>,
) -> String {
    destroy_constant_to_id(value, molecule_id_map)
}

fn destroy_constant_to_id(
    value: &str,
    molecule_id_map: &std::collections::HashMap<String, String>,
) -> String {
    let trimmed = value.trim();
    let name = trimmed
        .strip_prefix("DestroyMolecules.")
        .unwrap_or_else(|| panic!("expected DestroyMolecules reference, got '{trimmed}'"));
    molecule_id_map
        .get(name)
        .map(|id| format!("destroy:{id}"))
        .unwrap_or_else(|| format!("destroy:{}", name.to_ascii_lowercase()))
}

fn parse_molecule_id_map(source: &str) -> std::collections::HashMap<String, String> {
    let mut map = std::collections::HashMap::new();
    for block in split_reaction_blocks(source) {
        let text = block.text;
        let constant_name = reaction_constant_name(&text)
            .unwrap_or_else(|| panic!("missing molecule constant in block:\n{text}"));
        let cleaned = strip_call(&text, ".withResult(");
        let id = extract_single_call_arg(&cleaned, ".id(")
            .unwrap_or_else(|| panic!("missing molecule id in block:\n{cleaned}"));
        let id = unquote_java_string(&id);
        map.insert(constant_name, id);
    }
    map
}

fn reaction_constant_name(block: &str) -> Option<String> {
    for line in block.lines() {
        let trimmed = line.trim();
        if let Some(index) = trimmed.find("= builder()") {
            return Some(trimmed[..index].trim().to_string());
        }
    }
    None
}

fn unquote_java_string(value: &str) -> String {
    let trimmed = value.trim();
    if let Some(stripped) = trimmed.strip_prefix('"').and_then(|inner| inner.strip_suffix('"')) {
        stripped.to_string()
    } else {
        trimmed.to_string()
    }
}