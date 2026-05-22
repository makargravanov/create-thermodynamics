use std::env;
use std::fs;
use std::path::PathBuf;
use thermodynamics_datafile::{encode_database_file, ThermodynamicsDatabaseFile};
use thermodynamics_importer::import_database;

fn main() {
    if let Err(error) = run() {
        eprintln!("{error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.len() != 4 {
        return Err(
            "usage: thermodynamics-importer <phreeqc.dat> <cantera.yaml> <out.json> <out.ctdb.gz>"
                .to_owned(),
        );
    }

    let phreeqc_path = PathBuf::from(&args[0]);
    let cantera_path = PathBuf::from(&args[1]);
    let json_path = PathBuf::from(&args[2]);
    let packed_path = PathBuf::from(&args[3]);

    let phreeqc_bytes = fs::read(&phreeqc_path)
        .map_err(|error| format!("failed to read {}: {error}", phreeqc_path.display()))?;
    let cantera_bytes = fs::read(&cantera_path)
        .map_err(|error| format!("failed to read {}: {error}", cantera_path.display()))?;
    let phreeqc_text = String::from_utf8_lossy(&phreeqc_bytes);
    let cantera_yaml = String::from_utf8_lossy(&cantera_bytes);
    let database = import_database(&phreeqc_text, &cantera_yaml)
        .map_err(|error| format!("failed to import database: {error:?}"))?;

    let json = database_to_review_json(&database);
    fs::write(&json_path, json)
        .map_err(|error| format!("failed to write {}: {error}", json_path.display()))?;

    let packed = encode_database_file(&database)
        .map_err(|error| format!("failed to encode packed database: {error:?}"))?;
    fs::write(&packed_path, &packed)
        .map_err(|error| format!("failed to write {}: {error}", packed_path.display()))?;

    println!(
        "wrote {} species, packed size {} bytes",
        database.species.len(),
        packed.len()
    );
    Ok(())
}

fn database_to_review_json(database: &ThermodynamicsDatabaseFile) -> String {
    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!(
        "  \"name\": \"{}\",\n",
        escape_json(&database.metadata.name)
    ));
    json.push_str(&format!("  \"version\": {},\n", database.version));
    json.push_str("  \"elements\": [\n");
    for (index, element) in database.elements.iter().enumerate() {
        let comma = if index + 1 == database.elements.len() {
            ""
        } else {
            ","
        };
        json.push_str(&format!(
            "    {{\"id\": {}, \"atomic_number\": {}, \"symbol\": \"{}\"}}{}\n",
            element.id,
            element.atomic_number,
            escape_json(&element.symbol),
            comma
        ));
    }
    json.push_str("  ],\n");
    json.push_str("  \"species\": [\n");
    for (index, species) in database.species.iter().enumerate() {
        let comma = if index + 1 == database.species.len() {
            ""
        } else {
            ","
        };
        json.push_str(&format!(
            "    {{\"id\": {}, \"symbol\": \"{}\", \"charge_number\": {}, \"tags\": [{}]}}{}\n",
            species.id,
            escape_json(&species.symbol),
            species.charge_number,
            species
                .tags
                .iter()
                .map(|tag| format!("\"{}\"", escape_json(tag)))
                .collect::<Vec<_>>()
                .join(", "),
            comma
        ));
    }
    json.push_str("  ]\n");
    json.push_str("}\n");
    json
}

fn escape_json(value: &str) -> String {
    value.replace('\\', "\\\\").replace('"', "\\\"")
}
