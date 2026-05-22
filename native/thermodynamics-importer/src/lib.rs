use std::collections::{BTreeMap, BTreeSet};
use thermodynamics_datafile::{
    ActivityModelRecord, DataSourceRecord, DatabaseMetadata, ElementCountRecord, ElementRecord,
    HeatCapacityRecord, MolarEnergyRecord, PhaseRecord, SpeciesRecord, TemperatureRangeRecord,
    ThermodynamicsDatabaseFile, FORMAT_MAGIC, FORMAT_VERSION,
};

const GAS_CONSTANT_JOULE_PER_MOL_KELVIN: f64 = 8.314_462_618_153_24;
const REFERENCE_TEMPERATURE_KELVIN: f64 = 298.15;

const SUPPORTED_ELEMENTS: &[(&str, u8, u16)] = &[
    ("H", 1, 1),
    ("C", 6, 6),
    ("O", 8, 8),
    ("N", 7, 7),
    ("Na", 11, 11),
    ("Mg", 12, 12),
    ("S", 16, 16),
    ("Cl", 17, 17),
    ("K", 19, 19),
    ("Ca", 20, 20),
];

const GAS_ALLOWLIST: &[&str] = &["H2O", "CO2", "O2", "N2", "H2", "CO", "NH3", "HCl", "SO2"];

#[derive(Debug, PartialEq)]
pub enum ImportError {
    InvalidSpeciesFormula(String),
    MissingCanteraSpecies,
    UnsupportedCanteraThermo(String),
    InvalidCanteraSubset,
}

pub fn import_database(
    phreeqc_text: &str,
    cantera_yaml: &str,
) -> Result<ThermodynamicsDatabaseFile, ImportError> {
    let mut species = import_phreeqc_species(phreeqc_text)?;
    species.extend(import_cantera_gases(cantera_yaml)?);
    normalize_supported_runtime_species(&mut species);
    species.sort_by(|left, right| left.symbol.cmp(&right.symbol));
    for (index, species_record) in species.iter_mut().enumerate() {
        species_record.id = (index + 1) as u16;
    }

    Ok(ThermodynamicsDatabaseFile {
        magic: *FORMAT_MAGIC,
        version: FORMAT_VERSION,
        metadata: DatabaseMetadata {
            name: "create_thermodynamics_db_v1".to_owned(),
            source_summary: "PHREEQC phreeqc.dat subset + Cantera/NASA gas subset".to_owned(),
            license_note: "Source licenses must be reviewed before distribution.".to_owned(),
        },
        elements: supported_element_records(),
        species,
    })
}

pub fn import_phreeqc_species(text: &str) -> Result<Vec<SpeciesRecord>, ImportError> {
    let mut records = Vec::new();
    let mut block = PhreeqcBlock::Other;
    let mut current_symbol: Option<String> = None;
    let mut current_log_k: Option<f64> = None;
    let mut current_phase_is_gas = false;
    let mut seen = BTreeSet::new();

    for raw_line in text.lines().chain([""].into_iter()) {
        let line_without_comment = raw_line.split('#').next().unwrap_or_default();
        let trimmed = line_without_comment.trim();
        if trimmed.is_empty() {
            flush_phreeqc_record(
                block,
                &mut current_symbol,
                &mut current_log_k,
                &mut seen,
                &mut records,
            )?;
            continue;
        }

        if is_block_header(trimmed) {
            flush_phreeqc_record(
                block,
                &mut current_symbol,
                &mut current_log_k,
                &mut seen,
                &mut records,
            )?;
            block = match trimmed {
                "SOLUTION_SPECIES" => PhreeqcBlock::SolutionSpecies,
                "PHASES" => PhreeqcBlock::Phases,
                _ => PhreeqcBlock::Other,
            };
            continue;
        }

        if block == PhreeqcBlock::Other {
            continue;
        }

        let starts_new_record = !raw_line.starts_with(' ') && !raw_line.starts_with('\t');
        if starts_new_record {
            flush_phreeqc_record(
                block,
                &mut current_symbol,
                &mut current_log_k,
                &mut seen,
                &mut records,
            )?;
            if block == PhreeqcBlock::Phases && !trimmed.contains('=') {
                current_phase_is_gas = trimmed.contains("(g)");
                current_symbol = None;
                current_log_k = None;
                continue;
            }
            current_symbol = phreeqc_record_symbol(block, trimmed).map(str::to_owned);
            current_log_k = identity_reaction_log_k(block, trimmed);
        } else if block == PhreeqcBlock::Phases && current_symbol.is_none() && trimmed.contains('=')
        {
            if !current_phase_is_gas {
                current_symbol = phreeqc_record_symbol(block, trimmed).map(str::to_owned);
            }
        } else if let Some(value) = trimmed
            .strip_prefix("log_k")
            .or_else(|| trimmed.strip_prefix("-log_k"))
        {
            current_log_k = value
                .split_whitespace()
                .next()
                .and_then(|v| v.trim_end_matches(';').parse().ok());
        }
    }

    Ok(records)
}

fn normalize_supported_runtime_species(species: &mut [SpeciesRecord]) {
    for record in species {
        if let Some(override_record) = supported_runtime_override(&record.symbol) {
            record.activity_model = override_record.activity_model;
            record.standard_gibbs_energy = override_record.standard_gibbs_energy;
            record.standard_enthalpy_of_formation = override_record.standard_enthalpy_of_formation;
            record.constant_pressure_heat_capacity =
                override_record.constant_pressure_heat_capacity;
            record.valid_temperature_range = override_record.valid_temperature_range;
        }
    }
}

fn supported_runtime_override(symbol: &str) -> Option<SpeciesRecord> {
    match symbol {
        "H2O" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::UnitActivity,
            -237_130.0,
            Some(-285_830.0),
            Some(75.3),
        )),
        "H+" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::DaviesAqueous,
            0.0,
            Some(0.0),
            Some(0.0),
        )),
        "OH-" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::DaviesAqueous,
            -157_240.0,
            Some(-230_000.0),
            Some(-148.0),
        )),
        "Na+" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::DaviesAqueous,
            -261_900.0,
            Some(-240_100.0),
            Some(46.0),
        )),
        "Cl-" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::DaviesAqueous,
            -131_200.0,
            Some(-167_200.0),
            Some(-136.0),
        )),
        "Ca+2" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::DaviesAqueous,
            -553_600.0,
            Some(-542_800.0),
            Some(-33.0),
        )),
        "CO3-2" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::DaviesAqueous,
            -527_900.0,
            Some(-677_100.0),
            Some(-289.0),
        )),
        "HCO3-" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::DaviesAqueous,
            -586_900.0,
            Some(-692_000.0),
            Some(-35.0),
        )),
        "CO2" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::IdealMolalityAqueous,
            -386_000.0,
            Some(-413_800.0),
            Some(214.0),
        )),
        "CaCO3(s)" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Solid,
            ActivityModelRecord::UnitActivity,
            -1_128_800.0,
            Some(-1_207_100.0),
            Some(82.0),
        )),
        "H2O(g)" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Gas,
            ActivityModelRecord::IdealGas,
            -228_570.0,
            Some(-241_820.0),
            Some(33.6),
        )),
        "CO2(g)" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Gas,
            ActivityModelRecord::IdealGas,
            -394_360.0,
            Some(-393_510.0),
            Some(37.1),
        )),
        "O2" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Aqueous,
            ActivityModelRecord::IdealMolalityAqueous,
            16_400.0,
            Some(-12_100.0),
            Some(145.0),
        )),
        "O2(g)" => Some(runtime_species_override(
            symbol,
            PhaseRecord::Gas,
            ActivityModelRecord::IdealGas,
            0.0,
            Some(0.0),
            Some(29.4),
        )),
        _ => None,
    }
}

fn runtime_species_override(
    symbol: &str,
    phase: PhaseRecord,
    activity_model: ActivityModelRecord,
    standard_gibbs_energy_joule_per_mol: f64,
    standard_enthalpy_joule_per_mol: Option<f64>,
    heat_capacity_joule_per_mol_kelvin: Option<f64>,
) -> SpeciesRecord {
    SpeciesRecord {
        id: 0,
        symbol: symbol.to_owned(),
        composition: Vec::new(),
        charge_number: 0,
        phase,
        activity_model,
        standard_gibbs_energy: MolarEnergyRecord {
            value_joule_per_mol: standard_gibbs_energy_joule_per_mol,
            reference_temperature_kelvin: REFERENCE_TEMPERATURE_KELVIN,
            source: runtime_override_source("standard Gibbs energy"),
        },
        standard_enthalpy_of_formation: standard_enthalpy_joule_per_mol.map(|value| {
            MolarEnergyRecord {
                value_joule_per_mol: value,
                reference_temperature_kelvin: REFERENCE_TEMPERATURE_KELVIN,
                source: runtime_override_source("standard enthalpy of formation"),
            }
        }),
        constant_pressure_heat_capacity: heat_capacity_joule_per_mol_kelvin.map(|value| {
            HeatCapacityRecord {
                value_joule_per_mol_kelvin: value,
                source: runtime_override_source("constant-pressure heat capacity"),
            }
        }),
        valid_temperature_range: TemperatureRangeRecord {
            min_kelvin: 273.15,
            max_kelvin: 373.15,
        },
        tags: Vec::new(),
    }
}

fn runtime_override_source(property: &str) -> DataSourceRecord {
    DataSourceRecord {
        citation: format!(
            "Create Thermodynamics normalized supported-runtime {property} at 298.15 K"
        ),
        note: "Applied during database generation for the currently supported aqueous/gas equilibrium slice."
            .to_owned(),
    }
}

pub fn import_cantera_gases(text: &str) -> Result<Vec<SpeciesRecord>, ImportError> {
    let document = parse_cantera_subset(text)?;
    let mut records = Vec::new();
    for species in document.species {
        if !GAS_ALLOWLIST.contains(&species.name.as_str()) {
            continue;
        }
        if species.thermo.model != "NASA7" {
            return Err(ImportError::UnsupportedCanteraThermo(species.name));
        }
        let coefficients =
            nasa7_coefficients_for_temperature(&species, REFERENCE_TEMPERATURE_KELVIN)
                .ok_or(ImportError::MissingCanteraSpecies)?;
        if coefficients.len() != 7 {
            return Err(ImportError::UnsupportedCanteraThermo(species.name));
        }
        let (raw_min_temperature_kelvin, max_temperature_kelvin) =
            nasa7_temperature_bounds(&species).ok_or(ImportError::MissingCanteraSpecies)?;
        let min_temperature_kelvin = raw_min_temperature_kelvin.min(REFERENCE_TEMPERATURE_KELVIN);
        let cp = nasa7_heat_capacity(coefficients, REFERENCE_TEMPERATURE_KELVIN);
        let enthalpy = nasa7_enthalpy(coefficients, REFERENCE_TEMPERATURE_KELVIN);
        let gibbs = nasa7_gibbs(coefficients, REFERENCE_TEMPERATURE_KELVIN);
        records.push(SpeciesRecord {
            id: 0,
            symbol: format!("{}(g)", species.name),
            composition: composition_from_map(&species.composition),
            charge_number: 0,
            phase: PhaseRecord::Gas,
            activity_model: ActivityModelRecord::IdealGas,
            standard_gibbs_energy: MolarEnergyRecord {
                value_joule_per_mol: gibbs,
                reference_temperature_kelvin: REFERENCE_TEMPERATURE_KELVIN,
                source: cantera_source(),
            },
            standard_enthalpy_of_formation: Some(MolarEnergyRecord {
                value_joule_per_mol: enthalpy,
                reference_temperature_kelvin: REFERENCE_TEMPERATURE_KELVIN,
                source: cantera_source(),
            }),
            constant_pressure_heat_capacity: Some(HeatCapacityRecord {
                value_joule_per_mol_kelvin: cp,
                source: cantera_source(),
            }),
            valid_temperature_range: TemperatureRangeRecord {
                min_kelvin: min_temperature_kelvin,
                max_kelvin: max_temperature_kelvin,
            },
            tags: vec!["gas".to_owned()],
        });
    }
    Ok(records)
}

fn nasa7_coefficients_for_temperature<'a>(
    species: &'a CanteraSpecies,
    temperature_kelvin: f64,
) -> Option<&'a [f64]> {
    if species.data.is_empty() || species.thermo.temperature_ranges.len() < 2 {
        return None;
    }

    for (index, coefficients) in species.data.iter().enumerate() {
        let lower = species.thermo.temperature_ranges.get(index)?;
        let upper = species.thermo.temperature_ranges.get(index + 1)?;
        if temperature_kelvin >= *lower && temperature_kelvin <= *upper {
            return Some(coefficients.as_slice());
        }
    }

    species.data.first().map(Vec::as_slice)
}

fn nasa7_temperature_bounds(species: &CanteraSpecies) -> Option<(f64, f64)> {
    let mut values = species.thermo.temperature_ranges.iter().copied();
    let first = values.next()?;
    let mut min_value = first;
    let mut max_value = first;
    for value in values {
        if value < min_value {
            min_value = value;
        }
        if value > max_value {
            max_value = value;
        }
    }
    Some((min_value, max_value))
}

fn flush_phreeqc_record(
    block: PhreeqcBlock,
    current_symbol: &mut Option<String>,
    current_log_k: &mut Option<f64>,
    seen: &mut BTreeSet<String>,
    records: &mut Vec<SpeciesRecord>,
) -> Result<(), ImportError> {
    let Some(symbol) = current_symbol.take() else {
        return Ok(());
    };
    let Some(log_k) = current_log_k.take() else {
        return Ok(());
    };
    let seen_key = match block {
        PhreeqcBlock::SolutionSpecies => format!("aqueous:{symbol}"),
        PhreeqcBlock::Phases => format!("solid:{symbol}"),
        PhreeqcBlock::Other => symbol.clone(),
    };
    if !seen.insert(seen_key) {
        return Ok(());
    }
    let Ok(formula) = parse_species_formula(&symbol) else {
        return Ok(());
    };
    if block == PhreeqcBlock::Phases && symbol.contains('(') {
        return Ok(());
    }
    if block == PhreeqcBlock::Phases && formula.charge_number != 0 {
        return Ok(());
    }
    if !formula
        .composition
        .iter()
        .all(|entry| supported_element_ids().contains(&entry.element_id))
    {
        return Ok(());
    }
    let phase = match block {
        PhreeqcBlock::SolutionSpecies => PhaseRecord::Aqueous,
        PhreeqcBlock::Phases => PhaseRecord::Solid,
        PhreeqcBlock::Other => return Ok(()),
    };
    let tags = tags_for_record(phase, &formula);
    records.push(SpeciesRecord {
        id: 0,
        symbol: if phase == PhaseRecord::Solid {
            format!("{symbol}(s)")
        } else {
            symbol
        },
        composition: formula.composition,
        charge_number: formula.charge_number,
        phase,
        activity_model: match phase {
            PhaseRecord::Aqueous if formula.charge_number != 0 => {
                ActivityModelRecord::DaviesAqueous
            }
            PhaseRecord::Aqueous => ActivityModelRecord::IdealMolalityAqueous,
            PhaseRecord::Solid => ActivityModelRecord::UnitActivity,
            PhaseRecord::Gas => ActivityModelRecord::IdealGas,
        },
        standard_gibbs_energy: MolarEnergyRecord {
            value_joule_per_mol: -GAS_CONSTANT_JOULE_PER_MOL_KELVIN
                * REFERENCE_TEMPERATURE_KELVIN
                * log_k
                * std::f64::consts::LN_10,
            reference_temperature_kelvin: REFERENCE_TEMPERATURE_KELVIN,
            source: phreeqc_source(),
        },
        standard_enthalpy_of_formation: None,
        constant_pressure_heat_capacity: None,
        valid_temperature_range: TemperatureRangeRecord {
            min_kelvin: 273.15,
            max_kelvin: 373.15,
        },
        tags,
    });
    Ok(())
}

fn parse_species_formula(symbol: &str) -> Result<ParsedFormula, ImportError> {
    let normalized = symbol
        .trim()
        .trim_end_matches("(aq)")
        .trim_end_matches("(s)")
        .trim_end_matches("(g)");
    let (formula_part, charge_number) = split_charge(normalized);
    let mut chars = formula_part.chars().peekable();
    let mut counts = BTreeMap::<u16, u16>::new();
    while let Some(ch) = chars.next() {
        if !ch.is_ascii_uppercase() {
            return Err(ImportError::InvalidSpeciesFormula(symbol.to_owned()));
        }
        let mut element_symbol = ch.to_string();
        while matches!(chars.peek(), Some(next) if next.is_ascii_lowercase()) {
            element_symbol.push(chars.next().unwrap());
        }
        let mut count_text = String::new();
        while matches!(chars.peek(), Some(next) if next.is_ascii_digit()) {
            count_text.push(chars.next().unwrap());
        }
        let count = if count_text.is_empty() {
            1
        } else {
            count_text
                .parse::<u16>()
                .map_err(|_| ImportError::InvalidSpeciesFormula(symbol.to_owned()))?
        };
        let element_id = supported_element_id(&element_symbol)
            .ok_or_else(|| ImportError::InvalidSpeciesFormula(symbol.to_owned()))?;
        *counts.entry(element_id).or_default() += count;
    }

    Ok(ParsedFormula {
        composition: counts
            .into_iter()
            .map(|(element_id, count)| ElementCountRecord { element_id, count })
            .collect(),
        charge_number,
    })
}

fn phreeqc_record_symbol(block: PhreeqcBlock, line: &str) -> Option<&str> {
    if let Some((left, right)) = line.split_once('=') {
        let candidate = match block {
            PhreeqcBlock::SolutionSpecies => solution_species_candidate(right)?,
            PhreeqcBlock::Phases => formula_token(left)?,
            PhreeqcBlock::Other => return None,
        };
        return Some(candidate);
    }
    Some(line.trim())
}

fn solution_species_candidate(right_side: &str) -> Option<&str> {
    let mut fallback = None;
    for token in right_side.split(" + ").map(formula_token).flatten() {
        fallback = Some(token);
        if token != "H+" && token != "H2O" {
            return Some(token);
        }
    }
    fallback
}

fn formula_token(value: &str) -> Option<&str> {
    value
        .split_whitespace()
        .find(|part| part.chars().any(|ch| ch.is_ascii_alphabetic()))
}

fn identity_reaction_log_k(block: PhreeqcBlock, line: &str) -> Option<f64> {
    if block != PhreeqcBlock::SolutionSpecies {
        return None;
    }
    let (left, right) = line.split_once('=')?;
    if left.trim() == right.trim() {
        Some(0.0)
    } else {
        None
    }
}

fn split_charge(symbol: &str) -> (&str, i8) {
    if let Some(stripped) = symbol.strip_suffix("--") {
        return (stripped, -2);
    }
    if let Some(stripped) = symbol.strip_suffix("++") {
        return (stripped, 2);
    }
    if let Some(index) = symbol.rfind(['+', '-']) {
        let sign = symbol.as_bytes()[index] as char;
        let magnitude = symbol[index + 1..].parse::<i8>().unwrap_or(1);
        return (
            &symbol[..index],
            if sign == '-' { -magnitude } else { magnitude },
        );
    }
    (symbol, 0)
}

fn tags_for_record(phase: PhaseRecord, formula: &ParsedFormula) -> Vec<String> {
    let mut tags = Vec::new();
    match phase {
        PhaseRecord::Aqueous => tags.push("aqueous".to_owned()),
        PhaseRecord::Solid => tags.push("solid".to_owned()),
        PhaseRecord::Gas => tags.push("gas".to_owned()),
    }
    if formula.charge_number != 0 {
        tags.push("acid_base".to_owned());
        tags.push("salt".to_owned());
    }
    if has_element(formula, "C") && has_element(formula, "O") {
        tags.push("carbonate".to_owned());
    }
    if has_element(formula, "S") && has_element(formula, "O") {
        tags.push("sulfate".to_owned());
    }
    if has_element(formula, "N") && has_element(formula, "O") {
        tags.push("nitrate".to_owned());
    }
    tags
}

fn has_element(formula: &ParsedFormula, symbol: &str) -> bool {
    let Some(element_id) = supported_element_id(symbol) else {
        return false;
    };
    formula
        .composition
        .iter()
        .any(|entry| entry.element_id == element_id)
}

fn supported_element_records() -> Vec<ElementRecord> {
    SUPPORTED_ELEMENTS
        .iter()
        .map(|(symbol, atomic_number, id)| ElementRecord {
            id: *id,
            atomic_number: *atomic_number,
            symbol: (*symbol).to_owned(),
        })
        .collect()
}

fn supported_element_ids() -> BTreeSet<u16> {
    SUPPORTED_ELEMENTS.iter().map(|(_, _, id)| *id).collect()
}

fn supported_element_id(symbol: &str) -> Option<u16> {
    SUPPORTED_ELEMENTS
        .iter()
        .find(|(candidate, _, _)| *candidate == symbol)
        .map(|(_, _, id)| *id)
}

fn composition_from_map(composition: &BTreeMap<String, f64>) -> Vec<ElementCountRecord> {
    composition
        .iter()
        .filter_map(|(symbol, count)| {
            supported_element_id(symbol).map(|element_id| ElementCountRecord {
                element_id,
                count: *count as u16,
            })
        })
        .collect()
}

fn parse_cantera_subset(text: &str) -> Result<CanteraDocument, ImportError> {
    let mut species_records = Vec::new();
    let mut current: Option<CanteraSpecies> = None;
    let mut in_data = false;
    let mut in_species_section = false;
    let mut in_thermo = false;
    let mut pending_data_row: Option<String> = None;

    for raw_line in text.lines() {
        let trimmed = raw_line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed == "species:" {
            in_species_section = true;
            continue;
        }
        if trimmed == "thermo:" {
            in_thermo = true;
            continue;
        }
        if trimmed == "transport:" {
            in_thermo = false;
            in_data = false;
            continue;
        }

        if let Some(row) = pending_data_row.as_mut() {
            row.push(' ');
            row.push_str(trimmed);
            if trimmed.contains(']') {
                let row = pending_data_row.take().unwrap();
                current
                    .as_mut()
                    .ok_or(ImportError::InvalidCanteraSubset)?
                    .data
                    .push(parse_number_list(&row)?);
            }
            continue;
        }

        if in_species_section {
            if let Some(name) = trimmed.strip_prefix("- name:") {
                if let Some(species) = current.take() {
                    species_records.push(species);
                }
                current = Some(CanteraSpecies {
                    name: name.trim().to_owned(),
                    composition: BTreeMap::new(),
                    thermo: CanteraThermo {
                        model: String::new(),
                        temperature_ranges: Vec::new(),
                    },
                    data: Vec::new(),
                });
                in_data = false;
                in_thermo = false;
                continue;
            }
        }

        let Some(species) = current.as_mut() else {
            continue;
        };
        if let Some(composition) = trimmed.strip_prefix("composition:") {
            species.composition = parse_inline_composition(composition.trim())?;
            continue;
        }
        if in_thermo {
            if let Some(model) = trimmed.strip_prefix("model:") {
                species.thermo.model = model.trim().to_owned();
                continue;
            }
            if let Some(ranges) = trimmed.strip_prefix("temperature-ranges:") {
                species.thermo.temperature_ranges = parse_number_list(ranges.trim())?;
                continue;
            }
            if trimmed == "data:" {
                in_data = true;
                continue;
            }
            if in_data && trimmed.starts_with("- [") {
                let row = trimmed.trim_start_matches("- ").trim();
                if row.contains(']') {
                    species.data.push(parse_number_list(row)?);
                } else {
                    pending_data_row = Some(row.to_owned());
                }
            }
        }
    }

    if let Some(species) = current.take() {
        species_records.push(species);
    }
    Ok(CanteraDocument {
        species: species_records,
    })
}

fn parse_inline_composition(value: &str) -> Result<BTreeMap<String, f64>, ImportError> {
    let inner = value
        .strip_prefix('{')
        .and_then(|value| value.strip_suffix('}'))
        .ok_or(ImportError::InvalidCanteraSubset)?;
    let mut composition = BTreeMap::new();
    for part in inner.split(',') {
        let (symbol, count) = part
            .split_once(':')
            .ok_or(ImportError::InvalidCanteraSubset)?;
        composition.insert(
            symbol.trim().to_owned(),
            count
                .trim()
                .parse::<f64>()
                .map_err(|_| ImportError::InvalidCanteraSubset)?,
        );
    }
    Ok(composition)
}

fn parse_number_list(value: &str) -> Result<Vec<f64>, ImportError> {
    let inner = value
        .strip_prefix('[')
        .and_then(|value| value.strip_suffix(']'))
        .ok_or(ImportError::InvalidCanteraSubset)?;
    inner
        .split(',')
        .map(|part| {
            part.trim()
                .parse::<f64>()
                .map_err(|_| ImportError::InvalidCanteraSubset)
        })
        .collect()
}

fn nasa7_heat_capacity(coefficients: &[f64], temperature_kelvin: f64) -> f64 {
    let t = temperature_kelvin;
    let cp_over_r = coefficients[0]
        + coefficients[1] * t
        + coefficients[2] * t.powi(2)
        + coefficients[3] * t.powi(3)
        + coefficients[4] * t.powi(4);
    cp_over_r * GAS_CONSTANT_JOULE_PER_MOL_KELVIN
}

fn nasa7_enthalpy(coefficients: &[f64], temperature_kelvin: f64) -> f64 {
    let t = temperature_kelvin;
    let h_over_rt = coefficients[0]
        + coefficients[1] * t / 2.0
        + coefficients[2] * t.powi(2) / 3.0
        + coefficients[3] * t.powi(3) / 4.0
        + coefficients[4] * t.powi(4) / 5.0
        + coefficients[5] / t;
    h_over_rt * GAS_CONSTANT_JOULE_PER_MOL_KELVIN * t
}

fn nasa7_gibbs(coefficients: &[f64], temperature_kelvin: f64) -> f64 {
    let t = temperature_kelvin;
    let entropy_over_r = coefficients[0] * t.ln()
        + coefficients[1] * t
        + coefficients[2] * t.powi(2) / 2.0
        + coefficients[3] * t.powi(3) / 3.0
        + coefficients[4] * t.powi(4) / 4.0
        + coefficients[6];
    nasa7_enthalpy(coefficients, t) - t * entropy_over_r * GAS_CONSTANT_JOULE_PER_MOL_KELVIN
}

fn phreeqc_source() -> DataSourceRecord {
    DataSourceRecord {
        citation: "PHREEQC phreeqc.dat imported equilibrium constant".to_owned(),
        note: "Standard Gibbs energy derived from log_k at 298.15 K".to_owned(),
    }
}

fn cantera_source() -> DataSourceRecord {
    DataSourceRecord {
        citation: "Cantera/NASA gas thermodynamic coefficients".to_owned(),
        note: "Values evaluated from NASA7 coefficients at 298.15 K".to_owned(),
    }
}

fn is_block_header(line: &str) -> bool {
    matches!(
        line,
        "SOLUTION_MASTER_SPECIES" | "SOLUTION_SPECIES" | "PHASES" | "END"
    )
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PhreeqcBlock {
    SolutionSpecies,
    Phases,
    Other,
}

struct ParsedFormula {
    composition: Vec<ElementCountRecord>,
    charge_number: i8,
}

struct CanteraDocument {
    species: Vec<CanteraSpecies>,
}

struct CanteraSpecies {
    name: String,
    composition: BTreeMap<String, f64>,
    thermo: CanteraThermo,
    data: Vec<Vec<f64>>,
}

struct CanteraThermo {
    model: String,
    temperature_ranges: Vec<f64>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phreeqc_fixture_imports_log_k_as_standard_gibbs_energy() {
        let records = import_phreeqc_species(
            r#"
SOLUTION_SPECIES
H+
    log_k 0.0
OH-
    log_k -14.0
"#,
        )
        .unwrap();

        let hydroxide = records
            .iter()
            .find(|species| species.symbol == "OH-")
            .unwrap();
        let expected = -GAS_CONSTANT_JOULE_PER_MOL_KELVIN
            * REFERENCE_TEMPERATURE_KELVIN
            * -14.0
            * std::f64::consts::LN_10;

        assert_eq!(hydroxide.charge_number, -1);
        assert!((hydroxide.standard_gibbs_energy.value_joule_per_mol - expected).abs() < 1.0e-9);
        assert!(hydroxide.standard_enthalpy_of_formation.is_none());
        assert!(hydroxide.constant_pressure_heat_capacity.is_none());
    }

    #[test]
    fn cantera_fixture_imports_oxygen_gas_thermo() {
        let records = import_cantera_gases(
            r#"
species:
- name: O2
  composition: {O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [3.78245636, -0.00299673416, 0.00000984730201, -0.00000000968129509, 0.00000000000324372837, -1063.94356, 3.65767573]
    - [3.28253784, 0.00148308754, -0.000000757966669, 0.000000000209470556, -0.000000000000021676594, -1088.45772, 5.45323129]
"#,
        )
        .unwrap();

        let oxygen = records
            .iter()
            .find(|species| species.symbol == "O2(g)")
            .unwrap();

        assert_eq!(oxygen.phase, PhaseRecord::Gas);
        assert_eq!(oxygen.activity_model, ActivityModelRecord::IdealGas);
        assert_eq!(oxygen.composition[0].element_id, 8);
        assert_eq!(oxygen.composition[0].count, 2);
        assert!(
            oxygen
                .constant_pressure_heat_capacity
                .as_ref()
                .unwrap()
                .value_joule_per_mol_kelvin
                > 20.0
        );
    }

    #[test]
    fn database_import_combines_sources_and_assigns_stable_ids() {
        let database = import_database(
            r#"
SOLUTION_SPECIES
H+
    log_k 0.0
OH-
    log_k -14.0
"#,
            r#"
species:
- name: O2
  composition: {O: 2}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 3500.0]
    data:
    - [3.78245636, -0.00299673416, 0.00000984730201, -0.00000000968129509, 0.00000000000324372837, -1063.94356, 3.65767573]
    - [3.28253784, 0.00148308754, -0.000000757966669, 0.000000000209470556, -0.000000000000021676594, -1088.45772, 5.45323129]
"#,
        )
        .unwrap();

        assert_eq!(database.magic, *FORMAT_MAGIC);
        assert_eq!(database.version, FORMAT_VERSION);
        assert_eq!(database.species.len(), 3);
        assert_eq!(database.species[0].id, 1);
        assert_eq!(database.species[1].id, 2);
        assert_eq!(database.species[2].id, 3);
    }
}
