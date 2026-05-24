use std::collections::BTreeMap;
use thermodynamics_core::{
    ActivityModel, ConstantPressureHeatCapacity, DataSource, Element, ElementId, PhaseKind,
    Species, SpeciesId, SpeciesRegistry, SpeciesRegistryError, StandardEnthalpyOfFormation,
    StandardGibbsEnergy, StandardThermo, TemperatureRange,
};

pub const FORMAT_MAGIC: &[u8; 4] = b"CTDB";
pub const FORMAT_VERSION: u16 = 1;

#[derive(Debug, Clone, PartialEq)]
pub struct ThermodynamicsDatabaseFile {
    pub magic: [u8; 4],
    pub version: u16,
    pub metadata: DatabaseMetadata,
    pub elements: Vec<ElementRecord>,
    pub species: Vec<SpeciesRecord>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DatabaseMetadata {
    pub name: String,
    pub source_summary: String,
    pub license_note: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ElementRecord {
    pub id: u16,
    pub atomic_number: u8,
    pub symbol: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesRecord {
    pub id: u16,
    pub symbol: String,
    pub composition: Vec<ElementCountRecord>,
    pub charge_number: i8,
    pub phase: PhaseRecord,
    pub activity_model: ActivityModelRecord,
    pub standard_gibbs_energy: MolarEnergyRecord,
    pub standard_enthalpy_of_formation: Option<MolarEnergyRecord>,
    pub constant_pressure_heat_capacity: Option<HeatCapacityRecord>,
    pub valid_temperature_range: TemperatureRangeRecord,
    pub tags: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ElementCountRecord {
    pub element_id: u16,
    pub count: u16,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MolarEnergyRecord {
    pub value_joule_per_mol: f64,
    pub reference_temperature_kelvin: f64,
    pub source: DataSourceRecord,
}

#[derive(Debug, Clone, PartialEq)]
pub struct HeatCapacityRecord {
    pub value_joule_per_mol_kelvin: f64,
    pub source: DataSourceRecord,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TemperatureRangeRecord {
    pub min_kelvin: f64,
    pub max_kelvin: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataSourceRecord {
    pub citation: String,
    pub note: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseRecord {
    Aqueous,
    Solid,
    Gas,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActivityModelRecord {
    DaviesAqueous,
    IdealMolalityAqueous,
    IdealGas,
    UnitActivity,
}

#[derive(Debug, PartialEq)]
pub enum DatafileError {
    InvalidMagic,
    UnsupportedVersion(u16),
    Decode,
    Registry(SpeciesRegistryError),
}

pub fn encode_database_file(
    database: &ThermodynamicsDatabaseFile,
) -> Result<Vec<u8>, DatafileError> {
    let mut writer = BinaryWriter::default();
    writer.write_database(database)?;
    Ok(writer.bytes)
}

pub fn decode_database_file(bytes: &[u8]) -> Result<ThermodynamicsDatabaseFile, DatafileError> {
    let mut reader = BinaryReader { bytes, offset: 0 };
    let database = reader.read_database()?;
    if reader.offset != bytes.len() {
        return Err(DatafileError::Decode);
    }
    validate_header(&database)?;
    Ok(database)
}

pub fn registry_from_database_file(bytes: &[u8]) -> Result<SpeciesRegistry, DatafileError> {
    let database = decode_database_file(bytes)?;
    database_to_registry(database)
}

pub fn database_to_registry(
    database: ThermodynamicsDatabaseFile,
) -> Result<SpeciesRegistry, DatafileError> {
    validate_header(&database)?;
    let elements = database
        .elements
        .into_iter()
        .map(|element| Element {
            id: ElementId(element.id),
            atomic_number: element.atomic_number,
            symbol: element.symbol,
        })
        .collect();
    let species = database
        .species
        .into_iter()
        .map(species_record_to_species)
        .collect();
    SpeciesRegistry::new(elements, species).map_err(DatafileError::Registry)
}

pub fn registry_to_database(
    registry: &SpeciesRegistry,
    metadata: DatabaseMetadata,
) -> ThermodynamicsDatabaseFile {
    ThermodynamicsDatabaseFile {
        magic: *FORMAT_MAGIC,
        version: FORMAT_VERSION,
        metadata,
        elements: registry
            .elements()
            .map(|element| ElementRecord {
                id: element.id.0,
                atomic_number: element.atomic_number,
                symbol: element.symbol.clone(),
            })
            .collect(),
        species: registry
            .species_records()
            .map(species_to_species_record)
            .collect(),
    }
}

fn validate_header(database: &ThermodynamicsDatabaseFile) -> Result<(), DatafileError> {
    if database.magic != *FORMAT_MAGIC {
        return Err(DatafileError::InvalidMagic);
    }
    if database.version != FORMAT_VERSION {
        return Err(DatafileError::UnsupportedVersion(database.version));
    }
    Ok(())
}

fn species_record_to_species(record: SpeciesRecord) -> Species {
    Species {
        id: SpeciesId(record.id),
        symbol: record.symbol,
        composition: record
            .composition
            .into_iter()
            .map(|entry| (ElementId(entry.element_id), entry.count))
            .collect::<BTreeMap<_, _>>(),
        charge_number: record.charge_number,
        phase: phase_from_record(record.phase),
        activity_model: activity_model_from_record(record.activity_model),
        thermo: StandardThermo {
            standard_gibbs_energy: molar_energy_from_record(record.standard_gibbs_energy),
            standard_enthalpy_of_formation: record
                .standard_enthalpy_of_formation
                .map(molar_enthalpy_from_record),
            constant_pressure_heat_capacity: record
                .constant_pressure_heat_capacity
                .map(heat_capacity_from_record),
            valid_temperature_range: TemperatureRange {
                min_kelvin: record.valid_temperature_range.min_kelvin,
                max_kelvin: record.valid_temperature_range.max_kelvin,
            },
        },
    }
}

fn species_to_species_record(species: &Species) -> SpeciesRecord {
    SpeciesRecord {
        id: species.id.0,
        symbol: species.symbol.clone(),
        composition: species
            .composition
            .iter()
            .map(|(element_id, count)| ElementCountRecord {
                element_id: element_id.0,
                count: *count,
            })
            .collect(),
        charge_number: species.charge_number,
        phase: phase_to_record(species.phase),
        activity_model: activity_model_to_record(species.activity_model),
        standard_gibbs_energy: molar_energy_to_record(&species.thermo.standard_gibbs_energy),
        standard_enthalpy_of_formation: species
            .thermo
            .standard_enthalpy_of_formation
            .as_ref()
            .map(molar_enthalpy_to_record),
        constant_pressure_heat_capacity: species
            .thermo
            .constant_pressure_heat_capacity
            .as_ref()
            .map(heat_capacity_to_record),
        valid_temperature_range: TemperatureRangeRecord {
            min_kelvin: species.thermo.valid_temperature_range.min_kelvin,
            max_kelvin: species.thermo.valid_temperature_range.max_kelvin,
        },
        tags: Vec::new(),
    }
}

fn molar_energy_from_record(record: MolarEnergyRecord) -> StandardGibbsEnergy {
    StandardGibbsEnergy {
        value_joule_per_mol: record.value_joule_per_mol,
        reference_temperature_kelvin: record.reference_temperature_kelvin,
        source: source_from_record(record.source),
    }
}

fn molar_enthalpy_from_record(record: MolarEnergyRecord) -> StandardEnthalpyOfFormation {
    StandardEnthalpyOfFormation {
        value_joule_per_mol: record.value_joule_per_mol,
        reference_temperature_kelvin: record.reference_temperature_kelvin,
        source: source_from_record(record.source),
    }
}

fn heat_capacity_from_record(record: HeatCapacityRecord) -> ConstantPressureHeatCapacity {
    ConstantPressureHeatCapacity {
        value_joule_per_mol_kelvin: record.value_joule_per_mol_kelvin,
        source: source_from_record(record.source),
    }
}

fn molar_energy_to_record(energy: &StandardGibbsEnergy) -> MolarEnergyRecord {
    MolarEnergyRecord {
        value_joule_per_mol: energy.value_joule_per_mol,
        reference_temperature_kelvin: energy.reference_temperature_kelvin,
        source: source_to_record(&energy.source),
    }
}

fn molar_enthalpy_to_record(enthalpy: &StandardEnthalpyOfFormation) -> MolarEnergyRecord {
    MolarEnergyRecord {
        value_joule_per_mol: enthalpy.value_joule_per_mol,
        reference_temperature_kelvin: enthalpy.reference_temperature_kelvin,
        source: source_to_record(&enthalpy.source),
    }
}

fn heat_capacity_to_record(heat_capacity: &ConstantPressureHeatCapacity) -> HeatCapacityRecord {
    HeatCapacityRecord {
        value_joule_per_mol_kelvin: heat_capacity.value_joule_per_mol_kelvin,
        source: source_to_record(&heat_capacity.source),
    }
}

fn source_from_record(record: DataSourceRecord) -> DataSource {
    DataSource {
        citation: record.citation,
        note: record.note,
    }
}

fn source_to_record(source: &DataSource) -> DataSourceRecord {
    DataSourceRecord {
        citation: source.citation.clone(),
        note: source.note.clone(),
    }
}

fn phase_from_record(record: PhaseRecord) -> PhaseKind {
    match record {
        PhaseRecord::Aqueous => PhaseKind::Aqueous,
        PhaseRecord::Solid => PhaseKind::Solid,
        PhaseRecord::Gas => PhaseKind::Gas,
    }
}

fn phase_to_record(phase: PhaseKind) -> PhaseRecord {
    match phase {
        PhaseKind::Aqueous => PhaseRecord::Aqueous,
        PhaseKind::Solid => PhaseRecord::Solid,
        PhaseKind::Gas => PhaseRecord::Gas,
    }
}

fn activity_model_from_record(record: ActivityModelRecord) -> ActivityModel {
    match record {
        ActivityModelRecord::DaviesAqueous => ActivityModel::DaviesAqueous,
        ActivityModelRecord::IdealMolalityAqueous => ActivityModel::IdealMolalityAqueous,
        ActivityModelRecord::IdealGas => ActivityModel::IdealGas,
        ActivityModelRecord::UnitActivity => ActivityModel::UnitActivity,
    }
}

fn activity_model_to_record(activity_model: ActivityModel) -> ActivityModelRecord {
    match activity_model {
        ActivityModel::DaviesAqueous => ActivityModelRecord::DaviesAqueous,
        ActivityModel::IdealMolalityAqueous => ActivityModelRecord::IdealMolalityAqueous,
        ActivityModel::IdealGas => ActivityModelRecord::IdealGas,
        ActivityModel::UnitActivity => ActivityModelRecord::UnitActivity,
    }
}

#[derive(Default)]
struct BinaryWriter {
    bytes: Vec<u8>,
}

impl BinaryWriter {
    fn write_database(
        &mut self,
        database: &ThermodynamicsDatabaseFile,
    ) -> Result<(), DatafileError> {
        self.bytes.extend_from_slice(&database.magic);
        self.write_u16(database.version);
        self.write_string(&database.metadata.name)?;
        self.write_string(&database.metadata.source_summary)?;
        self.write_string(&database.metadata.license_note)?;
        self.write_len(database.elements.len())?;
        for element in &database.elements {
            self.write_u16(element.id);
            self.bytes.push(element.atomic_number);
            self.write_string(&element.symbol)?;
        }
        self.write_len(database.species.len())?;
        for species in &database.species {
            self.write_species(species)?;
        }
        Ok(())
    }

    fn write_species(&mut self, species: &SpeciesRecord) -> Result<(), DatafileError> {
        self.write_u16(species.id);
        self.write_string(&species.symbol)?;
        self.write_len(species.composition.len())?;
        for entry in &species.composition {
            self.write_u16(entry.element_id);
            self.write_u16(entry.count);
        }
        self.bytes.push(species.charge_number as u8);
        self.bytes.push(match species.phase {
            PhaseRecord::Aqueous => 0,
            PhaseRecord::Solid => 1,
            PhaseRecord::Gas => 2,
        });
        self.bytes.push(match species.activity_model {
            ActivityModelRecord::DaviesAqueous => 0,
            ActivityModelRecord::IdealMolalityAqueous => 1,
            ActivityModelRecord::IdealGas => 2,
            ActivityModelRecord::UnitActivity => 3,
        });
        self.write_molar_energy(&species.standard_gibbs_energy)?;
        self.write_optional_molar_energy(&species.standard_enthalpy_of_formation)?;
        self.write_optional_heat_capacity(&species.constant_pressure_heat_capacity)?;
        self.write_f64(species.valid_temperature_range.min_kelvin);
        self.write_f64(species.valid_temperature_range.max_kelvin);
        self.write_len(species.tags.len())?;
        for tag in &species.tags {
            self.write_string(tag)?;
        }
        Ok(())
    }

    fn write_optional_molar_energy(
        &mut self,
        value: &Option<MolarEnergyRecord>,
    ) -> Result<(), DatafileError> {
        match value {
            Some(record) => {
                self.bytes.push(1);
                self.write_molar_energy(record)
            }
            None => {
                self.bytes.push(0);
                Ok(())
            }
        }
    }

    fn write_optional_heat_capacity(
        &mut self,
        value: &Option<HeatCapacityRecord>,
    ) -> Result<(), DatafileError> {
        match value {
            Some(record) => {
                self.bytes.push(1);
                self.write_heat_capacity(record)
            }
            None => {
                self.bytes.push(0);
                Ok(())
            }
        }
    }

    fn write_molar_energy(&mut self, value: &MolarEnergyRecord) -> Result<(), DatafileError> {
        self.write_f64(value.value_joule_per_mol);
        self.write_f64(value.reference_temperature_kelvin);
        self.write_source(&value.source)
    }

    fn write_heat_capacity(&mut self, value: &HeatCapacityRecord) -> Result<(), DatafileError> {
        self.write_f64(value.value_joule_per_mol_kelvin);
        self.write_source(&value.source)
    }

    fn write_source(&mut self, source: &DataSourceRecord) -> Result<(), DatafileError> {
        self.write_string(&source.citation)?;
        self.write_string(&source.note)
    }

    fn write_string(&mut self, value: &str) -> Result<(), DatafileError> {
        self.write_len(value.len())?;
        self.bytes.extend_from_slice(value.as_bytes());
        Ok(())
    }

    fn write_len(&mut self, value: usize) -> Result<(), DatafileError> {
        let value = u32::try_from(value).map_err(|_| DatafileError::Decode)?;
        self.bytes.extend_from_slice(&value.to_le_bytes());
        Ok(())
    }

    fn write_u16(&mut self, value: u16) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn write_f64(&mut self, value: f64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }
}

struct BinaryReader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> BinaryReader<'a> {
    fn read_database(&mut self) -> Result<ThermodynamicsDatabaseFile, DatafileError> {
        let magic = self.read_magic()?;
        let version = self.read_u16()?;
        let metadata = DatabaseMetadata {
            name: self.read_string()?,
            source_summary: self.read_string()?,
            license_note: self.read_string()?,
        };
        let elements = (0..self.read_len()?)
            .map(|_| {
                Ok(ElementRecord {
                    id: self.read_u16()?,
                    atomic_number: self.read_u8()?,
                    symbol: self.read_string()?,
                })
            })
            .collect::<Result<Vec<_>, DatafileError>>()?;
        let species = (0..self.read_len()?)
            .map(|_| self.read_species())
            .collect::<Result<Vec<_>, DatafileError>>()?;
        Ok(ThermodynamicsDatabaseFile {
            magic,
            version,
            metadata,
            elements,
            species,
        })
    }

    fn read_species(&mut self) -> Result<SpeciesRecord, DatafileError> {
        let id = self.read_u16()?;
        let symbol = self.read_string()?;
        let composition = (0..self.read_len()?)
            .map(|_| {
                Ok(ElementCountRecord {
                    element_id: self.read_u16()?,
                    count: self.read_u16()?,
                })
            })
            .collect::<Result<Vec<_>, DatafileError>>()?;
        let charge_number = self.read_u8()? as i8;
        let phase = match self.read_u8()? {
            0 => PhaseRecord::Aqueous,
            1 => PhaseRecord::Solid,
            2 => PhaseRecord::Gas,
            _ => return Err(DatafileError::Decode),
        };
        let activity_model = match self.read_u8()? {
            0 => ActivityModelRecord::DaviesAqueous,
            1 => ActivityModelRecord::IdealMolalityAqueous,
            2 => ActivityModelRecord::IdealGas,
            3 => ActivityModelRecord::UnitActivity,
            _ => return Err(DatafileError::Decode),
        };
        let standard_gibbs_energy = self.read_molar_energy()?;
        let standard_enthalpy_of_formation = self.read_optional_molar_energy()?;
        let constant_pressure_heat_capacity = self.read_optional_heat_capacity()?;
        let valid_temperature_range = TemperatureRangeRecord {
            min_kelvin: self.read_f64()?,
            max_kelvin: self.read_f64()?,
        };
        let tags = (0..self.read_len()?)
            .map(|_| self.read_string())
            .collect::<Result<Vec<_>, DatafileError>>()?;
        Ok(SpeciesRecord {
            id,
            symbol,
            composition,
            charge_number,
            phase,
            activity_model,
            standard_gibbs_energy,
            standard_enthalpy_of_formation,
            constant_pressure_heat_capacity,
            valid_temperature_range,
            tags,
        })
    }

    fn read_optional_molar_energy(&mut self) -> Result<Option<MolarEnergyRecord>, DatafileError> {
        match self.read_u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.read_molar_energy()?)),
            _ => Err(DatafileError::Decode),
        }
    }

    fn read_optional_heat_capacity(&mut self) -> Result<Option<HeatCapacityRecord>, DatafileError> {
        match self.read_u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.read_heat_capacity()?)),
            _ => Err(DatafileError::Decode),
        }
    }

    fn read_molar_energy(&mut self) -> Result<MolarEnergyRecord, DatafileError> {
        Ok(MolarEnergyRecord {
            value_joule_per_mol: self.read_f64()?,
            reference_temperature_kelvin: self.read_f64()?,
            source: self.read_source()?,
        })
    }

    fn read_heat_capacity(&mut self) -> Result<HeatCapacityRecord, DatafileError> {
        Ok(HeatCapacityRecord {
            value_joule_per_mol_kelvin: self.read_f64()?,
            source: self.read_source()?,
        })
    }

    fn read_source(&mut self) -> Result<DataSourceRecord, DatafileError> {
        Ok(DataSourceRecord {
            citation: self.read_string()?,
            note: self.read_string()?,
        })
    }

    fn read_string(&mut self) -> Result<String, DatafileError> {
        let len = self.read_len()?;
        let end = self.offset.checked_add(len).ok_or(DatafileError::Decode)?;
        let bytes = self
            .bytes
            .get(self.offset..end)
            .ok_or(DatafileError::Decode)?;
        self.offset = end;
        String::from_utf8(bytes.to_vec()).map_err(|_| DatafileError::Decode)
    }

    fn read_magic(&mut self) -> Result<[u8; 4], DatafileError> {
        let bytes = self.take(4)?;
        Ok([bytes[0], bytes[1], bytes[2], bytes[3]])
    }

    fn read_len(&mut self) -> Result<usize, DatafileError> {
        Ok(self.read_u32()? as usize)
    }

    fn read_u8(&mut self) -> Result<u8, DatafileError> {
        Ok(*self.take(1)?.first().ok_or(DatafileError::Decode)?)
    }

    fn read_u16(&mut self) -> Result<u16, DatafileError> {
        let bytes = self.take(2)?;
        Ok(u16::from_le_bytes([bytes[0], bytes[1]]))
    }

    fn read_u32(&mut self) -> Result<u32, DatafileError> {
        let bytes = self.take(4)?;
        Ok(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
    }

    fn read_f64(&mut self) -> Result<f64, DatafileError> {
        let bytes = self.take(8)?;
        Ok(f64::from_le_bytes([
            bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
        ]))
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], DatafileError> {
        let end = self.offset.checked_add(len).ok_or(DatafileError::Decode)?;
        let bytes = self
            .bytes
            .get(self.offset..end)
            .ok_or(DatafileError::Decode)?;
        self.offset = end;
        Ok(bytes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source() -> DataSourceRecord {
        DataSourceRecord {
            citation: "fixture".to_owned(),
            note: "fixture".to_owned(),
        }
    }

    fn valid_database() -> ThermodynamicsDatabaseFile {
        ThermodynamicsDatabaseFile {
            magic: *FORMAT_MAGIC,
            version: FORMAT_VERSION,
            metadata: DatabaseMetadata {
                name: "fixture".to_owned(),
                source_summary: "fixture".to_owned(),
                license_note: "fixture".to_owned(),
            },
            elements: vec![ElementRecord {
                id: 1,
                atomic_number: 1,
                symbol: "H".to_owned(),
            }],
            species: vec![SpeciesRecord {
                id: 1,
                symbol: "H+".to_owned(),
                composition: vec![ElementCountRecord {
                    element_id: 1,
                    count: 1,
                }],
                charge_number: 1,
                phase: PhaseRecord::Aqueous,
                activity_model: ActivityModelRecord::DaviesAqueous,
                standard_gibbs_energy: MolarEnergyRecord {
                    value_joule_per_mol: 0.0,
                    reference_temperature_kelvin: 298.15,
                    source: source(),
                },
                standard_enthalpy_of_formation: None,
                constant_pressure_heat_capacity: None,
                valid_temperature_range: TemperatureRangeRecord {
                    min_kelvin: 273.15,
                    max_kelvin: 373.15,
                },
                tags: vec!["aqueous".to_owned(), "acid_base".to_owned()],
            }],
        }
    }

    #[test]
    fn valid_compressed_database_loads_registry() {
        let bytes = encode_database_file(&valid_database()).unwrap();
        let registry = registry_from_database_file(&bytes).unwrap();

        assert_eq!(registry.species(SpeciesId(1)).unwrap().symbol, "H+");
    }

    #[test]
    fn unsupported_version_is_rejected() {
        let mut database = valid_database();
        database.version = FORMAT_VERSION + 1;
        let bytes = encode_database_file(&database).unwrap();

        assert_eq!(
            decode_database_file(&bytes),
            Err(DatafileError::UnsupportedVersion(FORMAT_VERSION + 1))
        );
    }

    #[test]
    fn corrupt_bytes_are_rejected() {
        assert!(matches!(
            registry_from_database_file(b"not a database"),
            Err(DatafileError::Decode)
        ));
    }

    #[test]
    fn duplicate_species_are_rejected_by_registry_validation() {
        let mut database = valid_database();
        database.species.push(database.species[0].clone());

        assert!(matches!(
            database_to_registry(database),
            Err(DatafileError::Registry(
                SpeciesRegistryError::DuplicateSpecies(_)
            ))
        ));
    }
}
