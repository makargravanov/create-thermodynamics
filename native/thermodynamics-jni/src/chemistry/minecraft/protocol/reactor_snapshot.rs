use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::mixture::{
    CondensedPhaseAmount, Mixture, MixtureComponentPhaseSnapshot, MixturePhase, OrganicPhaseAmount,
    SolidComponentPhaseAmount,
};
use crate::chemistry::reactor::{
    Input, Output, PhaseEntry, Reactor, ReactorVolumeMode, ReactorZone, SubstanceEntry,
    TransitionMode, ZoneId, ZoneTransition,
};
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::SubstanceId;

use super::blob::{
    decode_sections, encode_sections, NativeBlob, NativeBlobKind, NativeBlobLimits,
    NativeBlobSection,
};

pub const REACTOR_SNAPSHOT_MODEL_VERSION: &str = "create-thermodynamics:reactor:1";

const SECTION_METADATA: u16 = 1;
const SECTION_ZONES: u16 = 2;
const SECTION_TRANSITIONS: u16 = 3;
const SECTION_INPUTS: u16 = 4;
const SECTION_OUTPUTS: u16 = 5;
const SECTION_OUTPUT_BUFFERS: u16 = 6;

const METADATA_VERSION: u16 = 1;
const ZONES_VERSION: u16 = 1;
const TRANSITIONS_VERSION: u16 = 1;
const INPUTS_VERSION: u16 = 1;
const OUTPUTS_VERSION: u16 = 1;
const OUTPUT_BUFFERS_VERSION: u16 = 1;

#[derive(Debug, Clone, PartialEq)]
pub struct ReactorSnapshotMetadata {
    pub catalog_version: String,
    pub content_version: u64,
    pub ambient_temperature_kelvin: Option<f64>,
    pub heat_transfer_coefficient_kw_per_kelvin: Option<f64>,
    pub last_ambient_energy_exchange_j: f64,
    pub vle_iterations: usize,
    pub vle_relaxation: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReactorZoneSnapshot {
    pub volume_cubic_meters: f64,
    pub volume_mode: ReactorVolumeMode,
    pub sealed: bool,
    pub elapsed_seconds: f64,
    pub mixture_temperature_kelvin: f64,
    pub mixture_gas_volume_cubic_meters: f64,
    pub mixture_components: Vec<MixtureComponentPhaseSnapshot>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OutputBufferSnapshot {
    pub output_index: usize,
    pub substances: Vec<OutputBufferSubstanceSnapshot>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OutputBufferSubstanceSnapshot {
    pub substance_id: SubstanceId,
    pub mol_per_bucket: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReactorSnapshot {
    pub metadata: ReactorSnapshotMetadata,
    pub zones: Vec<ReactorZoneSnapshot>,
    pub transitions: Vec<ZoneTransition>,
    pub inputs: Vec<Input>,
    pub outputs: Vec<Output>,
    pub output_buffers: Vec<OutputBufferSnapshot>,
}

impl ReactorSnapshot {
    pub fn from_reactor(
        registry: &ChemistryRegistry,
        reactor: &Reactor,
        catalog_version: impl Into<String>,
        content_version: u64,
        output_buffers: Vec<OutputBufferSnapshot>,
    ) -> ChemistryResult<Self> {
        let catalog_version = catalog_version.into();
        validate_non_empty("catalog version", &catalog_version)?;
        let zones = reactor
            .zones()
            .iter()
            .map(|zone| {
                Ok(ReactorZoneSnapshot {
                    volume_cubic_meters: zone.volume_cubic_meters(),
                    volume_mode: zone.volume_mode(),
                    sealed: zone.sealed(),
                    elapsed_seconds: zone.elapsed_seconds(),
                    mixture_temperature_kelvin: zone.mixture().temperature_kelvin(),
                    mixture_gas_volume_cubic_meters: zone.mixture().gas_volume_cubic_meters(),
                    mixture_components: zone.mixture().component_phase_snapshots(registry)?,
                })
            })
            .collect::<ChemistryResult<Vec<_>>>()?;
        let snapshot = Self {
            metadata: ReactorSnapshotMetadata {
                catalog_version,
                content_version,
                ambient_temperature_kelvin: reactor.ambient_temperature(),
                heat_transfer_coefficient_kw_per_kelvin: reactor.heat_transfer_coefficient(),
                last_ambient_energy_exchange_j: reactor.last_ambient_energy_exchange_j(),
                vle_iterations: reactor.vle_iterations(),
                vle_relaxation: reactor.vle_relaxation(),
            },
            zones,
            transitions: reactor.transitions().to_vec(),
            inputs: reactor.inputs().to_vec(),
            outputs: reactor.outputs().to_vec(),
            output_buffers,
        };
        snapshot.validate(registry)?;
        Ok(snapshot)
    }

    pub fn restore_reactor(&self, registry: &ChemistryRegistry) -> ChemistryResult<Reactor> {
        self.validate(registry)?;
        let zones = self
            .zones
            .iter()
            .map(|zone| {
                let mixture = Mixture::from_component_phase_snapshots(
                    registry,
                    zone.mixture_temperature_kelvin,
                    zone.mixture_gas_volume_cubic_meters,
                    &zone.mixture_components,
                )?;
                ReactorZone::from_parts(
                    mixture,
                    zone.volume_cubic_meters,
                    zone.volume_mode,
                    zone.sealed,
                    zone.elapsed_seconds,
                )
            })
            .collect::<ChemistryResult<Vec<_>>>()?;
        Reactor::from_parts(
            zones,
            self.transitions.clone(),
            self.inputs.clone(),
            self.outputs.clone(),
            self.metadata.ambient_temperature_kelvin,
            self.metadata.heat_transfer_coefficient_kw_per_kelvin,
            self.metadata.last_ambient_energy_exchange_j,
            self.metadata.vle_iterations,
            self.metadata.vle_relaxation,
        )
    }

    fn validate(&self, registry: &ChemistryRegistry) -> ChemistryResult<()> {
        validate_non_empty("catalog version", &self.metadata.catalog_version)?;
        validate_finite(
            "last ambient energy exchange",
            self.metadata.last_ambient_energy_exchange_j,
        )?;
        if self.metadata.vle_iterations == 0 {
            return Err(snapshot_error("vle iterations must be greater than zero"));
        }
        validate_range("vle relaxation", self.metadata.vle_relaxation, 0.01, 1.0)?;
        if let Some(value) = self.metadata.ambient_temperature_kelvin {
            validate_finite("ambient temperature", value)?;
        }
        if let Some(value) = self.metadata.heat_transfer_coefficient_kw_per_kelvin {
            validate_finite("heat transfer coefficient", value)?;
        }
        for zone in &self.zones {
            validate_finite_positive("zone volume", zone.volume_cubic_meters)?;
            validate_finite("zone elapsed seconds", zone.elapsed_seconds)?;
            validate_finite("mixture temperature", zone.mixture_temperature_kelvin)?;
            validate_finite_positive("mixture gas volume", zone.mixture_gas_volume_cubic_meters)?;
            for component in &zone.mixture_components {
                registry.substance(&component.substance_id)?;
            }
        }
        let zone_count = self.zones.len();
        for transition in &self.transitions {
            validate_zone("transition source", transition.from, zone_count)?;
            validate_zone("transition target", transition.to, zone_count)?;
            validate_transition_mode(&transition.mode, registry)?;
        }
        for input in &self.inputs {
            validate_zone("input target", input.to, zone_count)?;
            validate_transition_mode(&input.mode, registry)?;
        }
        for output in &self.outputs {
            validate_zone("output source", output.from, zone_count)?;
            validate_transition_mode(&output.mode, registry)?;
        }
        for buffer in &self.output_buffers {
            for substance in &buffer.substances {
                registry.substance(&substance.substance_id)?;
                validate_finite("output buffer amount", substance.mol_per_bucket)?;
                if substance.mol_per_bucket < 0.0 {
                    return Err(snapshot_error("output buffer amount must not be negative"));
                }
            }
        }
        Ok(())
    }
}

pub fn export_reactor_checkpoint(
    registry: &ChemistryRegistry,
    reactor: &Reactor,
    catalog_version: impl Into<String>,
    content_version: u64,
    output_buffers: Vec<OutputBufferSnapshot>,
    limits: &NativeBlobLimits,
) -> ChemistryResult<Vec<u8>> {
    let snapshot = ReactorSnapshot::from_reactor(
        registry,
        reactor,
        catalog_version,
        content_version,
        output_buffers,
    )?;
    let payload = encode_snapshot_sections(&snapshot)?;
    NativeBlob::encode(
        NativeBlobKind::ReactorSnapshot,
        REACTOR_SNAPSHOT_MODEL_VERSION,
        content_version,
        &payload,
        limits,
    )
}

pub fn read_reactor_checkpoint(
    encoded: &[u8],
    expected_catalog_version: &str,
    registry: &ChemistryRegistry,
    limits: &NativeBlobLimits,
) -> ChemistryResult<ReactorSnapshot> {
    validate_non_empty("expected catalog version", expected_catalog_version)?;
    let (blob, payload) =
        NativeBlob::decode_expected(encoded, NativeBlobKind::ReactorSnapshot, limits)?;
    if blob.model_version != REACTOR_SNAPSHOT_MODEL_VERSION {
        return Err(snapshot_error(format!(
            "reactor snapshot has model version '{}', expected '{}'",
            blob.model_version, REACTOR_SNAPSHOT_MODEL_VERSION
        )));
    }
    let snapshot = decode_snapshot_sections(&payload, limits)?;
    snapshot.validate(registry)?;
    if snapshot.metadata.content_version != blob.content_version {
        return Err(snapshot_error(format!(
            "reactor metadata content version {} does not match blob content version {}",
            snapshot.metadata.content_version, blob.content_version
        )));
    }
    if snapshot.metadata.catalog_version != expected_catalog_version {
        return Err(snapshot_error(format!(
            "reactor snapshot is for catalog '{}', expected '{}'",
            snapshot.metadata.catalog_version, expected_catalog_version
        )));
    }
    Ok(snapshot)
}

fn encode_snapshot_sections(snapshot: &ReactorSnapshot) -> ChemistryResult<Vec<u8>> {
    encode_sections(&[
        NativeBlobSection {
            section_kind: SECTION_METADATA,
            payload: encode_metadata(&snapshot.metadata)?,
        },
        NativeBlobSection {
            section_kind: SECTION_ZONES,
            payload: encode_zones(&snapshot.zones)?,
        },
        NativeBlobSection {
            section_kind: SECTION_TRANSITIONS,
            payload: encode_transitions(TRANSITIONS_VERSION, &snapshot.transitions)?,
        },
        NativeBlobSection {
            section_kind: SECTION_INPUTS,
            payload: encode_inputs(&snapshot.inputs)?,
        },
        NativeBlobSection {
            section_kind: SECTION_OUTPUTS,
            payload: encode_outputs(&snapshot.outputs)?,
        },
        NativeBlobSection {
            section_kind: SECTION_OUTPUT_BUFFERS,
            payload: encode_output_buffers(&snapshot.output_buffers)?,
        },
    ])
}

fn decode_snapshot_sections(
    payload: &[u8],
    limits: &NativeBlobLimits,
) -> ChemistryResult<ReactorSnapshot> {
    let mut metadata = None;
    let mut zones = None;
    let mut transitions = None;
    let mut inputs = None;
    let mut outputs = None;
    let mut output_buffers = None;
    for section in decode_sections(payload, limits)? {
        match section.section_kind {
            SECTION_METADATA => assign_section(
                &mut metadata,
                decode_metadata(&section.payload)?,
                "metadata",
            )?,
            SECTION_ZONES => assign_section(&mut zones, decode_zones(&section.payload)?, "zones")?,
            SECTION_TRANSITIONS => assign_section(
                &mut transitions,
                decode_transitions(&section.payload, TRANSITIONS_VERSION)?,
                "transitions",
            )?,
            SECTION_INPUTS => {
                assign_section(&mut inputs, decode_inputs(&section.payload)?, "inputs")?
            }
            SECTION_OUTPUTS => {
                assign_section(&mut outputs, decode_outputs(&section.payload)?, "outputs")?
            }
            SECTION_OUTPUT_BUFFERS => assign_section(
                &mut output_buffers,
                decode_output_buffers(&section.payload)?,
                "output buffers",
            )?,
            other => {
                return Err(snapshot_error(format!(
                    "unknown reactor snapshot section {other}"
                )))
            }
        }
    }
    Ok(ReactorSnapshot {
        metadata: metadata.ok_or_else(|| snapshot_error("reactor snapshot is missing metadata"))?,
        zones: zones.ok_or_else(|| snapshot_error("reactor snapshot is missing zones"))?,
        transitions: transitions
            .ok_or_else(|| snapshot_error("reactor snapshot is missing transitions"))?,
        inputs: inputs.ok_or_else(|| snapshot_error("reactor snapshot is missing inputs"))?,
        outputs: outputs.ok_or_else(|| snapshot_error("reactor snapshot is missing outputs"))?,
        output_buffers: output_buffers
            .ok_or_else(|| snapshot_error("reactor snapshot is missing output buffers"))?,
    })
}

fn encode_metadata(metadata: &ReactorSnapshotMetadata) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, METADATA_VERSION);
    write_string(&mut output, &metadata.catalog_version)?;
    write_u64(&mut output, metadata.content_version);
    write_optional_f64(&mut output, metadata.ambient_temperature_kelvin);
    write_optional_f64(
        &mut output,
        metadata.heat_transfer_coefficient_kw_per_kelvin,
    );
    write_f64(&mut output, metadata.last_ambient_energy_exchange_j);
    write_u64(&mut output, metadata.vle_iterations as u64);
    write_f64(&mut output, metadata.vle_relaxation);
    Ok(output)
}

fn decode_metadata(payload: &[u8]) -> ChemistryResult<ReactorSnapshotMetadata> {
    let mut cursor = Cursor::new(payload);
    expect_version(cursor.read_u16()?, METADATA_VERSION, "metadata")?;
    let metadata = ReactorSnapshotMetadata {
        catalog_version: cursor.read_string()?,
        content_version: cursor.read_u64()?,
        ambient_temperature_kelvin: cursor.read_optional_f64()?,
        heat_transfer_coefficient_kw_per_kelvin: cursor.read_optional_f64()?,
        last_ambient_energy_exchange_j: cursor.read_f64()?,
        vle_iterations: cursor.read_len()?,
        vle_relaxation: cursor.read_f64()?,
    };
    cursor.finish()?;
    Ok(metadata)
}

fn encode_zones(zones: &[ReactorZoneSnapshot]) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, ZONES_VERSION);
    write_u64(&mut output, zones.len() as u64);
    for zone in zones {
        write_f64(&mut output, zone.volume_cubic_meters);
        write_u8(&mut output, volume_mode_to_wire(zone.volume_mode));
        write_bool(&mut output, zone.sealed);
        write_f64(&mut output, zone.elapsed_seconds);
        write_f64(&mut output, zone.mixture_temperature_kelvin);
        write_f64(&mut output, zone.mixture_gas_volume_cubic_meters);
        write_components(&mut output, &zone.mixture_components)?;
    }
    Ok(output)
}

fn decode_zones(payload: &[u8]) -> ChemistryResult<Vec<ReactorZoneSnapshot>> {
    let mut cursor = Cursor::new(payload);
    expect_version(cursor.read_u16()?, ZONES_VERSION, "zones")?;
    let count = cursor.read_len()?;
    let mut zones = Vec::with_capacity(count);
    for _ in 0..count {
        zones.push(ReactorZoneSnapshot {
            volume_cubic_meters: cursor.read_f64()?,
            volume_mode: volume_mode_from_wire(cursor.read_u8()?)?,
            sealed: cursor.read_bool()?,
            elapsed_seconds: cursor.read_f64()?,
            mixture_temperature_kelvin: cursor.read_f64()?,
            mixture_gas_volume_cubic_meters: cursor.read_f64()?,
            mixture_components: read_components(&mut cursor)?,
        });
    }
    cursor.finish()?;
    Ok(zones)
}

fn write_components(
    output: &mut Vec<u8>,
    components: &[MixtureComponentPhaseSnapshot],
) -> ChemistryResult<()> {
    write_u64(output, components.len() as u64);
    for component in components {
        write_substance_id(output, &component.substance_id)?;
        write_f64(output, component.aqueous_mol_per_bucket);
        write_u64(
            output,
            component.organic_mol_per_bucket_by_solvent.len() as u64,
        );
        for amount in &component.organic_mol_per_bucket_by_solvent {
            write_substance_id(output, &amount.solvent_id)?;
            write_f64(output, amount.concentration_mol_per_bucket);
        }
        write_u64(
            output,
            component.molten_mol_per_bucket_by_phase.len() as u64,
        );
        for amount in &component.molten_mol_per_bucket_by_phase {
            write_u8(output, phase_to_wire(amount.coarse_phase));
            write_substance_id(output, &amount.anchor_substance_id)?;
            write_f64(output, amount.concentration_mol_per_bucket);
        }
        write_f64(output, component.gas_mol_per_bucket);
        write_f64(output, component.supercritical_mol_per_bucket);
        write_u64(output, component.solid_mol_per_bucket_by_phase.len() as u64);
        for amount in &component.solid_mol_per_bucket_by_phase {
            write_substance_id(output, &amount.anchor_substance_id)?;
            write_f64(output, amount.concentration_mol_per_bucket);
        }
    }
    Ok(())
}

fn read_components(cursor: &mut Cursor<'_>) -> ChemistryResult<Vec<MixtureComponentPhaseSnapshot>> {
    let count = cursor.read_len()?;
    let mut components = Vec::with_capacity(count);
    for _ in 0..count {
        let substance_id = cursor.read_substance_id()?;
        let aqueous_mol_per_bucket = cursor.read_f64()?;
        let organic_count = cursor.read_len()?;
        let mut organic_mol_per_bucket_by_solvent = Vec::with_capacity(organic_count);
        for _ in 0..organic_count {
            organic_mol_per_bucket_by_solvent.push(OrganicPhaseAmount {
                solvent_id: cursor.read_substance_id()?,
                concentration_mol_per_bucket: cursor.read_f64()?,
            });
        }
        let molten_count = cursor.read_len()?;
        let mut molten_mol_per_bucket_by_phase = Vec::with_capacity(molten_count);
        for _ in 0..molten_count {
            molten_mol_per_bucket_by_phase.push(CondensedPhaseAmount {
                coarse_phase: phase_from_wire(cursor.read_u8()?)?,
                anchor_substance_id: cursor.read_substance_id()?,
                concentration_mol_per_bucket: cursor.read_f64()?,
            });
        }
        let gas_mol_per_bucket = cursor.read_f64()?;
        let supercritical_mol_per_bucket = cursor.read_f64()?;
        let solid_count = cursor.read_len()?;
        let mut solid_mol_per_bucket_by_phase = Vec::with_capacity(solid_count);
        for _ in 0..solid_count {
            solid_mol_per_bucket_by_phase.push(SolidComponentPhaseAmount {
                anchor_substance_id: cursor.read_substance_id()?,
                concentration_mol_per_bucket: cursor.read_f64()?,
            });
        }
        components.push(MixtureComponentPhaseSnapshot {
            substance_id,
            aqueous_mol_per_bucket,
            organic_mol_per_bucket_by_solvent,
            molten_mol_per_bucket_by_phase,
            gas_mol_per_bucket,
            supercritical_mol_per_bucket,
            solid_mol_per_bucket_by_phase,
        });
    }
    Ok(components)
}

fn encode_inputs(inputs: &[Input]) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, INPUTS_VERSION);
    write_u64(&mut output, inputs.len() as u64);
    for input in inputs {
        write_u64(&mut output, input.to.0 as u64);
        write_bool(&mut output, input.enabled);
        write_transition_mode(&mut output, &input.mode)?;
    }
    Ok(output)
}

fn decode_inputs(payload: &[u8]) -> ChemistryResult<Vec<Input>> {
    let mut cursor = Cursor::new(payload);
    expect_version(cursor.read_u16()?, INPUTS_VERSION, "inputs")?;
    let count = cursor.read_len()?;
    let mut inputs = Vec::with_capacity(count);
    for _ in 0..count {
        inputs.push(Input {
            to: ZoneId(cursor.read_len()?),
            enabled: cursor.read_bool()?,
            mode: read_transition_mode(&mut cursor)?,
        });
    }
    cursor.finish()?;
    Ok(inputs)
}

fn encode_outputs(outputs: &[Output]) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, OUTPUTS_VERSION);
    write_u64(&mut output, outputs.len() as u64);
    for output_entry in outputs {
        write_u64(&mut output, output_entry.from.0 as u64);
        write_bool(&mut output, output_entry.enabled);
        write_transition_mode(&mut output, &output_entry.mode)?;
    }
    Ok(output)
}

fn decode_outputs(payload: &[u8]) -> ChemistryResult<Vec<Output>> {
    let mut cursor = Cursor::new(payload);
    expect_version(cursor.read_u16()?, OUTPUTS_VERSION, "outputs")?;
    let count = cursor.read_len()?;
    let mut outputs = Vec::with_capacity(count);
    for _ in 0..count {
        outputs.push(Output {
            from: ZoneId(cursor.read_len()?),
            enabled: cursor.read_bool()?,
            mode: read_transition_mode(&mut cursor)?,
        });
    }
    cursor.finish()?;
    Ok(outputs)
}

fn encode_transitions(version: u16, transitions: &[ZoneTransition]) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, version);
    write_u64(&mut output, transitions.len() as u64);
    for transition in transitions {
        write_u64(&mut output, transition.from.0 as u64);
        write_u64(&mut output, transition.to.0 as u64);
        write_bool(&mut output, transition.enabled);
        write_transition_mode(&mut output, &transition.mode)?;
    }
    Ok(output)
}

fn decode_transitions(payload: &[u8], version: u16) -> ChemistryResult<Vec<ZoneTransition>> {
    let mut cursor = Cursor::new(payload);
    expect_version(cursor.read_u16()?, version, "transitions")?;
    let count = cursor.read_len()?;
    let mut transitions = Vec::with_capacity(count);
    for _ in 0..count {
        transitions.push(ZoneTransition {
            from: ZoneId(cursor.read_len()?),
            to: ZoneId(cursor.read_len()?),
            enabled: cursor.read_bool()?,
            mode: read_transition_mode(&mut cursor)?,
        });
    }
    cursor.finish()?;
    Ok(transitions)
}

fn write_transition_mode(output: &mut Vec<u8>, mode: &TransitionMode) -> ChemistryResult<()> {
    match mode {
        TransitionMode::Substances { entries } => {
            write_u8(output, 1);
            write_substance_entries(output, entries)?;
        }
        TransitionMode::Phases { entries } => {
            write_u8(output, 2);
            write_u64(output, entries.len() as u64);
            for entry in entries {
                write_u8(output, phase_to_wire(entry.phase));
                write_f64(output, entry.rate_mol_per_second);
            }
        }
        TransitionMode::All {
            rate_mol_per_second,
        } => {
            write_u8(output, 3);
            write_f64(output, *rate_mol_per_second);
        }
        TransitionMode::SubstancesThreshold {
            entries,
            threshold_mol_per_bucket,
        } => {
            write_u8(output, 4);
            write_f64(output, *threshold_mol_per_bucket);
            write_substance_entries(output, entries)?;
        }
    }
    Ok(())
}

fn read_transition_mode(cursor: &mut Cursor<'_>) -> ChemistryResult<TransitionMode> {
    match cursor.read_u8()? {
        1 => Ok(TransitionMode::Substances {
            entries: read_substance_entries(cursor)?,
        }),
        2 => {
            let count = cursor.read_len()?;
            let mut entries = Vec::with_capacity(count);
            for _ in 0..count {
                entries.push(PhaseEntry {
                    phase: phase_from_wire(cursor.read_u8()?)?,
                    rate_mol_per_second: cursor.read_f64()?,
                });
            }
            Ok(TransitionMode::Phases { entries })
        }
        3 => Ok(TransitionMode::All {
            rate_mol_per_second: cursor.read_f64()?,
        }),
        4 => {
            let threshold_mol_per_bucket = cursor.read_f64()?;
            Ok(TransitionMode::SubstancesThreshold {
                entries: read_substance_entries(cursor)?,
                threshold_mol_per_bucket,
            })
        }
        other => Err(snapshot_error(format!("unknown transition mode {other}"))),
    }
}

fn write_substance_entries(
    output: &mut Vec<u8>,
    entries: &[SubstanceEntry],
) -> ChemistryResult<()> {
    write_u64(output, entries.len() as u64);
    for entry in entries {
        write_substance_id(output, &entry.id)?;
        write_f64(output, entry.rate_mol_per_second);
    }
    Ok(())
}

fn read_substance_entries(cursor: &mut Cursor<'_>) -> ChemistryResult<Vec<SubstanceEntry>> {
    let count = cursor.read_len()?;
    let mut entries = Vec::with_capacity(count);
    for _ in 0..count {
        entries.push(SubstanceEntry {
            id: cursor.read_substance_id()?,
            rate_mol_per_second: cursor.read_f64()?,
        });
    }
    Ok(entries)
}

fn encode_output_buffers(buffers: &[OutputBufferSnapshot]) -> ChemistryResult<Vec<u8>> {
    let mut output = Vec::new();
    write_u16(&mut output, OUTPUT_BUFFERS_VERSION);
    write_u64(&mut output, buffers.len() as u64);
    for buffer in buffers {
        write_u64(&mut output, buffer.output_index as u64);
        write_u64(&mut output, buffer.substances.len() as u64);
        for substance in &buffer.substances {
            write_substance_id(&mut output, &substance.substance_id)?;
            write_f64(&mut output, substance.mol_per_bucket);
        }
    }
    Ok(output)
}

fn decode_output_buffers(payload: &[u8]) -> ChemistryResult<Vec<OutputBufferSnapshot>> {
    let mut cursor = Cursor::new(payload);
    expect_version(cursor.read_u16()?, OUTPUT_BUFFERS_VERSION, "output buffers")?;
    let count = cursor.read_len()?;
    let mut buffers = Vec::with_capacity(count);
    for _ in 0..count {
        let output_index = cursor.read_len()?;
        let substance_count = cursor.read_len()?;
        let mut substances = Vec::with_capacity(substance_count);
        for _ in 0..substance_count {
            substances.push(OutputBufferSubstanceSnapshot {
                substance_id: cursor.read_substance_id()?,
                mol_per_bucket: cursor.read_f64()?,
            });
        }
        buffers.push(OutputBufferSnapshot {
            output_index,
            substances,
        });
    }
    cursor.finish()?;
    Ok(buffers)
}

fn validate_transition_mode(
    mode: &TransitionMode,
    registry: &ChemistryRegistry,
) -> ChemistryResult<()> {
    match mode {
        TransitionMode::Substances { entries }
        | TransitionMode::SubstancesThreshold { entries, .. } => {
            for entry in entries {
                registry.substance(&entry.id)?;
                validate_non_negative_finite("substance transfer rate", entry.rate_mol_per_second)?;
            }
            if let TransitionMode::SubstancesThreshold {
                threshold_mol_per_bucket,
                ..
            } = mode
            {
                validate_non_negative_finite("substance threshold", *threshold_mol_per_bucket)?;
            }
        }
        TransitionMode::Phases { entries } => {
            for entry in entries {
                validate_non_negative_finite("phase transfer rate", entry.rate_mol_per_second)?;
            }
        }
        TransitionMode::All {
            rate_mol_per_second,
        } => validate_non_negative_finite("all transfer rate", *rate_mol_per_second)?,
    }
    Ok(())
}

fn validate_zone(name: &str, zone_id: ZoneId, zone_count: usize) -> ChemistryResult<()> {
    if zone_id.0 < zone_count {
        Ok(())
    } else {
        Err(snapshot_error(format!(
            "{name} zone index {} out of range (have {zone_count})",
            zone_id.0
        )))
    }
}

fn validate_non_empty(name: &str, value: &str) -> ChemistryResult<()> {
    if value.trim().is_empty() {
        Err(snapshot_error(format!("{name} must not be empty")))
    } else {
        Ok(())
    }
}

fn validate_finite(name: &str, value: f64) -> ChemistryResult<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(snapshot_error(format!(
            "{name} must be finite, got {value}"
        )))
    }
}

fn validate_finite_positive(name: &str, value: f64) -> ChemistryResult<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(snapshot_error(format!(
            "{name} must be positive and finite, got {value}"
        )))
    }
}

fn validate_non_negative_finite(name: &str, value: f64) -> ChemistryResult<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(snapshot_error(format!(
            "{name} must be non-negative and finite, got {value}"
        )))
    }
}

fn validate_range(name: &str, value: f64, min: f64, max: f64) -> ChemistryResult<()> {
    if value.is_finite() && value >= min && value <= max {
        Ok(())
    } else {
        Err(snapshot_error(format!(
            "{name} must be finite and in {min}..={max}, got {value}"
        )))
    }
}

fn assign_section<T>(target: &mut Option<T>, value: T, name: &str) -> ChemistryResult<()> {
    if target.is_some() {
        return Err(snapshot_error(format!(
            "reactor snapshot has duplicate {name} section"
        )));
    }
    *target = Some(value);
    Ok(())
}

fn expect_version(actual: u16, expected: u16, name: &str) -> ChemistryResult<()> {
    if actual == expected {
        Ok(())
    } else {
        Err(snapshot_error(format!(
            "unsupported reactor {name} version {actual}"
        )))
    }
}

fn write_u8(output: &mut Vec<u8>, value: u8) {
    output.push(value);
}

fn write_bool(output: &mut Vec<u8>, value: bool) {
    output.push(u8::from(value));
}

fn write_u16(output: &mut Vec<u8>, value: u16) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_u64(output: &mut Vec<u8>, value: u64) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_f64(output: &mut Vec<u8>, value: f64) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_optional_f64(output: &mut Vec<u8>, value: Option<f64>) {
    match value {
        Some(value) => {
            write_u8(output, 1);
            write_f64(output, value);
        }
        None => write_u8(output, 0),
    }
}

fn write_string(output: &mut Vec<u8>, value: &str) -> ChemistryResult<()> {
    validate_non_empty("string", value)?;
    write_u64(output, value.len() as u64);
    output.extend_from_slice(value.as_bytes());
    Ok(())
}

fn write_substance_id(output: &mut Vec<u8>, value: &SubstanceId) -> ChemistryResult<()> {
    write_string(output, value.as_str())
}

fn phase_to_wire(phase: MixturePhase) -> u8 {
    match phase {
        MixturePhase::Aqueous => 1,
        MixturePhase::Organic => 2,
        MixturePhase::MoltenMetal => 3,
        MixturePhase::MoltenSlag => 4,
        MixturePhase::Gas => 5,
        MixturePhase::SupercriticalFluid => 6,
        MixturePhase::Solid => 7,
    }
}

fn phase_from_wire(value: u8) -> ChemistryResult<MixturePhase> {
    match value {
        1 => Ok(MixturePhase::Aqueous),
        2 => Ok(MixturePhase::Organic),
        3 => Ok(MixturePhase::MoltenMetal),
        4 => Ok(MixturePhase::MoltenSlag),
        5 => Ok(MixturePhase::Gas),
        6 => Ok(MixturePhase::SupercriticalFluid),
        7 => Ok(MixturePhase::Solid),
        _ => Err(snapshot_error(format!("unknown mixture phase {value}"))),
    }
}

fn volume_mode_to_wire(mode: ReactorVolumeMode) -> u8 {
    match mode {
        ReactorVolumeMode::HeadspaceRequired => 1,
        ReactorVolumeMode::LiquidFilled => 2,
    }
}

fn volume_mode_from_wire(value: u8) -> ChemistryResult<ReactorVolumeMode> {
    match value {
        1 => Ok(ReactorVolumeMode::HeadspaceRequired),
        2 => Ok(ReactorVolumeMode::LiquidFilled),
        _ => Err(snapshot_error(format!(
            "unknown reactor volume mode {value}"
        ))),
    }
}

fn snapshot_error(reason: impl Into<String>) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!("invalid reactor snapshot: {}", reason.into()))
}

struct Cursor<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Cursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    fn read_u8(&mut self) -> ChemistryResult<u8> {
        Ok(self.read_slice(1)?[0])
    }

    fn read_bool(&mut self) -> ChemistryResult<bool> {
        match self.read_u8()? {
            0 => Ok(false),
            1 => Ok(true),
            other => Err(snapshot_error(format!("invalid bool marker {other}"))),
        }
    }

    fn read_u16(&mut self) -> ChemistryResult<u16> {
        Ok(u16::from_le_bytes(self.read_array::<2>()?))
    }

    fn read_u64(&mut self) -> ChemistryResult<u64> {
        Ok(u64::from_le_bytes(self.read_array::<8>()?))
    }

    fn read_len(&mut self) -> ChemistryResult<usize> {
        usize::try_from(self.read_u64()?)
            .map_err(|_| snapshot_error("length does not fit in usize"))
    }

    fn read_f64(&mut self) -> ChemistryResult<f64> {
        Ok(f64::from_le_bytes(self.read_array::<8>()?))
    }

    fn read_optional_f64(&mut self) -> ChemistryResult<Option<f64>> {
        match self.read_u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.read_f64()?)),
            other => Err(snapshot_error(format!(
                "invalid optional f64 marker {other}"
            ))),
        }
    }

    fn read_string(&mut self) -> ChemistryResult<String> {
        let len = self.read_len()?;
        let bytes = self.read_slice(len)?;
        let value = String::from_utf8(bytes.to_vec())
            .map_err(|error| snapshot_error(format!("invalid UTF-8 string: {error}")))?;
        validate_non_empty("string", &value)?;
        Ok(value)
    }

    fn read_substance_id(&mut self) -> ChemistryResult<SubstanceId> {
        let value = self.read_string()?;
        Ok(SubstanceId::from(value.as_str()))
    }

    fn read_array<const N: usize>(&mut self) -> ChemistryResult<[u8; N]> {
        let slice = self.read_slice(N)?;
        let mut output = [0_u8; N];
        output.copy_from_slice(slice);
        Ok(output)
    }

    fn read_slice(&mut self, len: usize) -> ChemistryResult<&'a [u8]> {
        let end = self
            .offset
            .checked_add(len)
            .ok_or_else(|| snapshot_error("offset overflow"))?;
        if end > self.bytes.len() {
            return Err(snapshot_error(format!(
                "snapshot ended early at byte {}, needed {end}",
                self.bytes.len()
            )));
        }
        let slice = &self.bytes[self.offset..end];
        self.offset = end;
        Ok(slice)
    }

    fn finish(self) -> ChemistryResult<()> {
        if self.offset != self.bytes.len() {
            Err(snapshot_error(format!(
                "snapshot has {} trailing bytes",
                self.bytes.len() - self.offset
            )))
        } else {
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::catalog;

    fn limits() -> NativeBlobLimits {
        NativeBlobLimits::new(1024 * 1024, 8 * 1024 * 1024)
    }

    fn registry() -> ChemistryRegistry {
        catalog::destroy_substances_registry_builder()
            .unwrap()
            .build()
            .unwrap()
    }

    #[test]
    fn reactor_snapshot_roundtrip_preserves_reactor_and_output_buffers() {
        let registry = registry();
        let water = SubstanceId::from("destroy:water");
        let oxygen = SubstanceId::from("destroy:oxygen");
        let mut mixture = Mixture::new(420.0).unwrap();
        mixture.set_gas_volume_cubic_meters(0.0004).unwrap();
        mixture
            .add_substance(&registry, water.clone(), 1.0)
            .unwrap();
        mixture
            .set_gaseous_fraction(&registry, water.clone(), 0.25)
            .unwrap();
        mixture
            .add_substance(&registry, oxygen.clone(), 0.2)
            .unwrap();
        let zone = ReactorZone::from_parts(
            mixture,
            0.001,
            ReactorVolumeMode::HeadspaceRequired,
            true,
            12.5,
        )
        .unwrap();
        let mut reactor = Reactor::new();
        let zone_id = reactor.add_zone(zone);
        reactor.set_ambient_temperature(300.0);
        reactor.set_heat_transfer_coefficient(0.05);
        reactor.set_vle_iterations(3);
        reactor.set_vle_relaxation(0.5);
        reactor.add_input(Input::new(
            zone_id,
            TransitionMode::Substances {
                entries: vec![SubstanceEntry {
                    id: water.clone(),
                    rate_mol_per_second: 0.1,
                }],
            },
        ));
        reactor.add_output(Output::new(
            zone_id,
            TransitionMode::Phases {
                entries: vec![PhaseEntry {
                    phase: MixturePhase::Gas,
                    rate_mol_per_second: 0.05,
                }],
            },
        ));
        let expected_water_gas = reactor
            .zone(&ZoneId(0))
            .unwrap()
            .mixture()
            .concentration_in_phase(&water, MixturePhase::Gas);
        let output_buffers = vec![OutputBufferSnapshot {
            output_index: 0,
            substances: vec![OutputBufferSubstanceSnapshot {
                substance_id: water.clone(),
                mol_per_bucket: 0.125,
            }],
        }];

        let encoded = export_reactor_checkpoint(
            &registry,
            &reactor,
            "test-catalog",
            7,
            output_buffers.clone(),
            &limits(),
        )
        .unwrap();
        let decoded =
            read_reactor_checkpoint(&encoded, "test-catalog", &registry, &limits()).unwrap();
        let restored = decoded.restore_reactor(&registry).unwrap();

        assert_eq!(decoded.output_buffers, output_buffers);
        assert_eq!(restored.zone_count(), 1);
        assert_eq!(restored.input_count(), 1);
        assert_eq!(restored.output_count(), 1);
        assert_eq!(restored.vle_iterations(), 3);
        assert_eq!(restored.vle_relaxation(), 0.5);
        let restored_zone = restored.zone(&ZoneId(0)).unwrap();
        assert!(restored_zone.sealed());
        assert_eq!(
            restored_zone.volume_mode(),
            ReactorVolumeMode::HeadspaceRequired
        );
        assert!((restored_zone.elapsed_seconds() - 12.5).abs() < 1.0e-12);
        assert!((restored_zone.mixture().gas_volume_cubic_meters() - 0.0004).abs() < 1.0e-12);
        assert!((restored_zone.concentration_of(&water) - 1.0).abs() < 1.0e-12);
        assert!(
            (restored_zone
                .mixture()
                .concentration_in_phase(&water, MixturePhase::Gas)
                - expected_water_gas)
                .abs()
                < 1.0e-12
        );
    }

    #[test]
    fn reactor_snapshot_rejects_wrong_catalog_version() {
        let registry = registry();
        let reactor = Reactor::new();
        let encoded =
            export_reactor_checkpoint(&registry, &reactor, "catalog-a", 1, Vec::new(), &limits())
                .unwrap();
        let error =
            read_reactor_checkpoint(&encoded, "catalog-b", &registry, &limits()).unwrap_err();

        assert!(error.to_string().contains("expected 'catalog-b'"));
    }
}
