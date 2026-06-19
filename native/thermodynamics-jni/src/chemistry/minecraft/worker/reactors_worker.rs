use std::collections::BTreeMap;

use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::minecraft::mregistry::item_to_substance::MinecraftId;
use crate::chemistry::minecraft::mregistry::mregistry::MinecraftChemicalRegistry;
use crate::chemistry::minecraft::protocol::blob::NativeBlobLimits;
use crate::chemistry::minecraft::protocol::reactor_snapshot::{
    export_reactor_checkpoint, read_reactor_checkpoint, OutputBufferSnapshot,
    OutputBufferSubstanceSnapshot,
};
use crate::chemistry::reactor::io;
use crate::chemistry::reactor::{Reactor, TransitionMode};
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::SubstanceId;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReactorInstanceId(pub u64);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OutputPortId {
    pub reactor_id: ReactorInstanceId,
    pub output_index: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ItemInputReport {
    pub reactor_id: ReactorInstanceId,
    pub input_index: usize,
    pub item_id: MinecraftId,
    pub substance_id: SubstanceId,
    pub items_consumed: u32,
    pub mol_inserted: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OutputDrainReport {
    pub reactor_id: ReactorInstanceId,
    pub output_index: usize,
    pub extracted: Vec<OutputSubstanceAmount>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OutputSubstanceAmount {
    pub substance_id: SubstanceId,
    pub mol_per_bucket: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ItemOutputReport {
    pub reactor_id: ReactorInstanceId,
    pub output_index: usize,
    pub item_id: MinecraftId,
    pub substance_id: SubstanceId,
    pub items_extracted: u32,
    pub mol_consumed: f64,
    pub mol_remaining_in_buffer: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReactorTickReport {
    pub reactor_id: ReactorInstanceId,
    pub dt_seconds: f64,
    pub zone_count: usize,
    pub total_electrical_draw_w: f64,
}

struct ReactorRuntime {
    reactor: Reactor,
    output_buffers_mol_per_bucket: BTreeMap<usize, BTreeMap<SubstanceId, f64>>,
}

#[derive(Default)]
pub struct ReactorsWorker {
    reactors: BTreeMap<ReactorInstanceId, ReactorRuntime>,
    next_reactor_id: u64,
}

impl ReactorsWorker {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register_reactor(&mut self, reactor: Reactor) -> ReactorInstanceId {
        let id = ReactorInstanceId(self.next_reactor_id);
        self.next_reactor_id += 1;
        self.register_runtime(
            id,
            ReactorRuntime {
                reactor,
                output_buffers_mol_per_bucket: BTreeMap::new(),
            },
        );
        id
    }

    pub fn register_reactor_from_checkpoint(
        &mut self,
        registry: &ChemistryRegistry,
        encoded: &[u8],
        expected_catalog_version: &str,
        limits: &NativeBlobLimits,
    ) -> ChemistryResult<ReactorInstanceId> {
        let snapshot =
            read_reactor_checkpoint(encoded, expected_catalog_version, registry, limits)?;
        let reactor = snapshot.restore_reactor(registry)?;
        let output_buffers_mol_per_bucket =
            output_buffers_from_snapshots(registry, &snapshot.output_buffers)?;
        let id = ReactorInstanceId(self.next_reactor_id);
        self.next_reactor_id += 1;
        self.register_runtime(
            id,
            ReactorRuntime {
                reactor,
                output_buffers_mol_per_bucket,
            },
        );
        Ok(id)
    }

    pub fn remove_reactor(&mut self, reactor_id: ReactorInstanceId) -> ChemistryResult<Reactor> {
        self.reactors
            .remove(&reactor_id)
            .map(|runtime| runtime.reactor)
            .ok_or_else(|| unknown_reactor_error(reactor_id))
    }

    pub fn reactor(&self, reactor_id: ReactorInstanceId) -> ChemistryResult<&Reactor> {
        Ok(&self.runtime(reactor_id)?.reactor)
    }

    pub fn reactor_mut(&mut self, reactor_id: ReactorInstanceId) -> ChemistryResult<&mut Reactor> {
        Ok(&mut self.runtime_mut(reactor_id)?.reactor)
    }

    pub fn reactor_count(&self) -> usize {
        self.reactors.len()
    }

    pub fn reactor_ids(&self) -> impl Iterator<Item = ReactorInstanceId> + '_ {
        self.reactors.keys().copied()
    }

    pub fn tick_reactor(
        &mut self,
        registry: &ChemistryRegistry,
        reactor_id: ReactorInstanceId,
        dt_seconds: f64,
    ) -> ChemistryResult<ReactorTickReport> {
        validate_positive_finite("dt_seconds", dt_seconds)?;
        let runtime = self.runtime_mut(reactor_id)?;
        runtime.reactor.tick(registry, dt_seconds)?;
        Ok(ReactorTickReport {
            reactor_id,
            dt_seconds,
            zone_count: runtime.reactor.zone_count(),
            total_electrical_draw_w: runtime.reactor.total_electrical_draw_w(),
        })
    }

    pub fn tick_all(
        &mut self,
        registry: &ChemistryRegistry,
        dt_seconds: f64,
    ) -> ChemistryResult<Vec<ReactorTickReport>> {
        validate_positive_finite("dt_seconds", dt_seconds)?;
        let reactor_ids: Vec<_> = self.reactor_ids().collect();
        let mut reports = Vec::with_capacity(reactor_ids.len());
        for reactor_id in reactor_ids {
            reports.push(self.tick_reactor(registry, reactor_id, dt_seconds)?);
        }
        Ok(reports)
    }

    pub fn insert_item_stack_to_input(
        &mut self,
        registry: &ChemistryRegistry,
        minecraft_registry: &MinecraftChemicalRegistry,
        reactor_id: ReactorInstanceId,
        input_index: usize,
        item_id: &str,
        item_count: u32,
    ) -> ChemistryResult<ItemInputReport> {
        if item_count == 0 {
            return Err(ChemistryError::InvalidMixtureState(
                "item_count must be positive".to_string(),
            ));
        }
        let mapping = minecraft_registry
            .lookup_by_item(item_id)
            .ok_or_else(|| unknown_item_error(item_id))?;
        let mol_inserted = mapping.mol_per_item * f64::from(item_count);
        validate_positive_finite("mol_inserted", mol_inserted)?;

        let runtime = self.runtime_mut(reactor_id)?;
        let input = runtime
            .reactor
            .input(input_index)
            .ok_or_else(|| input_index_error(reactor_id, input_index))?
            .clone();
        if !input.enabled {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "input {input_index} on reactor {} is disabled",
                reactor_id.0
            )));
        }
        if !input_accepts_substance(&input.mode, &mapping.substance_id) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "input {input_index} on reactor {} does not accept substance '{}'",
                reactor_id.0, mapping.substance_id
            )));
        }
        let zone = runtime
            .reactor
            .zone_mut(&input.to)
            .ok_or_else(|| zone_index_error(reactor_id, input.to.0))?;
        io::insert_substance(zone, registry, &mapping.substance_id, mol_inserted)?;

        Ok(ItemInputReport {
            reactor_id,
            input_index,
            item_id: MinecraftId::from(item_id),
            substance_id: mapping.substance_id.clone(),
            items_consumed: item_count,
            mol_inserted,
        })
    }

    pub fn drain_output_to_buffer(
        &mut self,
        registry: &ChemistryRegistry,
        reactor_id: ReactorInstanceId,
        output_index: usize,
        dt_seconds: f64,
    ) -> ChemistryResult<OutputDrainReport> {
        validate_positive_finite("dt_seconds", dt_seconds)?;
        let runtime = self.runtime_mut(reactor_id)?;
        if runtime.reactor.output(output_index).is_none() {
            return Err(output_index_error(reactor_id, output_index));
        }
        let extracted = runtime
            .reactor
            .apply_output(registry, output_index, dt_seconds)?
            .into_iter()
            .filter(|(_, amount)| *amount > 0.0)
            .map(|(substance_id, mol_per_bucket)| {
                let buffer = runtime
                    .output_buffers_mol_per_bucket
                    .entry(output_index)
                    .or_default();
                *buffer.entry(substance_id.clone()).or_insert(0.0) += mol_per_bucket;
                OutputSubstanceAmount {
                    substance_id,
                    mol_per_bucket,
                }
            })
            .collect();
        Ok(OutputDrainReport {
            reactor_id,
            output_index,
            extracted,
        })
    }

    pub fn take_output_items(
        &mut self,
        minecraft_registry: &MinecraftChemicalRegistry,
        reactor_id: ReactorInstanceId,
        output_index: usize,
        item_id: &str,
        max_items: u32,
    ) -> ChemistryResult<ItemOutputReport> {
        if max_items == 0 {
            return Err(ChemistryError::InvalidMixtureState(
                "max_items must be positive".to_string(),
            ));
        }
        let mapping = minecraft_registry
            .lookup_by_item(item_id)
            .ok_or_else(|| unknown_item_error(item_id))?;
        validate_positive_finite("mol_per_item", mapping.mol_per_item)?;

        let runtime = self.runtime_mut(reactor_id)?;
        if runtime.reactor.output(output_index).is_none() {
            return Err(output_index_error(reactor_id, output_index));
        }
        let buffer = runtime
            .output_buffers_mol_per_bucket
            .entry(output_index)
            .or_default();
        let available = *buffer.get(&mapping.substance_id).unwrap_or(&0.0);
        let possible_items = (available / mapping.mol_per_item).floor();
        let items_extracted = if possible_items >= f64::from(max_items) {
            max_items
        } else {
            possible_items as u32
        };
        let mol_consumed = mapping.mol_per_item * f64::from(items_extracted);
        if mol_consumed > 0.0 {
            let remaining = available - mol_consumed;
            if remaining > 1.0e-12 {
                buffer.insert(mapping.substance_id.clone(), remaining);
            } else {
                buffer.remove(&mapping.substance_id);
            }
        }
        Ok(ItemOutputReport {
            reactor_id,
            output_index,
            item_id: MinecraftId::from(item_id),
            substance_id: mapping.substance_id.clone(),
            items_extracted,
            mol_consumed,
            mol_remaining_in_buffer: buffer.get(&mapping.substance_id).copied().unwrap_or(0.0),
        })
    }

    pub fn buffered_mol(
        &self,
        reactor_id: ReactorInstanceId,
        output_index: usize,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<f64> {
        let runtime = self.runtime(reactor_id)?;
        Ok(runtime
            .output_buffers_mol_per_bucket
            .get(&output_index)
            .and_then(|buffer| buffer.get(substance_id))
            .copied()
            .unwrap_or(0.0))
    }

    pub fn export_reactor_checkpoint(
        &self,
        registry: &ChemistryRegistry,
        reactor_id: ReactorInstanceId,
        catalog_version: &str,
        content_version: u64,
        limits: &NativeBlobLimits,
    ) -> ChemistryResult<Vec<u8>> {
        let runtime = self.runtime(reactor_id)?;
        export_reactor_checkpoint(
            registry,
            &runtime.reactor,
            catalog_version,
            content_version,
            output_buffer_snapshots(&runtime.output_buffers_mol_per_bucket),
            limits,
        )
    }

    fn runtime(&self, reactor_id: ReactorInstanceId) -> ChemistryResult<&ReactorRuntime> {
        self.reactors
            .get(&reactor_id)
            .ok_or_else(|| unknown_reactor_error(reactor_id))
    }

    fn runtime_mut(
        &mut self,
        reactor_id: ReactorInstanceId,
    ) -> ChemistryResult<&mut ReactorRuntime> {
        self.reactors
            .get_mut(&reactor_id)
            .ok_or_else(|| unknown_reactor_error(reactor_id))
    }

    fn register_runtime(&mut self, id: ReactorInstanceId, runtime: ReactorRuntime) {
        let previous = self.reactors.insert(id, runtime);
        debug_assert!(
            previous.is_none(),
            "reactor instance id {} was allocated twice",
            id.0
        );
    }
}

fn output_buffer_snapshots(
    buffers: &BTreeMap<usize, BTreeMap<SubstanceId, f64>>,
) -> Vec<OutputBufferSnapshot> {
    buffers
        .iter()
        .map(|(output_index, substances)| OutputBufferSnapshot {
            output_index: *output_index,
            substances: substances
                .iter()
                .map(
                    |(substance_id, mol_per_bucket)| OutputBufferSubstanceSnapshot {
                        substance_id: substance_id.clone(),
                        mol_per_bucket: *mol_per_bucket,
                    },
                )
                .collect(),
        })
        .collect()
}

fn output_buffers_from_snapshots(
    registry: &ChemistryRegistry,
    snapshots: &[OutputBufferSnapshot],
) -> ChemistryResult<BTreeMap<usize, BTreeMap<SubstanceId, f64>>> {
    let mut buffers = BTreeMap::new();
    for snapshot in snapshots {
        if buffers.contains_key(&snapshot.output_index) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "reactor checkpoint contains duplicate output buffer {}",
                snapshot.output_index
            )));
        }
        let buffer = buffers
            .entry(snapshot.output_index)
            .or_insert_with(BTreeMap::new);
        for substance in &snapshot.substances {
            registry.substance(&substance.substance_id)?;
            validate_non_negative_finite("output buffer amount", substance.mol_per_bucket)?;
            let previous = buffer.insert(substance.substance_id.clone(), substance.mol_per_bucket);
            if previous.is_some() {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "output buffer {} contains duplicate substance '{}'",
                    snapshot.output_index, substance.substance_id
                )));
            }
        }
    }
    Ok(buffers)
}

fn input_accepts_substance(mode: &TransitionMode, substance_id: &SubstanceId) -> bool {
    match mode {
        TransitionMode::Substances { entries }
        | TransitionMode::SubstancesThreshold { entries, .. } => {
            entries.iter().any(|entry| &entry.id == substance_id)
        }
        TransitionMode::All { .. } => true,
        TransitionMode::Phases { .. } => false,
    }
}

fn validate_positive_finite(name: &str, value: f64) -> ChemistryResult<()> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(ChemistryError::InvalidMixtureState(format!(
            "{name} must be positive and finite, got {value}"
        )))
    }
}

fn validate_non_negative_finite(name: &str, value: f64) -> ChemistryResult<()> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(ChemistryError::InvalidMixtureState(format!(
            "{name} must be non-negative and finite, got {value}"
        )))
    }
}

fn unknown_reactor_error(reactor_id: ReactorInstanceId) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!("unknown reactor instance {}", reactor_id.0))
}

fn unknown_item_error(item_id: &str) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!("unknown minecraft chemical item '{item_id}'"))
}

fn input_index_error(reactor_id: ReactorInstanceId, input_index: usize) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!(
        "input index {input_index} out of range for reactor {}",
        reactor_id.0
    ))
}

fn output_index_error(reactor_id: ReactorInstanceId, output_index: usize) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!(
        "output index {output_index} out of range for reactor {}",
        reactor_id.0
    ))
}

fn zone_index_error(reactor_id: ReactorInstanceId, zone_index: usize) -> ChemistryError {
    ChemistryError::InvalidMixtureState(format!(
        "zone index {zone_index} out of range for reactor {}",
        reactor_id.0
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::catalog;
    use crate::chemistry::minecraft::mregistry::item_to_substance::MinecraftId;
    use crate::chemistry::reactor::{
        Input, Output, ReactorZone, SubstanceEntry, TransitionMode, ZoneId,
    };

    fn registry() -> ChemistryRegistry {
        catalog::destroy_substances_registry_builder()
            .unwrap()
            .build()
            .unwrap()
    }

    fn minecraft_registry(registry: &ChemistryRegistry) -> MinecraftChemicalRegistry {
        let mut minecraft_registry = MinecraftChemicalRegistry::new();
        minecraft_registry
            .register(
                MinecraftId::from("minecraft:water_bucket"),
                water_id(),
                1.0,
                registry,
            )
            .unwrap();
        minecraft_registry
            .register(
                MinecraftId::from("minecraft:small_water_cell"),
                water_id(),
                0.25,
                registry,
            )
            .unwrap();
        minecraft_registry
    }

    fn water_id() -> SubstanceId {
        SubstanceId::from("destroy:water")
    }

    fn reactor_with_water_input_output() -> Reactor {
        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(ReactorZone::new(0.001).unwrap());
        reactor.add_input(Input::new(
            zone,
            TransitionMode::Substances {
                entries: vec![SubstanceEntry {
                    id: water_id(),
                    rate_mol_per_second: 1.0,
                }],
            },
        ));
        reactor.add_output(Output::new(
            zone,
            TransitionMode::Substances {
                entries: vec![SubstanceEntry {
                    id: water_id(),
                    rate_mol_per_second: 0.6,
                }],
            },
        ));
        reactor
    }

    #[test]
    fn worker_stores_reactors_by_stable_ids() {
        let mut worker = ReactorsWorker::new();
        let first = worker.register_reactor(Reactor::new());
        let second = worker.register_reactor(Reactor::new());

        assert_eq!(first, ReactorInstanceId(0));
        assert_eq!(second, ReactorInstanceId(1));
        assert_eq!(worker.reactor_count(), 2);

        worker.remove_reactor(first).unwrap();
        assert_eq!(worker.reactor_count(), 1);
        assert!(worker.reactor(first).is_err());
        assert!(worker.reactor(second).is_ok());
    }

    #[test]
    fn item_input_inserts_mapped_substance_into_reactor_zone() {
        let registry = registry();
        let minecraft_registry = minecraft_registry(&registry);
        let mut worker = ReactorsWorker::new();
        let reactor_id = worker.register_reactor(reactor_with_water_input_output());

        let report = worker
            .insert_item_stack_to_input(
                &registry,
                &minecraft_registry,
                reactor_id,
                0,
                "minecraft:water_bucket",
                2,
            )
            .unwrap();

        assert_eq!(report.items_consumed, 2);
        assert_eq!(report.substance_id, water_id());
        assert_eq!(report.mol_inserted, 2.0);
        let concentration = worker
            .reactor(reactor_id)
            .unwrap()
            .zone(&ZoneId(0))
            .unwrap()
            .concentration_of(&water_id());
        assert!((concentration - 2.0).abs() < 1.0e-9);
    }

    #[test]
    fn output_drain_buffers_fractional_items_without_loss() {
        let registry = registry();
        let minecraft_registry = minecraft_registry(&registry);
        let mut worker = ReactorsWorker::new();
        let reactor_id = worker.register_reactor(reactor_with_water_input_output());
        worker
            .insert_item_stack_to_input(
                &registry,
                &minecraft_registry,
                reactor_id,
                0,
                "minecraft:water_bucket",
                1,
            )
            .unwrap();

        let drain = worker
            .drain_output_to_buffer(&registry, reactor_id, 0, 1.0)
            .unwrap();
        assert_eq!(drain.extracted.len(), 1);
        assert!((worker.buffered_mol(reactor_id, 0, &water_id()).unwrap() - 0.6).abs() < 1.0e-9);

        let whole_bucket = worker
            .take_output_items(
                &minecraft_registry,
                reactor_id,
                0,
                "minecraft:water_bucket",
                1,
            )
            .unwrap();
        assert_eq!(whole_bucket.items_extracted, 0);
        assert!((whole_bucket.mol_remaining_in_buffer - 0.6).abs() < 1.0e-9);

        let small_cells = worker
            .take_output_items(
                &minecraft_registry,
                reactor_id,
                0,
                "minecraft:small_water_cell",
                10,
            )
            .unwrap();
        assert_eq!(small_cells.items_extracted, 2);
        assert!((small_cells.mol_consumed - 0.5).abs() < 1.0e-9);
        assert!((small_cells.mol_remaining_in_buffer - 0.1).abs() < 1.0e-9);
    }

    #[test]
    fn reactor_checkpoint_restores_reactor_and_output_buffers() {
        let registry = registry();
        let minecraft_registry = minecraft_registry(&registry);
        let limits = NativeBlobLimits::default();
        let mut worker = ReactorsWorker::new();
        let reactor_id = worker.register_reactor(reactor_with_water_input_output());
        worker
            .insert_item_stack_to_input(
                &registry,
                &minecraft_registry,
                reactor_id,
                0,
                "minecraft:water_bucket",
                1,
            )
            .unwrap();
        worker
            .drain_output_to_buffer(&registry, reactor_id, 0, 0.5)
            .unwrap();

        let encoded = worker
            .export_reactor_checkpoint(&registry, reactor_id, "test-catalog", 7, &limits)
            .unwrap();
        let restored_id = worker
            .register_reactor_from_checkpoint(&registry, &encoded, "test-catalog", &limits)
            .unwrap();

        let concentration = worker
            .reactor(restored_id)
            .unwrap()
            .zone(&ZoneId(0))
            .unwrap()
            .concentration_of(&water_id());
        assert!((concentration - 0.7).abs() < 1.0e-9);
        assert!((worker.buffered_mol(restored_id, 0, &water_id()).unwrap() - 0.3).abs() < 1.0e-9);
    }

    #[test]
    fn tick_all_ticks_registered_reactors() {
        let registry = registry();
        let mut worker = ReactorsWorker::new();
        let first = worker.register_reactor(Reactor::new());
        let second = worker.register_reactor(Reactor::new());

        let reports = worker.tick_all(&registry, 0.5).unwrap();

        assert_eq!(reports.len(), 2);
        assert_eq!(reports[0].reactor_id, first);
        assert_eq!(reports[1].reactor_id, second);
        assert_eq!(reports[0].dt_seconds, 0.5);
    }

    #[test]
    fn input_rejects_unmapped_item() {
        let registry = registry();
        let minecraft_registry = minecraft_registry(&registry);
        let mut worker = ReactorsWorker::new();
        let reactor_id = worker.register_reactor(reactor_with_water_input_output());

        let error = worker
            .insert_item_stack_to_input(
                &registry,
                &minecraft_registry,
                reactor_id,
                0,
                "minecraft:stone",
                1,
            )
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidMixtureState(_)));
    }
}
