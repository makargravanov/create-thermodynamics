use std::sync::{Mutex, OnceLock};

use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::minecraft::mregistry::item_to_substance::MinecraftId;
use crate::chemistry::minecraft::mregistry::mregistry::{
    MinecraftChemicalRegistry, RegistrationError,
};
use crate::chemistry::minecraft::worker::reactors_worker::{ReactorInstanceId, ReactorsWorker};
use crate::chemistry::reactor::{Input, Output, Reactor, ReactorZone, TransitionMode};
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::SubstanceId;

const DEFAULT_MULTIBLOCK_OUTPUT_RATE_MOL_PER_SECOND: f64 = 1.0;

#[derive(Debug, Clone, PartialEq)]
pub struct ItemChemicalBinding {
    pub item_id: String,
    pub substance_id: SubstanceId,
    pub mol_per_item: f64,
}

impl ItemChemicalBinding {
    pub fn new(
        item_id: impl Into<String>,
        substance_id: impl Into<SubstanceId>,
        mol_per_item: f64,
    ) -> Self {
        Self {
            item_id: item_id.into(),
            substance_id: substance_id.into(),
            mol_per_item,
        }
    }
}

pub struct MinecraftChemistryState {
    chemistry_registry: ChemistryRegistry,
    item_bindings: MinecraftChemicalRegistry,
    reactors_worker: ReactorsWorker,
}

impl MinecraftChemistryState {
    fn new() -> ChemistryResult<Self> {
        Ok(Self {
            chemistry_registry: crate::chemistry::destroy_registry_builder()?.build()?,
            item_bindings: MinecraftChemicalRegistry::new(),
            reactors_worker: ReactorsWorker::new(),
        })
    }

    pub fn chemistry_registry(&self) -> &ChemistryRegistry {
        &self.chemistry_registry
    }

    pub fn item_bindings(&self) -> &MinecraftChemicalRegistry {
        &self.item_bindings
    }

    pub fn reactors_worker(&self) -> &ReactorsWorker {
        &self.reactors_worker
    }
}

static MINECRAFT_CHEMISTRY_STATE: OnceLock<ChemistryResult<Mutex<MinecraftChemistryState>>> =
    OnceLock::new();

pub fn with_minecraft_chemistry_state<T>(
    f: impl FnOnce(&MinecraftChemistryState) -> ChemistryResult<T>,
) -> ChemistryResult<T> {
    let state = minecraft_chemistry_state()?;
    let guard = state.lock().map_err(|_| {
        ChemistryError::InvalidMixtureState(
            "minecraft chemistry state mutex is poisoned".to_string(),
        )
    })?;
    f(&guard)
}

pub fn with_minecraft_chemistry_state_mut<T>(
    f: impl FnOnce(&mut MinecraftChemistryState) -> ChemistryResult<T>,
) -> ChemistryResult<T> {
    let state = minecraft_chemistry_state()?;
    let mut guard = state.lock().map_err(|_| {
        ChemistryError::InvalidMixtureState(
            "minecraft chemistry state mutex is poisoned".to_string(),
        )
    })?;
    f(&mut guard)
}

pub fn replace_item_chemical_bindings(bindings: Vec<ItemChemicalBinding>) -> ChemistryResult<()> {
    with_minecraft_chemistry_state_mut(|state| {
        let mut next_registry = MinecraftChemicalRegistry::new();
        for binding in bindings {
            next_registry
                .register(
                    MinecraftId::new(binding.item_id),
                    binding.substance_id,
                    binding.mol_per_item,
                    &state.chemistry_registry,
                )
                .map_err(registration_error_to_chemistry_error)?;
        }
        state.item_bindings = next_registry;
        Ok(())
    })
}

pub fn clear_item_chemical_bindings() -> ChemistryResult<()> {
    replace_item_chemical_bindings(Vec::new())
}

pub fn item_chemical_binding_count() -> ChemistryResult<usize> {
    with_minecraft_chemistry_state(|state| Ok(state.item_bindings.item_count()))
}

pub fn has_item_chemical_binding(item_id: &str) -> ChemistryResult<bool> {
    with_minecraft_chemistry_state(|state| Ok(state.item_bindings.contains_item(item_id)))
}

pub fn static_substance_ids() -> ChemistryResult<Vec<String>> {
    with_minecraft_chemistry_state(|state| {
        Ok(state
            .chemistry_registry
            .substances()
            .map(|substance| substance.id.as_str().to_string())
            .collect())
    })
}

pub fn register_single_zone_reactor(
    volume_cubic_meters: f64,
    item_input_count: usize,
    item_output_count: usize,
    fluid_input_count: usize,
    fluid_output_count: usize,
) -> ChemistryResult<ReactorInstanceId> {
    if !volume_cubic_meters.is_finite() || volume_cubic_meters <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "reactor volume must be positive and finite, got {volume_cubic_meters}"
        )));
    }
    with_minecraft_chemistry_state_mut(|state| {
        let mut reactor = Reactor::new();
        let zone = reactor.add_zone(ReactorZone::new(volume_cubic_meters)?);
        for _ in 0..item_input_count + fluid_input_count {
            reactor.add_input(Input::new(
                zone,
                TransitionMode::All {
                    rate_mol_per_second: f64::MAX,
                },
            ));
        }
        for _ in 0..item_output_count + fluid_output_count {
            reactor.add_output(Output::new(
                zone,
                TransitionMode::All {
                    rate_mol_per_second: DEFAULT_MULTIBLOCK_OUTPUT_RATE_MOL_PER_SECOND,
                },
            ));
        }
        Ok(state.reactors_worker.register_reactor(reactor))
    })
}

pub fn remove_reactor(reactor_id: ReactorInstanceId) -> ChemistryResult<()> {
    with_minecraft_chemistry_state_mut(|state| {
        state.reactors_worker.remove_reactor(reactor_id)?;
        Ok(())
    })
}

pub fn reactor_count() -> ChemistryResult<usize> {
    with_minecraft_chemistry_state(|state| Ok(state.reactors_worker.reactor_count()))
}

pub fn tick_reactor(reactor_id: ReactorInstanceId, dt_seconds: f64) -> ChemistryResult<()> {
    with_minecraft_chemistry_state_mut(|state| {
        state
            .reactors_worker
            .tick_reactor(&state.chemistry_registry, reactor_id, dt_seconds)?;
        Ok(())
    })
}

pub fn insert_item_stack_to_reactor_input(
    reactor_id: ReactorInstanceId,
    input_index: usize,
    item_id: &str,
    item_count: u32,
) -> ChemistryResult<f64> {
    with_minecraft_chemistry_state_mut(|state| {
        let report = state.reactors_worker.insert_item_stack_to_input(
            &state.chemistry_registry,
            &state.item_bindings,
            reactor_id,
            input_index,
            item_id,
            item_count,
        )?;
        Ok(report.mol_inserted)
    })
}

fn minecraft_chemistry_state() -> ChemistryResult<&'static Mutex<MinecraftChemistryState>> {
    MINECRAFT_CHEMISTRY_STATE
        .get_or_init(|| MinecraftChemistryState::new().map(Mutex::new))
        .as_ref()
        .map_err(Clone::clone)
}

fn registration_error_to_chemistry_error(error: RegistrationError) -> ChemistryError {
    match error {
        RegistrationError::DuplicateItem {
            item_id,
            existing_substance_id,
            new_substance_id,
        } => ChemistryError::InvalidMixtureState(format!(
            "minecraft item '{}' is already bound to '{}' and cannot also bind to '{}'",
            item_id.as_str(),
            existing_substance_id,
            new_substance_id
        )),
        RegistrationError::UnknownSubstance { substance_id } => {
            ChemistryError::InvalidMixtureState(format!(
                "minecraft item binding refers to unknown substance '{}'",
                substance_id
            ))
        }
        RegistrationError::InvalidAmount {
            item_id,
            mol_per_item,
        } => ChemistryError::InvalidMixtureState(format!(
            "minecraft item '{}' has invalid mol_per_item {}",
            item_id.as_str(),
            mol_per_item
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn item_bindings_are_replaced_atomically() {
        replace_item_chemical_bindings(vec![ItemChemicalBinding::new(
            "minecraft:water_bucket",
            "destroy:water",
            1.0,
        )])
        .unwrap();
        assert_eq!(item_chemical_binding_count().unwrap(), 1);
        assert!(has_item_chemical_binding("minecraft:water_bucket").unwrap());

        let error = replace_item_chemical_bindings(vec![
            ItemChemicalBinding::new("minecraft:ethanol_bucket", "destroy:ethanol", 1.0),
            ItemChemicalBinding::new("minecraft:bad_bucket", "destroy:not_real", 1.0),
        ])
        .unwrap_err();
        assert!(matches!(error, ChemistryError::InvalidMixtureState(_)));

        assert_eq!(item_chemical_binding_count().unwrap(), 1);
        assert!(has_item_chemical_binding("minecraft:water_bucket").unwrap());
        assert!(!has_item_chemical_binding("minecraft:ethanol_bucket").unwrap());
    }

    #[test]
    fn clearing_item_bindings_removes_all_entries() {
        replace_item_chemical_bindings(vec![ItemChemicalBinding::new(
            "minecraft:water_bucket",
            "destroy:water",
            1.0,
        )])
        .unwrap();
        clear_item_chemical_bindings().unwrap();

        assert_eq!(item_chemical_binding_count().unwrap(), 0);
        assert!(!has_item_chemical_binding("minecraft:water_bucket").unwrap());
    }

    #[test]
    fn static_substance_ids_expose_destroy_catalog() {
        let ids = static_substance_ids().unwrap();

        assert!(ids.contains(&"destroy:water".to_string()));
        assert!(ids.contains(&"destroy:ethanol".to_string()));
        assert!(ids.len() >= 152);
    }
}
