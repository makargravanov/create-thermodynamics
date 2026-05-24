use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::registry::{ChemistryRegistry, SubstanceIndex};
use super::substance::{SubstanceId, SubstanceTagId};

pub const DEFAULT_TEMPERATURE_KELVIN: f64 = 298.0;
pub const TRACE_CONCENTRATION_MOL_PER_BUCKET: f64 = 1.0 / 512.0 / 512.0;
const THERMAL_STATE_EPSILON: f64 = 1.0e-12;

#[derive(Debug, Clone)]
pub struct Mixture {
    temperature_kelvin: f64,
    components: Vec<MixtureComponent>,
    positions_by_substance: Vec<Option<usize>>,
}

#[derive(Debug, Clone)]
struct MixtureComponent {
    substance: SubstanceIndex,
    substance_id: SubstanceId,
    concentration_mol_per_bucket: f64,
    gaseous_fraction: f64,
}

#[derive(Debug, Clone)]
pub(crate) struct MixtureCheckpoint {
    temperature_kelvin: f64,
    components: Vec<ComponentCheckpoint>,
}

#[derive(Debug, Clone)]
struct ComponentCheckpoint {
    substance: SubstanceIndex,
    previous: Option<MixtureComponent>,
}

impl Mixture {
    pub fn new(temperature_kelvin: f64) -> ChemistryResult<Self> {
        if !temperature_kelvin.is_finite() || temperature_kelvin < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "temperature must be non-negative and finite".to_string(),
            ));
        }
        Ok(Self {
            temperature_kelvin,
            components: Vec::new(),
            positions_by_substance: Vec::new(),
        })
    }

    pub fn empty() -> Self {
        Self::new(DEFAULT_TEMPERATURE_KELVIN).expect("default temperature must be valid")
    }

    pub fn temperature_kelvin(&self) -> f64 {
        self.temperature_kelvin
    }

    pub fn concentration_of(&self, substance_id: &SubstanceId) -> f64 {
        self.components
            .iter()
            .find(|component| &component.substance_id == substance_id)
            .map(|component| component.concentration_mol_per_bucket)
            .unwrap_or(0.0)
    }

    pub fn gaseous_fraction_of(&self, substance_id: &SubstanceId) -> f64 {
        self.components
            .iter()
            .find(|component| &component.substance_id == substance_id)
            .map(|component| component.gaseous_fraction)
            .unwrap_or(0.0)
    }

    pub fn substances(&self) -> impl Iterator<Item = &SubstanceId> {
        self.components
            .iter()
            .map(|component| &component.substance_id)
    }

    pub(crate) fn component_indices(&self) -> impl Iterator<Item = SubstanceIndex> + '_ {
        self.components.iter().map(|component| component.substance)
    }

    pub(crate) fn concentration_of_index(&self, substance: SubstanceIndex) -> f64 {
        self.position_of_substance(substance)
            .map(|position| self.components[position].concentration_mol_per_bucket)
            .unwrap_or(0.0)
    }

    pub(crate) fn gaseous_fraction_of_index(&self, substance: SubstanceIndex) -> f64 {
        self.position_of_substance(substance)
            .map(|position| self.components[position].gaseous_fraction)
            .unwrap_or(0.0)
    }

    pub fn add_substance(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: impl Into<SubstanceId>,
        concentration_mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        let substance_id = substance_id.into();
        let substance_index = registry.substance_index(&substance_id).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!("unknown substance '{substance_id}'"))
        })?;
        let substance = registry.substance(&substance_id)?;
        if substance
            .tags
            .iter()
            .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
        {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "hypothetical substance '{substance_id}' cannot be added to a mixture"
            )));
        }
        validate_concentration(concentration_mol_per_bucket)?;
        if concentration_mol_per_bucket == 0.0 {
            return Ok(());
        }
        self.ensure_position_capacity(registry);
        let initial_gas = if substance.boiling_point_kelvin < self.temperature_kelvin {
            1.0
        } else {
            0.0
        };
        self.change_concentration_by_index_unchecked(
            registry,
            substance_index,
            substance_id,
            concentration_mol_per_bucket,
            initial_gas,
        )?;
        self.validate(registry)
    }

    pub fn set_gaseous_fraction(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: impl Into<SubstanceId>,
        gaseous_fraction: f64,
    ) -> ChemistryResult<()> {
        let substance_id = substance_id.into();
        let substance_index = registry.substance_index(&substance_id).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!("unknown substance '{substance_id}'"))
        })?;
        validate_gaseous_fraction(gaseous_fraction)?;
        if let Some(position) = self.position_of_substance(substance_index) {
            self.components[position].gaseous_fraction = gaseous_fraction;
        }
        self.validate(registry)
    }

    pub fn change_concentration(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
        delta_mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        let substance_index = registry.substance_index(substance_id).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!("unknown substance '{substance_id}'"))
        })?;
        self.change_concentration_by_index(registry, substance_index, delta_mol_per_bucket)
    }

    pub(crate) fn change_concentration_by_index(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        delta_mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        let substance_id = registry.substance_by_index(substance)?.id.clone();
        if !delta_mol_per_bucket.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(
                "concentration change must be finite".to_string(),
            ));
        }
        self.ensure_position_capacity(registry);
        let current = self.concentration_of_index(substance);
        let next = current + delta_mol_per_bucket;
        if next < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "substance '{}' would become negative: {next}",
                substance_id
            )));
        }
        if next <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            self.remove_component(substance);
        } else {
            let initial_gas = if registry.substance_by_index(substance)?.boiling_point_kelvin
                < self.temperature_kelvin
            {
                1.0
            } else {
                0.0
            };
            self.change_concentration_by_index_unchecked(
                registry,
                substance,
                substance_id,
                delta_mol_per_bucket,
                initial_gas,
            )?;
        }
        self.validate(registry)
    }

    pub(crate) fn apply_concentration_deltas_by_index(
        &mut self,
        registry: &ChemistryRegistry,
        deltas: &[(SubstanceIndex, f64)],
    ) -> ChemistryResult<f64> {
        self.ensure_position_capacity(registry);
        let mut next = Vec::new();
        let mut max_delta = 0.0_f64;
        for (substance, delta) in deltas {
            if !delta.is_finite() {
                return Err(ChemistryError::InvalidMixtureState(
                    "concentration change must be finite".to_string(),
                ));
            }
            max_delta = max_delta.max(delta.abs());
            let current = self.concentration_of_index(*substance);
            let value = current + delta;
            if value < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
                let substance_id = &registry.substance_by_index(*substance)?.id;
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "substance '{substance_id}' would become negative: {value}"
                )));
            }
            next.push((*substance, value));
        }
        for (substance, value) in next {
            if value <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                self.remove_component(substance);
            } else if let Some(position) = self.position_of_substance(substance) {
                self.components[position].concentration_mol_per_bucket = value;
            } else {
                let substance_data = registry.substance_by_index(substance)?;
                let initial_gas = if substance_data.boiling_point_kelvin < self.temperature_kelvin {
                    1.0
                } else {
                    0.0
                };
                self.insert_component(substance, substance_data.id.clone(), value, initial_gas);
            }
        }
        self.validate(registry)?;
        Ok(max_delta)
    }

    pub fn heat(
        &mut self,
        registry: &ChemistryRegistry,
        mut energy_j_per_bucket: f64,
    ) -> ChemistryResult<()> {
        if !energy_j_per_bucket.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(
                "heat energy must be finite".to_string(),
            ));
        }
        if energy_j_per_bucket == 0.0 || self.components.is_empty() {
            return self.validate(registry);
        }

        let mut guard = 0usize;
        while energy_j_per_bucket.abs() > 1.0e-12 {
            guard += 1;
            if guard > 10_000 {
                return Err(ChemistryError::InvalidMixtureState(
                    "heating did not converge".to_string(),
                ));
            }

            let heat_capacity = self.volumetric_heat_capacity_j_per_bucket_kelvin(registry)?;
            if heat_capacity == 0.0 {
                return self.validate(registry);
            }
            let temperature_change = energy_j_per_bucket / heat_capacity;
            if temperature_change > 0.0 {
                if let Some((substance_index, boiling_point)) =
                    self.next_higher_boiling_point(registry)?
                {
                    if self.temperature_kelvin + temperature_change >= boiling_point {
                        let energy_to_boiling =
                            (boiling_point - self.temperature_kelvin) * heat_capacity;
                        self.temperature_kelvin = boiling_point;
                        energy_j_per_bucket -= energy_to_boiling;

                        let concentration = self.concentration_of_index(substance_index);
                        let current_gas = self.gaseous_fraction_of_index(substance_index);
                        let liquid_concentration = concentration * (1.0 - current_gas);
                        let latent_heat = registry.substance_properties().latent_heat_j_per_mol
                            [substance_index.as_usize()];
                        let energy_to_fully_boil = liquid_concentration * latent_heat;
                        if energy_to_fully_boil <= 0.0 {
                            self.set_gaseous_fraction_by_index(substance_index, 1.0);
                            continue;
                        }
                        if energy_j_per_bucket >= energy_to_fully_boil {
                            self.set_gaseous_fraction_by_index(substance_index, 1.0);
                            energy_j_per_bucket -= energy_to_fully_boil;
                        } else {
                            let boiled_fraction =
                                energy_j_per_bucket / (latent_heat * concentration);
                            self.set_gaseous_fraction_by_index(
                                substance_index,
                                current_gas + boiled_fraction,
                            );
                            energy_j_per_bucket = 0.0;
                        }
                        continue;
                    }
                }
                self.temperature_kelvin += temperature_change;
                energy_j_per_bucket = 0.0;
            } else {
                if let Some((substance_index, boiling_point)) =
                    self.next_lower_boiling_point(registry)?
                {
                    if self.temperature_kelvin + temperature_change < boiling_point {
                        let energy_to_boiling =
                            (boiling_point - self.temperature_kelvin) * heat_capacity;
                        self.temperature_kelvin = boiling_point;
                        energy_j_per_bucket -= energy_to_boiling;

                        let concentration = self.concentration_of_index(substance_index);
                        let current_gas = self.gaseous_fraction_of_index(substance_index);
                        let gas_concentration = concentration * current_gas;
                        let latent_heat = registry.substance_properties().latent_heat_j_per_mol
                            [substance_index.as_usize()];
                        let released_by_full_condensation = gas_concentration * latent_heat;
                        if released_by_full_condensation <= 0.0 {
                            self.set_gaseous_fraction_by_index(substance_index, 0.0);
                            continue;
                        }
                        if -energy_j_per_bucket >= released_by_full_condensation {
                            self.set_gaseous_fraction_by_index(substance_index, 0.0);
                            energy_j_per_bucket += released_by_full_condensation;
                        } else {
                            let condensed_fraction =
                                -energy_j_per_bucket / (latent_heat * concentration);
                            self.set_gaseous_fraction_by_index(
                                substance_index,
                                current_gas - condensed_fraction,
                            );
                            energy_j_per_bucket = 0.0;
                        }
                        continue;
                    }
                }
                self.temperature_kelvin =
                    (self.temperature_kelvin + temperature_change).max(0.0001);
                energy_j_per_bucket = 0.0;
            }
        }
        self.validate(registry)
    }

    pub(crate) fn checkpoint_for_reaction(
        &self,
        deltas: &[(SubstanceIndex, f64)],
    ) -> MixtureCheckpoint {
        let mut substances = self
            .components
            .iter()
            .map(|component| component.substance)
            .collect::<Vec<_>>();
        for (substance, _) in deltas {
            insert_sorted_unique(&mut substances, *substance);
        }
        MixtureCheckpoint {
            temperature_kelvin: self.temperature_kelvin,
            components: substances
                .into_iter()
                .map(|substance| ComponentCheckpoint {
                    substance,
                    previous: self
                        .position_of_substance(substance)
                        .map(|position| self.components[position].clone()),
                })
                .collect(),
        }
    }

    pub(crate) fn restore_checkpoint(&mut self, checkpoint: MixtureCheckpoint) {
        self.temperature_kelvin = checkpoint.temperature_kelvin;
        for component in checkpoint.components {
            match component.previous {
                Some(previous) => {
                    if let Some(position) = self.position_of_substance(component.substance) {
                        self.components[position] = previous;
                    } else {
                        self.components.push(previous);
                    }
                }
                None => self.remove_component(component.substance),
            }
        }
        self.rebuild_positions();
    }

    pub(crate) fn changed_since_checkpoint(&self, checkpoint: &MixtureCheckpoint) -> bool {
        if (self.temperature_kelvin - checkpoint.temperature_kelvin).abs() > THERMAL_STATE_EPSILON {
            return true;
        }
        checkpoint.components.iter().any(|component| {
            let current = self
                .position_of_substance(component.substance)
                .map(|position| &self.components[position]);
            match (&component.previous, current) {
                (Some(previous), Some(current)) => {
                    (current.concentration_mol_per_bucket - previous.concentration_mol_per_bucket)
                        .abs()
                        > TRACE_CONCENTRATION_MOL_PER_BUCKET
                        || (current.gaseous_fraction - previous.gaseous_fraction).abs()
                            > THERMAL_STATE_EPSILON
                }
                (None, Some(current)) => {
                    current.concentration_mol_per_bucket > TRACE_CONCENTRATION_MOL_PER_BUCKET
                }
                (Some(previous), None) => {
                    previous.concentration_mol_per_bucket > TRACE_CONCENTRATION_MOL_PER_BUCKET
                }
                (None, None) => false,
            }
        })
    }

    pub fn volumetric_heat_capacity_j_per_bucket_kelvin(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<f64> {
        self.validate_registry_shape(registry)?;
        let properties = registry.substance_properties();
        Ok(self.components.iter().fold(0.0, |acc, component| {
            acc + properties.molar_heat_capacity_j_per_mol_kelvin[component.substance.as_usize()]
                * component.concentration_mol_per_bucket
        }))
    }

    pub fn recalculate_volume_millibuckets(
        &self,
        registry: &ChemistryRegistry,
        initial_millibuckets: u32,
    ) -> ChemistryResult<u32> {
        self.validate(registry)?;
        let mut liquid_buckets = 0.0;
        let properties = registry.substance_properties();
        for component in &self.components {
            let liquid_fraction = 1.0 - component.gaseous_fraction;
            let index = component.substance.as_usize();
            liquid_buckets += component.concentration_mol_per_bucket
                * properties.molar_mass_grams[index]
                * liquid_fraction
                / properties.liquid_density_grams_per_bucket[index];
        }
        let calculated = (liquid_buckets * 1000.0).round();
        if calculated.is_finite() && calculated > 0.0 && calculated <= u32::MAX as f64 {
            Ok(calculated as u32)
        } else if calculated.is_finite() && calculated <= 0.0 {
            Ok(initial_millibuckets)
        } else {
            Err(ChemistryError::InvalidMixtureState(format!(
                "calculated volume must fit into u32 millibuckets: {calculated}"
            )))
        }
    }

    pub fn mix(
        registry: &ChemistryRegistry,
        mixtures: &[(Mixture, f64)],
    ) -> ChemistryResult<Mixture> {
        if mixtures.is_empty() {
            return Ok(Mixture::empty());
        }
        let total_amount = mixtures.iter().map(|(_, amount)| *amount).sum::<f64>();
        if !total_amount.is_finite() || total_amount <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "total mixed amount must be positive and finite".to_string(),
            ));
        }

        let mut moles_by_substance: BTreeMap<SubstanceId, f64> = BTreeMap::new();
        let mut total_energy = 0.0;
        for (mixture, amount) in mixtures {
            if !amount.is_finite() || *amount < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "mixed amount must be non-negative and finite".to_string(),
                ));
            }
            mixture.validate(registry)?;
            for component in &mixture.components {
                let moles = component.concentration_mol_per_bucket * amount;
                *moles_by_substance
                    .entry(component.substance_id.clone())
                    .or_insert(0.0) += moles;
                let substance = registry.substance_by_index(component.substance)?;
                total_energy += substance.molar_heat_capacity_j_per_mol_kelvin
                    * component.concentration_mol_per_bucket
                    * mixture.temperature_kelvin
                    * amount;
                total_energy += substance.latent_heat_j_per_mol
                    * component.concentration_mol_per_bucket
                    * component.gaseous_fraction
                    * amount;
            }
        }

        let mut result = Mixture::new(0.0)?;
        for (substance_id, moles) in moles_by_substance {
            result.add_substance(registry, substance_id.clone(), moles / total_amount)?;
            let substance_index = registry.substance_index(&substance_id).ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!("unknown substance '{substance_id}'"))
            })?;
            result.set_gaseous_fraction_by_index(substance_index, 0.0);
        }
        result.heat(registry, total_energy / total_amount)?;
        Ok(result)
    }

    pub fn validate(&self, registry: &ChemistryRegistry) -> ChemistryResult<()> {
        if !self.temperature_kelvin.is_finite() || self.temperature_kelvin < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "temperature must be non-negative and finite".to_string(),
            ));
        }
        self.validate_registry_shape(registry)?;
        let mut seen = BTreeSet::new();
        for (position, component) in self.components.iter().enumerate() {
            let substance = registry.substance_by_index(component.substance)?;
            if substance.id != component.substance_id {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "mixture component '{}' does not match registry index {}",
                    component.substance_id,
                    component.substance.as_usize()
                )));
            }
            if !seen.insert(component.substance) {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "duplicate mixture component '{}'",
                    component.substance_id
                )));
            }
            validate_concentration(component.concentration_mol_per_bucket)?;
            validate_gaseous_fraction(component.gaseous_fraction)?;
            if self
                .positions_by_substance
                .get(component.substance.as_usize())
                .copied()
                .flatten()
                != Some(position)
            {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "mixture position index is inconsistent for '{}'",
                    component.substance_id
                )));
            }
        }
        Ok(())
    }

    fn next_higher_boiling_point(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Option<(SubstanceIndex, f64)>> {
        let mut candidates = BTreeSet::new();
        let properties = registry.substance_properties();
        for component in &self.components {
            let boiling_point = properties.boiling_point_kelvin[component.substance.as_usize()];
            if boiling_point >= self.temperature_kelvin && component.gaseous_fraction < 1.0 {
                candidates.insert((ordered_f64(boiling_point), component.substance));
            }
        }
        Ok(candidates.into_iter().next().map(|(bp, id)| (id, bp.0)))
    }

    fn next_lower_boiling_point(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Option<(SubstanceIndex, f64)>> {
        let mut best: Option<(SubstanceIndex, f64)> = None;
        let properties = registry.substance_properties();
        for component in &self.components {
            let boiling_point = properties.boiling_point_kelvin[component.substance.as_usize()];
            if boiling_point <= self.temperature_kelvin
                && component.gaseous_fraction > 0.0
                && best
                    .as_ref()
                    .map(|(_, bp)| boiling_point > *bp)
                    .unwrap_or(true)
            {
                best = Some((component.substance, boiling_point));
            }
        }
        Ok(best)
    }

    fn position_of_substance(&self, substance: SubstanceIndex) -> Option<usize> {
        self.positions_by_substance
            .get(substance.as_usize())
            .copied()
            .flatten()
    }

    fn ensure_position_capacity(&mut self, registry: &ChemistryRegistry) {
        if self.positions_by_substance.len() < registry.substance_count() {
            self.positions_by_substance
                .resize(registry.substance_count(), None);
        }
    }

    fn insert_component(
        &mut self,
        substance: SubstanceIndex,
        substance_id: SubstanceId,
        concentration_mol_per_bucket: f64,
        gaseous_fraction: f64,
    ) {
        let position = self.components.len();
        self.components.push(MixtureComponent {
            substance,
            substance_id,
            concentration_mol_per_bucket,
            gaseous_fraction,
        });
        if self.positions_by_substance.len() <= substance.as_usize() {
            self.positions_by_substance
                .resize(substance.as_usize() + 1, None);
        }
        self.positions_by_substance[substance.as_usize()] = Some(position);
    }

    fn remove_component(&mut self, substance: SubstanceIndex) {
        let Some(position) = self.position_of_substance(substance) else {
            return;
        };
        self.components.remove(position);
        self.rebuild_positions();
    }

    fn rebuild_positions(&mut self) {
        self.positions_by_substance.fill(None);
        for (position, component) in self.components.iter().enumerate() {
            if self.positions_by_substance.len() <= component.substance.as_usize() {
                self.positions_by_substance
                    .resize(component.substance.as_usize() + 1, None);
            }
            self.positions_by_substance[component.substance.as_usize()] = Some(position);
        }
    }

    fn set_gaseous_fraction_by_index(&mut self, substance: SubstanceIndex, value: f64) {
        if let Some(position) = self.position_of_substance(substance) {
            self.components[position].gaseous_fraction = value;
        }
    }

    fn change_concentration_by_index_unchecked(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        substance_id: SubstanceId,
        delta_mol_per_bucket: f64,
        initial_gas: f64,
    ) -> ChemistryResult<()> {
        registry.substance_by_index(substance)?;
        if let Some(position) = self.position_of_substance(substance) {
            self.components[position].concentration_mol_per_bucket += delta_mol_per_bucket;
        } else {
            self.insert_component(substance, substance_id, delta_mol_per_bucket, initial_gas);
        }
        Ok(())
    }

    fn validate_registry_shape(&self, registry: &ChemistryRegistry) -> ChemistryResult<()> {
        for component in &self.components {
            if component.substance.as_usize() >= registry.substance_count() {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "mixture component '{}' uses substance index {} outside registry",
                    component.substance_id,
                    component.substance.as_usize()
                )));
            }
        }
        Ok(())
    }
}

fn validate_concentration(concentration: f64) -> ChemistryResult<()> {
    if !concentration.is_finite() || concentration < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "concentration must be non-negative and finite".to_string(),
        ));
    }
    Ok(())
}

fn validate_gaseous_fraction(gaseous_fraction: f64) -> ChemistryResult<()> {
    if !gaseous_fraction.is_finite() || !(0.0..=1.0).contains(&gaseous_fraction) {
        return Err(ChemistryError::InvalidMixtureState(
            "gas fraction must be within 0.0..=1.0".to_string(),
        ));
    }
    Ok(())
}

#[derive(Debug, Copy, Clone, PartialEq)]
struct OrderedF64(f64);

impl Eq for OrderedF64 {}

impl PartialOrd for OrderedF64 {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for OrderedF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&other.0)
    }
}

fn ordered_f64(value: f64) -> OrderedF64 {
    OrderedF64(value)
}

fn insert_sorted_unique<T: Ord + Copy>(values: &mut Vec<T>, value: T) {
    match values.binary_search(&value) {
        Ok(_) => {}
        Err(index) => values.insert(index, value),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::registry::ChemistryRegistryBuilder;
    use crate::chemistry::substance::Substance;

    fn test_registry() -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:water",
                0,
                18.0,
                18_000.0,
                373.0,
                75.0,
                40_650.0,
            ))
            .build()
            .unwrap()
    }

    #[test]
    fn recalculate_volume_rejects_invalid_mixture_state() {
        let registry = test_registry();
        let water: SubstanceId = "destroy:water".into();
        let water_index = registry.substance_index(&water).unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.components.push(MixtureComponent {
            substance: water_index,
            substance_id: water,
            concentration_mol_per_bucket: 1.0,
            gaseous_fraction: 1.5,
        });
        mixture.positions_by_substance = vec![Some(0)];

        let error = mixture
            .recalculate_volume_millibuckets(&registry, 1000)
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidMixtureState(_)));
    }

    #[test]
    fn public_substance_accessors_match_index_storage() {
        let registry = test_registry();
        let water: SubstanceId = "destroy:water".into();
        let water_index = registry.substance_index(&water).unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, water.clone(), 2.0)
            .unwrap();
        mixture
            .set_gaseous_fraction(&registry, water.clone(), 0.25)
            .unwrap();

        assert_eq!(mixture.concentration_of(&water), 2.0);
        assert_eq!(mixture.concentration_of_index(water_index), 2.0);
        assert_eq!(mixture.gaseous_fraction_of(&water), 0.25);
        assert_eq!(mixture.gaseous_fraction_of_index(water_index), 0.25);
        assert_eq!(
            mixture.substances().cloned().collect::<Vec<_>>(),
            vec![water]
        );
    }
}
