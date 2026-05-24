use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::registry::ChemistryRegistry;
use super::substance::{SubstanceId, SubstanceTagId};

pub const DEFAULT_TEMPERATURE_KELVIN: f64 = 298.0;
pub const TRACE_CONCENTRATION_MOL_PER_BUCKET: f64 = 1.0 / 512.0 / 512.0;

#[derive(Debug, Clone)]
pub struct Mixture {
    temperature_kelvin: f64,
    concentrations_mol_per_bucket: BTreeMap<SubstanceId, f64>,
    gaseous_fractions: BTreeMap<SubstanceId, f64>,
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
            concentrations_mol_per_bucket: BTreeMap::new(),
            gaseous_fractions: BTreeMap::new(),
        })
    }

    pub fn empty() -> Self {
        Self::new(DEFAULT_TEMPERATURE_KELVIN).expect("default temperature must be valid")
    }

    pub fn temperature_kelvin(&self) -> f64 {
        self.temperature_kelvin
    }

    pub fn concentration_of(&self, substance_id: &SubstanceId) -> f64 {
        self.concentrations_mol_per_bucket
            .get(substance_id)
            .copied()
            .unwrap_or(0.0)
    }

    pub fn gaseous_fraction_of(&self, substance_id: &SubstanceId) -> f64 {
        self.gaseous_fractions
            .get(substance_id)
            .copied()
            .unwrap_or(0.0)
    }

    pub fn substances(&self) -> impl Iterator<Item = &SubstanceId> {
        self.concentrations_mol_per_bucket.keys()
    }

    pub fn add_substance(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: impl Into<SubstanceId>,
        concentration_mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        let substance_id = substance_id.into();
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
        let new_concentration = self.concentration_of(&substance_id) + concentration_mol_per_bucket;
        self.concentrations_mol_per_bucket
            .insert(substance_id.clone(), new_concentration);
        self.gaseous_fractions
            .entry(substance_id.clone())
            .or_insert_with(|| {
                if registry
                    .substance(&substance_id)
                    .expect("already validated")
                    .boiling_point_kelvin
                    < self.temperature_kelvin
                {
                    1.0
                } else {
                    0.0
                }
            });
        self.validate(registry)
    }

    pub fn set_gaseous_fraction(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: impl Into<SubstanceId>,
        gaseous_fraction: f64,
    ) -> ChemistryResult<()> {
        let substance_id = substance_id.into();
        registry.substance(&substance_id)?;
        validate_gaseous_fraction(gaseous_fraction)?;
        if self.concentration_of(&substance_id) > 0.0 {
            self.gaseous_fractions
                .insert(substance_id, gaseous_fraction);
        }
        self.validate(registry)
    }

    pub fn change_concentration(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
        delta_mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        registry.substance(substance_id)?;
        if !delta_mol_per_bucket.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(
                "concentration change must be finite".to_string(),
            ));
        }
        let current = self.concentration_of(substance_id);
        let next = current + delta_mol_per_bucket;
        if next < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "substance '{substance_id}' would become negative: {next}"
            )));
        }
        if next <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            self.concentrations_mol_per_bucket.remove(substance_id);
            self.gaseous_fractions.remove(substance_id);
        } else {
            self.concentrations_mol_per_bucket
                .insert(substance_id.clone(), next);
            self.gaseous_fractions
                .entry(substance_id.clone())
                .or_insert_with(|| {
                    if registry
                        .substance(substance_id)
                        .expect("already validated")
                        .boiling_point_kelvin
                        < self.temperature_kelvin
                    {
                        1.0
                    } else {
                        0.0
                    }
                });
        }
        self.validate(registry)
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
        if energy_j_per_bucket == 0.0 || self.concentrations_mol_per_bucket.is_empty() {
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
                if let Some((substance_id, boiling_point)) =
                    self.next_higher_boiling_point(registry)?
                {
                    if self.temperature_kelvin + temperature_change >= boiling_point {
                        let energy_to_boiling =
                            (boiling_point - self.temperature_kelvin) * heat_capacity;
                        self.temperature_kelvin = boiling_point;
                        energy_j_per_bucket -= energy_to_boiling;

                        let concentration = self.concentration_of(&substance_id);
                        let current_gas = self.gaseous_fraction_of(&substance_id);
                        let liquid_concentration = concentration * (1.0 - current_gas);
                        let latent_heat = registry.substance(&substance_id)?.latent_heat_j_per_mol;
                        let energy_to_fully_boil = liquid_concentration * latent_heat;
                        if energy_to_fully_boil <= 0.0 {
                            self.gaseous_fractions.insert(substance_id, 1.0);
                            continue;
                        }
                        if energy_j_per_bucket >= energy_to_fully_boil {
                            self.gaseous_fractions.insert(substance_id, 1.0);
                            energy_j_per_bucket -= energy_to_fully_boil;
                        } else {
                            let boiled_fraction =
                                energy_j_per_bucket / (latent_heat * concentration);
                            self.gaseous_fractions
                                .insert(substance_id, current_gas + boiled_fraction);
                            energy_j_per_bucket = 0.0;
                        }
                        continue;
                    }
                }
                self.temperature_kelvin += temperature_change;
                energy_j_per_bucket = 0.0;
            } else {
                if let Some((substance_id, boiling_point)) =
                    self.next_lower_boiling_point(registry)?
                {
                    if self.temperature_kelvin + temperature_change < boiling_point {
                        let energy_to_boiling =
                            (boiling_point - self.temperature_kelvin) * heat_capacity;
                        self.temperature_kelvin = boiling_point;
                        energy_j_per_bucket -= energy_to_boiling;

                        let concentration = self.concentration_of(&substance_id);
                        let current_gas = self.gaseous_fraction_of(&substance_id);
                        let gas_concentration = concentration * current_gas;
                        let latent_heat = registry.substance(&substance_id)?.latent_heat_j_per_mol;
                        let released_by_full_condensation = gas_concentration * latent_heat;
                        if released_by_full_condensation <= 0.0 {
                            self.gaseous_fractions.insert(substance_id, 0.0);
                            continue;
                        }
                        if -energy_j_per_bucket >= released_by_full_condensation {
                            self.gaseous_fractions.insert(substance_id, 0.0);
                            energy_j_per_bucket += released_by_full_condensation;
                        } else {
                            let condensed_fraction =
                                -energy_j_per_bucket / (latent_heat * concentration);
                            self.gaseous_fractions
                                .insert(substance_id, current_gas - condensed_fraction);
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

    pub fn volumetric_heat_capacity_j_per_bucket_kelvin(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<f64> {
        self.concentrations_mol_per_bucket.iter().try_fold(
            0.0,
            |acc, (substance_id, concentration)| {
                Ok(acc
                    + registry
                        .substance(substance_id)?
                        .molar_heat_capacity_j_per_mol_kelvin
                        * concentration)
            },
        )
    }

    pub fn recalculate_volume_millibuckets(
        &self,
        registry: &ChemistryRegistry,
        initial_millibuckets: u32,
    ) -> ChemistryResult<u32> {
        self.validate(registry)?;
        let mut liquid_buckets = 0.0;
        for (substance_id, concentration) in &self.concentrations_mol_per_bucket {
            let substance = registry.substance(substance_id)?;
            let liquid_fraction = 1.0 - self.gaseous_fraction_of(substance_id);
            liquid_buckets += concentration * substance.molar_mass_grams * liquid_fraction
                / substance.liquid_density_grams_per_bucket;
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
            for (substance_id, concentration) in &mixture.concentrations_mol_per_bucket {
                let moles = concentration * amount;
                *moles_by_substance
                    .entry(substance_id.clone())
                    .or_insert(0.0) += moles;
                let substance = registry.substance(substance_id)?;
                total_energy += substance.molar_heat_capacity_j_per_mol_kelvin
                    * concentration
                    * mixture.temperature_kelvin
                    * amount;
                total_energy += substance.latent_heat_j_per_mol
                    * concentration
                    * mixture.gaseous_fraction_of(substance_id)
                    * amount;
            }
        }

        let mut result = Mixture::new(0.0)?;
        for (substance_id, moles) in moles_by_substance {
            result.add_substance(registry, substance_id.clone(), moles / total_amount)?;
            result.gaseous_fractions.insert(substance_id, 0.0);
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
        for (substance_id, concentration) in &self.concentrations_mol_per_bucket {
            registry.substance(substance_id)?;
            validate_concentration(*concentration)?;
            let gas = self.gaseous_fraction_of(substance_id);
            validate_gaseous_fraction(gas)?;
        }
        for substance_id in self.gaseous_fractions.keys() {
            if !self
                .concentrations_mol_per_bucket
                .contains_key(substance_id)
            {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "substance '{substance_id}' has gas fraction but no concentration"
                )));
            }
        }
        Ok(())
    }

    fn next_higher_boiling_point(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Option<(SubstanceId, f64)>> {
        let mut candidates = BTreeSet::new();
        for substance_id in self.substances() {
            let substance = registry.substance(substance_id)?;
            if substance.boiling_point_kelvin >= self.temperature_kelvin
                && self.gaseous_fraction_of(substance_id) < 1.0
            {
                candidates.insert((
                    ordered_f64(substance.boiling_point_kelvin),
                    substance_id.clone(),
                ));
            }
        }
        Ok(candidates.into_iter().next().map(|(bp, id)| (id, bp.0)))
    }

    fn next_lower_boiling_point(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Option<(SubstanceId, f64)>> {
        let mut best: Option<(SubstanceId, f64)> = None;
        for substance_id in self.substances() {
            let substance = registry.substance(substance_id)?;
            if substance.boiling_point_kelvin <= self.temperature_kelvin
                && self.gaseous_fraction_of(substance_id) > 0.0
                && best
                    .as_ref()
                    .map(|(_, bp)| substance.boiling_point_kelvin > *bp)
                    .unwrap_or(true)
            {
                best = Some((substance_id.clone(), substance.boiling_point_kelvin));
            }
        }
        Ok(best)
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
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .concentrations_mol_per_bucket
            .insert(water.clone(), 1.0);
        mixture.gaseous_fractions.insert(water, 1.5);

        let error = mixture
            .recalculate_volume_millibuckets(&registry, 1000)
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidMixtureState(_)));
    }
}
