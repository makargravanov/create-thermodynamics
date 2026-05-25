use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::reaction::GAS_CONSTANT_J_PER_MOL_KELVIN;
use super::registry::{ChemistryRegistry, GasSolubilityModel, SolventMiscibility, SubstanceIndex};
use super::solution::{self, SolutionState};
use super::substance::{
    LiquidPhasePreference, Substance, SubstanceAggregateState, SubstanceId, SubstanceTagId,
};

pub const DEFAULT_TEMPERATURE_KELVIN: f64 = 298.0;
pub const TRACE_CONCENTRATION_MOL_PER_BUCKET: f64 = 1.0 / 512.0 / 512.0;
const BUCKET_VOLUME_CUBIC_METERS: f64 = 0.001;

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MixturePhase {
    Aqueous,
    Organic,
    Gas,
    Solid,
}

impl MixturePhase {
    pub const ALL: [MixturePhase; 4] = [
        MixturePhase::Aqueous,
        MixturePhase::Organic,
        MixturePhase::Gas,
        MixturePhase::Solid,
    ];

    fn is_liquid(self) -> bool {
        matches!(self, MixturePhase::Aqueous | MixturePhase::Organic)
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PhaseAmount {
    pub phase: MixturePhase,
    pub concentration_mol_per_bucket: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OrganicPhaseAmount {
    pub solvent_id: SubstanceId,
    pub concentration_mol_per_bucket: f64,
}

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
    aqueous_mol_per_bucket: f64,
    organic_mol_per_bucket_by_solvent: BTreeMap<SubstanceIndex, f64>,
    gas_mol_per_bucket: f64,
    solid_mol_per_bucket: f64,
}

#[derive(Debug, Clone, Default)]
struct ComponentPhaseAmounts {
    aqueous_mol_per_bucket: f64,
    organic_mol_per_bucket_by_solvent: BTreeMap<SubstanceIndex, f64>,
    gas_mol_per_bucket: f64,
    solid_mol_per_bucket: f64,
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
            .map(MixtureComponent::total_concentration)
            .unwrap_or(0.0)
    }

    pub fn gaseous_fraction_of(&self, substance_id: &SubstanceId) -> f64 {
        self.components
            .iter()
            .find(|component| &component.substance_id == substance_id)
            .map(MixtureComponent::gaseous_fraction)
            .unwrap_or(0.0)
    }

    pub fn concentration_in_phase(&self, substance_id: &SubstanceId, phase: MixturePhase) -> f64 {
        self.components
            .iter()
            .find(|component| &component.substance_id == substance_id)
            .map(|component| component.amount_in_phase(phase))
            .unwrap_or(0.0)
    }

    pub fn phase_amounts_of(&self, substance_id: &SubstanceId) -> Vec<PhaseAmount> {
        self.components
            .iter()
            .find(|component| &component.substance_id == substance_id)
            .map(|component| {
                MixturePhase::ALL
                    .iter()
                    .copied()
                    .map(|phase| PhaseAmount {
                        phase,
                        concentration_mol_per_bucket: component.amount_in_phase(phase),
                    })
                    .collect()
            })
            .unwrap_or_default()
    }

    pub fn organic_phase_amounts_of(
        &self,
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<Vec<OrganicPhaseAmount>> {
        let Some(component) = self
            .components
            .iter()
            .find(|component| &component.substance_id == substance_id)
        else {
            return Ok(Vec::new());
        };
        component
            .organic_mol_per_bucket_by_solvent
            .iter()
            .map(|(solvent, concentration)| {
                Ok(OrganicPhaseAmount {
                    solvent_id: registry.substance_by_index(*solvent)?.id.clone(),
                    concentration_mol_per_bucket: *concentration,
                })
            })
            .collect()
    }

    pub fn concentration_in_organic_solvent(
        &self,
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
        solvent_id: &SubstanceId,
    ) -> ChemistryResult<f64> {
        let solvent = registry.substance_index(solvent_id).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!("unknown solvent substance '{solvent_id}'"))
        })?;
        Ok(self
            .components
            .iter()
            .find(|component| &component.substance_id == substance_id)
            .map(|component| component.amount_in_organic_solvent(solvent))
            .unwrap_or(0.0))
    }

    pub fn gas_pressure_pascal(&self) -> f64 {
        let gas_mol_per_bucket = self
            .components
            .iter()
            .map(|component| component.amount_in_phase(MixturePhase::Gas))
            .sum::<f64>();
        gas_mol_per_bucket * GAS_CONSTANT_J_PER_MOL_KELVIN * self.temperature_kelvin
            / BUCKET_VOLUME_CUBIC_METERS
    }

    pub fn aqueous_ionic_strength(&self, registry: &ChemistryRegistry) -> ChemistryResult<f64> {
        self.validate_registry_shape(registry)?;
        let mut ionic_strength = 0.0;
        let properties = registry.substance_properties();
        for component in &self.components {
            let charge = properties.charge[component.substance.as_usize()] as f64;
            if charge != 0.0 {
                ionic_strength +=
                    0.5 * component.amount_in_phase(MixturePhase::Aqueous) * charge * charge;
            }
        }
        if !ionic_strength.is_finite() || ionic_strength < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "ionic strength must be non-negative and finite: {ionic_strength}"
            )));
        }
        Ok(ionic_strength)
    }

    pub fn activity_of(
        &self,
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
        phase: MixturePhase,
    ) -> ChemistryResult<f64> {
        solution::activity_of(registry, self, substance_id, phase)
    }

    pub fn ph(&self, registry: &ChemistryRegistry) -> ChemistryResult<Option<f64>> {
        Ok(self.solution_state(registry)?.ph)
    }

    pub fn solution_state(&self, registry: &ChemistryRegistry) -> ChemistryResult<SolutionState> {
        solution::solution_state(registry, self)
    }

    pub fn transfer_gases_toward_solubility_equilibrium(
        &mut self,
        registry: &ChemistryRegistry,
        ticks: f64,
    ) -> ChemistryResult<f64> {
        if !ticks.is_finite() || ticks < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "gas transfer time must be non-negative and finite".to_string(),
            ));
        }
        if ticks == 0.0 {
            return self.validate(registry).map(|_| 0.0);
        }
        self.validate_registry_shape(registry)?;
        let gas_pressure_pascal = self.gas_pressure_pascal();
        let ionic_strength = self.aqueous_ionic_strength(registry)?;
        let mut max_delta = 0.0_f64;
        let gas_substances = self
            .components
            .iter()
            .filter_map(|component| {
                registry
                    .gas_solubility(component.substance)
                    .map(|model| (component.substance, model.clone()))
            })
            .collect::<Vec<_>>();
        for (substance, model) in gas_substances {
            let target_dissolved = gas_dissolved_limit(
                Some(&model),
                gas_pressure_pascal,
                ionic_strength,
                self.temperature_kelvin,
            )?;
            let coefficient = gas_transfer_coefficient_per_tick(&model);
            let fraction = 1.0 - (-coefficient * ticks).exp();
            if fraction <= 0.0 {
                continue;
            }
            let Some(position) = self.position_of_substance(substance) else {
                continue;
            };
            let current_dissolved = self.components[position]
                .amount_in_phases(&[MixturePhase::Aqueous, MixturePhase::Organic]);
            let difference = target_dissolved - current_dissolved;
            if difference.abs() <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                continue;
            }
            if difference > 0.0 {
                let available_gas = self.components[position].amount_in_phase(MixturePhase::Gas);
                let moved = (difference * fraction).min(available_gas);
                if moved <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                    continue;
                }
                self.components[position].remove_from_phase(MixturePhase::Gas, moved)?;
                let preferred = preferred_phase(
                    registry
                        .substance_by_index(substance)?
                        .phase_properties
                        .preferred_liquid_phase,
                );
                self.add_to_phase_for_substance(registry, substance, preferred, moved)?;
                max_delta = max_delta.max(moved);
            } else {
                let moved = (-difference * fraction).min(current_dissolved);
                if moved <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                    continue;
                }
                self.components[position]
                    .remove_from_phases(&[MixturePhase::Aqueous, MixturePhase::Organic], moved)?;
                self.components[position].add_to_phase(MixturePhase::Gas, moved);
                max_delta = max_delta.max(moved);
            }
        }
        self.equilibrate_phases(registry)?;
        Ok(max_delta)
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
            .map(|position| self.components[position].total_concentration())
            .unwrap_or(0.0)
    }

    pub(crate) fn concentration_of_index_in_phases(
        &self,
        substance: SubstanceIndex,
        phases: &[MixturePhase],
    ) -> f64 {
        self.position_of_substance(substance)
            .map(|position| self.components[position].amount_in_phases(phases))
            .unwrap_or(0.0)
    }

    pub fn move_between_phases(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: impl Into<SubstanceId>,
        from: MixturePhase,
        to: MixturePhase,
        concentration_mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        let substance_id = substance_id.into();
        let substance = registry.substance_index(&substance_id).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!("unknown substance '{substance_id}'"))
        })?;
        validate_concentration(concentration_mol_per_bucket)?;
        if concentration_mol_per_bucket == 0.0 {
            return self.validate(registry);
        }
        let position = self.position_of_substance(substance).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!(
                "substance '{substance_id}' is not present in the mixture"
            ))
        })?;
        self.components[position].remove_from_phase(from, concentration_mol_per_bucket)?;
        self.components[position].add_to_phase(to, concentration_mol_per_bucket);
        self.remove_trace_component(substance);
        self.validate(registry)
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
        let initial_phase = initial_phase_for_substance(substance, self.temperature_kelvin)?;
        self.change_concentration_by_index_unchecked(
            registry,
            substance_index,
            substance_id,
            concentration_mol_per_bucket,
            initial_phase,
        )?;
        self.equilibrate_phases(registry)?;
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
            let total = self.components[position].total_concentration();
            let gas = total * gaseous_fraction;
            let liquid = total - gas;
            let phase = preferred_phase(
                registry
                    .substance_by_index(substance_index)?
                    .phase_properties
                    .preferred_liquid_phase,
            );
            self.components[position].clear_phase_amounts();
            self.components[position].add_to_phase(MixturePhase::Gas, gas);
            self.add_to_phase_for_substance(registry, substance_index, phase, liquid)?;
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
            let substance_data = registry.substance_by_index(substance)?;
            let initial_phase =
                initial_phase_for_substance(substance_data, self.temperature_kelvin)?;
            self.change_concentration_by_index_unchecked(
                registry,
                substance,
                substance_id,
                delta_mol_per_bucket,
                initial_phase,
            )?;
        }
        self.equilibrate_phases(registry)?;
        self.validate(registry)
    }

    pub(crate) fn apply_reaction_phase_deltas_by_index(
        &mut self,
        registry: &ChemistryRegistry,
        reactants: &[(SubstanceIndex, u32, Vec<MixturePhase>)],
        products: &[(SubstanceIndex, f64, MixturePhase)],
        moles_per_bucket: f64,
    ) -> ChemistryResult<f64> {
        if !moles_per_bucket.is_finite() || moles_per_bucket < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "reaction amount must be non-negative and finite".to_string(),
            ));
        }
        self.ensure_position_capacity(registry);
        let mut max_delta = 0.0_f64;
        for (substance, coefficient, phases) in reactants {
            let delta = *coefficient as f64 * moles_per_bucket;
            max_delta = max_delta.max(delta);
            let available = self.concentration_of_index_in_phases(*substance, phases);
            if available - delta < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
                let substance_id = &registry.substance_by_index(*substance)?.id;
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "substance '{substance_id}' would become negative in reaction-accessible phases: {}",
                    available - delta
                )));
            }
        }
        for (substance, coefficient, phases) in reactants {
            let amount = *coefficient as f64 * moles_per_bucket;
            if let Some(position) = self.position_of_substance(*substance) {
                self.components[position].remove_from_phases(phases, amount)?;
                self.remove_trace_component(*substance);
            }
        }
        for (substance, coefficient, phase) in products {
            if !coefficient.is_finite() || *coefficient <= 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "product coefficient must be positive and finite".to_string(),
                ));
            }
            let amount = *coefficient * moles_per_bucket;
            max_delta = max_delta.max(amount);
            let substance_data = registry.substance_by_index(*substance)?;
            if let Some(position) = self.position_of_substance(*substance) {
                self.components[position].add_to_phase(*phase, amount);
            } else {
                self.insert_component(*substance, substance_data.id.clone(), amount, *phase);
            }
        }
        self.equilibrate_phases(registry)?;
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
                if let Some((substance_index, melting_point)) =
                    self.next_higher_melting_point(registry)?
                {
                    if self.temperature_kelvin + temperature_change >= melting_point {
                        let energy_to_melting =
                            (melting_point - self.temperature_kelvin) * heat_capacity;
                        self.temperature_kelvin = melting_point;
                        energy_j_per_bucket -= energy_to_melting;

                        let solid_concentration = self.concentration_of_index_in_phases(
                            substance_index,
                            &[MixturePhase::Solid],
                        );
                        let fusion_heat = registry.substance_properties().fusion_heat_j_per_mol
                            [substance_index.as_usize()];
                        let energy_to_fully_melt = solid_concentration * fusion_heat;
                        if energy_to_fully_melt <= 0.0 {
                            self.move_all_solid_to_preferred_liquid(registry, substance_index)?;
                            continue;
                        }
                        if energy_j_per_bucket >= energy_to_fully_melt {
                            self.move_all_solid_to_preferred_liquid(registry, substance_index)?;
                            energy_j_per_bucket -= energy_to_fully_melt;
                        } else {
                            let melted_amount = energy_j_per_bucket / fusion_heat;
                            self.move_solid_to_preferred_liquid(
                                registry,
                                substance_index,
                                melted_amount,
                            )?;
                            energy_j_per_bucket = 0.0;
                        }
                        continue;
                    }
                }
                if let Some((substance_index, boiling_point)) =
                    self.next_higher_boiling_point(registry)?
                {
                    if self.temperature_kelvin + temperature_change >= boiling_point {
                        let energy_to_boiling =
                            (boiling_point - self.temperature_kelvin) * heat_capacity;
                        self.temperature_kelvin = boiling_point;
                        energy_j_per_bucket -= energy_to_boiling;

                        let liquid_concentration =
                            self.liquid_concentration_of_index(substance_index);
                        let latent_heat = registry.substance_properties().latent_heat_j_per_mol
                            [substance_index.as_usize()];
                        let energy_to_fully_boil = liquid_concentration * latent_heat;
                        if energy_to_fully_boil <= 0.0 {
                            self.move_all_liquid_to_gas(substance_index);
                            continue;
                        }
                        if energy_j_per_bucket >= energy_to_fully_boil {
                            self.move_all_liquid_to_gas(substance_index);
                            energy_j_per_bucket -= energy_to_fully_boil;
                        } else {
                            let boiled_amount = energy_j_per_bucket / latent_heat;
                            self.move_liquid_to_gas(substance_index, boiled_amount)?;
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

                        let gas_concentration = self.concentration_of_index_in_phases(
                            substance_index,
                            &[MixturePhase::Gas],
                        );
                        let latent_heat = registry.substance_properties().latent_heat_j_per_mol
                            [substance_index.as_usize()];
                        let released_by_full_condensation = gas_concentration * latent_heat;
                        if released_by_full_condensation <= 0.0 {
                            self.move_all_gas_to_preferred_liquid(registry, substance_index)?;
                            continue;
                        }
                        if -energy_j_per_bucket >= released_by_full_condensation {
                            self.move_all_gas_to_preferred_liquid(registry, substance_index)?;
                            energy_j_per_bucket += released_by_full_condensation;
                        } else {
                            let condensed_amount = -energy_j_per_bucket / latent_heat;
                            self.move_gas_to_preferred_liquid(
                                registry,
                                substance_index,
                                condensed_amount,
                            )?;
                            energy_j_per_bucket = 0.0;
                        }
                        continue;
                    }
                }
                if let Some((substance_index, melting_point)) =
                    self.next_lower_melting_point(registry)?
                {
                    if self.temperature_kelvin + temperature_change < melting_point {
                        let energy_to_melting =
                            (melting_point - self.temperature_kelvin) * heat_capacity;
                        self.temperature_kelvin = melting_point;
                        energy_j_per_bucket -= energy_to_melting;

                        let liquid_concentration =
                            self.liquid_concentration_of_index(substance_index);
                        let fusion_heat = registry.substance_properties().fusion_heat_j_per_mol
                            [substance_index.as_usize()];
                        let released_by_full_freezing = liquid_concentration * fusion_heat;
                        if released_by_full_freezing <= 0.0 {
                            self.move_all_liquid_to_solid(substance_index);
                            continue;
                        }
                        if -energy_j_per_bucket >= released_by_full_freezing {
                            self.move_all_liquid_to_solid(substance_index);
                            energy_j_per_bucket += released_by_full_freezing;
                        } else {
                            let frozen_amount = -energy_j_per_bucket / fusion_heat;
                            self.move_liquid_to_solid(substance_index, frozen_amount)?;
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

    pub fn volumetric_heat_capacity_j_per_bucket_kelvin(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<f64> {
        self.validate_registry_shape(registry)?;
        let properties = registry.substance_properties();
        Ok(self.components.iter().fold(0.0, |acc, component| {
            acc + properties.molar_heat_capacity_j_per_mol_kelvin[component.substance.as_usize()]
                * component.total_concentration()
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
            let index = component.substance.as_usize();
            liquid_buckets += component
                .amount_in_phases(&[MixturePhase::Aqueous, MixturePhase::Organic])
                * properties.molar_mass_grams[index]
                / properties.liquid_density_grams_per_bucket[index];
            liquid_buckets += component.amount_in_phase(MixturePhase::Solid)
                * properties.molar_mass_grams[index]
                / properties.solid_density_grams_per_bucket[index];
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

        let mut moles_by_substance: BTreeMap<SubstanceId, ComponentPhaseAmounts> = BTreeMap::new();
        let mut total_energy = 0.0;
        for (mixture, amount) in mixtures {
            if !amount.is_finite() || *amount < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "mixed amount must be non-negative and finite".to_string(),
                ));
            }
            mixture.validate(registry)?;
            for component in &mixture.components {
                let entry = moles_by_substance
                    .entry(component.substance_id.clone())
                    .or_default();
                entry.aqueous_mol_per_bucket +=
                    component.amount_in_phase(MixturePhase::Aqueous) * amount;
                entry.gas_mol_per_bucket += component.amount_in_phase(MixturePhase::Gas) * amount;
                entry.solid_mol_per_bucket +=
                    component.amount_in_phase(MixturePhase::Solid) * amount;
                for (solvent, concentration) in &component.organic_mol_per_bucket_by_solvent {
                    *entry
                        .organic_mol_per_bucket_by_solvent
                        .entry(*solvent)
                        .or_insert(0.0) += concentration * amount;
                }
                let substance = registry.substance_by_index(component.substance)?;
                total_energy += substance.molar_heat_capacity_j_per_mol_kelvin
                    * component.total_concentration()
                    * mixture.temperature_kelvin
                    * amount;
                total_energy += substance.latent_heat_j_per_mol
                    * component.amount_in_phase(MixturePhase::Gas)
                    * amount;
            }
        }

        let mut result = Mixture::new(0.0)?;
        result.ensure_position_capacity(registry);
        for (substance_id, phase_amounts) in moles_by_substance {
            let substance_index = registry.substance_index(&substance_id).ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!("unknown substance '{substance_id}'"))
            })?;
            result.insert_component_with_phase_amounts(
                substance_index,
                substance_id,
                normalize_phase_amounts(phase_amounts, total_amount),
            );
        }
        result.heat(registry, total_energy / total_amount)?;
        result.equilibrate_phases(registry)?;
        Ok(result)
    }

    pub fn equilibrate_phases(&mut self, registry: &ChemistryRegistry) -> ChemistryResult<()> {
        self.validate_registry_shape(registry)?;
        let totals_before = self
            .components
            .iter()
            .map(|component| (component.substance, component.total_concentration()))
            .collect::<Vec<_>>();
        let solvent_clusters = self.solvent_clusters(registry)?;
        for position in 0..self.components.len() {
            let substance = registry.substance_by_index(self.components[position].substance)?;
            let total = self.components[position].total_concentration();
            let mut next = ComponentPhaseAmounts::default();
            match substance.aggregate_state_at(self.temperature_kelvin)? {
                SubstanceAggregateState::Gas => {
                    let current_gas = self.components[position].amount_in_phase(MixturePhase::Gas);
                    let dissolved = (total - current_gas).max(0.0);
                    next.gas_mol_per_bucket = total - dissolved;
                    if dissolved > 0.0 {
                        distribute_condensed_amount(
                            substance.phase_properties.preferred_liquid_phase,
                            substance.phase_properties.aqueous_solubility_mol_per_bucket,
                            substance.phase_properties.organic_solubility_mol_per_bucket,
                            substance.phase_properties.can_precipitate,
                            self.organic_solvent_for_component(
                                registry,
                                position,
                                &solvent_clusters,
                            )?,
                            dissolved,
                            &substance.id,
                            &mut next,
                        )?;
                    }
                }
                SubstanceAggregateState::Solid => {
                    let dissolved = condensed_solubility_capacity(substance)
                        .map(|limit| total.min(limit))
                        .unwrap_or(0.0);
                    if dissolved > 0.0 {
                        distribute_condensed_amount(
                            substance.phase_properties.preferred_liquid_phase,
                            substance.phase_properties.aqueous_solubility_mol_per_bucket,
                            substance.phase_properties.organic_solubility_mol_per_bucket,
                            substance.phase_properties.can_precipitate,
                            self.organic_solvent_for_component(
                                registry,
                                position,
                                &solvent_clusters,
                            )?,
                            dissolved,
                            &substance.id,
                            &mut next,
                        )?;
                    }
                    let remaining = total - dissolved;
                    if remaining > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                        if substance.phase_properties.can_precipitate {
                            next.solid_mol_per_bucket += remaining;
                        } else {
                            return Err(ChemistryError::InvalidMixtureState(format!(
                                "substance '{}' is solid at current temperature but cannot precipitate",
                                substance.id
                            )));
                        }
                    }
                }
                SubstanceAggregateState::Liquid => {
                    distribute_condensed_amount(
                        substance.phase_properties.preferred_liquid_phase,
                        substance.phase_properties.aqueous_solubility_mol_per_bucket,
                        substance.phase_properties.organic_solubility_mol_per_bucket,
                        substance.phase_properties.can_precipitate,
                        self.organic_solvent_for_component(registry, position, &solvent_clusters)?,
                        total,
                        &substance.id,
                        &mut next,
                    )?;
                }
            }
            self.components[position].set_phase_amounts(next);
        }
        for (substance, before) in totals_before {
            let after = self.concentration_of_index(substance);
            if (before - after).abs() > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                let substance_id = &registry.substance_by_index(substance)?.id;
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "phase equilibration changed amount of '{substance_id}': before {before}, after {after}"
                )));
            }
        }
        self.validate(registry)
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
            for amount in component.phase_amounts() {
                validate_concentration(amount)?;
            }
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
            if boiling_point >= self.temperature_kelvin && component.condensed_concentration() > 0.0
            {
                candidates.insert((ordered_f64(boiling_point), component.substance));
            }
        }
        Ok(candidates.into_iter().next().map(|(bp, id)| (id, bp.0)))
    }

    fn next_higher_melting_point(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Option<(SubstanceIndex, f64)>> {
        let mut candidates = BTreeSet::new();
        let properties = registry.substance_properties();
        for component in &self.components {
            let melting_point = properties.melting_point_kelvin[component.substance.as_usize()];
            if melting_point >= self.temperature_kelvin
                && component.amount_in_phase(MixturePhase::Solid) > 0.0
            {
                candidates.insert((ordered_f64(melting_point), component.substance));
            }
        }
        Ok(candidates.into_iter().next().map(|(mp, id)| (id, mp.0)))
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
                && component.amount_in_phase(MixturePhase::Gas) > 0.0
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

    fn next_lower_melting_point(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Option<(SubstanceIndex, f64)>> {
        let mut best: Option<(SubstanceIndex, f64)> = None;
        let properties = registry.substance_properties();
        for component in &self.components {
            let melting_point = properties.melting_point_kelvin[component.substance.as_usize()];
            if melting_point <= self.temperature_kelvin
                && self.liquid_concentration_of_index(component.substance) > 0.0
                && best
                    .as_ref()
                    .map(|(_, mp)| melting_point > *mp)
                    .unwrap_or(true)
            {
                best = Some((component.substance, melting_point));
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
        phase: MixturePhase,
    ) {
        let position = self.components.len();
        let mut component = MixtureComponent {
            substance,
            substance_id,
            aqueous_mol_per_bucket: 0.0,
            organic_mol_per_bucket_by_solvent: BTreeMap::new(),
            gas_mol_per_bucket: 0.0,
            solid_mol_per_bucket: 0.0,
        };
        component.add_to_phase_for_solvent(phase, substance, concentration_mol_per_bucket);
        self.components.push(component);
        if self.positions_by_substance.len() <= substance.as_usize() {
            self.positions_by_substance
                .resize(substance.as_usize() + 1, None);
        }
        self.positions_by_substance[substance.as_usize()] = Some(position);
    }

    fn insert_component_with_phase_amounts(
        &mut self,
        substance: SubstanceIndex,
        substance_id: SubstanceId,
        phase_amounts: ComponentPhaseAmounts,
    ) {
        let position = self.components.len();
        self.components.push(MixtureComponent {
            substance,
            substance_id,
            aqueous_mol_per_bucket: phase_amounts.aqueous_mol_per_bucket,
            organic_mol_per_bucket_by_solvent: phase_amounts.organic_mol_per_bucket_by_solvent,
            gas_mol_per_bucket: phase_amounts.gas_mol_per_bucket,
            solid_mol_per_bucket: phase_amounts.solid_mol_per_bucket,
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

    fn change_concentration_by_index_unchecked(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        substance_id: SubstanceId,
        delta_mol_per_bucket: f64,
        initial_phase: MixturePhase,
    ) -> ChemistryResult<()> {
        registry.substance_by_index(substance)?;
        if let Some(position) = self.position_of_substance(substance) {
            self.components[position].add_to_phase(initial_phase, delta_mol_per_bucket);
        } else {
            self.insert_component(substance, substance_id, delta_mol_per_bucket, initial_phase);
        }
        Ok(())
    }

    fn liquid_concentration_of_index(&self, substance: SubstanceIndex) -> f64 {
        self.concentration_of_index_in_phases(
            substance,
            &[MixturePhase::Aqueous, MixturePhase::Organic],
        )
    }

    fn move_liquid_to_gas(
        &mut self,
        substance: SubstanceIndex,
        amount: f64,
    ) -> ChemistryResult<()> {
        if let Some(position) = self.position_of_substance(substance) {
            self.components[position]
                .remove_from_phases(&[MixturePhase::Aqueous, MixturePhase::Organic], amount)?;
            self.components[position].add_to_phase(MixturePhase::Gas, amount);
        }
        Ok(())
    }

    fn move_all_liquid_to_gas(&mut self, substance: SubstanceIndex) {
        if let Some(position) = self.position_of_substance(substance) {
            let amount = self.components[position]
                .amount_in_phases(&[MixturePhase::Aqueous, MixturePhase::Organic]);
            self.components[position].aqueous_mol_per_bucket = 0.0;
            self.components[position]
                .organic_mol_per_bucket_by_solvent
                .clear();
            self.components[position].add_to_phase(MixturePhase::Gas, amount);
        }
    }

    fn move_gas_to_preferred_liquid(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        amount: f64,
    ) -> ChemistryResult<()> {
        if let Some(position) = self.position_of_substance(substance) {
            let preferred = preferred_phase(
                registry
                    .substance_by_index(substance)?
                    .phase_properties
                    .preferred_liquid_phase,
            );
            self.components[position].remove_from_phase(MixturePhase::Gas, amount)?;
            self.add_to_phase_for_substance(registry, substance, preferred, amount)?;
        }
        Ok(())
    }

    fn move_solid_to_preferred_liquid(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        amount: f64,
    ) -> ChemistryResult<()> {
        if let Some(position) = self.position_of_substance(substance) {
            let preferred = preferred_phase(
                registry
                    .substance_by_index(substance)?
                    .phase_properties
                    .preferred_liquid_phase,
            );
            self.components[position].remove_from_phase(MixturePhase::Solid, amount)?;
            self.add_to_phase_for_substance(registry, substance, preferred, amount)?;
        }
        Ok(())
    }

    fn move_all_solid_to_preferred_liquid(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
    ) -> ChemistryResult<()> {
        let amount = self.concentration_of_index_in_phases(substance, &[MixturePhase::Solid]);
        self.move_solid_to_preferred_liquid(registry, substance, amount)
    }

    fn move_liquid_to_solid(
        &mut self,
        substance: SubstanceIndex,
        amount: f64,
    ) -> ChemistryResult<()> {
        if let Some(position) = self.position_of_substance(substance) {
            self.components[position]
                .remove_from_phases(&[MixturePhase::Aqueous, MixturePhase::Organic], amount)?;
            self.components[position].add_to_phase(MixturePhase::Solid, amount);
        }
        Ok(())
    }

    fn move_all_liquid_to_solid(&mut self, substance: SubstanceIndex) {
        if let Some(position) = self.position_of_substance(substance) {
            let amount = self.components[position]
                .amount_in_phases(&[MixturePhase::Aqueous, MixturePhase::Organic]);
            self.components[position].aqueous_mol_per_bucket = 0.0;
            self.components[position]
                .organic_mol_per_bucket_by_solvent
                .clear();
            self.components[position].add_to_phase(MixturePhase::Solid, amount);
        }
    }

    fn move_all_gas_to_preferred_liquid(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
    ) -> ChemistryResult<()> {
        let amount = self.concentration_of_index_in_phases(substance, &[MixturePhase::Gas]);
        self.move_gas_to_preferred_liquid(registry, substance, amount)
    }

    fn remove_trace_component(&mut self, substance: SubstanceIndex) {
        if self.concentration_of_index(substance) <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            self.remove_component(substance);
        }
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
            for solvent in component.organic_mol_per_bucket_by_solvent.keys() {
                if solvent.as_usize() >= registry.substance_count() {
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "mixture component '{}' uses solvent index {} outside registry",
                        component.substance_id,
                        solvent.as_usize()
                    )));
                }
            }
        }
        Ok(())
    }

    fn add_to_phase_for_substance(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        phase: MixturePhase,
        amount: f64,
    ) -> ChemistryResult<()> {
        if amount == 0.0 {
            return Ok(());
        }
        let Some(position) = self.position_of_substance(substance) else {
            return Ok(());
        };
        if phase == MixturePhase::Organic {
            let substance_id = registry.substance_by_index(substance)?.id.clone();
            let solvent = self.organic_solvent_for_substance(registry, substance)?.ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!(
                    "substance '{}' cannot enter an organic phase because no concrete organic solvent is available",
                    substance_id
                ))
            })?;
            self.components[position].add_to_phase_for_solvent(phase, solvent, amount);
        } else {
            self.components[position].add_to_phase(phase, amount);
        }
        Ok(())
    }

    fn organic_solvent_for_component(
        &self,
        registry: &ChemistryRegistry,
        position: usize,
        solvent_clusters: &[BTreeSet<SubstanceIndex>],
    ) -> ChemistryResult<Option<SubstanceIndex>> {
        let component = &self.components[position];
        if registry
            .substance_by_index(component.substance)?
            .phase_properties
            .can_form_liquid_phase
        {
            if let Some(cluster) = solvent_clusters
                .iter()
                .find(|cluster| cluster.contains(&component.substance))
            {
                return Ok(cluster.iter().next().copied());
            }
            return Ok(Some(component.substance));
        }
        if let Some(solvent) =
            component
                .organic_mol_per_bucket_by_solvent
                .iter()
                .find_map(|(solvent, amount)| {
                    (*amount > TRACE_CONCENTRATION_MOL_PER_BUCKET).then_some(*solvent)
                })
        {
            return Ok(Some(solvent));
        }
        self.first_available_organic_solvent(registry, solvent_clusters)
    }

    fn organic_solvent_for_substance(
        &self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
    ) -> ChemistryResult<Option<SubstanceIndex>> {
        if registry
            .substance_by_index(substance)?
            .phase_properties
            .can_form_liquid_phase
        {
            Ok(Some(substance))
        } else {
            self.first_available_organic_solvent(registry, &self.solvent_clusters(registry)?)
        }
    }

    fn first_available_organic_solvent(
        &self,
        registry: &ChemistryRegistry,
        solvent_clusters: &[BTreeSet<SubstanceIndex>],
    ) -> ChemistryResult<Option<SubstanceIndex>> {
        for cluster in solvent_clusters {
            if let Some(solvent) = cluster.iter().next().copied() {
                return Ok(Some(solvent));
            }
        }
        for component in &self.components {
            if registry
                .substance_by_index(component.substance)?
                .phase_properties
                .can_form_liquid_phase
                && component.total_concentration() > TRACE_CONCENTRATION_MOL_PER_BUCKET
            {
                return Ok(Some(component.substance));
            }
        }
        Ok(None)
    }

    fn solvent_clusters(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Vec<BTreeSet<SubstanceIndex>>> {
        let solvents = self
            .components
            .iter()
            .filter_map(|component| {
                registry
                    .substance_by_index(component.substance)
                    .ok()
                    .filter(|substance| substance.phase_properties.can_form_liquid_phase)
                    .and_then(|_| {
                        (component.total_concentration() > TRACE_CONCENTRATION_MOL_PER_BUCKET)
                            .then_some(component.substance)
                    })
            })
            .collect::<Vec<_>>();
        let mut clusters: Vec<BTreeSet<SubstanceIndex>> = Vec::new();
        'solvent: for solvent in solvents {
            for cluster in &mut clusters {
                let miscible = cluster.iter().any(|existing| {
                    let existing_amount = self.concentration_of_index(*existing);
                    let solvent_amount = self.concentration_of_index(solvent);
                    match registry.solvent_miscibility(*existing, solvent) {
                        SolventMiscibility::FullyMiscible => true,
                        SolventMiscibility::PartiallyMiscible {
                            limit_mol_per_bucket,
                        } => existing_amount.min(solvent_amount) <= limit_mol_per_bucket,
                        SolventMiscibility::Immiscible => false,
                    }
                });
                if miscible {
                    cluster.insert(solvent);
                    continue 'solvent;
                }
            }
            clusters.push(BTreeSet::from([solvent]));
        }
        Ok(clusters)
    }

    pub(crate) fn total_in_phase(&self, phase: MixturePhase) -> f64 {
        self.components
            .iter()
            .map(|component| component.amount_in_phase(phase))
            .sum()
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

impl MixtureComponent {
    fn total_concentration(&self) -> f64 {
        self.aqueous_mol_per_bucket
            + self.organic_mol_per_bucket_by_solvent.values().sum::<f64>()
            + self.gas_mol_per_bucket
            + self.solid_mol_per_bucket
    }

    fn condensed_concentration(&self) -> f64 {
        self.amount_in_phase(MixturePhase::Aqueous)
            + self.amount_in_phase(MixturePhase::Organic)
            + self.amount_in_phase(MixturePhase::Solid)
    }

    fn gaseous_fraction(&self) -> f64 {
        let total = self.total_concentration();
        if total <= 0.0 {
            0.0
        } else {
            self.amount_in_phase(MixturePhase::Gas) / total
        }
    }

    fn amount_in_phase(&self, phase: MixturePhase) -> f64 {
        match phase {
            MixturePhase::Aqueous => self.aqueous_mol_per_bucket,
            MixturePhase::Organic => self.organic_mol_per_bucket_by_solvent.values().sum(),
            MixturePhase::Gas => self.gas_mol_per_bucket,
            MixturePhase::Solid => self.solid_mol_per_bucket,
        }
    }

    fn amount_in_organic_solvent(&self, solvent: SubstanceIndex) -> f64 {
        self.organic_mol_per_bucket_by_solvent
            .get(&solvent)
            .copied()
            .unwrap_or(0.0)
    }

    fn phase_amounts(&self) -> Vec<f64> {
        let mut amounts = vec![
            self.aqueous_mol_per_bucket,
            self.gas_mol_per_bucket,
            self.solid_mol_per_bucket,
        ];
        amounts.extend(self.organic_mol_per_bucket_by_solvent.values().copied());
        amounts
    }

    fn amount_in_phases(&self, phases: &[MixturePhase]) -> f64 {
        phases
            .iter()
            .map(|phase| self.amount_in_phase(*phase))
            .sum()
    }

    fn add_to_phase(&mut self, phase: MixturePhase, amount: f64) {
        match phase {
            MixturePhase::Aqueous => self.aqueous_mol_per_bucket += amount,
            MixturePhase::Organic => self.add_to_phase_for_solvent(phase, self.substance, amount),
            MixturePhase::Gas => self.gas_mol_per_bucket += amount,
            MixturePhase::Solid => self.solid_mol_per_bucket += amount,
        }
    }

    fn add_to_phase_for_solvent(
        &mut self,
        phase: MixturePhase,
        solvent: SubstanceIndex,
        amount: f64,
    ) {
        if phase == MixturePhase::Organic {
            *self
                .organic_mol_per_bucket_by_solvent
                .entry(solvent)
                .or_insert(0.0) += amount;
        } else {
            self.add_to_phase(phase, amount);
        }
    }

    fn clear_phase_amounts(&mut self) {
        self.aqueous_mol_per_bucket = 0.0;
        self.organic_mol_per_bucket_by_solvent.clear();
        self.gas_mol_per_bucket = 0.0;
        self.solid_mol_per_bucket = 0.0;
    }

    fn set_phase_amounts(&mut self, phase_amounts: ComponentPhaseAmounts) {
        self.aqueous_mol_per_bucket = phase_amounts.aqueous_mol_per_bucket;
        self.organic_mol_per_bucket_by_solvent = phase_amounts.organic_mol_per_bucket_by_solvent;
        self.gas_mol_per_bucket = phase_amounts.gas_mol_per_bucket;
        self.solid_mol_per_bucket = phase_amounts.solid_mol_per_bucket;
    }

    fn remove_from_phase(&mut self, phase: MixturePhase, amount: f64) -> ChemistryResult<()> {
        self.remove_from_phases(&[phase], amount)
    }

    fn remove_from_phases(
        &mut self,
        phases: &[MixturePhase],
        mut amount: f64,
    ) -> ChemistryResult<()> {
        if !amount.is_finite() || amount < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "removed phase amount must be non-negative and finite".to_string(),
            ));
        }
        for phase in phases {
            let available = self.amount_in_phase(*phase);
            let removed = available.min(amount);
            self.remove_exact_from_phase(*phase, removed);
            amount -= removed;
            if amount <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                return Ok(());
            }
        }
        Err(ChemistryError::InvalidMixtureState(format!(
            "not enough substance '{}' in requested phases",
            self.substance_id
        )))
    }

    fn remove_exact_from_phase(&mut self, phase: MixturePhase, amount: f64) {
        match phase {
            MixturePhase::Aqueous => self.aqueous_mol_per_bucket -= amount,
            MixturePhase::Organic => {
                let mut remaining = amount;
                let solvents = self
                    .organic_mol_per_bucket_by_solvent
                    .keys()
                    .copied()
                    .collect::<Vec<_>>();
                for solvent in solvents {
                    if remaining <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                        break;
                    }
                    let Some(current) = self.organic_mol_per_bucket_by_solvent.get_mut(&solvent)
                    else {
                        continue;
                    };
                    let removed = (*current).min(remaining);
                    *current -= removed;
                    remaining -= removed;
                }
                self.organic_mol_per_bucket_by_solvent
                    .retain(|_, value| *value > TRACE_CONCENTRATION_MOL_PER_BUCKET);
            }
            MixturePhase::Gas => self.gas_mol_per_bucket -= amount,
            MixturePhase::Solid => self.solid_mol_per_bucket -= amount,
        }
    }
}

fn preferred_phase(preference: LiquidPhasePreference) -> MixturePhase {
    match preference {
        LiquidPhasePreference::Aqueous => MixturePhase::Aqueous,
        LiquidPhasePreference::Organic => MixturePhase::Organic,
    }
}

fn initial_phase_for_substance(
    substance: &Substance,
    temperature_kelvin: f64,
) -> ChemistryResult<MixturePhase> {
    match substance.aggregate_state_at(temperature_kelvin)? {
        SubstanceAggregateState::Gas => Ok(MixturePhase::Gas),
        SubstanceAggregateState::Solid => Ok(MixturePhase::Solid),
        SubstanceAggregateState::Liquid => Ok(preferred_phase(
            substance.phase_properties.preferred_liquid_phase,
        )),
    }
}

fn condensed_solubility_capacity(substance: &Substance) -> Option<f64> {
    let aqueous = substance.phase_properties.aqueous_solubility_mol_per_bucket;
    let organic = substance.phase_properties.organic_solubility_mol_per_bucket;
    match (aqueous, organic) {
        (None, _) | (_, None) => None,
        (Some(aqueous), Some(organic)) => Some(aqueous + organic),
    }
}

fn gas_dissolved_limit(
    model: Option<&GasSolubilityModel>,
    pressure_pascal: f64,
    ionic_strength: f64,
    temperature_kelvin: f64,
) -> ChemistryResult<f64> {
    let Some(model) = model else {
        return Ok(0.0);
    };
    match model {
        GasSolubilityModel::Henry {
            henry_mol_per_bucket_pascal,
            temperature_kelvin: reference_temperature,
            salting_out_coefficient,
            transfer_coefficient_per_tick: _,
            estimated: _,
        } => {
            if !pressure_pascal.is_finite() || pressure_pascal < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "gas pressure must be non-negative and finite".to_string(),
                ));
            }
            if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "temperature must be positive and finite for gas solubility".to_string(),
                ));
            }
            let temperature_factor = (*reference_temperature / temperature_kelvin).sqrt();
            let activity_factor = (-salting_out_coefficient * ionic_strength).exp();
            let dissolved = pressure_pascal
                * henry_mol_per_bucket_pascal
                * temperature_factor
                * activity_factor;
            if !dissolved.is_finite() || dissolved < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "gas dissolved amount must be non-negative and finite".to_string(),
                ));
            }
            Ok(dissolved)
        }
    }
}

fn gas_transfer_coefficient_per_tick(model: &GasSolubilityModel) -> f64 {
    match model {
        GasSolubilityModel::Henry {
            transfer_coefficient_per_tick,
            ..
        } => *transfer_coefficient_per_tick,
    }
}

fn distribute_condensed_amount(
    preference: LiquidPhasePreference,
    aqueous_limit: Option<f64>,
    organic_limit: Option<f64>,
    can_precipitate: bool,
    organic_solvent: Option<SubstanceIndex>,
    amount: f64,
    substance_id: &SubstanceId,
    target: &mut ComponentPhaseAmounts,
) -> ChemistryResult<()> {
    if amount <= 0.0 {
        return Ok(());
    }
    let mut remaining = amount;
    let first = preferred_phase(preference);
    let second = match first {
        MixturePhase::Aqueous => MixturePhase::Organic,
        MixturePhase::Organic => MixturePhase::Aqueous,
        MixturePhase::Gas | MixturePhase::Solid => unreachable!("preferred phase must be liquid"),
    };
    fill_liquid_phase(
        first,
        liquid_limit(first, aqueous_limit, organic_limit),
        organic_solvent,
        &mut remaining,
        target,
    );
    fill_liquid_phase(
        second,
        liquid_limit(second, aqueous_limit, organic_limit),
        organic_solvent,
        &mut remaining,
        target,
    );
    if remaining <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
        return Ok(());
    }
    if can_precipitate {
        target.solid_mol_per_bucket += remaining;
        return Ok(());
    }
    Err(ChemistryError::InvalidMixtureState(format!(
        "substance '{substance_id}' has {remaining} mol/bucket that cannot fit any liquid phase and cannot precipitate"
    )))
}

fn liquid_limit(
    phase: MixturePhase,
    aqueous_limit: Option<f64>,
    organic_limit: Option<f64>,
) -> Option<f64> {
    match phase {
        MixturePhase::Aqueous => aqueous_limit,
        MixturePhase::Organic => organic_limit,
        MixturePhase::Gas | MixturePhase::Solid => None,
    }
}

fn fill_liquid_phase(
    phase: MixturePhase,
    limit: Option<f64>,
    organic_solvent: Option<SubstanceIndex>,
    remaining: &mut f64,
    target: &mut ComponentPhaseAmounts,
) {
    if !phase.is_liquid() || *remaining <= 0.0 {
        return;
    }
    let existing = match phase {
        MixturePhase::Aqueous => target.aqueous_mol_per_bucket,
        MixturePhase::Organic => organic_solvent
            .and_then(|solvent| {
                target
                    .organic_mol_per_bucket_by_solvent
                    .get(&solvent)
                    .copied()
            })
            .unwrap_or(0.0),
        MixturePhase::Gas | MixturePhase::Solid => 0.0,
    };
    let capacity = limit
        .map(|limit| (limit - existing).max(0.0))
        .unwrap_or(*remaining);
    let moved = capacity.min(*remaining);
    match phase {
        MixturePhase::Aqueous => target.aqueous_mol_per_bucket += moved,
        MixturePhase::Organic => {
            let Some(solvent) = organic_solvent else {
                return;
            };
            *target
                .organic_mol_per_bucket_by_solvent
                .entry(solvent)
                .or_insert(0.0) += moved;
        }
        MixturePhase::Gas | MixturePhase::Solid => {}
    }
    *remaining -= moved;
}

fn normalize_phase_amounts(
    mut phase_amounts: ComponentPhaseAmounts,
    total_amount: f64,
) -> ComponentPhaseAmounts {
    phase_amounts.aqueous_mol_per_bucket /= total_amount;
    phase_amounts.gas_mol_per_bucket /= total_amount;
    phase_amounts.solid_mol_per_bucket /= total_amount;
    for concentration in phase_amounts.organic_mol_per_bucket_by_solvent.values_mut() {
        *concentration /= total_amount;
    }
    phase_amounts
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
    use crate::chemistry::registry::{
        ChemistryRegistryBuilder, GasSolubilityModel, SolventMiscibility,
    };
    use crate::chemistry::substance::{LiquidPhasePreference, Substance, SubstancePhaseProperties};

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

    fn gas_registry() -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_unlimited()),
            )
            .substance(
                Substance::new("destroy:oxygen", 0, 32.0, 1_140.0, 90.0, 29.4, 6_820.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: None,
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: false,
                        can_form_liquid_phase: false,
                    }),
            )
            .substance(
                Substance::new("destroy:unknown_gas", 0, 10.0, 1_000.0, 100.0, 20.0, 100.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: None,
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: false,
                        can_form_liquid_phase: false,
                    }),
            )
            .substance(
                Substance::new("destroy:sodium", 1, 23.0, 23_000.0, 10_000.0, 10.0, 0.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_unlimited()),
            )
            .gas_solubility(
                "destroy:oxygen",
                GasSolubilityModel::Henry {
                    henry_mol_per_bucket_pascal: 1.0e-8,
                    temperature_kelvin: 298.0,
                    salting_out_coefficient: 0.5,
                    transfer_coefficient_per_tick: 0.25,
                    estimated: false,
                },
            )
            .build()
            .unwrap()
    }

    fn melting_registry() -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:wax", 0, 100.0, 100_000.0, 600.0, 10.0, 1_000.0)
                    .with_melting_point_kelvin(300.0)
                    .with_fusion_heat_j_per_mol(100.0)
                    .with_solid_density_grams_per_bucket(120_000.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Organic,
                        aqueous_solubility_mol_per_bucket: Some(0.0),
                        organic_solubility_mol_per_bucket: None,
                        can_precipitate: true,
                        can_form_liquid_phase: true,
                    }),
            )
            .build()
            .unwrap()
    }

    fn phase_registry() -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_unlimited()),
            )
            .substance(
                Substance::new("destroy:salt", 0, 58.0, 58_000.0, 1_000.0, 80.0, 20_000.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: Some(0.5),
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: true,
                        can_form_liquid_phase: false,
                    }),
            )
            .substance(
                Substance::new("destroy:oil", 0, 100.0, 80_000.0, 450.0, 120.0, 20_000.0)
                    .with_phase_properties(SubstancePhaseProperties::organic_unlimited(0.05)),
            )
            .substance(
                Substance::new(
                    "destroy:chloroform",
                    0,
                    119.0,
                    119_000.0,
                    334.0,
                    114.0,
                    20_000.0,
                )
                .with_phase_properties(SubstancePhaseProperties::organic_unlimited(0.01)),
            )
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
            aqueous_mol_per_bucket: 1.0,
            organic_mol_per_bucket_by_solvent: BTreeMap::new(),
            gas_mol_per_bucket: f64::INFINITY,
            solid_mol_per_bucket: 0.0,
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
        assert_eq!(
            mixture.substances().cloned().collect::<Vec<_>>(),
            vec![water]
        );
    }

    #[test]
    fn phase_amounts_sum_to_public_concentration() {
        let registry = phase_registry();
        let water: SubstanceId = "destroy:water".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, water.clone(), 1.0)
            .unwrap();
        mixture
            .move_between_phases(
                &registry,
                water.clone(),
                MixturePhase::Aqueous,
                MixturePhase::Gas,
                0.25,
            )
            .unwrap();

        assert_eq!(mixture.concentration_of(&water), 1.0);
        assert_eq!(
            mixture.concentration_in_phase(&water, MixturePhase::Aqueous),
            0.75
        );
        assert_eq!(
            mixture.concentration_in_phase(&water, MixturePhase::Gas),
            0.25
        );
        assert_eq!(mixture.gaseous_fraction_of(&water), 0.25);
    }

    #[test]
    fn solubility_excess_precipitates_and_can_redissolve() {
        let registry = phase_registry();
        let salt: SubstanceId = "destroy:salt".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture.add_substance(&registry, salt.clone(), 1.0).unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&salt, MixturePhase::Aqueous),
            0.5
        );
        assert_eq!(
            mixture.concentration_in_phase(&salt, MixturePhase::Solid),
            0.5
        );

        mixture
            .change_concentration(&registry, &salt, -0.5)
            .unwrap();

        assert_eq!(mixture.concentration_of(&salt), 0.5);
        assert_eq!(
            mixture.concentration_in_phase(&salt, MixturePhase::Aqueous),
            0.5
        );
        assert_eq!(
            mixture.concentration_in_phase(&salt, MixturePhase::Solid),
            0.0
        );
    }

    #[test]
    fn neutral_organic_substance_prefers_organic_phase() {
        let registry = phase_registry();
        let oil: SubstanceId = "destroy:oil".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture.add_substance(&registry, oil.clone(), 1.0).unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&oil, MixturePhase::Organic),
            1.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&oil, MixturePhase::Aqueous),
            0.0
        );
    }

    #[test]
    fn organic_solvents_create_distinct_liquid_phases() {
        let registry = phase_registry();
        let oil: SubstanceId = "destroy:oil".into();
        let chloroform: SubstanceId = "destroy:chloroform".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture.add_substance(&registry, oil.clone(), 1.0).unwrap();
        mixture
            .add_substance(&registry, chloroform.clone(), 2.0)
            .unwrap();

        assert_eq!(
            mixture
                .concentration_in_organic_solvent(&registry, &oil, &oil)
                .unwrap(),
            1.0
        );
        assert_eq!(
            mixture
                .concentration_in_organic_solvent(&registry, &oil, &chloroform)
                .unwrap(),
            0.0
        );
        assert_eq!(
            mixture
                .concentration_in_organic_solvent(&registry, &chloroform, &chloroform)
                .unwrap(),
            2.0
        );
        assert_eq!(
            mixture
                .organic_phase_amounts_of(&registry, &oil)
                .unwrap()
                .len(),
            1
        );
        assert_eq!(
            mixture
                .organic_phase_amounts_of(&registry, &chloroform)
                .unwrap()
                .len(),
            1
        );
    }

    #[test]
    fn aggregate_state_controls_initial_phase() {
        let registry = melting_registry();
        let wax: SubstanceId = "destroy:wax".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture.add_substance(&registry, wax.clone(), 1.0).unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&wax, MixturePhase::Solid),
            1.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&wax, MixturePhase::Organic),
            0.0
        );
    }

    #[test]
    fn heating_solid_consumes_fusion_heat_before_liquid_phase() {
        let registry = melting_registry();
        let wax: SubstanceId = "destroy:wax".into();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, wax.clone(), 1.0).unwrap();

        mixture.heat(&registry, 70.0).unwrap();

        assert_eq!(mixture.temperature_kelvin(), 300.0);
        assert!((mixture.concentration_in_phase(&wax, MixturePhase::Solid) - 0.5).abs() < 1.0e-9);
        assert!((mixture.concentration_in_phase(&wax, MixturePhase::Organic) - 0.5).abs() < 1.0e-9);
    }

    #[test]
    fn gas_solubility_uses_pressure_and_activity() {
        let registry = gas_registry();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let sodium: SubstanceId = "destroy:sodium".into();
        let mut pure = Mixture::new(298.0).unwrap();
        pure.add_substance(&registry, oxygen.clone(), 1.0).unwrap();
        pure.transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();
        let pure_dissolved = pure.concentration_in_phase(&oxygen, MixturePhase::Aqueous);

        let mut salty = Mixture::new(298.0).unwrap();
        salty.add_substance(&registry, sodium, 1.0).unwrap();
        salty.add_substance(&registry, oxygen.clone(), 1.0).unwrap();
        salty
            .transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();
        let salty_dissolved = salty.concentration_in_phase(&oxygen, MixturePhase::Aqueous);

        assert!(pure.gas_pressure_pascal() > 0.0);
        assert!(pure_dissolved > 0.0);
        assert!(salty_dissolved < pure_dissolved);
    }

    #[test]
    fn gas_without_solubility_data_does_not_enter_solution() {
        let registry = gas_registry();
        let gas: SubstanceId = "destroy:unknown_gas".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture.add_substance(&registry, gas.clone(), 1.0).unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&gas, MixturePhase::Aqueous),
            0.0
        );
        assert_eq!(mixture.concentration_in_phase(&gas, MixturePhase::Gas), 1.0);
    }

    #[test]
    fn gas_transfer_moves_toward_henry_limit_gradually() {
        let registry = gas_registry();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, oxygen.clone(), 1.0)
            .unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&oxygen, MixturePhase::Aqueous),
            0.0
        );

        let first_delta = mixture
            .transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();
        let after_first = mixture.concentration_in_phase(&oxygen, MixturePhase::Aqueous);
        let second_delta = mixture
            .transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();
        let after_second = mixture.concentration_in_phase(&oxygen, MixturePhase::Aqueous);

        assert!(first_delta > 0.0);
        assert!(second_delta > 0.0);
        assert!(after_first > 0.0);
        assert!(after_second > after_first);
        assert!(mixture.concentration_in_phase(&oxygen, MixturePhase::Gas) > 0.0);
    }

    #[test]
    fn fully_miscible_solvents_share_one_liquid_phase() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_unlimited()),
            )
            .substance(
                Substance::new("destroy:ethanol", 0, 46.0, 46_000.0, 351.0, 110.0, 20_000.0)
                    .with_phase_properties(SubstancePhaseProperties::organic_unlimited(1.0)),
            )
            .solvent_miscibility(
                "destroy:water",
                "destroy:ethanol",
                SolventMiscibility::FullyMiscible,
            )
            .build()
            .unwrap();
        let ethanol: SubstanceId = "destroy:ethanol".into();
        let water: SubstanceId = "destroy:water".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, water.clone(), 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, ethanol.clone(), 1.0)
            .unwrap();

        let organic_amounts = mixture
            .organic_phase_amounts_of(&registry, &ethanol)
            .unwrap();
        assert_eq!(organic_amounts.len(), 1);
        assert_eq!(organic_amounts[0].concentration_mol_per_bucket, 1.0);
        assert!(organic_amounts[0].solvent_id == water || organic_amounts[0].solvent_id == ethanol);
    }

    #[test]
    fn partial_miscibility_keeps_large_solvent_amounts_separate() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new(
                    "destroy:solvent_a",
                    0,
                    80.0,
                    80_000.0,
                    420.0,
                    80.0,
                    20_000.0,
                )
                .with_phase_properties(SubstancePhaseProperties::organic_unlimited(0.0)),
            )
            .substance(
                Substance::new(
                    "destroy:solvent_b",
                    0,
                    90.0,
                    90_000.0,
                    430.0,
                    90.0,
                    20_000.0,
                )
                .with_phase_properties(SubstancePhaseProperties::organic_unlimited(0.0)),
            )
            .solvent_miscibility(
                "destroy:solvent_a",
                "destroy:solvent_b",
                SolventMiscibility::PartiallyMiscible {
                    limit_mol_per_bucket: 0.2,
                },
            )
            .build()
            .unwrap();
        let solvent_a: SubstanceId = "destroy:solvent_a".into();
        let solvent_b: SubstanceId = "destroy:solvent_b".into();

        let mut small = Mixture::new(298.0).unwrap();
        small
            .add_substance(&registry, solvent_a.clone(), 1.0)
            .unwrap();
        small
            .add_substance(&registry, solvent_b.clone(), 0.1)
            .unwrap();
        assert_eq!(
            small
                .organic_phase_amounts_of(&registry, &solvent_b)
                .unwrap()
                .len(),
            1
        );

        let mut large = Mixture::new(298.0).unwrap();
        large
            .add_substance(&registry, solvent_a.clone(), 1.0)
            .unwrap();
        large
            .add_substance(&registry, solvent_b.clone(), 0.5)
            .unwrap();
        assert_eq!(
            large
                .concentration_in_organic_solvent(&registry, &solvent_b, &solvent_b)
                .unwrap(),
            0.5
        );
    }
}
