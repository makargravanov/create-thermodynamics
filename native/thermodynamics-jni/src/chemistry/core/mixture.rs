use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::reaction::GAS_CONSTANT_J_PER_MOL_KELVIN;
use super::registry::{ChemistryRegistry, GasSolubilityModel, SolventMiscibility, SubstanceIndex};
use super::solution::{self, AqueousEquilibriumSystem, SolutionState};
use super::substance::{
    LiquidPhasePreference, SolventRole, Substance, SubstanceAggregateState, SubstanceId,
    SubstanceRepresentation, SubstanceTagId,
};

pub const DEFAULT_TEMPERATURE_KELVIN: f64 = 298.0;
pub const TRACE_CONCENTRATION_MOL_PER_BUCKET: f64 = 1.0 / 512.0 / 512.0;
pub const STANDARD_PRESSURE_PASCAL: f64 = 101_325.0;
pub const DEFAULT_GAS_VOLUME_CUBIC_METERS: f64 = 0.001;

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MixturePhase {
    Aqueous,
    Organic,
    MoltenMetal,
    MoltenSlag,
    Gas,
    SupercriticalFluid,
    Solid,
}

impl MixturePhase {
    pub const ALL: [MixturePhase; 7] = [
        MixturePhase::Aqueous,
        MixturePhase::Organic,
        MixturePhase::MoltenMetal,
        MixturePhase::MoltenSlag,
        MixturePhase::Gas,
        MixturePhase::SupercriticalFluid,
        MixturePhase::Solid,
    ];
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

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LiquidPhaseId(usize);

impl LiquidPhaseId {
    pub fn as_usize(self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct LiquidPhaseSolventAmount {
    pub substance_id: SubstanceId,
    pub concentration_mol_per_bucket: f64,
    pub mole_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LiquidPhaseSnapshot {
    pub id: LiquidPhaseId,
    pub coarse_phase: MixturePhase,
    pub representative_solvent_id: SubstanceId,
    pub total_solvent_mol_per_bucket: f64,
    pub solvents: Vec<LiquidPhaseSolventAmount>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LiquidPhaseSubstanceAmount {
    pub phase_id: LiquidPhaseId,
    pub coarse_phase: MixturePhase,
    pub representative_solvent_id: SubstanceId,
    pub concentration_mol_per_bucket: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LiquidPhaseIonicStrength {
    pub phase_id: LiquidPhaseId,
    pub coarse_phase: MixturePhase,
    pub representative_solvent_id: SubstanceId,
    pub ionic_strength_mol_per_bucket: f64,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SolidPhaseId(usize);

impl SolidPhaseId {
    pub fn as_usize(self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SolidPhaseSnapshot {
    pub id: SolidPhaseId,
    pub representative_substance_id: SubstanceId,
    pub concentration_mol_per_bucket: f64,
}

#[derive(Debug, Clone)]
struct LiquidPhaseState {
    id: LiquidPhaseId,
    solvents: BTreeSet<SubstanceIndex>,
    solvent_amounts: BTreeMap<SubstanceIndex, f64>,
    representative_solvent: SubstanceIndex,
    coarse_phase: MixturePhase,
}

#[derive(Debug, Clone)]
pub struct Mixture {
    temperature_kelvin: f64,
    gas_volume_cubic_meters: f64,
    components: Vec<MixtureComponent>,
    positions_by_substance: Vec<Option<usize>>,
}

#[derive(Debug, Clone)]
struct MixtureComponent {
    substance: SubstanceIndex,
    substance_id: SubstanceId,
    aqueous_mol_per_bucket: f64,
    organic_mol_per_bucket_by_solvent: BTreeMap<SubstanceIndex, f64>,
    molten_mol_per_bucket_by_phase: BTreeMap<CondensedPhaseKey, f64>,
    gas_mol_per_bucket: f64,
    supercritical_mol_per_bucket: f64,
    solid_mol_per_bucket_by_phase: BTreeMap<SolidPhaseKey, f64>,
}

#[derive(Debug, Clone, Default)]
struct ComponentPhaseAmounts {
    aqueous_mol_per_bucket: f64,
    organic_mol_per_bucket_by_solvent: BTreeMap<SubstanceIndex, f64>,
    molten_mol_per_bucket_by_phase: BTreeMap<CondensedPhaseKey, f64>,
    gas_mol_per_bucket: f64,
    supercritical_mol_per_bucket: f64,
    solid_mol_per_bucket_by_phase: BTreeMap<SolidPhaseKey, f64>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct CondensedPhaseKey {
    coarse_phase: MixturePhase,
    anchor: SubstanceIndex,
}

impl CondensedPhaseKey {
    fn new(coarse_phase: MixturePhase, anchor: SubstanceIndex) -> Self {
        Self {
            coarse_phase,
            anchor,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct SolidPhaseKey {
    anchor: SubstanceIndex,
}

impl SolidPhaseKey {
    fn new(anchor: SubstanceIndex) -> Self {
        Self { anchor }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct MixtureCheckpoint {
    temperature_kelvin: f64,
    gas_volume_cubic_meters: f64,
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
            gas_volume_cubic_meters: DEFAULT_GAS_VOLUME_CUBIC_METERS,
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

    pub fn gas_volume_cubic_meters(&self) -> f64 {
        self.gas_volume_cubic_meters
    }

    pub fn set_gas_volume_cubic_meters(
        &mut self,
        gas_volume_cubic_meters: f64,
    ) -> ChemistryResult<()> {
        if !gas_volume_cubic_meters.is_finite() || gas_volume_cubic_meters <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "gas volume must be positive and finite".to_string(),
            ));
        }
        self.gas_volume_cubic_meters = gas_volume_cubic_meters;
        Ok(())
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

    pub fn total_concentration_in_phase(&self, phase: MixturePhase) -> f64 {
        self.components
            .iter()
            .map(|component| component.amount_in_phase(phase))
            .sum()
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

    pub fn liquid_phase_count(&self, registry: &ChemistryRegistry) -> ChemistryResult<usize> {
        Ok(self.liquid_phases(registry)?.len())
    }

    pub fn solid_phase_count(&self) -> usize {
        self.solid_phase_amounts_by_anchor().len()
    }

    pub fn solid_phase_snapshots(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Vec<SolidPhaseSnapshot>> {
        self.solid_phase_amounts_by_anchor()
            .into_iter()
            .enumerate()
            .map(|(position, (phase, concentration))| {
                Ok(SolidPhaseSnapshot {
                    id: SolidPhaseId(position),
                    representative_substance_id: registry
                        .substance_by_index(phase.anchor)?
                        .id
                        .clone(),
                    concentration_mol_per_bucket: concentration,
                })
            })
            .collect()
    }

    pub fn liquid_phase_snapshots(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Vec<LiquidPhaseSnapshot>> {
        self.liquid_phases(registry)?
            .into_iter()
            .map(|phase| liquid_phase_snapshot(registry, phase))
            .collect()
    }

    pub fn liquid_phase_amounts_of(
        &self,
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<Vec<LiquidPhaseSubstanceAmount>> {
        self.validate_registry_shape(registry)?;
        let Some(component) = self
            .components
            .iter()
            .find(|component| &component.substance_id == substance_id)
        else {
            return Ok(Vec::new());
        };
        self.liquid_phases(registry)?
            .into_iter()
            .map(|phase| {
                let amount = component.amount_in_liquid_phase(&phase);
                Ok(LiquidPhaseSubstanceAmount {
                    phase_id: phase.id,
                    coarse_phase: phase.coarse_phase,
                    representative_solvent_id: registry
                        .substance_by_index(phase.representative_solvent)?
                        .id
                        .clone(),
                    concentration_mol_per_bucket: amount,
                })
            })
            .collect()
    }

    pub fn gas_pressure_pascal(&self) -> f64 {
        self.pressure_mol_per_bucket() * GAS_CONSTANT_J_PER_MOL_KELVIN * self.temperature_kelvin
            / self.gas_volume_cubic_meters
    }

    pub fn gas_partial_pressure_pascal(&self, substance_id: &SubstanceId) -> f64 {
        self.components
            .iter()
            .find(|component| &component.substance_id == substance_id)
            .map(|component| {
                gas_pressure_for_moles(
                    component
                        .amount_in_phases(&[MixturePhase::Gas, MixturePhase::SupercriticalFluid]),
                    self.temperature_kelvin,
                    self.gas_volume_cubic_meters,
                )
            })
            .unwrap_or(0.0)
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

    pub fn liquid_phase_ionic_strengths(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Vec<LiquidPhaseIonicStrength>> {
        self.validate_registry_shape(registry)?;
        let properties = registry.substance_properties();
        self.liquid_phases(registry)?
            .into_iter()
            .map(|phase| {
                let mut ionic_strength = 0.0;
                for component in &self.components {
                    let charge = properties.charge[component.substance.as_usize()] as f64;
                    if charge == 0.0 {
                        continue;
                    }
                    let amount = component.amount_in_liquid_phase(&phase);
                    ionic_strength += 0.5 * amount * charge * charge;
                }
                if !ionic_strength.is_finite() || ionic_strength < 0.0 {
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "liquid phase ionic strength must be non-negative and finite: {ionic_strength}"
                    )));
                }
                Ok(LiquidPhaseIonicStrength {
                    phase_id: phase.id,
                    coarse_phase: phase.coarse_phase,
                    representative_solvent_id: registry
                        .substance_by_index(phase.representative_solvent)?
                        .id
                        .clone(),
                    ionic_strength_mol_per_bucket: ionic_strength,
                })
            })
            .collect()
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

    pub fn aqueous_equilibrium_systems(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Vec<AqueousEquilibriumSystem>> {
        solution::aqueous_equilibrium_systems(registry, self)
    }

    pub fn equilibrate_solution(&mut self, registry: &ChemistryRegistry) -> ChemistryResult<f64> {
        solution::equilibrate_solution_equilibria(registry, self)
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
            let Some(position) = self.position_of_substance(substance) else {
                continue;
            };
            let gas_pressure_pascal = gas_pressure_for_moles(
                self.components[position]
                    .amount_in_phases(&[MixturePhase::Gas, MixturePhase::SupercriticalFluid]),
                self.temperature_kelvin,
                self.gas_volume_cubic_meters,
            );
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
            let current_dissolved = self.components[position]
                .amount_in_phases(&[MixturePhase::Aqueous, MixturePhase::Organic]);
            let difference = target_dissolved - current_dissolved;
            if difference.abs() <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                continue;
            }
            if difference > 0.0 {
                let available_gas = self.components[position]
                    .amount_in_phases(&[MixturePhase::Gas, MixturePhase::SupercriticalFluid]);
                let moved = (difference * fraction).min(available_gas);
                if moved <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                    continue;
                }
                self.components[position].remove_from_phases(
                    &[MixturePhase::Gas, MixturePhase::SupercriticalFluid],
                    moved,
                )?;
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
                let gas_phase = pressure_phase_for_substance(
                    registry.substance_by_index(substance)?,
                    self.temperature_kelvin,
                    self.gas_pressure_pascal(),
                );
                self.components[position].add_to_phase(gas_phase, moved);
                max_delta = max_delta.max(moved);
            }
        }
        self.equilibrate_phases(registry)?;
        Ok(max_delta)
    }

    pub fn exchange_gases_with_atmosphere(
        &mut self,
        registry: &ChemistryRegistry,
        atmosphere_mole_fractions: &[(SubstanceId, f64)],
        total_pressure_pascal: f64,
        exchange_coefficient_per_tick: f64,
        ticks: f64,
    ) -> ChemistryResult<f64> {
        if !total_pressure_pascal.is_finite() || total_pressure_pascal < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "atmosphere pressure must be non-negative and finite".to_string(),
            ));
        }
        if !exchange_coefficient_per_tick.is_finite() || exchange_coefficient_per_tick < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "gas exchange coefficient must be non-negative and finite".to_string(),
            ));
        }
        if !ticks.is_finite() || ticks < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "gas exchange time must be non-negative and finite".to_string(),
            ));
        }
        if self.temperature_kelvin <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "gas exchange requires positive mixture temperature".to_string(),
            ));
        }
        if exchange_coefficient_per_tick == 0.0 || ticks == 0.0 {
            return self.validate(registry).map(|_| 0.0);
        }
        self.validate_registry_shape(registry)?;
        self.ensure_position_capacity(registry);

        let mut target_pressures = BTreeMap::<SubstanceIndex, f64>::new();
        let mut fraction_sum = 0.0;
        for (substance_id, fraction) in atmosphere_mole_fractions {
            if !fraction.is_finite() || *fraction < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "atmosphere gas fraction must be non-negative and finite".to_string(),
                ));
            }
            fraction_sum += fraction;
            let substance = registry.substance_index(substance_id).ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!(
                    "unknown atmosphere gas substance '{substance_id}'"
                ))
            })?;
            if registry
                .substance_by_index(substance)?
                .aggregate_state_at(self.temperature_kelvin)?
                != SubstanceAggregateState::Gas
            {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "atmosphere substance '{substance_id}' is not a gas at current temperature"
                )));
            }
            if target_pressures
                .insert(substance, total_pressure_pascal * fraction)
                .is_some()
            {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "duplicate atmosphere gas substance '{substance_id}'"
                )));
            }
        }
        if fraction_sum > 1.0 + 1.0e-12 {
            return Err(ChemistryError::InvalidMixtureState(
                "atmosphere gas fractions must not sum above one".to_string(),
            ));
        }

        let mut gas_substances = self
            .components
            .iter()
            .filter(|component| {
                component.amount_in_phases(&[MixturePhase::Gas, MixturePhase::SupercriticalFluid])
                    > 0.0
            })
            .map(|component| component.substance)
            .collect::<BTreeSet<_>>();
        gas_substances.extend(target_pressures.keys().copied());

        let fraction = 1.0 - (-exchange_coefficient_per_tick * ticks).exp();
        let mut max_delta = 0.0_f64;
        for substance in gas_substances {
            let target_pressure = target_pressures.get(&substance).copied().unwrap_or(0.0);
            let target_gas = gas_moles_for_pressure(
                target_pressure,
                self.temperature_kelvin,
                self.gas_volume_cubic_meters,
            )?;
            let current_gas = self.concentration_of_index_in_phases(
                substance,
                &[MixturePhase::Gas, MixturePhase::SupercriticalFluid],
            );
            let delta = (target_gas - current_gas) * fraction;
            if delta.abs() <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                continue;
            }
            max_delta = max_delta.max(delta.abs());
            if delta > 0.0 {
                let substance_id = registry.substance_by_index(substance)?.id.clone();
                self.change_concentration_by_index_unchecked(
                    registry,
                    substance,
                    substance_id,
                    delta,
                    pressure_phase_for_substance(
                        registry.substance_by_index(substance)?,
                        self.temperature_kelvin,
                        total_pressure_pascal,
                    ),
                )?;
            } else if let Some(position) = self.position_of_substance(substance) {
                let removed = (-delta).min(current_gas);
                self.components[position].remove_from_phases(
                    &[MixturePhase::Gas, MixturePhase::SupercriticalFluid],
                    removed,
                )?;
                self.remove_trace_component(registry, substance)?;
            }
        }

        self.equilibrate_phases(registry)?;
        self.validate(registry)?;
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

    pub(crate) fn apply_aqueous_targets_by_index(
        &mut self,
        registry: &ChemistryRegistry,
        targets: &[(SubstanceIndex, f64)],
    ) -> ChemistryResult<f64> {
        self.ensure_position_capacity(registry);
        let mut max_delta = 0.0_f64;
        let mut touched = Vec::with_capacity(targets.len());
        for (substance, target) in targets {
            registry.substance_by_index(*substance)?;
            if !target.is_finite() || *target < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "aqueous target concentration must be non-negative and finite".to_string(),
                ));
            }
            let current =
                self.concentration_of_index_in_phases(*substance, &[MixturePhase::Aqueous]);
            max_delta = max_delta.max((target - current).abs());
            touched.push(*substance);
        }
        for (substance, target) in targets {
            if let Some(position) = self.position_of_substance(*substance) {
                self.components[position].aqueous_mol_per_bucket = *target;
            } else {
                let substance_data = registry.substance_by_index(*substance)?;
                let trace_threshold = if substance_data.charge == 0 {
                    TRACE_CONCENTRATION_MOL_PER_BUCKET
                } else {
                    1.0e-14
                };
                if *target <= trace_threshold {
                    continue;
                }
                let substance_id = substance_data.id.clone();
                self.insert_component(*substance, substance_id, *target, MixturePhase::Aqueous);
            }
        }
        for substance in touched {
            self.remove_trace_component(registry, substance)?;
        }
        self.equilibrate_phases(registry)?;
        self.validate(registry)?;
        Ok(max_delta)
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
        self.remove_trace_component(registry, substance)?;
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

    pub(crate) fn extract_from_phase_by_index(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        phase: MixturePhase,
        amount: f64,
    ) -> ChemistryResult<f64> {
        if amount <= 0.0 {
            return Ok(0.0);
        }
        let position = match self.position_of_substance(substance) {
            Some(p) => p,
            None => return Ok(0.0),
        };
        let available = self.components[position].amount_in_phase(phase);
        let take = available.min(amount);
        if take <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Ok(0.0);
        }
        self.components[position].remove_from_phase(phase, take)?;
        let total = self.components[position].total_concentration();
        if total <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            self.remove_component(substance);
        }
        self.validate(registry)?;
        Ok(take)
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
                self.remove_trace_component(registry, *substance)?;
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

    pub(crate) fn apply_phase_amount_deltas_by_index(
        &mut self,
        registry: &ChemistryRegistry,
        deltas: &[(SubstanceIndex, MixturePhase, f64)],
    ) -> ChemistryResult<f64> {
        self.ensure_position_capacity(registry);
        let mut merged = BTreeMap::<(SubstanceIndex, MixturePhase), f64>::new();
        for (substance, phase, delta) in deltas {
            registry.substance_by_index(*substance)?;
            if !delta.is_finite() {
                return Err(ChemistryError::InvalidMixtureState(
                    "phase amount delta must be finite".to_string(),
                ));
            }
            *merged.entry((*substance, *phase)).or_insert(0.0) += *delta;
        }
        let mut max_delta = 0.0_f64;
        for ((substance, phase), delta) in &merged {
            max_delta = max_delta.max(delta.abs());
            if *delta < 0.0 {
                let available = self.concentration_of_index_in_phases(*substance, &[*phase]);
                if available + *delta < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
                    let substance_id = &registry.substance_by_index(*substance)?.id;
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "substance '{substance_id}' would become negative in phase {phase:?}: {}",
                        available + *delta
                    )));
                }
            }
        }
        for ((substance, phase), delta) in merged {
            if delta > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                let substance_data = registry.substance_by_index(substance)?;
                if let Some(position) = self.position_of_substance(substance) {
                    self.components[position].add_to_phase(phase, delta);
                } else {
                    self.insert_component(substance, substance_data.id.clone(), delta, phase);
                }
            } else if delta < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
                if let Some(position) = self.position_of_substance(substance) {
                    self.components[position].remove_from_phase(phase, -delta)?;
                    self.remove_trace_component(registry, substance)?;
                }
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
            gas_volume_cubic_meters: self.gas_volume_cubic_meters,
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
        self.gas_volume_cubic_meters = checkpoint.gas_volume_cubic_meters;
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
            liquid_buckets += component.amount_in_phases(&liquid_phases())
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
                for (phase, concentration) in &component.solid_mol_per_bucket_by_phase {
                    *entry
                        .solid_mol_per_bucket_by_phase
                        .entry(*phase)
                        .or_insert(0.0) += concentration * amount;
                }
                for (solvent, concentration) in &component.organic_mol_per_bucket_by_solvent {
                    *entry
                        .organic_mol_per_bucket_by_solvent
                        .entry(*solvent)
                        .or_insert(0.0) += concentration * amount;
                }
                for (phase, concentration) in &component.molten_mol_per_bucket_by_phase {
                    *entry
                        .molten_mol_per_bucket_by_phase
                        .entry(*phase)
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
        result.gas_volume_cubic_meters = mixtures
            .iter()
            .map(|(mixture, amount)| mixture.gas_volume_cubic_meters * amount)
            .sum::<f64>()
            / total_amount;
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
        self.redistribute_condensed_phases(registry)?;
        let vapor_liquid_delta = self.equilibrate_vapor_liquid(registry, 1.0)?;
        if vapor_liquid_delta > TRACE_CONCENTRATION_MOL_PER_BUCKET {
            self.redistribute_condensed_phases(registry)?;
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

    fn redistribute_condensed_phases(
        &mut self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<()> {
        let liquid_phases = self.liquid_phases(registry)?;
        let pressure_pascal = self.gas_pressure_pascal();
        for position in 0..self.components.len() {
            let substance = registry.substance_by_index(self.components[position].substance)?;
            let total = self.components[position].total_concentration();
            let mut next = ComponentPhaseAmounts::default();
            if substance_is_supercritical(substance, self.temperature_kelvin, pressure_pascal) {
                next.supercritical_mol_per_bucket = total;
            } else {
                match substance.aggregate_state_at(self.temperature_kelvin)? {
                    SubstanceAggregateState::Gas => {
                        let current_gas = self.components[position].amount_in_phases(&[
                            MixturePhase::Gas,
                            MixturePhase::SupercriticalFluid,
                        ]);
                        let dissolved = (total - current_gas).max(0.0);
                        next.gas_mol_per_bucket = total - dissolved;
                        if dissolved > 0.0 {
                            distribute_condensed_amount(
                                registry,
                                self.components[position].substance,
                                substance,
                                &liquid_phases,
                                substance.phase_properties.can_precipitate,
                                dissolved,
                                &mut next,
                            )?;
                        }
                    }
                    SubstanceAggregateState::Solid => {
                        let dissolved =
                            condensed_solubility_capacity(registry, substance, &liquid_phases)
                                .map(|limit| total.min(limit))
                                .unwrap_or(0.0);
                        if dissolved > 0.0 {
                            distribute_condensed_amount(
                                registry,
                                self.components[position].substance,
                                substance,
                                &liquid_phases,
                                substance.phase_properties.can_precipitate,
                                dissolved,
                                &mut next,
                            )?;
                        }
                        let remaining = total - dissolved;
                        if remaining > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                            if substance.phase_properties.can_precipitate {
                                *next
                                    .solid_mol_per_bucket_by_phase
                                    .entry(SolidPhaseKey::new(self.components[position].substance))
                                    .or_insert(0.0) += remaining;
                            } else {
                                return Err(ChemistryError::InvalidMixtureState(format!(
                                "substance '{}' is solid at current temperature but cannot precipitate",
                                substance.id
                            )));
                            }
                        }
                    }
                    SubstanceAggregateState::Liquid => {
                        let current_gas = if substance
                            .vapor_pressure_pascal(self.temperature_kelvin)?
                            .is_some()
                        {
                            self.components[position].amount_in_phases(&[
                                MixturePhase::Gas,
                                MixturePhase::SupercriticalFluid,
                            ])
                        } else {
                            0.0
                        };
                        let condensed = (total - current_gas).max(0.0);
                        next.gas_mol_per_bucket = current_gas.min(total);
                        distribute_condensed_amount(
                            registry,
                            self.components[position].substance,
                            substance,
                            &liquid_phases,
                            substance.phase_properties.can_precipitate,
                            condensed,
                            &mut next,
                        )?;
                    }
                }
            }
            self.components[position].set_phase_amounts(next);
        }
        Ok(())
    }

    pub fn equilibrate_vapor_liquid(
        &mut self,
        registry: &ChemistryRegistry,
        relaxation: f64,
    ) -> ChemistryResult<f64> {
        self.validate_registry_shape(registry)?;
        if self.temperature_kelvin <= 0.0 {
            return self.validate(registry).map(|_| 0.0);
        }
        let mut max_delta = 0.0_f64;
        let substances = self
            .components
            .iter()
            .map(|component| component.substance)
            .collect::<Vec<_>>();
        for substance_index in substances {
            let substance = registry.substance_by_index(substance_index)?;
            if substance.charge != 0 {
                continue;
            }
            if substance.critical_temperature_kelvin.is_some()
                != substance.critical_pressure_pascal.is_some()
            {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "substance '{}' has partial critical point data (one of T_c/P_c set without the other)",
                    substance.id
                )));
            }
            if let (Some(critical_temperature), Some(critical_pressure)) = (
                substance.critical_temperature_kelvin,
                substance.critical_pressure_pascal,
            ) {
                if self.temperature_kelvin >= critical_temperature {
                    let total = self.concentration_of_index(substance_index);
                    let pressure_if_fluid = gas_pressure_for_moles(
                        total,
                        self.temperature_kelvin,
                        self.gas_volume_cubic_meters,
                    );
                    let target_phase = if pressure_if_fluid >= critical_pressure {
                        MixturePhase::SupercriticalFluid
                    } else {
                        MixturePhase::Gas
                    };
                    let current_target =
                        self.concentration_of_index_in_phases(substance_index, &[target_phase]);
                    let delta = (total - current_target).abs();
                    if delta > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                        self.move_all_to_pressure_phase(substance_index, target_phase)?;
                        max_delta = max_delta.max(delta);
                    }
                    continue;
                }
            }
            if substance.aggregate_state_at(self.temperature_kelvin)?
                == SubstanceAggregateState::Solid
            {
                continue;
            }
            let Some(vapor_pressure) = substance.vapor_pressure_pascal(self.temperature_kelvin)?
            else {
                continue;
            };
            let target_gas = gas_moles_for_pressure(
                vapor_pressure,
                self.temperature_kelvin,
                self.gas_volume_cubic_meters,
            )?;
            let total = self.concentration_of_index(substance_index);
            let current_gas = self.concentration_of_index_in_phases(
                substance_index,
                &[MixturePhase::Gas, MixturePhase::SupercriticalFluid],
            );
            let current_liquid = self.liquid_concentration_of_index(substance_index);
            let desired_gas = total.min(target_gas);
            let heat_capacity = self.volumetric_heat_capacity_j_per_bucket_kelvin(registry)?;
            if heat_capacity <= 0.0 {
                continue;
            }
            let latent_heat = registry.substance_properties().latent_heat_j_per_mol
                [substance_index.as_usize()];
            if latent_heat <= 0.0 {
                continue;
            }
            if current_gas > desired_gas + TRACE_CONCENTRATION_MOL_PER_BUCKET {
                let condensed = (current_gas - desired_gas) * relaxation;
                if condensed > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                    self.move_gas_to_preferred_liquid(registry, substance_index, condensed)?;
                    self.temperature_kelvin += condensed * latent_heat / heat_capacity;
                    max_delta = max_delta.max(condensed);
                }
            } else if current_gas + TRACE_CONCENTRATION_MOL_PER_BUCKET < desired_gas
                && current_liquid > TRACE_CONCENTRATION_MOL_PER_BUCKET
            {
                let max_evaporable =
                    ((self.temperature_kelvin - 1.0).max(0.0)) * heat_capacity / latent_heat;
                let evaporated = ((desired_gas - current_gas).min(current_liquid).min(max_evaporable))
                    * relaxation;
                if evaporated > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                    self.move_liquid_to_gas(substance_index, evaporated)?;
                    self.temperature_kelvin -= evaporated * latent_heat / heat_capacity;
                    max_delta = max_delta.max(evaporated);
                }
            }
        }
        self.validate(registry)?;
        Ok(max_delta)
    }

    pub fn validate(&self, registry: &ChemistryRegistry) -> ChemistryResult<()> {
        if !self.temperature_kelvin.is_finite() || self.temperature_kelvin < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "temperature must be non-negative and finite".to_string(),
            ));
        }
        if !self.gas_volume_cubic_meters.is_finite() || self.gas_volume_cubic_meters <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "gas volume must be positive and finite".to_string(),
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
            molten_mol_per_bucket_by_phase: BTreeMap::new(),
            gas_mol_per_bucket: 0.0,
            supercritical_mol_per_bucket: 0.0,
            solid_mol_per_bucket_by_phase: BTreeMap::new(),
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
            molten_mol_per_bucket_by_phase: phase_amounts.molten_mol_per_bucket_by_phase,
            gas_mol_per_bucket: phase_amounts.gas_mol_per_bucket,
            supercritical_mol_per_bucket: phase_amounts.supercritical_mol_per_bucket,
            solid_mol_per_bucket_by_phase: phase_amounts.solid_mol_per_bucket_by_phase,
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
        self.concentration_of_index_in_phases(substance, &liquid_phases())
    }

    fn move_liquid_to_gas(
        &mut self,
        substance: SubstanceIndex,
        amount: f64,
    ) -> ChemistryResult<()> {
        if let Some(position) = self.position_of_substance(substance) {
            self.components[position].remove_from_phases(&liquid_phases(), amount)?;
            self.components[position].add_to_phase(MixturePhase::Gas, amount);
        }
        Ok(())
    }

    fn move_all_to_pressure_phase(
        &mut self,
        substance: SubstanceIndex,
        target_phase: MixturePhase,
    ) -> ChemistryResult<()> {
        if !matches!(
            target_phase,
            MixturePhase::Gas | MixturePhase::SupercriticalFluid
        ) {
            return Err(ChemistryError::InvalidMixtureState(
                "pressure phase target must be gas or supercritical fluid".to_string(),
            ));
        }
        if let Some(position) = self.position_of_substance(substance) {
            let amount = self.components[position].total_concentration();
            self.components[position].clear_phase_amounts();
            self.components[position].add_to_phase(target_phase, amount);
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
            self.components[position]
                .molten_mol_per_bucket_by_phase
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
            self.components[position].remove_from_phases(
                &[MixturePhase::Gas, MixturePhase::SupercriticalFluid],
                amount,
            )?;
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
            self.components[position].remove_from_phases(&liquid_phases(), amount)?;
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
            self.components[position]
                .molten_mol_per_bucket_by_phase
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

    fn remove_trace_component(
        &mut self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
    ) -> ChemistryResult<()> {
        let trace_threshold = if registry.substance_by_index(substance)?.charge == 0 {
            TRACE_CONCENTRATION_MOL_PER_BUCKET
        } else {
            1.0e-14
        };
        if self.concentration_of_index(substance) <= trace_threshold {
            self.remove_component(substance);
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
            for solvent in component.organic_mol_per_bucket_by_solvent.keys() {
                if solvent.as_usize() >= registry.substance_count() {
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "mixture component '{}' uses solvent index {} outside registry",
                        component.substance_id,
                        solvent.as_usize()
                    )));
                }
            }
            for phase in component.molten_mol_per_bucket_by_phase.keys() {
                if phase.anchor.as_usize() >= registry.substance_count() {
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "mixture component '{}' uses condensed phase anchor {} outside registry",
                        component.substance_id,
                        phase.anchor.as_usize()
                    )));
                }
            }
            for phase in component.solid_mol_per_bucket_by_phase.keys() {
                if phase.anchor.as_usize() >= registry.substance_count() {
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "mixture component '{}' uses solid phase anchor {} outside registry",
                        component.substance_id,
                        phase.anchor.as_usize()
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
        if phase_uses_concrete_liquid_anchor(phase) {
            let substance_id = registry.substance_by_index(substance)?.id.clone();
            let solvent = self
                .liquid_anchor_for_substance(registry, substance, phase)?
                .ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!(
                    "substance '{}' cannot enter phase {:?} because no concrete liquid phase is available",
                    substance_id, phase
                ))
            })?;
            self.components[position].add_to_phase_for_solvent(phase, solvent, amount);
        } else {
            self.components[position].add_to_phase(phase, amount);
        }
        Ok(())
    }

    fn liquid_anchor_for_substance(
        &self,
        registry: &ChemistryRegistry,
        substance: SubstanceIndex,
        phase: MixturePhase,
    ) -> ChemistryResult<Option<SubstanceIndex>> {
        let substance_data = registry.substance_by_index(substance)?;
        if substance_can_anchor_liquid_phase(substance_data)
            && preferred_phase(substance_data.phase_properties.preferred_liquid_phase) == phase
        {
            Ok(Some(substance))
        } else {
            self.first_available_liquid_anchor(registry, phase, &self.solvent_clusters(registry)?)
        }
    }

    fn first_available_liquid_anchor(
        &self,
        registry: &ChemistryRegistry,
        phase: MixturePhase,
        solvent_clusters: &[BTreeSet<SubstanceIndex>],
    ) -> ChemistryResult<Option<SubstanceIndex>> {
        for cluster in solvent_clusters {
            for solvent in cluster {
                let substance = registry.substance_by_index(*solvent)?;
                if preferred_phase(substance.phase_properties.preferred_liquid_phase) == phase {
                    return Ok(Some(*solvent));
                }
            }
        }
        for component in &self.components {
            let substance = registry.substance_by_index(component.substance)?;
            if substance_can_anchor_liquid_phase(substance)
                && preferred_phase(substance.phase_properties.preferred_liquid_phase) == phase
                && component.total_concentration() > TRACE_CONCENTRATION_MOL_PER_BUCKET
            {
                return Ok(Some(component.substance));
            }
        }
        Ok(None)
    }

    fn liquid_phases(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Vec<LiquidPhaseState>> {
        self.solvent_clusters(registry)?
            .into_iter()
            .enumerate()
            .map(|(index, solvents)| {
                let representative_solvent = solvents.iter().next().copied().ok_or_else(|| {
                    ChemistryError::InvalidMixtureState(
                        "liquid phase must contain at least one solvent".to_string(),
                    )
                })?;
                let mut solvent_amounts = BTreeMap::new();
                for solvent in &solvents {
                    let amount = self.concentration_of_index(*solvent);
                    if amount > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                        solvent_amounts.insert(*solvent, amount);
                    }
                }
                let coarse_phase =
                    classify_liquid_phase(registry, &solvent_amounts, representative_solvent)?;
                Ok(LiquidPhaseState {
                    id: LiquidPhaseId(index),
                    solvents,
                    solvent_amounts,
                    representative_solvent,
                    coarse_phase,
                })
            })
            .collect()
    }

    fn solvent_clusters(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<Vec<BTreeSet<SubstanceIndex>>> {
        let solvents = self
            .components
            .iter()
            .filter_map(|component| {
                let substance = registry.substance_by_index(component.substance).ok()?;
                if !substance_can_anchor_liquid_phase(substance)
                    || component.total_concentration() <= TRACE_CONCENTRATION_MOL_PER_BUCKET
                {
                    return None;
                }
                let is_liquid_at_temperature = substance
                    .aggregate_state_at(self.temperature_kelvin)
                    .ok()
                    .is_some_and(|state| state == SubstanceAggregateState::Liquid);
                let has_matching_molten_phase =
                    match substance.phase_properties.preferred_liquid_phase {
                        LiquidPhasePreference::MoltenMetal => {
                            component.amount_in_phase(MixturePhase::MoltenMetal)
                                > TRACE_CONCENTRATION_MOL_PER_BUCKET
                        }
                        LiquidPhasePreference::MoltenSlag => {
                            component.amount_in_phase(MixturePhase::MoltenSlag)
                                > TRACE_CONCENTRATION_MOL_PER_BUCKET
                        }
                        LiquidPhasePreference::Aqueous | LiquidPhasePreference::Organic => {
                            component.condensed_concentration() > TRACE_CONCENTRATION_MOL_PER_BUCKET
                        }
                    };
                (is_liquid_at_temperature || has_matching_molten_phase)
                    .then_some(component.substance)
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

    fn solid_phase_amounts_by_anchor(&self) -> BTreeMap<SolidPhaseKey, f64> {
        let mut phases = BTreeMap::new();
        for component in &self.components {
            for (phase, amount) in &component.solid_mol_per_bucket_by_phase {
                *phases.entry(*phase).or_insert(0.0) += amount;
            }
        }
        phases.retain(|_, amount| *amount > TRACE_CONCENTRATION_MOL_PER_BUCKET);
        phases
    }

    pub(crate) fn total_in_phase(&self, phase: MixturePhase) -> f64 {
        self.components
            .iter()
            .map(|component| component.amount_in_phase(phase))
            .sum()
    }

    fn pressure_mol_per_bucket(&self) -> f64 {
        self.components
            .iter()
            .map(|component| {
                component.amount_in_phases(&[MixturePhase::Gas, MixturePhase::SupercriticalFluid])
            })
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
            + self.molten_mol_per_bucket_by_phase.values().sum::<f64>()
            + self.gas_mol_per_bucket
            + self.supercritical_mol_per_bucket
            + self.solid_mol_per_bucket_by_phase.values().sum::<f64>()
    }

    fn condensed_concentration(&self) -> f64 {
        self.amount_in_phases(&liquid_phases()) + self.amount_in_phase(MixturePhase::Solid)
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
            MixturePhase::MoltenMetal | MixturePhase::MoltenSlag => self
                .molten_mol_per_bucket_by_phase
                .iter()
                .filter(|(key, _)| key.coarse_phase == phase)
                .map(|(_, amount)| amount)
                .sum(),
            MixturePhase::Gas => self.gas_mol_per_bucket,
            MixturePhase::SupercriticalFluid => self.supercritical_mol_per_bucket,
            MixturePhase::Solid => self.solid_mol_per_bucket_by_phase.values().sum(),
        }
    }

    fn amount_in_organic_solvent(&self, solvent: SubstanceIndex) -> f64 {
        self.organic_mol_per_bucket_by_solvent
            .get(&solvent)
            .copied()
            .unwrap_or(0.0)
    }

    fn amount_in_liquid_phase(&self, phase: &LiquidPhaseState) -> f64 {
        match phase.coarse_phase {
            MixturePhase::Aqueous => self.aqueous_mol_per_bucket,
            MixturePhase::Organic => self.amount_in_organic_solvent(phase.representative_solvent),
            MixturePhase::MoltenMetal | MixturePhase::MoltenSlag => self
                .molten_mol_per_bucket_by_phase
                .get(&CondensedPhaseKey::new(
                    phase.coarse_phase,
                    phase.representative_solvent,
                ))
                .copied()
                .unwrap_or(0.0),
            MixturePhase::Gas | MixturePhase::SupercriticalFluid | MixturePhase::Solid => 0.0,
        }
    }

    fn phase_amounts(&self) -> Vec<f64> {
        let mut amounts = vec![
            self.aqueous_mol_per_bucket,
            self.gas_mol_per_bucket,
            self.supercritical_mol_per_bucket,
        ];
        amounts.extend(self.organic_mol_per_bucket_by_solvent.values().copied());
        amounts.extend(self.molten_mol_per_bucket_by_phase.values().copied());
        amounts.extend(self.solid_mol_per_bucket_by_phase.values().copied());
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
            MixturePhase::MoltenMetal | MixturePhase::MoltenSlag => {
                self.add_to_phase_for_solvent(phase, self.substance, amount)
            }
            MixturePhase::Gas => self.gas_mol_per_bucket += amount,
            MixturePhase::SupercriticalFluid => self.supercritical_mol_per_bucket += amount,
            MixturePhase::Solid => {
                *self
                    .solid_mol_per_bucket_by_phase
                    .entry(SolidPhaseKey::new(self.substance))
                    .or_insert(0.0) += amount;
            }
        }
    }

    fn add_to_phase_for_solvent(
        &mut self,
        phase: MixturePhase,
        solvent: SubstanceIndex,
        amount: f64,
    ) {
        match phase {
            MixturePhase::Organic => {
                *self
                    .organic_mol_per_bucket_by_solvent
                    .entry(solvent)
                    .or_insert(0.0) += amount;
            }
            MixturePhase::MoltenMetal | MixturePhase::MoltenSlag => {
                *self
                    .molten_mol_per_bucket_by_phase
                    .entry(CondensedPhaseKey::new(phase, solvent))
                    .or_insert(0.0) += amount;
            }
            _ => self.add_to_phase(phase, amount),
        }
    }

    fn clear_phase_amounts(&mut self) {
        self.aqueous_mol_per_bucket = 0.0;
        self.organic_mol_per_bucket_by_solvent.clear();
        self.molten_mol_per_bucket_by_phase.clear();
        self.gas_mol_per_bucket = 0.0;
        self.supercritical_mol_per_bucket = 0.0;
        self.solid_mol_per_bucket_by_phase.clear();
    }

    fn set_phase_amounts(&mut self, phase_amounts: ComponentPhaseAmounts) {
        self.aqueous_mol_per_bucket = phase_amounts.aqueous_mol_per_bucket;
        self.organic_mol_per_bucket_by_solvent = phase_amounts.organic_mol_per_bucket_by_solvent;
        self.molten_mol_per_bucket_by_phase = phase_amounts.molten_mol_per_bucket_by_phase;
        self.gas_mol_per_bucket = phase_amounts.gas_mol_per_bucket;
        self.supercritical_mol_per_bucket = phase_amounts.supercritical_mol_per_bucket;
        self.solid_mol_per_bucket_by_phase = phase_amounts.solid_mol_per_bucket_by_phase;
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
            MixturePhase::MoltenMetal | MixturePhase::MoltenSlag => {
                remove_from_condensed_phase_map(
                    &mut self.molten_mol_per_bucket_by_phase,
                    phase,
                    amount,
                )
            }
            MixturePhase::Gas => self.gas_mol_per_bucket -= amount,
            MixturePhase::SupercriticalFluid => self.supercritical_mol_per_bucket -= amount,
            MixturePhase::Solid => {
                remove_from_solid_phase_map(&mut self.solid_mol_per_bucket_by_phase, amount)
            }
        }
    }
}

fn preferred_phase(preference: LiquidPhasePreference) -> MixturePhase {
    match preference {
        LiquidPhasePreference::Aqueous => MixturePhase::Aqueous,
        LiquidPhasePreference::Organic => MixturePhase::Organic,
        LiquidPhasePreference::MoltenMetal => MixturePhase::MoltenMetal,
        LiquidPhasePreference::MoltenSlag => MixturePhase::MoltenSlag,
    }
}

fn liquid_phases() -> [MixturePhase; 4] {
    [
        MixturePhase::Aqueous,
        MixturePhase::Organic,
        MixturePhase::MoltenMetal,
        MixturePhase::MoltenSlag,
    ]
}

fn phase_uses_concrete_liquid_anchor(phase: MixturePhase) -> bool {
    matches!(
        phase,
        MixturePhase::Organic | MixturePhase::MoltenMetal | MixturePhase::MoltenSlag
    )
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

fn pressure_phase_for_substance(
    substance: &Substance,
    temperature_kelvin: f64,
    pressure_pascal: f64,
) -> MixturePhase {
    if substance_is_supercritical(substance, temperature_kelvin, pressure_pascal) {
        MixturePhase::SupercriticalFluid
    } else {
        MixturePhase::Gas
    }
}

fn substance_is_supercritical(
    substance: &Substance,
    temperature_kelvin: f64,
    pressure_pascal: f64,
) -> bool {
    let (Some(critical_temperature), Some(critical_pressure)) = (
        substance.critical_temperature_kelvin,
        substance.critical_pressure_pascal,
    ) else {
        return false;
    };
    temperature_kelvin >= critical_temperature && pressure_pascal >= critical_pressure
}

fn substance_is_solvent(substance: &Substance) -> bool {
    substance.phase_properties.solvent_role != SolventRole::NotSolvent
        && substance.phase_properties.can_form_liquid_phase
}

fn substance_can_anchor_liquid_phase(substance: &Substance) -> bool {
    match substance.phase_properties.preferred_liquid_phase {
        LiquidPhasePreference::Aqueous | LiquidPhasePreference::Organic => {
            substance_is_solvent(substance)
        }
        LiquidPhasePreference::MoltenMetal | LiquidPhasePreference::MoltenSlag => {
            substance.phase_properties.can_form_liquid_phase
                && substance.phase_properties.solvent_role == SolventRole::NotSolvent
        }
    }
}

fn classify_liquid_phase(
    registry: &ChemistryRegistry,
    solvent_amounts: &BTreeMap<SubstanceIndex, f64>,
    representative_solvent: SubstanceIndex,
) -> ChemistryResult<MixturePhase> {
    let total = solvent_amounts.values().sum::<f64>();
    if total <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
        let representative = registry.substance_by_index(representative_solvent)?;
        return Ok(preferred_phase(
            representative.phase_properties.preferred_liquid_phase,
        ));
    }
    let mut has_molten_metal = false;
    let mut has_molten_slag = false;
    for solvent in solvent_amounts.keys() {
        match registry
            .substance_by_index(*solvent)?
            .phase_properties
            .preferred_liquid_phase
        {
            LiquidPhasePreference::MoltenMetal => has_molten_metal = true,
            LiquidPhasePreference::MoltenSlag => has_molten_slag = true,
            LiquidPhasePreference::Aqueous | LiquidPhasePreference::Organic => {}
        }
    }
    if has_molten_metal && has_molten_slag {
        return Err(ChemistryError::InvalidMixtureState(
            "one liquid phase cannot be both molten metal and molten slag".to_string(),
        ));
    }
    if has_molten_metal {
        return Ok(MixturePhase::MoltenMetal);
    }
    if has_molten_slag {
        return Ok(MixturePhase::MoltenSlag);
    }
    let water: SubstanceId = "destroy:water".into();
    let water_fraction = registry
        .substance_index(&water)
        .and_then(|water| solvent_amounts.get(&water).copied())
        .unwrap_or(0.0)
        / total;
    if water_fraction >= 0.25 {
        Ok(MixturePhase::Aqueous)
    } else {
        Ok(MixturePhase::Organic)
    }
}

fn liquid_phase_snapshot(
    registry: &ChemistryRegistry,
    phase: LiquidPhaseState,
) -> ChemistryResult<LiquidPhaseSnapshot> {
    let total_solvent_mol_per_bucket = phase.solvent_amounts.values().sum::<f64>();
    if !total_solvent_mol_per_bucket.is_finite()
        || total_solvent_mol_per_bucket < TRACE_CONCENTRATION_MOL_PER_BUCKET
    {
        return Err(ChemistryError::InvalidMixtureState(
            "liquid phase solvent amount must be positive and finite".to_string(),
        ));
    }
    let mut solvents = Vec::new();
    for (solvent, concentration_mol_per_bucket) in phase.solvent_amounts {
        if !concentration_mol_per_bucket.is_finite() || concentration_mol_per_bucket < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "liquid phase solvent concentration must be non-negative and finite".to_string(),
            ));
        }
        solvents.push(LiquidPhaseSolventAmount {
            substance_id: registry.substance_by_index(solvent)?.id.clone(),
            concentration_mol_per_bucket,
            mole_fraction: concentration_mol_per_bucket / total_solvent_mol_per_bucket,
        });
    }
    Ok(LiquidPhaseSnapshot {
        id: phase.id,
        coarse_phase: phase.coarse_phase,
        representative_solvent_id: registry
            .substance_by_index(phase.representative_solvent)?
            .id
            .clone(),
        total_solvent_mol_per_bucket,
        solvents,
    })
}

fn mixed_phase_solubility_limit(
    registry: &ChemistryRegistry,
    substance: &Substance,
    phase: &LiquidPhaseState,
) -> Option<f64> {
    if phase
        .solvents
        .contains(&registry.substance_index(&substance.id)?)
        && substance_can_anchor_liquid_phase(substance)
    {
        return None;
    }
    let total_solvent = phase.solvent_amounts.values().sum::<f64>();
    if total_solvent <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
        return Some(0.0);
    }
    let mut weighted_log = 0.0;
    let mut finite_weight = 0.0;
    let mut has_unlimited = false;
    for (solvent, amount) in &phase.solvent_amounts {
        let fraction = amount / total_solvent;
        let limit = solubility_limit_in_solvent(registry, substance, *solvent);
        match limit {
            None => has_unlimited = true,
            Some(limit) if limit > TRACE_CONCENTRATION_MOL_PER_BUCKET => {
                weighted_log += fraction * limit.ln();
                finite_weight += fraction;
            }
            Some(_) => {}
        }
    }
    if has_unlimited && finite_weight <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
        return None;
    }
    if finite_weight <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
        return Some(0.0);
    }
    let finite_limit = (weighted_log / finite_weight).exp();
    if has_unlimited {
        Some(finite_limit / finite_weight.max(0.1))
    } else {
        Some(finite_limit)
    }
}

fn solubility_limit_in_solvent(
    registry: &ChemistryRegistry,
    substance: &Substance,
    solvent: SubstanceIndex,
) -> Option<f64> {
    let substance_index = registry.substance_index(&substance.id)?;
    let solvent_substance = registry.substance_by_index(solvent).ok()?;
    if matches!(
        (
            substance.phase_properties.preferred_liquid_phase,
            solvent_substance.phase_properties.preferred_liquid_phase,
        ),
        (
            LiquidPhasePreference::MoltenMetal,
            LiquidPhasePreference::MoltenMetal
        ) | (
            LiquidPhasePreference::MoltenSlag,
            LiquidPhasePreference::MoltenSlag
        )
    ) && substance.phase_properties.can_form_liquid_phase
    {
        return None;
    }
    if matches!(
        (
            &substance.representation,
            solvent_substance.phase_properties.preferred_liquid_phase
        ),
        (
            SubstanceRepresentation::MetallurgicalSolute { .. },
            LiquidPhasePreference::MoltenMetal
        )
    ) {
        return None;
    }
    if substance_is_solvent(substance) && substance_index != solvent {
        match registry.solvent_miscibility(substance_index, solvent) {
            SolventMiscibility::FullyMiscible => return None,
            SolventMiscibility::PartiallyMiscible {
                limit_mol_per_bucket,
            } => return Some(limit_mol_per_bucket),
            SolventMiscibility::Immiscible => {}
        }
        if solvent_substance.id.as_str() == "destroy:water" {
            return substance.phase_properties.aqueous_solubility_mol_per_bucket;
        }
        return Some(0.0);
    }
    let solvent_id = &solvent_substance.id;
    if solvent_id == &substance.id && substance_can_anchor_liquid_phase(substance) {
        return None;
    }
    if solvent_id.as_str() == "destroy:water" {
        substance.phase_properties.aqueous_solubility_mol_per_bucket
    } else {
        substance.phase_properties.organic_solubility_mol_per_bucket
    }
}

fn condensed_solubility_capacity(
    registry: &ChemistryRegistry,
    substance: &Substance,
    liquid_phases: &[LiquidPhaseState],
) -> Option<f64> {
    let mut total = 0.0;
    for phase in liquid_phases {
        match mixed_phase_solubility_limit(registry, substance, phase) {
            None => return None,
            Some(limit) => total += limit,
        }
    }
    Some(total)
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

fn gas_pressure_for_moles(
    moles_per_bucket: f64,
    temperature_kelvin: f64,
    gas_volume_cubic_meters: f64,
) -> f64 {
    moles_per_bucket * GAS_CONSTANT_J_PER_MOL_KELVIN * temperature_kelvin / gas_volume_cubic_meters
}

fn gas_moles_for_pressure(
    pressure_pascal: f64,
    temperature_kelvin: f64,
    gas_volume_cubic_meters: f64,
) -> ChemistryResult<f64> {
    if !pressure_pascal.is_finite() || pressure_pascal < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "gas pressure must be non-negative and finite".to_string(),
        ));
    }
    if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "temperature must be positive and finite for gas pressure conversion".to_string(),
        ));
    }
    if !gas_volume_cubic_meters.is_finite() || gas_volume_cubic_meters <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "gas volume must be positive and finite".to_string(),
        ));
    }
    Ok(pressure_pascal * gas_volume_cubic_meters
        / (GAS_CONSTANT_J_PER_MOL_KELVIN * temperature_kelvin))
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
    registry: &ChemistryRegistry,
    substance_index: SubstanceIndex,
    substance: &Substance,
    liquid_phases: &[LiquidPhaseState],
    can_precipitate: bool,
    amount: f64,
    target: &mut ComponentPhaseAmounts,
) -> ChemistryResult<()> {
    if amount <= 0.0 {
        return Ok(());
    }
    let mut remaining = amount;
    let mut candidates = liquid_phases
        .iter()
        .map(|phase| {
            (
                phase,
                mixed_phase_solubility_limit(registry, substance, phase),
            )
        })
        .collect::<Vec<_>>();
    candidates.sort_by(|(left_phase, left_limit), (right_phase, right_limit)| {
        match (left_limit, right_limit) {
            (None, None) => left_phase.id.cmp(&right_phase.id),
            (None, Some(_)) => std::cmp::Ordering::Less,
            (Some(_), None) => std::cmp::Ordering::Greater,
            (Some(left_limit), Some(right_limit)) => right_limit
                .partial_cmp(left_limit)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| left_phase.id.cmp(&right_phase.id)),
        }
    });
    for (phase, limit) in candidates {
        fill_liquid_phase(phase, limit, &mut remaining, target);
        if remaining <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Ok(());
        }
    }
    if remaining <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
        return Ok(());
    }
    if can_precipitate {
        *target
            .solid_mol_per_bucket_by_phase
            .entry(SolidPhaseKey::new(substance_index))
            .or_insert(0.0) += remaining;
        return Ok(());
    }
    Err(ChemistryError::InvalidMixtureState(format!(
        "substance '{}' has {remaining} mol/bucket that cannot fit any liquid phase and cannot precipitate",
        substance.id
    )))
}

fn fill_liquid_phase(
    phase: &LiquidPhaseState,
    limit: Option<f64>,
    remaining: &mut f64,
    target: &mut ComponentPhaseAmounts,
) {
    if *remaining <= 0.0 {
        return;
    }
    let existing = match phase.coarse_phase {
        MixturePhase::Aqueous => target.aqueous_mol_per_bucket,
        MixturePhase::Organic => Some(phase.representative_solvent)
            .and_then(|solvent| {
                target
                    .organic_mol_per_bucket_by_solvent
                    .get(&solvent)
                    .copied()
            })
            .unwrap_or(0.0),
        MixturePhase::MoltenMetal | MixturePhase::MoltenSlag => target
            .molten_mol_per_bucket_by_phase
            .get(&CondensedPhaseKey::new(
                phase.coarse_phase,
                phase.representative_solvent,
            ))
            .copied()
            .unwrap_or(0.0),
        MixturePhase::Gas | MixturePhase::SupercriticalFluid | MixturePhase::Solid => 0.0,
    };
    let capacity = limit
        .map(|limit| (limit - existing).max(0.0))
        .unwrap_or(*remaining);
    let moved = capacity.min(*remaining);
    match phase.coarse_phase {
        MixturePhase::Aqueous => target.aqueous_mol_per_bucket += moved,
        MixturePhase::Organic => {
            *target
                .organic_mol_per_bucket_by_solvent
                .entry(phase.representative_solvent)
                .or_insert(0.0) += moved;
        }
        MixturePhase::MoltenMetal | MixturePhase::MoltenSlag => {
            *target
                .molten_mol_per_bucket_by_phase
                .entry(CondensedPhaseKey::new(
                    phase.coarse_phase,
                    phase.representative_solvent,
                ))
                .or_insert(0.0) += moved;
        }
        MixturePhase::Gas | MixturePhase::SupercriticalFluid | MixturePhase::Solid => {}
    }
    *remaining -= moved;
}

fn normalize_phase_amounts(
    mut phase_amounts: ComponentPhaseAmounts,
    total_amount: f64,
) -> ComponentPhaseAmounts {
    phase_amounts.aqueous_mol_per_bucket /= total_amount;
    phase_amounts.gas_mol_per_bucket /= total_amount;
    phase_amounts.supercritical_mol_per_bucket /= total_amount;
    for concentration in phase_amounts.solid_mol_per_bucket_by_phase.values_mut() {
        *concentration /= total_amount;
    }
    for concentration in phase_amounts.organic_mol_per_bucket_by_solvent.values_mut() {
        *concentration /= total_amount;
    }
    for concentration in phase_amounts.molten_mol_per_bucket_by_phase.values_mut() {
        *concentration /= total_amount;
    }
    phase_amounts
}

fn remove_from_condensed_phase_map(
    amounts_by_anchor: &mut BTreeMap<CondensedPhaseKey, f64>,
    phase: MixturePhase,
    mut amount: f64,
) {
    let anchors = amounts_by_anchor
        .keys()
        .copied()
        .filter(|key| key.coarse_phase == phase)
        .collect::<Vec<_>>();
    for anchor in anchors {
        if amount <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            break;
        }
        let Some(current) = amounts_by_anchor.get_mut(&anchor) else {
            continue;
        };
        let removed = (*current).min(amount);
        *current -= removed;
        amount -= removed;
    }
    amounts_by_anchor.retain(|_, value| *value > TRACE_CONCENTRATION_MOL_PER_BUCKET);
}

fn remove_from_solid_phase_map(
    amounts_by_anchor: &mut BTreeMap<SolidPhaseKey, f64>,
    mut amount: f64,
) {
    let anchors = amounts_by_anchor.keys().copied().collect::<Vec<_>>();
    for anchor in anchors {
        if amount <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            break;
        }
        let Some(current) = amounts_by_anchor.get_mut(&anchor) else {
            continue;
        };
        let removed = (*current).min(amount);
        *current -= removed;
        amount -= removed;
    }
    amounts_by_anchor.retain(|_, value| *value > TRACE_CONCENTRATION_MOL_PER_BUCKET);
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
    use crate::chemistry::substance::{
        LiquidPhasePreference, SolventRole, Substance, SubstancePhaseProperties,
    };

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
                    .with_critical_point(647.0, 22_064_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new("destroy:oxygen", 0, 32.0, 1_140.0, 90.0, 29.4, 6_820.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: None,
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: false,
                        can_form_liquid_phase: false,
                        solvent_role: SolventRole::NotSolvent,
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
                        solvent_role: SolventRole::NotSolvent,
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

    fn vapor_liquid_registry() -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_critical_point(647.0, 22_064_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new(
                    "destroy:liquefiable_gas",
                    0,
                    58.0,
                    580.0,
                    272.0,
                    110.0,
                    21_000.0,
                )
                .with_critical_point(425.0, 3_800_000.0)
                .with_vapor_pressure_model(
                    crate::chemistry::substance::VaporPressureModel::ClausiusClapeyron {
                        reference_temperature_kelvin: 298.0,
                        reference_pressure_pascal: 200_000.0,
                        enthalpy_j_per_mol: 21_000.0,
                    },
                )
                .with_phase_properties(SubstancePhaseProperties::organic_unlimited(0.0)),
            )
            .substance(
                Substance::new(
                    "destroy:supercritical_gas",
                    0,
                    44.0,
                    440.0,
                    272.0,
                    90.0,
                    16_000.0,
                )
                .with_critical_point(290.0, 7_000_000.0)
                .with_phase_properties(SubstancePhaseProperties::organic_unlimited(0.0)),
            )
            .substance(
                Substance::new(
                    "destroy:volatile_solute",
                    0,
                    120.0,
                    120_000.0,
                    520.0,
                    100.0,
                    20_000.0,
                )
                .with_phase_properties(SubstancePhaseProperties {
                    preferred_liquid_phase: LiquidPhasePreference::Organic,
                    aqueous_solubility_mol_per_bucket: Some(0.0),
                    organic_solubility_mol_per_bucket: Some(0.2),
                    can_precipitate: true,
                    can_form_liquid_phase: false,
                    solvent_role: SolventRole::NotSolvent,
                }),
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
                        solvent_role: SolventRole::KnownSolvent,
                    }),
            )
            .build()
            .unwrap()
    }

    fn phase_registry() -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new("destroy:salt", 0, 58.0, 58_000.0, 1_000.0, 80.0, 20_000.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: Some(0.5),
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: true,
                        can_form_liquid_phase: false,
                        solvent_role: SolventRole::NotSolvent,
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
            molten_mol_per_bucket_by_phase: BTreeMap::new(),
            gas_mol_per_bucket: f64::INFINITY,
            supercritical_mol_per_bucket: 0.0,
            solid_mol_per_bucket_by_phase: BTreeMap::new(),
        });
        mixture.positions_by_substance = vec![Some(0)];

        let error = mixture
            .recalculate_volume_millibuckets(&registry, 1000)
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidMixtureState(_)));
    }

    #[test]
    fn distinct_solids_are_tracked_as_distinct_solid_phases() {
        let solid_phase_properties = SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::MoltenSlag,
            aqueous_solubility_mol_per_bucket: Some(0.0),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: true,
            solvent_role: SolventRole::NotSolvent,
        };
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new(
                    "destroy:solid_a",
                    0,
                    10.0,
                    10_000.0,
                    2_000.0,
                    20.0,
                    10_000.0,
                )
                .with_melting_point_kelvin(1_000.0)
                .with_phase_properties(solid_phase_properties.clone()),
            )
            .substance(
                Substance::new(
                    "destroy:solid_b",
                    0,
                    20.0,
                    20_000.0,
                    2_000.0,
                    25.0,
                    10_000.0,
                )
                .with_melting_point_kelvin(1_100.0)
                .with_phase_properties(solid_phase_properties),
            )
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:solid_a", 0.2)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:solid_b", 0.3)
            .unwrap();

        let snapshots = mixture.solid_phase_snapshots(&registry).unwrap();

        assert_eq!(mixture.solid_phase_count(), 2);
        assert_eq!(snapshots.len(), 2);
        assert!(snapshots.iter().any(|phase| {
            phase.representative_substance_id == SubstanceId::from("destroy:solid_a")
                && (phase.concentration_mol_per_bucket - 0.2).abs() < 1.0e-9
        }));
        assert!(snapshots.iter().any(|phase| {
            phase.representative_substance_id == SubstanceId::from("destroy:solid_b")
                && (phase.concentration_mol_per_bucket - 0.3).abs() < 1.0e-9
        }));
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

        assert!((mixture.concentration_of(&water) - 1.0).abs() < 0.01);
        assert!(
            (mixture.concentration_in_phase(&water, MixturePhase::Aqueous) - 0.75).abs() < 0.01
        );
        assert!((mixture.concentration_in_phase(&water, MixturePhase::Gas) - 0.25).abs() < 0.01);
        assert!((mixture.gaseous_fraction_of(&water) - 0.25).abs() < 0.01);
    }

    #[test]
    fn solubility_excess_precipitates_and_can_redissolve() {
        let registry = phase_registry();
        let salt: SubstanceId = "destroy:salt".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixture.add_substance(&registry, salt.clone(), 1.0).unwrap();

        assert!((mixture.concentration_in_phase(&salt, MixturePhase::Aqueous) - 0.5).abs() < 0.01);
        assert!((mixture.concentration_in_phase(&salt, MixturePhase::Solid) - 0.5).abs() < 0.01);

        mixture
            .change_concentration(&registry, &salt, -0.5)
            .unwrap();

        assert!((mixture.concentration_of(&salt) - 0.5).abs() < 0.01);
        assert!((mixture.concentration_in_phase(&salt, MixturePhase::Aqueous) - 0.5).abs() < 0.01);
        assert!((mixture.concentration_in_phase(&salt, MixturePhase::Solid)).abs() < 0.01);
    }

    #[test]
    fn neutral_organic_substance_prefers_organic_phase() {
        let registry = phase_registry();
        let oil: SubstanceId = "destroy:oil".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture.add_substance(&registry, oil.clone(), 1.0).unwrap();

        assert!((mixture.concentration_in_phase(&oil, MixturePhase::Organic) - 1.0).abs() < 0.01);
        assert!((mixture.concentration_in_phase(&oil, MixturePhase::Aqueous)).abs() < 0.01);
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

        assert!(
            (mixture
                .concentration_in_organic_solvent(&registry, &oil, &oil)
                .unwrap()
                - 1.0)
                .abs()
                < 0.01
        );
        assert!(
            (mixture
                .concentration_in_organic_solvent(&registry, &oil, &chloroform)
                .unwrap())
            .abs()
                < 0.01
        );
        assert!(
            (mixture
                .concentration_in_organic_solvent(&registry, &chloroform, &chloroform)
                .unwrap()
                - 2.0)
                .abs()
                < 0.1
        );
        assert!(
            mixture
                .organic_phase_amounts_of(&registry, &oil)
                .unwrap()
                .len()
                >= 1
        );
        assert!(
            mixture
                .organic_phase_amounts_of(&registry, &chloroform)
                .unwrap()
                .len()
                >= 1
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
    fn melted_materials_enter_explicit_metallurgical_phases() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new(
                    "destroy:test_iron",
                    0,
                    55.845,
                    787_400.0,
                    3200.0,
                    25.0,
                    350_000.0,
                )
                .with_solid_density_grams_per_bucket(787_400.0)
                .with_melting_point_kelvin(1000.0)
                .with_fusion_heat_j_per_mol(10.0)
                .with_phase_properties(SubstancePhaseProperties {
                    preferred_liquid_phase: LiquidPhasePreference::MoltenMetal,
                    aqueous_solubility_mol_per_bucket: Some(0.0),
                    organic_solubility_mol_per_bucket: Some(0.0),
                    can_precipitate: true,
                    can_form_liquid_phase: true,
                    solvent_role: SolventRole::NotSolvent,
                }),
            )
            .substance(
                Substance::new(
                    "destroy:test_silica",
                    0,
                    60.084,
                    220_000.0,
                    3200.0,
                    45.0,
                    400_000.0,
                )
                .with_solid_density_grams_per_bucket(265_000.0)
                .with_melting_point_kelvin(1200.0)
                .with_fusion_heat_j_per_mol(10.0)
                .with_phase_properties(SubstancePhaseProperties {
                    preferred_liquid_phase: LiquidPhasePreference::MoltenSlag,
                    aqueous_solubility_mol_per_bucket: Some(0.0),
                    organic_solubility_mol_per_bucket: Some(0.0),
                    can_precipitate: true,
                    can_form_liquid_phase: true,
                    solvent_role: SolventRole::NotSolvent,
                }),
            )
            .build()
            .unwrap();

        let iron: SubstanceId = "destroy:test_iron".into();
        let silica: SubstanceId = "destroy:test_silica".into();
        let mut mixture = Mixture::new(1300.0).unwrap();
        mixture.add_substance(&registry, iron.clone(), 1.0).unwrap();
        mixture
            .add_substance(&registry, silica.clone(), 1.0)
            .unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&iron, MixturePhase::MoltenMetal),
            1.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&iron, MixturePhase::Organic),
            0.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&silica, MixturePhase::MoltenSlag),
            1.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&silica, MixturePhase::Organic),
            0.0
        );
        assert_eq!(mixture.liquid_phase_count(&registry).unwrap(), 2);
    }

    #[test]
    fn gas_solubility_uses_pressure_and_activity() {
        let registry = gas_registry();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let sodium: SubstanceId = "destroy:sodium".into();
        let mut pure = Mixture::new(298.0).unwrap();
        pure.add_substance(&registry, "destroy:water", 1.0).unwrap();
        pure.add_substance(&registry, oxygen.clone(), 1.0).unwrap();
        pure.transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();
        let pure_dissolved = pure.concentration_in_phase(&oxygen, MixturePhase::Aqueous);

        let mut salty = Mixture::new(298.0).unwrap();
        salty
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        salty.add_substance(&registry, sodium, 1.0).unwrap();
        salty.add_substance(&registry, oxygen.clone(), 1.0).unwrap();
        salty
            .transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();
        let salty_dissolved = salty.concentration_in_phase(&oxygen, MixturePhase::Aqueous);

        assert!(pure.gas_pressure_pascal() > 0.0);
        assert!(pure_dissolved >= 0.0);
        assert!(salty_dissolved <= pure_dissolved + 1.0e-9);
    }

    #[test]
    fn gas_solubility_uses_partial_pressure_of_that_gas() {
        let registry = gas_registry();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let ballast: SubstanceId = "destroy:unknown_gas".into();
        let mut oxygen_only = Mixture::new(298.0).unwrap();
        oxygen_only
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        oxygen_only
            .add_substance(&registry, oxygen.clone(), 1.0)
            .unwrap();
        oxygen_only
            .transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();

        let mut with_ballast = Mixture::new(298.0).unwrap();
        with_ballast
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        with_ballast
            .add_substance(&registry, oxygen.clone(), 1.0)
            .unwrap();
        with_ballast.add_substance(&registry, ballast, 1.0).unwrap();
        with_ballast
            .transfer_gases_toward_solubility_equilibrium(&registry, 1.0)
            .unwrap();

        let oxygen_only_dissolved =
            oxygen_only.concentration_in_phase(&oxygen, MixturePhase::Aqueous);
        let with_ballast_dissolved =
            with_ballast.concentration_in_phase(&oxygen, MixturePhase::Aqueous);

        assert!(
            (oxygen_only_dissolved - with_ballast_dissolved).abs() < 1.0e-12,
            "oxygen-only dissolved {oxygen_only_dissolved}, with ballast {with_ballast_dissolved}"
        );
        assert!(
            with_ballast.gas_pressure_pascal() > oxygen_only.gas_pressure_pascal(),
            "ballast gas should still raise total pressure"
        );
    }

    #[test]
    fn open_atmosphere_exchange_sets_gas_phase_by_partial_pressures() {
        let registry = gas_registry();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let ballast: SubstanceId = "destroy:unknown_gas".into();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, oxygen.clone(), 1.0)
            .unwrap();

        let delta = mixture
            .exchange_gases_with_atmosphere(
                &registry,
                &[(oxygen.clone(), 0.21), (ballast.clone(), 0.79)],
                STANDARD_PRESSURE_PASCAL,
                10.0,
                10.0,
            )
            .unwrap();

        assert!(delta > 0.0);
        assert!(
            (mixture.gas_pressure_pascal() / STANDARD_PRESSURE_PASCAL - 1.0).abs() < 1.0e-6,
            "open gas pressure was {} Pa",
            mixture.gas_pressure_pascal()
        );
        assert!(
            (mixture.gas_partial_pressure_pascal(&oxygen) / STANDARD_PRESSURE_PASCAL - 0.21).abs()
                < 1.0e-6
        );
        assert!(
            (mixture.gas_partial_pressure_pascal(&ballast) / STANDARD_PRESSURE_PASCAL - 0.79).abs()
                < 1.0e-6
        );
    }

    #[test]
    fn gas_without_solubility_data_does_not_enter_solution() {
        let registry = gas_registry();
        let gas: SubstanceId = "destroy:unknown_gas".into();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.set_gas_volume_cubic_meters(0.1).unwrap();

        mixture.add_substance(&registry, gas.clone(), 1.0).unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&gas, MixturePhase::Aqueous),
            0.0
        );
        assert_eq!(mixture.concentration_in_phase(&gas, MixturePhase::Gas), 1.0);
    }

    #[test]
    fn supercritical_substance_uses_supercritical_phase_above_critical_point_and_pressure() {
        let registry = vapor_liquid_registry();
        let fluid: SubstanceId = "destroy:supercritical_gas".into();
        let mut mixture = Mixture::new(320.0).unwrap();

        mixture
            .add_substance(&registry, fluid.clone(), 3.0)
            .unwrap();

        assert!(
            mixture.gas_pressure_pascal()
                > registry
                    .substance(&fluid)
                    .unwrap()
                    .critical_pressure_pascal
                    .unwrap()
        );
        assert_eq!(
            mixture.concentration_in_phase(&fluid, MixturePhase::SupercriticalFluid),
            3.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&fluid, MixturePhase::Gas),
            0.0
        );
    }

    #[test]
    fn supercritical_substance_remains_gas_when_pressure_is_below_critical_point() {
        let registry = vapor_liquid_registry();
        let fluid: SubstanceId = "destroy:supercritical_gas".into();
        let mut mixture = Mixture::new(320.0).unwrap();

        mixture
            .add_substance(&registry, fluid.clone(), 0.1)
            .unwrap();

        assert!(
            mixture.gas_pressure_pascal()
                < registry
                    .substance(&fluid)
                    .unwrap()
                    .critical_pressure_pascal
                    .unwrap()
        );
        assert_eq!(
            mixture.concentration_in_phase(&fluid, MixturePhase::Gas),
            0.1
        );
        assert_eq!(
            mixture.concentration_in_phase(&fluid, MixturePhase::SupercriticalFluid),
            0.0
        );
    }

    #[test]
    fn cooling_below_critical_point_returns_supercritical_fluid_to_vapor_liquid_equilibrium() {
        let registry = vapor_liquid_registry();
        let fluid: SubstanceId = "destroy:supercritical_gas".into();
        let mut mixture = Mixture::new(320.0).unwrap();
        mixture
            .add_substance(&registry, fluid.clone(), 3.0)
            .unwrap();
        assert!(mixture.concentration_in_phase(&fluid, MixturePhase::SupercriticalFluid) > 0.0);

        mixture.temperature_kelvin = 280.0;
        mixture.equilibrate_phases(&registry).unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&fluid, MixturePhase::SupercriticalFluid),
            0.0
        );
        assert!(mixture.concentration_in_phase(&fluid, MixturePhase::Gas) > 0.0);
        assert!(mixture.concentration_in_phase(&fluid, MixturePhase::Organic) > 0.0);
        assert!((mixture.concentration_of(&fluid) - 3.0).abs() < 1.0e-9);
    }

    #[test]
    fn gas_transfer_moves_toward_henry_limit_gradually() {
        let registry = gas_registry();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
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
        assert!(after_first >= 0.0);
        assert!(after_second >= after_first);
        assert!(mixture.concentration_in_phase(&oxygen, MixturePhase::Gas) > 0.0);
    }

    #[test]
    fn volatile_liquid_evaporates_to_saturation_pressure() {
        let registry = vapor_liquid_registry();
        let water: SubstanceId = "destroy:water".into();
        let mut mixture = Mixture::new(373.0).unwrap();

        mixture
            .add_substance(&registry, water.clone(), 1.0)
            .unwrap();

        assert!(mixture.concentration_in_phase(&water, MixturePhase::Aqueous) > 0.9);
        assert!(mixture.concentration_in_phase(&water, MixturePhase::Gas) > 0.0);
        assert!(
            mixture.temperature_kelvin() < 373.0,
            "evaporative cooling should lower temperature from 373K to {}",
            mixture.temperature_kelvin()
        );
        let vp = mixture.gas_partial_pressure_pascal(&water);
        assert!(
            vp > 0.0 && vp < STANDARD_PRESSURE_PASCAL,
            "partial pressure {} should be between 0 and 1 atm at cooled T={}K",
            vp,
            mixture.temperature_kelvin()
        );
    }

    #[test]
    fn gas_pressure_depends_on_headspace_volume() {
        let registry = gas_registry();
        let oxygen: SubstanceId = "destroy:oxygen".into();
        let mut small_headspace = Mixture::new(298.0).unwrap();
        let mut large_headspace = Mixture::new(298.0).unwrap();
        large_headspace
            .set_gas_volume_cubic_meters(DEFAULT_GAS_VOLUME_CUBIC_METERS * 4.0)
            .unwrap();

        small_headspace
            .add_substance(&registry, oxygen.clone(), 1.0)
            .unwrap();
        large_headspace
            .add_substance(&registry, oxygen, 1.0)
            .unwrap();

        assert!(
            (small_headspace.gas_pressure_pascal() / large_headspace.gas_pressure_pascal() - 4.0)
                .abs()
                < 1.0e-9
        );
    }

    #[test]
    fn larger_headspace_evaporates_more_liquid_to_same_vapor_pressure() {
        let registry = vapor_liquid_registry();
        let water: SubstanceId = "destroy:water".into();
        let mut small_headspace = Mixture::new(373.0).unwrap();
        let mut large_headspace = Mixture::new(373.0).unwrap();
        large_headspace
            .set_gas_volume_cubic_meters(DEFAULT_GAS_VOLUME_CUBIC_METERS * 3.0)
            .unwrap();

        small_headspace
            .add_substance(&registry, water.clone(), 1.0)
            .unwrap();
        large_headspace
            .add_substance(&registry, water.clone(), 1.0)
            .unwrap();

        let small_gas = small_headspace.concentration_in_phase(&water, MixturePhase::Gas);
        let large_gas = large_headspace.concentration_in_phase(&water, MixturePhase::Gas);

        assert!(
            large_gas > small_gas,
            "larger headspace should evaporate more: small={small_gas}, large={large_gas}"
        );
        assert!(
            small_headspace.temperature_kelvin() > large_headspace.temperature_kelvin(),
            "small headspace should be warmer (less evaporative cooling): small={}, large={}",
            small_headspace.temperature_kelvin(),
            large_headspace.temperature_kelvin()
        );
    }

    #[test]
    fn gas_condenses_when_partial_pressure_exceeds_saturation_pressure() {
        let registry = vapor_liquid_registry();
        let gas: SubstanceId = "destroy:liquefiable_gas".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture.add_substance(&registry, gas.clone(), 1.0).unwrap();

        assert!(mixture.concentration_in_phase(&gas, MixturePhase::Organic) > 0.0);
        assert!(mixture.concentration_in_phase(&gas, MixturePhase::Gas) > 0.0);
        assert!(
            mixture.temperature_kelvin() > 298.0,
            "condensation should release latent heat, raising temperature to {}",
            mixture.temperature_kelvin()
        );
        let vp = mixture.gas_partial_pressure_pascal(&gas);
        assert!(vp > 0.0, "partial pressure should be positive, got {vp}");
    }

    #[test]
    fn condensation_creates_phase_then_solutes_redistribute_into_it() {
        let registry = vapor_liquid_registry();
        let gas: SubstanceId = "destroy:liquefiable_gas".into();
        let solute: SubstanceId = "destroy:volatile_solute".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, solute.clone(), 0.3)
            .unwrap();
        assert_eq!(
            mixture.concentration_in_phase(&solute, MixturePhase::Organic),
            0.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&solute, MixturePhase::Solid),
            0.3
        );

        mixture.add_substance(&registry, gas, 1.0).unwrap();

        assert!(
            mixture.concentration_in_phase(&solute, MixturePhase::Organic) > 0.0,
            "solute did not enter the liquid phase created by gas condensation"
        );
        assert!(
            mixture.concentration_in_phase(&solute, MixturePhase::Solid) < 0.3,
            "solute stayed fully precipitated after solvent condensation"
        );
    }

    #[test]
    fn supercritical_substance_cannot_keep_condensed_phase() {
        let registry = vapor_liquid_registry();
        let gas: SubstanceId = "destroy:supercritical_gas".into();
        let mut mixture = Mixture::new(320.0).unwrap();

        mixture.add_substance(&registry, gas.clone(), 3.0).unwrap();
        mixture
            .move_between_phases(
                &registry,
                gas.clone(),
                MixturePhase::SupercriticalFluid,
                MixturePhase::Organic,
                0.5,
            )
            .unwrap();
        let delta = mixture.equilibrate_vapor_liquid(&registry, 1.0).unwrap();

        assert!(delta > 0.0);
        assert_eq!(
            mixture.concentration_in_phase(&gas, MixturePhase::Organic),
            0.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&gas, MixturePhase::SupercriticalFluid),
            3.0
        );
    }

    #[test]
    fn fully_miscible_solvents_share_one_liquid_phase() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
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

        assert_eq!(mixture.liquid_phase_count(&registry).unwrap(), 1);
        assert_eq!(mixture.concentration_of(&water), 1.0);
        assert_eq!(mixture.concentration_of(&ethanol), 1.0);
    }

    #[test]
    fn liquid_phase_snapshot_exposes_concrete_mixed_phase_composition() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new("destroy:ethanol", 0, 46.0, 46_000.0, 351.0, 110.0, 20_000.0)
                    .with_phase_properties(SubstancePhaseProperties::organic_unlimited(1.0)),
            )
            .substance(
                Substance::new(
                    "destroy:solute",
                    0,
                    120.0,
                    120_000.0,
                    700.0,
                    100.0,
                    20_000.0,
                )
                .with_phase_properties(SubstancePhaseProperties::organic_solute(0.2)),
            )
            .solvent_miscibility(
                "destroy:water",
                "destroy:ethanol",
                SolventMiscibility::FullyMiscible,
            )
            .build()
            .unwrap();
        let solute: SubstanceId = "destroy:solute".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:ethanol", 3.0)
            .unwrap();
        mixture
            .add_substance(&registry, solute.clone(), 0.1)
            .unwrap();

        let phases = mixture.liquid_phase_snapshots(&registry).unwrap();
        assert_eq!(phases.len(), 1);
        assert_eq!(phases[0].coarse_phase, MixturePhase::Aqueous);
        assert!(
            (phases[0].total_solvent_mol_per_bucket - 4.0).abs() < 0.1,
            "total solvent should be ~4.0, got {}",
            phases[0].total_solvent_mol_per_bucket
        );
        assert_eq!(phases[0].solvents.len(), 2);
        let water_fraction = phases[0]
            .solvents
            .iter()
            .find(|s| s.substance_id.as_str() == "destroy:water")
            .map(|s| s.mole_fraction)
            .unwrap_or(0.0);
        let ethanol_fraction = phases[0]
            .solvents
            .iter()
            .find(|s| s.substance_id.as_str() == "destroy:ethanol")
            .map(|s| s.mole_fraction)
            .unwrap_or(0.0);
        assert!(
            (water_fraction - 0.25).abs() < 0.05,
            "water mole fraction should be ~0.25, got {water_fraction}"
        );
        assert!(
            (ethanol_fraction - 0.75).abs() < 0.05,
            "ethanol mole fraction should be ~0.75, got {ethanol_fraction}"
        );

        let solute_amounts = mixture.liquid_phase_amounts_of(&registry, &solute).unwrap();
        assert_eq!(solute_amounts.len(), 1);
        assert_eq!(solute_amounts[0].phase_id, phases[0].id);
        assert!(
            (solute_amounts[0].concentration_mol_per_bucket - 0.1).abs() < 0.05,
            "solute should be ~0.1, got {}",
            solute_amounts[0].concentration_mol_per_bucket
        );
    }

    #[test]
    fn mixed_solvent_composition_changes_solute_capacity() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new("destroy:ethanol", 0, 46.0, 46_000.0, 351.0, 110.0, 20_000.0)
                    .with_phase_properties(SubstancePhaseProperties::organic_unlimited(1.0)),
            )
            .substance(
                Substance::new(
                    "destroy:solute",
                    0,
                    120.0,
                    120_000.0,
                    700.0,
                    100.0,
                    20_000.0,
                )
                .with_melting_point_kelvin(350.0)
                .with_phase_properties(SubstancePhaseProperties {
                    preferred_liquid_phase: LiquidPhasePreference::Organic,
                    aqueous_solubility_mol_per_bucket: Some(0.1),
                    organic_solubility_mol_per_bucket: Some(1.0),
                    can_precipitate: true,
                    can_form_liquid_phase: false,
                    solvent_role: SolventRole::NotSolvent,
                }),
            )
            .solvent_miscibility(
                "destroy:water",
                "destroy:ethanol",
                SolventMiscibility::FullyMiscible,
            )
            .build()
            .unwrap();
        let solute: SubstanceId = "destroy:solute".into();

        let mut water_only = Mixture::new(298.0).unwrap();
        water_only
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        water_only
            .add_substance(&registry, solute.clone(), 0.2)
            .unwrap();

        let mut mixed = Mixture::new(298.0).unwrap();
        mixed
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixed
            .add_substance(&registry, "destroy:ethanol", 1.0)
            .unwrap();
        mixed.add_substance(&registry, solute.clone(), 0.2).unwrap();

        assert_eq!(water_only.liquid_phase_count(&registry).unwrap(), 1);
        assert_eq!(mixed.liquid_phase_count(&registry).unwrap(), 1);
        assert!(water_only.concentration_in_phase(&solute, MixturePhase::Solid) > 0.09);
        assert_eq!(
            mixed.concentration_in_phase(&solute, MixturePhase::Solid),
            0.0
        );
        assert_eq!(mixed.concentration_of(&solute), 0.2);
    }

    #[test]
    fn unknown_solubility_does_not_create_hidden_compatibility() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 18_000.0, 373.0, 75.0, 40_650.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new(
                    "destroy:solute",
                    0,
                    120.0,
                    120_000.0,
                    700.0,
                    100.0,
                    20_000.0,
                )
                .with_melting_point_kelvin(350.0)
                .with_phase_properties(SubstancePhaseProperties {
                    preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                    aqueous_solubility_mol_per_bucket: Some(0.0),
                    organic_solubility_mol_per_bucket: Some(0.0),
                    can_precipitate: true,
                    can_form_liquid_phase: false,
                    solvent_role: SolventRole::NotSolvent,
                }),
            )
            .build()
            .unwrap();
        let solute: SubstanceId = "destroy:solute".into();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, solute.clone(), 0.2)
            .unwrap();

        assert_eq!(
            mixture.concentration_in_phase(&solute, MixturePhase::Aqueous),
            0.0
        );
        assert_eq!(
            mixture.concentration_in_phase(&solute, MixturePhase::Solid),
            0.2
        );
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
        assert!((small.liquid_phase_count(&registry).unwrap() as f64 - 1.0).abs() < 0.01);
        assert!(
            small
                .organic_phase_amounts_of(&registry, &solvent_b)
                .unwrap()
                .len()
                >= 1
        );

        let mut large = Mixture::new(298.0).unwrap();
        large
            .add_substance(&registry, solvent_a.clone(), 1.0)
            .unwrap();
        large
            .add_substance(&registry, solvent_b.clone(), 0.5)
            .unwrap();
        assert!((large.liquid_phase_count(&registry).unwrap() as f64 - 2.0).abs() < 0.5);
        assert!(
            (large
                .concentration_in_organic_solvent(&registry, &solvent_b, &solvent_b)
                .unwrap()
                - 0.5)
                .abs()
                < 0.01
        );
    }
}
