use crate::chemistry::error::ChemistryResult;
use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::SubstanceId;

use super::io;
use super::zone::ReactorZone;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ZoneId(pub usize);

#[derive(Debug, Clone)]
pub struct SubstanceEntry {
    pub id: SubstanceId,
    pub rate_mol_per_second: f64,
}

#[derive(Debug, Clone)]
pub struct PhaseEntry {
    pub phase: MixturePhase,
    pub rate_mol_per_second: f64,
}

#[derive(Debug, Clone)]
pub enum TransitionMode {
    Substances {
        entries: Vec<SubstanceEntry>,
    },
    Phases {
        entries: Vec<PhaseEntry>,
    },
    All {
        rate_mol_per_second: f64,
    },
    SubstancesThreshold {
        entries: Vec<SubstanceEntry>,
        threshold_mol_per_bucket: f64,
    },
}



#[derive(Debug, Clone)]
pub struct ZoneTransition {
    pub from: ZoneId,
    pub to: ZoneId,
    pub mode: TransitionMode,
    pub enabled: bool,
}

impl ZoneTransition {
    pub fn new(from: ZoneId, to: ZoneId, mode: TransitionMode) -> Self {
        Self {
            from,
            to,
            mode,
            enabled: true,
        }
    }

    pub fn with_enabled(mut self, enabled: bool) -> Self {
        self.enabled = enabled;
        self
    }
}

#[derive(Debug, Clone)]
pub struct Input {
    pub to: ZoneId,
    pub mode: TransitionMode,
    pub enabled: bool,
}

impl Input {
    pub fn new(to: ZoneId, mode: TransitionMode) -> Self {
        Self {
            to,
            mode,
            enabled: true,
        }
    }

    pub fn with_enabled(mut self, enabled: bool) -> Self {
        self.enabled = enabled;
        self
    }
}

#[derive(Debug, Clone)]
pub struct Output {
    pub from: ZoneId,
    pub mode: TransitionMode,
    pub enabled: bool,
}

impl Output {
    pub fn new(from: ZoneId, mode: TransitionMode) -> Self {
        Self {
            from,
            mode,
            enabled: true,
        }
    }

    pub fn with_enabled(mut self, enabled: bool) -> Self {
        self.enabled = enabled;
        self
    }
}

pub struct Reactor {
    zones: Vec<ReactorZone>,
    transitions: Vec<ZoneTransition>,
    inputs: Vec<Input>,
    outputs: Vec<Output>,
    ambient_temperature_kelvin: Option<f64>,
    heat_transfer_coefficient_kw_per_kelvin: Option<f64>,
    last_ambient_energy_exchange_j: f64,
    vle_iterations: usize,
}

impl Reactor {
    pub fn new() -> Self {
        Self {
            zones: Vec::new(),
            transitions: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
            ambient_temperature_kelvin: None,
            heat_transfer_coefficient_kw_per_kelvin: None,
            last_ambient_energy_exchange_j: 0.0,
            vle_iterations: 1,
        }
    }

    pub fn set_vle_iterations(&mut self, iterations: usize) {
        self.vle_iterations = iterations.max(1);
    }

    pub fn vle_iterations(&self) -> usize {
        self.vle_iterations
    }

    pub fn add_input(&mut self, input: Input) -> usize {
        self.inputs.push(input);
        self.inputs.len() - 1
    }

    pub fn input(&self, index: usize) -> Option<&Input> {
        self.inputs.get(index)
    }

    pub fn input_mut(&mut self, index: usize) -> Option<&mut Input> {
        self.inputs.get_mut(index)
    }

    pub fn input_count(&self) -> usize {
        self.inputs.len()
    }

    pub fn add_output(&mut self, output: Output) -> usize {
        self.outputs.push(output);
        self.outputs.len() - 1
    }

    pub fn output(&self, index: usize) -> Option<&Output> {
        self.outputs.get(index)
    }

    pub fn output_mut(&mut self, index: usize) -> Option<&mut Output> {
        self.outputs.get_mut(index)
    }

    pub fn output_count(&self) -> usize {
        self.outputs.len()
    }

    pub fn add_zone(&mut self, zone: ReactorZone) -> ZoneId {
        let id = ZoneId(self.zones.len());
        self.zones.push(zone);
        id
    }

    pub fn zone(&self, id: &ZoneId) -> Option<&ReactorZone> {
        self.zones.get(id.0)
    }

    pub fn zone_mut(&mut self, id: &ZoneId) -> Option<&mut ReactorZone> {
        self.zones.get_mut(id.0)
    }

    pub fn zone_count(&self) -> usize {
        self.zones.len()
    }

    pub fn add_transition(&mut self, transition: ZoneTransition) -> usize {
        self.transitions.push(transition);
        self.transitions.len() - 1
    }

    pub fn transition(&self, index: usize) -> Option<&ZoneTransition> {
        self.transitions.get(index)
    }

    pub fn transition_mut(&mut self, index: usize) -> Option<&mut ZoneTransition> {
        self.transitions.get_mut(index)
    }

    pub fn transition_count(&self) -> usize {
        self.transitions.len()
    }

    pub fn transitions(&self) -> &[ZoneTransition] {
        &self.transitions
    }

    pub fn inputs(&self) -> &[Input] {
        &self.inputs
    }

    pub fn outputs(&self) -> &[Output] {
        &self.outputs
    }

    pub fn set_ambient_temperature(&mut self, temperature_kelvin: f64) {
        self.ambient_temperature_kelvin = Some(temperature_kelvin);
    }

    pub fn ambient_temperature(&self) -> Option<f64> {
        self.ambient_temperature_kelvin
    }

    pub fn set_heat_transfer_coefficient(&mut self, u_kw_per_kelvin: f64) {
        self.heat_transfer_coefficient_kw_per_kelvin = Some(u_kw_per_kelvin);
    }

    pub fn heat_transfer_coefficient(&self) -> Option<f64> {
        self.heat_transfer_coefficient_kw_per_kelvin
    }

    pub fn last_ambient_energy_exchange_j(&self) -> f64 {
        self.last_ambient_energy_exchange_j
    }

    pub fn total_electrical_draw_w(&self) -> f64 {
        self.zones.iter().map(|z| z.total_electrical_draw_w()).sum()
    }

    pub fn tick(&mut self, registry: &ChemistryRegistry, dt_seconds: f64) -> ChemistryResult<()> {
        let transitions: Vec<ZoneTransition> = self.transitions.clone();
        for transition in &transitions {
            if !transition.enabled {
                continue;
            }
            self.apply_transition_mode(
                registry,
                &transition.from,
                &transition.to,
                &transition.mode,
                dt_seconds,
            )?;
        }
        for zone in &mut self.zones {
            zone.tick(registry, dt_seconds);
        }
        for _ in 0..self.vle_iterations {
            for zone in &mut self.zones {
                zone.mixture_mut().equilibrate_vapor_liquid(registry)?;
            }
        }
        self.apply_ambient_heat_exchange(registry, dt_seconds)?;
        Ok(())
    }

    fn apply_ambient_heat_exchange(
        &mut self,
        registry: &ChemistryRegistry,
        dt_seconds: f64,
    ) -> ChemistryResult<()> {
        let (ambient_temp, u_kw) = match (
            self.ambient_temperature_kelvin,
            self.heat_transfer_coefficient_kw_per_kelvin,
        ) {
            (Some(t), Some(u)) => (t, u),
            _ => {
                self.last_ambient_energy_exchange_j = 0.0;
                return Ok(());
            }
        };

        let u_w_per_k = u_kw * 1000.0;
        let mut total_energy_j = 0.0;

        for zone in &mut self.zones {
            let zone_temp = zone.temperature_kelvin();
            let delta_t = ambient_temp - zone_temp;
            if delta_t.abs() < 0.01 {
                continue;
            }

            let heat_capacity = match zone
                .mixture()
                .volumetric_heat_capacity_j_per_bucket_kelvin(registry)
            {
                Ok(hc) if hc > 0.0 => hc,
                _ => continue,
            };

            let max_energy = u_w_per_k * delta_t.abs() * dt_seconds;
            let energy_to_equilibrium = delta_t * heat_capacity;
            let energy = energy_to_equilibrium.clamp(-max_energy, max_energy);

            if energy.abs() > 1.0e-12 {
                let _ = zone.mixture_mut().heat(registry, energy);
                total_energy_j += energy;
            }
        }

        self.last_ambient_energy_exchange_j = total_energy_j;
        Ok(())
    }

    pub fn apply_input(
        &mut self,
        registry: &ChemistryRegistry,
        input_index: usize,
        dt_seconds: f64,
    ) -> ChemistryResult<()> {
        let input = self
            .inputs
            .get(input_index)
            .ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!(
                    "input index {} out of range (have {})",
                    input_index,
                    self.inputs.len()
                ))
            })?
            .clone();
        if !input.enabled {
            return Ok(());
        }
        let zone_id = input.to;
        match &input.mode {
            TransitionMode::Substances { entries } => {
                for entry in entries {
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    if max_amount > 0.0 {
                        io::insert_substance(
                            &mut self.zones[zone_id.0],
                            registry,
                            &entry.id,
                            max_amount,
                        )?;
                    }
                }
            }
            TransitionMode::Phases { .. }
            | TransitionMode::All { .. }
            | TransitionMode::SubstancesThreshold { .. } => {
                return Err(ChemistryError::InvalidMixtureState(
                    "input mode must be Substances — Phases/All/Threshold are not valid for external input"
                        .to_string(),
                ));
            }
        }
        Ok(())
    }

    pub fn apply_output(
        &mut self,
        registry: &ChemistryRegistry,
        output_index: usize,
        dt_seconds: f64,
    ) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
        let output = self
            .outputs
            .get(output_index)
            .ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!(
                    "output index {} out of range (have {})",
                    output_index,
                    self.outputs.len()
                ))
            })?
            .clone();
        if !output.enabled {
            return Ok(Vec::new());
        }
        let zone_id = output.from;
        self.apply_output_mode(registry, &zone_id, &output.mode, dt_seconds)
    }

    fn apply_transition_mode(
        &mut self,
        registry: &ChemistryRegistry,
        from: &ZoneId,
        to: &ZoneId,
        mode: &TransitionMode,
        dt_seconds: f64,
    ) -> ChemistryResult<()> {
        match mode {
            TransitionMode::Substances { entries } => {
                for entry in entries {
                    let available = self.zones[from.0].concentration_of(&entry.id);
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    let take = available.min(max_amount);
                    if take > 0.0 {
                        io::extract_substance(&mut self.zones[from.0], registry, &entry.id, take)?;
                        io::insert_substance(&mut self.zones[to.0], registry, &entry.id, take)?;
                    }
                }
            }
            TransitionMode::Phases { entries } => {
                for entry in entries {
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    if max_amount <= 0.0 {
                        continue;
                    }
                    io::extract_from_phase(
                        &mut self.zones[from.0],
                        registry,
                        entry.phase,
                        max_amount,
                    )?
                    .into_iter()
                    .try_for_each(|(id, amount)| -> ChemistryResult<()> {
                        io::insert_substance(&mut self.zones[to.0], registry, &id, amount)?;
                        Ok(())
                    })?;
                }
            }
            TransitionMode::All {
                rate_mol_per_second,
            } => {
                let max_amount = rate_mol_per_second * dt_seconds;
                if max_amount <= 0.0 {
                    return Ok(());
                }
                let snapshot = io::mixture_snapshot(&self.zones[from.0]);
                let total: f64 = snapshot
                    .substances
                    .iter()
                    .map(|s| s.total_mol_per_bucket)
                    .sum();
                if total <= 0.0 {
                    return Ok(());
                }
                let scale = (max_amount / total).min(1.0);
                for component in &snapshot.substances {
                    let take = component.total_mol_per_bucket * scale;
                    if take > 0.0 {
                        let amount = io::extract_substance(
                            &mut self.zones[from.0],
                            registry,
                            &component.id,
                            take,
                        )?;
                        if amount > 0.0 {
                            io::insert_substance(
                                &mut self.zones[to.0],
                                registry,
                                &component.id,
                                amount,
                            )?;
                        }
                    }
                }
            }
            TransitionMode::SubstancesThreshold {
                entries,
                threshold_mol_per_bucket,
            } => {
                for entry in entries {
                    let concentration = self.zones[from.0].concentration_of(&entry.id);
                    if concentration <= *threshold_mol_per_bucket {
                        continue;
                    }
                    let excess = concentration - threshold_mol_per_bucket;
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    let take = excess.min(max_amount);
                    if take > 0.0 {
                        io::extract_substance(&mut self.zones[from.0], registry, &entry.id, take)?;
                        io::insert_substance(&mut self.zones[to.0], registry, &entry.id, take)?;
                    }
                }
            }
        }
        Ok(())
    }

    fn apply_output_mode(
        &mut self,
        registry: &ChemistryRegistry,
        zone_id: &ZoneId,
        mode: &TransitionMode,
        dt_seconds: f64,
    ) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
        let mut extracted = Vec::new();
        match mode {
            TransitionMode::Substances { entries } => {
                for entry in entries {
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    let amount = io::extract_substance(
                        &mut self.zones[zone_id.0],
                        registry,
                        &entry.id,
                        max_amount,
                    )?;
                    if amount > 0.0 {
                        extracted.push((entry.id.clone(), amount));
                    }
                }
            }
            TransitionMode::Phases { entries } => {
                for entry in entries {
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    if max_amount <= 0.0 {
                        continue;
                    }
                    let phase_extracted = io::extract_from_phase(
                        &mut self.zones[zone_id.0],
                        registry,
                        entry.phase,
                        max_amount,
                    )?;
                    extracted.extend(phase_extracted);
                }
            }
            TransitionMode::All {
                rate_mol_per_second,
            } => {
                let max_amount = rate_mol_per_second * dt_seconds;
                if max_amount <= 0.0 {
                    return Ok(extracted);
                }
                let snapshot = io::mixture_snapshot(&self.zones[zone_id.0]);
                let total: f64 = snapshot
                    .substances
                    .iter()
                    .map(|s| s.total_mol_per_bucket)
                    .sum();
                if total <= 0.0 {
                    return Ok(extracted);
                }
                let scale = (max_amount / total).min(1.0);
                for component in &snapshot.substances {
                    let take = component.total_mol_per_bucket * scale;
                    if take > 0.0 {
                        let amount = io::extract_substance(
                            &mut self.zones[zone_id.0],
                            registry,
                            &component.id,
                            take,
                        )?;
                        if amount > 0.0 {
                            extracted.push((component.id.clone(), amount));
                        }
                    }
                }
            }
            TransitionMode::SubstancesThreshold {
                entries,
                threshold_mol_per_bucket,
            } => {
                for entry in entries {
                    let concentration = self.zones[zone_id.0].concentration_of(&entry.id);
                    if concentration <= *threshold_mol_per_bucket {
                        continue;
                    }
                    let excess = concentration - threshold_mol_per_bucket;
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    let take = excess.min(max_amount);
                    if take > 0.0 {
                        let amount = io::extract_substance(
                            &mut self.zones[zone_id.0],
                            registry,
                            &entry.id,
                            take,
                        )?;
                        if amount > 0.0 {
                            extracted.push((entry.id.clone(), amount));
                        }
                    }
                }
            }
        }
        Ok(extracted)
    }
}

use crate::chemistry::error::ChemistryError;
