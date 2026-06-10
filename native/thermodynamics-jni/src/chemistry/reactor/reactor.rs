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
}

impl Reactor {
    pub fn new() -> Self {
        Self {
            zones: Vec::new(),
            transitions: Vec::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
        }
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
        for zone in &mut self.zones {
            zone.mixture_mut().equilibrate_vapor_liquid(registry)?;
        }
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
            TransitionMode::Phases { entries } => {
                for entry in entries {
                    let max_amount = entry.rate_mol_per_second * dt_seconds;
                    if max_amount > 0.0 {
                        io::insert_from_phase(
                            &mut self.zones[zone_id.0],
                            registry,
                            entry.phase,
                            max_amount,
                        )?;
                    }
                }
            }
            TransitionMode::All {
                rate_mol_per_second,
            } => {
                let max_amount = rate_mol_per_second * dt_seconds;
                if max_amount > 0.0 {
                    io::insert_all(&mut self.zones[zone_id.0], registry, max_amount)?;
                }
            }
            TransitionMode::SubstancesThreshold {
                entries,
                threshold_mol_per_bucket,
            } => {
                for entry in entries {
                    let concentration = self.zones[zone_id.0].concentration_of(&entry.id);
                    let excess = concentration - threshold_mol_per_bucket;
                    if excess > 0.0 {
                        let max_amount = entry.rate_mol_per_second * dt_seconds;
                        let take = excess.min(max_amount);
                        io::insert_substance(
                            &mut self.zones[zone_id.0],
                            registry,
                            &entry.id,
                            take,
                        )?;
                    }
                }
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
                let (left, right) = self.zones.split_at_mut(from.0 + 1);
                let from_zone = &mut left[from.0];
                let to_zone = &mut right[to.0 - from.0 - 1];
                let extracted = io::extract_all(from_zone, registry)?;
                let total: f64 = extracted.iter().map(|(_, a)| a).sum();
                if total > 0.0 {
                    let scale = (max_amount / total).min(1.0);
                    for (id, amount) in &extracted {
                        let take = amount * scale;
                        if take > 0.0 {
                            io::insert_substance(to_zone, registry, id, take)?;
                        }
                    }
                    let remain = total - total * scale;
                    if remain > 0.0 {
                        for (id, amount) in &extracted {
                            let leftover = amount * (1.0 - scale);
                            if leftover > 0.0 {
                                io::insert_substance(from_zone, registry, id, leftover)?;
                            }
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
                let all_extracted = io::extract_all(&mut self.zones[zone_id.0], registry)?;
                let total: f64 = all_extracted.iter().map(|(_, a)| a).sum();
                if total > 0.0 && max_amount > 0.0 {
                    let scale = (max_amount / total).min(1.0);
                    for (id, amount) in &all_extracted {
                        let take = amount * scale;
                        if take > 0.0 {
                            extracted.push((id.clone(), take));
                        }
                    }
                    let remain: f64 = all_extracted
                        .iter()
                        .map(|(id, a)| {
                            let taken = extracted
                                .iter()
                                .find(|(eid, _)| eid == id)
                                .map(|(_, t)| *t)
                                .unwrap_or(0.0);
                            a - taken
                        })
                        .sum();
                    if remain > 0.0 {
                        for (id, amount) in &all_extracted {
                            let leftover = amount
                                - extracted
                                    .iter()
                                    .find(|(eid, _)| eid == id)
                                    .map(|(_, t)| *t)
                                    .unwrap_or(0.0);
                            if leftover > 0.0 {
                                io::insert_substance(
                                    &mut self.zones[zone_id.0],
                                    registry,
                                    id,
                                    leftover,
                                )?;
                            }
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
