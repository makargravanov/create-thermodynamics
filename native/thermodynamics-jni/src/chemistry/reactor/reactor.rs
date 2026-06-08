use crate::chemistry::error::ChemistryResult;
use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::SubstanceId;

use super::io;
use super::zone::ReactorZone;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ZoneId(pub usize);

#[derive(Debug, Clone)]
pub struct ZoneTransition {
    pub from: ZoneId,
    pub to: ZoneId,
    pub mode: TransitionMode,
}

#[derive(Debug, Clone)]
pub enum TransitionMode {
    /// Fixed rate transfer of specific substances (mol per second each)
    Substances { ids: Vec<SubstanceId>, rate_mol_per_second: f64 },
    /// Fixed rate transfer of entire phases (mol per second each)
    Phases { phases: Vec<MixturePhase>, rate_mol_per_second: f64 },
    /// Transfer everything from source zone
    All { rate_mol_per_second: f64 },
    /// Transfer substances above a concentration threshold
    SubstancesThreshold { ids: Vec<SubstanceId>, threshold_mol_per_bucket: f64, rate_mol_per_second: f64 },
}

impl TransitionMode {
    pub fn rate(&self) -> f64 {
        match self {
            TransitionMode::Substances { rate_mol_per_second, .. }
            | TransitionMode::Phases { rate_mol_per_second, .. }
            | TransitionMode::All { rate_mol_per_second }
            | TransitionMode::SubstancesThreshold { rate_mol_per_second, .. } => *rate_mol_per_second,
        }
    }
}

pub struct Reactor {
    zones: Vec<ReactorZone>,
    transitions: Vec<ZoneTransition>,
}

impl Reactor {
    pub fn new() -> Self {
        Self {
            zones: Vec::new(),
            transitions: Vec::new(),
        }
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

    pub fn add_transition(&mut self, transition: ZoneTransition) {
        self.transitions.push(transition);
    }

    pub fn transitions(&self) -> &[ZoneTransition] {
        &self.transitions
    }

    pub fn tick(&mut self, registry: &ChemistryRegistry, dt_seconds: f64) -> ChemistryResult<()> {
        let transitions: Vec<ZoneTransition> = self.transitions.clone();
        for transition in &transitions {
            self.apply_transition(registry, transition, dt_seconds)?;
        }
        for zone in &mut self.zones {
            zone.tick(registry, dt_seconds);
        }
        Ok(())
    }

    fn apply_transition(
        &mut self,
        registry: &ChemistryRegistry,
        transition: &ZoneTransition,
        dt_seconds: f64,
    ) -> ChemistryResult<()> {
        let max_amount = transition.mode.rate() * dt_seconds;
        if max_amount <= 0.0 {
            return Ok(());
        }

        match &transition.mode {
            TransitionMode::Substances { ids, .. } => {
                for id in ids {
                    let available = self.zones[transition.from.0].concentration_of(id);
                    let take = available.min(max_amount);
                    if take > 0.0 {
                        io::extract_substance(
                            &mut self.zones[transition.from.0], registry, id, take,
                        )?;
                        io::insert_substance(
                            &mut self.zones[transition.to.0], registry, id, take,
                        )?;
                    }
                }
            }
            TransitionMode::Phases { phases, .. } => {
                for &phase in phases {
                    io::extract_from_phase(
                        &mut self.zones[transition.from.0], registry, phase, max_amount,
                    )?
                    .into_iter()
                    .try_for_each(|(id, amount)| -> ChemistryResult<()> {
                        io::insert_substance(
                            &mut self.zones[transition.to.0], registry, &id, amount,
                        )?;
                        Ok(())
                    })?;
                }
            }
            TransitionMode::All { .. } => {
                let (left, right) = self.zones.split_at_mut(transition.from.0 + 1);
                let from = &mut left[transition.from.0];
                let to = &mut right[transition.to.0 - transition.from.0 - 1];
                io::transfer_all(from, to, registry)?;
            }
            TransitionMode::SubstancesThreshold { ids, threshold_mol_per_bucket, .. } => {
                for id in ids {
                    let concentration = self.zones[transition.from.0].concentration_of(id);
                    if concentration <= *threshold_mol_per_bucket {
                        continue;
                    }
                    let excess = concentration - threshold_mol_per_bucket;
                    let take = excess.min(max_amount);
                    if take > 0.0 {
                        io::extract_substance(
                            &mut self.zones[transition.from.0], registry, id, take,
                        )?;
                        io::insert_substance(
                            &mut self.zones[transition.to.0], registry, id, take,
                        )?;
                    }
                }
            }
        }
        Ok(())
    }
}
