use crate::chemistry::mixture::{Mixture, DEFAULT_TEMPERATURE_KELVIN};
use crate::chemistry::substance::SubstanceId;

use super::peripheral::Peripheral;

#[derive(Debug, Clone)]
pub struct ReactorZone {
    mixture: Mixture,
    volume_cubic_meters: f64,
    sealed: bool,
    elapsed_seconds: f64,
    peripherals: Vec<Peripheral>,
}

impl ReactorZone {
    pub fn new(volume_cubic_meters: f64) -> Self {
        Self {
            mixture: Mixture::new(DEFAULT_TEMPERATURE_KELVIN).unwrap(),
            volume_cubic_meters,
            sealed: false,
            elapsed_seconds: 0.0,
            peripherals: Vec::new(),
        }
    }

    pub fn with_peripheral(mut self, peripheral: Peripheral) -> Self {
        self.peripherals.push(peripheral);
        self
    }

    pub fn add_peripheral(&mut self, peripheral: Peripheral) {
        self.peripherals.push(peripheral);
    }

    pub fn remove_peripheral(&mut self, name: &str) -> bool {
        let before = self.peripherals.len();
        self.peripherals.retain(|p| p.name() != name);
        self.peripherals.len() < before
    }

    pub fn peripherals(&self) -> &[Peripheral] {
        &self.peripherals
    }

    pub fn peripherals_mut(&mut self) -> &mut Vec<Peripheral> {
        &mut self.peripherals
    }

    pub fn total_uv_intensity(&self) -> f64 {
        self.peripherals.iter().map(|p| p.uv_intensity()).sum()
    }

    pub fn mixture(&self) -> &Mixture {
        &self.mixture
    }

    pub fn mixture_mut(&mut self) -> &mut Mixture {
        &mut self.mixture
    }

    pub fn volume_cubic_meters(&self) -> f64 {
        self.volume_cubic_meters
    }

    pub fn set_volume_cubic_meters(&mut self, volume: f64) {
        self.volume_cubic_meters = volume;
        self.mixture.set_gas_volume_cubic_meters(volume);
    }

    pub fn sealed(&self) -> bool {
        self.sealed
    }

    pub fn set_sealed(&mut self, sealed: bool) {
        self.sealed = sealed;
    }

    pub fn elapsed_seconds(&self) -> f64 {
        self.elapsed_seconds
    }

    pub fn temperature_kelvin(&self) -> f64 {
        self.mixture.temperature_kelvin()
    }

    pub fn pressure_pascal(&self) -> f64 {
        self.mixture.gas_pressure_pascal()
    }

    pub fn concentration_of(&self, substance_id: &SubstanceId) -> f64 {
        self.mixture.concentration_of(substance_id)
    }

    pub fn tick(&mut self, registry: &crate::chemistry::registry::ChemistryRegistry, dt_seconds: f64) {
        self.elapsed_seconds += dt_seconds;
        let mut peripherals = std::mem::take(&mut self.peripherals);
        for peripheral in &mut peripherals {
            peripheral.apply(self, registry, dt_seconds);
        }
        self.peripherals = peripherals;
    }
}
