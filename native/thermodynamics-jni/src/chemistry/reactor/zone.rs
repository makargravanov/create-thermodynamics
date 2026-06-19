use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::mixture::{
    Mixture, MixturePhase, DEFAULT_TEMPERATURE_KELVIN, TRACE_CONCENTRATION_MOL_PER_BUCKET,
};
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::SubstanceId;

use super::peripheral::Peripheral;

const MIN_HEADSPACE_CUBIC_METERS: f64 = 1.0e-9;
const VOLUME_TOLERANCE_CUBIC_METERS: f64 = 1.0e-12;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactorVolumeMode {
    HeadspaceRequired,
    LiquidFilled,
}

#[derive(Debug, Clone)]
pub struct ReactorZone {
    mixture: Mixture,
    volume_cubic_meters: f64,
    volume_mode: ReactorVolumeMode,
    sealed: bool,
    elapsed_seconds: f64,
    peripherals: Vec<Peripheral>,
}

impl ReactorZone {
    pub fn new(volume_cubic_meters: f64) -> ChemistryResult<Self> {
        let mut mixture = Mixture::new(DEFAULT_TEMPERATURE_KELVIN)?;
        mixture.set_gas_volume_cubic_meters(volume_cubic_meters)?;
        Ok(Self {
            mixture,
            volume_cubic_meters,
            volume_mode: ReactorVolumeMode::HeadspaceRequired,
            sealed: false,
            elapsed_seconds: 0.0,
            peripherals: Vec::new(),
        })
    }

    pub fn from_parts(
        mixture: Mixture,
        volume_cubic_meters: f64,
        volume_mode: ReactorVolumeMode,
        sealed: bool,
        elapsed_seconds: f64,
    ) -> ChemistryResult<Self> {
        if !volume_cubic_meters.is_finite() || volume_cubic_meters <= MIN_HEADSPACE_CUBIC_METERS {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "reactor zone volume must be finite and greater than {MIN_HEADSPACE_CUBIC_METERS}, got {volume_cubic_meters}"
            )));
        }
        if !elapsed_seconds.is_finite() || elapsed_seconds < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "reactor zone elapsed time must be non-negative and finite, got {elapsed_seconds}"
            )));
        }
        Ok(Self {
            mixture,
            volume_cubic_meters,
            volume_mode,
            sealed,
            elapsed_seconds,
            peripherals: Vec::new(),
        })
    }

    pub fn with_volume_mode(mut self, volume_mode: ReactorVolumeMode) -> Self {
        self.volume_mode = volume_mode;
        self
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

    pub fn total_electrical_draw_w(&self) -> f64 {
        self.peripherals.iter().map(|p| p.electrical_draw_w()).sum()
    }

    pub fn electrical_peripherals(&self) -> impl Iterator<Item = &Peripheral> {
        self.peripherals.iter().filter(|p| p.is_electrical())
    }

    pub fn electrical_peripherals_mut(&mut self) -> impl Iterator<Item = &mut Peripheral> {
        self.peripherals.iter_mut().filter(|p| p.is_electrical())
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

    pub fn volume_mode(&self) -> ReactorVolumeMode {
        self.volume_mode
    }

    pub fn set_volume_mode(&mut self, volume_mode: ReactorVolumeMode) {
        self.volume_mode = volume_mode;
    }

    pub fn set_volume_cubic_meters(
        &mut self,
        registry: &ChemistryRegistry,
        volume: f64,
    ) -> ChemistryResult<()> {
        if !volume.is_finite() || volume <= MIN_HEADSPACE_CUBIC_METERS {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "reactor zone volume must be finite and greater than {MIN_HEADSPACE_CUBIC_METERS}, got {volume}"
            )));
        }
        let old_volume = self.volume_cubic_meters;
        self.volume_cubic_meters = volume;
        let result = (|| {
            let condensed = self.condensed_volume_cubic_meters(registry)?;
            let headspace =
                self.validated_headspace_for_mixture(&self.mixture, condensed, "reactor zone")?;
            self.mixture.set_gas_volume_cubic_meters(headspace)
        })();
        if result.is_err() {
            self.volume_cubic_meters = old_volume;
        }
        result
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
        if self.non_condensed_mol_per_bucket() <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            0.0
        } else {
            self.mixture.gas_pressure_pascal()
        }
    }

    pub fn concentration_of(&self, substance_id: &SubstanceId) -> f64 {
        self.mixture.concentration_of(substance_id)
    }

    pub fn condensed_volume_cubic_meters(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<f64> {
        self.mixture.condensed_volume_cubic_meters(registry)
    }

    pub fn headspace_volume_cubic_meters(
        &self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<f64> {
        let condensed = self.condensed_volume_cubic_meters(registry)?;
        let headspace = self.volume_cubic_meters - condensed;
        if !headspace.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "reactor zone headspace is not finite: volume={}, condensed={condensed}",
                self.volume_cubic_meters
            )));
        }
        Ok(headspace.max(0.0))
    }

    pub fn refresh_headspace_volume(
        &mut self,
        registry: &ChemistryRegistry,
    ) -> ChemistryResult<()> {
        let condensed = self.condensed_volume_cubic_meters(registry)?;
        let headspace =
            self.validated_headspace_for_mixture(&self.mixture, condensed, "reactor zone")?;
        self.mixture.set_gas_volume_cubic_meters(headspace)
    }

    pub fn can_accept_substance(
        &self,
        registry: &ChemistryRegistry,
        substance_id: &SubstanceId,
        mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        if !mol_per_bucket.is_finite() || mol_per_bucket <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "inserted amount must be positive and finite, got {mol_per_bucket}"
            )));
        }
        let mut predicted = self.mixture.clone();
        match self.volume_mode {
            ReactorVolumeMode::HeadspaceRequired => {
                predicted.add_substance(registry, substance_id.clone(), mol_per_bucket)?;
            }
            ReactorVolumeMode::LiquidFilled => {
                predicted.add_substance_without_vapor_liquid_equilibrium(
                    registry,
                    substance_id.clone(),
                    mol_per_bucket,
                )?;
            }
        }
        let condensed_before_vle = predicted.condensed_volume_cubic_meters(registry)?;
        let headspace_before_vle = self.validated_headspace_for_mixture(
            &predicted,
            condensed_before_vle,
            "predicted reactor zone",
        )?;
        if headspace_before_vle > MIN_HEADSPACE_CUBIC_METERS {
            predicted.set_gas_volume_cubic_meters(headspace_before_vle)?;
            predicted.equilibrate_vapor_liquid(registry, 1.0)?;
        }
        let condensed = predicted.condensed_volume_cubic_meters(registry)?;
        self.validated_headspace_for_mixture(&predicted, condensed, "predicted reactor zone")
            .map_err(|error| {
                ChemistryError::InvalidMixtureState(format!(
                    "reactor zone cannot accept {mol_per_bucket} mol/bucket of '{substance_id}': {error}"
                ))
            })?;
        Ok(())
    }

    pub fn add_substance_checked(
        &mut self,
        registry: &ChemistryRegistry,
        substance_id: SubstanceId,
        mol_per_bucket: f64,
    ) -> ChemistryResult<()> {
        self.can_accept_substance(registry, &substance_id, mol_per_bucket)?;
        match self.volume_mode {
            ReactorVolumeMode::HeadspaceRequired => {
                self.mixture
                    .add_substance(registry, substance_id, mol_per_bucket)?;
            }
            ReactorVolumeMode::LiquidFilled => {
                self.mixture
                    .add_substance_without_vapor_liquid_equilibrium(
                        registry,
                        substance_id,
                        mol_per_bucket,
                    )?;
            }
        }
        self.refresh_headspace_volume(registry)
    }

    pub fn tick(&mut self, registry: &ChemistryRegistry, dt_seconds: f64) -> ChemistryResult<()> {
        self.elapsed_seconds += dt_seconds;
        let mut peripherals = std::mem::take(&mut self.peripherals);
        let mut result = Ok(());
        for peripheral in &mut peripherals {
            if let Err(error) = peripheral.apply(self, registry, dt_seconds) {
                result = Err(error);
                break;
            }
        }
        self.peripherals = peripherals;
        result?;
        self.refresh_headspace_volume(registry)
    }

    fn non_condensed_mol_per_bucket(&self) -> f64 {
        non_condensed_mol_per_bucket(&self.mixture)
    }

    fn validated_headspace_for_mixture(
        &self,
        mixture: &Mixture,
        condensed: f64,
        context: &str,
    ) -> ChemistryResult<f64> {
        let headspace = self.volume_cubic_meters - condensed;
        if !headspace.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "{context} headspace is not finite: volume={} m^3, condensed={condensed} m^3",
                self.volume_cubic_meters
            )));
        }
        if headspace > MIN_HEADSPACE_CUBIC_METERS {
            return Ok(headspace);
        }
        if headspace < -VOLUME_TOLERANCE_CUBIC_METERS {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "{context} is overfilled: volume={} m^3, condensed={condensed} m^3",
                self.volume_cubic_meters
            )));
        }
        match self.volume_mode {
            ReactorVolumeMode::HeadspaceRequired => Err(ChemistryError::InvalidMixtureState(
                format!(
                    "{context} requires positive headspace>{MIN_HEADSPACE_CUBIC_METERS} m^3: volume={} m^3, condensed={condensed} m^3",
                    self.volume_cubic_meters
                ),
            )),
            ReactorVolumeMode::LiquidFilled => {
                let gas = non_condensed_mol_per_bucket(mixture);
                if gas > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "{context} is liquid-filled but contains {gas} mol/bucket of gas or supercritical fluid"
                    )));
                }
                Ok(MIN_HEADSPACE_CUBIC_METERS)
            }
        }
    }
}

fn non_condensed_mol_per_bucket(mixture: &Mixture) -> f64 {
    mixture.total_concentration_in_phase(MixturePhase::Gas)
        + mixture.total_concentration_in_phase(MixturePhase::SupercriticalFluid)
}
