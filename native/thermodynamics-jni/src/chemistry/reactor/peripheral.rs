use crate::chemistry::error::ChemistryResult;
use crate::chemistry::registry::ChemistryRegistry;

use super::zone::ReactorZone;

#[derive(Debug, Clone)]
pub enum Peripheral {
    Heater {
        power_kw: f64,
    },
    SmartHeater(SmartHeaterState),
    HeatExchanger {
        coolant_temperature_kelvin: f64,
        u_kw_per_kelvin: f64,
    },
    UVLamp {
        intensity: f64,
    },
    Electrode(ElectrodeState),
}

#[derive(Debug, Clone)]
pub struct SmartHeaterState {
    pub target_kelvin: f64,
    pub resistance_at_ref_ohm: f64,
    pub reference_temperature_kelvin: f64,
    pub temperature_coefficient_per_kelvin: f64,
    pub melting_temperature_kelvin: f64,
    pub voltage: f64,
    pub mode: ControlMode,
    pub integral_error: f64,
    pub prev_error: f64,
    pub last_electrical_draw_w: f64,
    pub last_heating_power_w: f64,
    pub last_resistance_ohm: f64,
    pub failed: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ControlMode {
    P,
    PD,
    PID,
}

#[derive(Debug, Clone)]
pub struct ElectrodeState {
    pub voltage: f64,
    pub max_current_a: f64,
    pub current_a: f64,
    pub electrode_potential_v: f64,
    pub last_electrical_draw_w: f64,
    pub last_energy_delivered_j: f64,
    pub failed: bool,
}

impl SmartHeaterState {
    pub fn new(
        target_kelvin: f64,
        resistance_at_ref_ohm: f64,
        reference_temperature_kelvin: f64,
        temperature_coefficient_per_kelvin: f64,
        melting_temperature_kelvin: f64,
        voltage: f64,
        mode: ControlMode,
    ) -> Self {
        Self {
            target_kelvin,
            resistance_at_ref_ohm,
            reference_temperature_kelvin,
            temperature_coefficient_per_kelvin,
            melting_temperature_kelvin,
            voltage,
            mode,
            integral_error: 0.0,
            prev_error: 0.0,
            last_electrical_draw_w: 0.0,
            last_heating_power_w: 0.0,
            last_resistance_ohm: resistance_at_ref_ohm,
            failed: false,
        }
    }

    pub fn resistance_at(&self, temperature_kelvin: f64) -> f64 {
        let delta_t = temperature_kelvin - self.reference_temperature_kelvin;
        self.resistance_at_ref_ohm * (1.0 + self.temperature_coefficient_per_kelvin * delta_t)
    }

    pub fn power_at(&self, temperature_kelvin: f64) -> f64 {
        let r = self.resistance_at(temperature_kelvin);
        if r <= 0.0 {
            return 0.0;
        }
        self.voltage * self.voltage / r
    }

    pub fn set_voltage(&mut self, voltage: f64) {
        self.voltage = voltage.max(0.0);
    }

    pub fn voltage(&self) -> f64 {
        self.voltage
    }

    pub fn electrical_draw_w(&self) -> f64 {
        self.last_electrical_draw_w
    }

    pub fn heating_power_w(&self) -> f64 {
        self.last_heating_power_w
    }

    pub fn last_resistance_ohm(&self) -> f64 {
        self.last_resistance_ohm
    }
}

impl ElectrodeState {
    pub fn new(voltage: f64, max_current_a: f64, electrode_potential_v: f64) -> Self {
        Self {
            voltage,
            max_current_a,
            current_a: 0.0,
            electrode_potential_v,
            last_electrical_draw_w: 0.0,
            last_energy_delivered_j: 0.0,
            failed: false,
        }
    }

    pub fn with_current(mut self, current_a: f64) -> Self {
        self.current_a = current_a.clamp(0.0, self.max_current_a);
        self
    }

    pub fn set_voltage(&mut self, voltage: f64) {
        self.voltage = voltage.max(0.0);
    }

    pub fn voltage(&self) -> f64 {
        self.voltage
    }

    pub fn set_current(&mut self, current_a: f64) {
        self.current_a = current_a.clamp(0.0, self.max_current_a);
    }

    pub fn current_a(&self) -> f64 {
        self.current_a
    }

    pub fn max_current_a(&self) -> f64 {
        self.max_current_a
    }

    pub fn electrical_draw_w(&self) -> f64 {
        self.last_electrical_draw_w
    }

    pub fn energy_delivered_j(&self) -> f64 {
        self.last_energy_delivered_j
    }

    /// Electrical power drawn from the grid: P = V × I
    pub fn electrical_power_w(&self) -> f64 {
        self.voltage * self.current_a
    }

    /// Faraday constant for electrochemistry
    pub fn coulombs_this_tick(&self, dt_seconds: f64) -> f64 {
        self.current_a * dt_seconds
    }

    /// Moles of electrons transferred this tick
    pub fn moles_of_electrons(&self, dt_seconds: f64) -> f64 {
        const FARADAY: f64 = 96485.0;
        self.coulombs_this_tick(dt_seconds) / FARADAY
    }
}

impl Peripheral {
    pub fn name(&self) -> &str {
        match self {
            Peripheral::Heater { .. } => "heater",
            Peripheral::SmartHeater { .. } => "smart_heater",
            Peripheral::HeatExchanger { .. } => "heat_exchanger",
            Peripheral::UVLamp { .. } => "uv_lamp",
            Peripheral::Electrode { .. } => "electrode",
        }
    }

    pub fn apply(
        &mut self,
        zone: &mut ReactorZone,
        registry: &ChemistryRegistry,
        dt_seconds: f64,
    ) -> ChemistryResult<()> {
        match self {
            Peripheral::Heater { power_kw } => {
                let energy = *power_kw * 1000.0 * dt_seconds;
                zone.mixture_mut().heat(registry, energy)?;
            }
            Peripheral::SmartHeater(state) => {
                apply_smart_heater(zone, registry, state, dt_seconds)?;
            }
            Peripheral::HeatExchanger {
                coolant_temperature_kelvin,
                u_kw_per_kelvin,
            } => {
                apply_heat_exchanger(
                    zone,
                    registry,
                    *coolant_temperature_kelvin,
                    *u_kw_per_kelvin,
                    dt_seconds,
                )?;
            }
            Peripheral::UVLamp { intensity: _ } => {}
            Peripheral::Electrode(state) => {
                apply_electrode(zone, registry, state, dt_seconds)?;
            }
        }
        Ok(())
    }

    pub fn uv_intensity(&self) -> f64 {
        match self {
            Peripheral::UVLamp { intensity } => *intensity,
            _ => 0.0,
        }
    }

    pub fn electrical_draw_w(&self) -> f64 {
        match self {
            Peripheral::SmartHeater(state) => state.electrical_draw_w(),
            Peripheral::Electrode(state) => state.electrical_draw_w(),
            _ => 0.0,
        }
    }

    pub fn is_electrical(&self) -> bool {
        matches!(self, Peripheral::SmartHeater(_) | Peripheral::Electrode(_))
    }
}

fn apply_heat_exchanger(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    coolant_temperature_kelvin: f64,
    u_kw_per_kelvin: f64,
    dt_seconds: f64,
) -> ChemistryResult<()> {
    let mixture_temperature = zone.temperature_kelvin();
    let delta_t = coolant_temperature_kelvin - mixture_temperature;

    if delta_t.abs() < 0.01 {
        return Ok(());
    }

    let heat_capacity = match zone
        .mixture()
        .volumetric_heat_capacity_j_per_bucket_kelvin(registry)
    {
        Ok(hc) if hc > 0.0 => hc,
        _ => return Ok(()),
    };

    let max_energy = u_kw_per_kelvin * 1000.0 * delta_t.abs() * dt_seconds;
    let energy_to_equilibrium = delta_t * heat_capacity;
    let energy = energy_to_equilibrium.clamp(-max_energy, max_energy);

    zone.mixture_mut().heat(registry, energy)
}

fn apply_smart_heater(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    state: &mut SmartHeaterState,
    dt_seconds: f64,
) -> ChemistryResult<()> {
    if state.failed {
        state.last_electrical_draw_w = 0.0;
        state.last_heating_power_w = 0.0;
        return Ok(());
    }

    let mixture_temperature = zone.temperature_kelvin();

    if mixture_temperature >= state.melting_temperature_kelvin {
        state.failed = true;
        state.last_electrical_draw_w = 0.0;
        state.last_heating_power_w = 0.0;
        return Ok(());
    }

    let resistance = state.resistance_at(mixture_temperature);
    state.last_resistance_ohm = resistance;

    if resistance <= 0.0 {
        state.failed = true;
        state.last_electrical_draw_w = 0.0;
        state.last_heating_power_w = 0.0;
        return Ok(());
    }

    let electrical_power_w = state.voltage * state.voltage / resistance;

    let error = state.target_kelvin - mixture_temperature;

    if error.abs() < 0.01 {
        state.integral_error = 0.0;
        state.prev_error = error;
        state.last_electrical_draw_w = 0.0;
        state.last_heating_power_w = 0.0;
        return Ok(());
    }

    let _heat_capacity = match zone
        .mixture()
        .volumetric_heat_capacity_j_per_bucket_kelvin(registry)
    {
        Ok(hc) if hc > 0.0 => hc,
        _ => {
            state.last_electrical_draw_w = 0.0;
            state.last_heating_power_w = 0.0;
            return Ok(());
        }
    };

    let kp = 0.01;
    let kd = 0.005;
    let ki = 0.001;

    let proportional = kp * error;

    let derivative = match state.mode {
        ControlMode::P => 0.0,
        ControlMode::PD | ControlMode::PID => {
            let d = (error - state.prev_error) / dt_seconds.max(0.001);
            state.prev_error = error;
            kd * d
        }
    };

    let integral = match state.mode {
        ControlMode::P | ControlMode::PD => {
            state.integral_error = 0.0;
            0.0
        }
        ControlMode::PID => {
            state.integral_error += error * dt_seconds;
            state.integral_error = state.integral_error.clamp(-100.0, 100.0);
            ki * state.integral_error
        }
    };

    let duty_cycle = (proportional + derivative + integral).clamp(0.0, 1.0);

    let actual_electrical_w = electrical_power_w * duty_cycle;

    state.last_electrical_draw_w = actual_electrical_w;
    state.last_heating_power_w = actual_electrical_w;

    let energy = actual_electrical_w * dt_seconds;
    zone.mixture_mut().heat(registry, energy)
}

fn apply_electrode(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    state: &mut ElectrodeState,
    dt_seconds: f64,
) -> ChemistryResult<()> {
    if state.failed || state.current_a <= 0.0 {
        state.last_electrical_draw_w = 0.0;
        state.last_energy_delivered_j = 0.0;
        return Ok(());
    }

    let current = state.current_a.min(state.max_current_a);
    let electrical_w = state.voltage * current;
    let energy_j = electrical_w * dt_seconds;

    state.last_electrical_draw_w = electrical_w;
    state.last_energy_delivered_j = energy_j;

    zone.mixture_mut().heat(registry, energy_j)
}
