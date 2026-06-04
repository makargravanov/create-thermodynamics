use std::collections::{BTreeMap, BTreeSet};

use crate::chemistry::alloy::AlloyPhaseSnapshot;
use crate::chemistry::error::{ChemistryError, ChemistryResult};

use super::constants::{GAS_CONSTANT_J_PER_MOL_KELVIN, TRACE_COMPONENT_FRACTION};
use super::validation::{
    validate_kinetic_model, validate_non_negative_finite, validate_phase_boundary_point,
    validate_phase_model, validate_positive_finite,
};

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MetallurgicalComponentId(String);

impl MetallurgicalComponentId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(ChemistryError::InvalidMixtureState(
                "metallurgical component id must not be empty".to_string(),
            ));
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl From<&str> for MetallurgicalComponentId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

impl From<String> for MetallurgicalComponentId {
    fn from(value: String) -> Self {
        Self(value)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalComposition {
    pub components: BTreeMap<MetallurgicalComponentId, f64>,
}

impl MetallurgicalComposition {
    pub fn from_alloy_phase(alloy: &AlloyPhaseSnapshot) -> ChemistryResult<Self> {
        let mut components = BTreeMap::new();
        for constituent in &alloy.constituents {
            let component =
                MetallurgicalComponentId::new(constituent.metallurgical_component_id.clone())?;
            if constituent.mole_fraction > TRACE_COMPONENT_FRACTION {
                components.insert(component, constituent.mole_fraction);
            }
        }
        let composition = Self { components };
        composition.validate()?;
        Ok(composition)
    }

    pub fn fraction_of(&self, component: &MetallurgicalComponentId) -> f64 {
        self.components.get(component).copied().unwrap_or(0.0)
    }

    pub fn validate(&self) -> ChemistryResult<()> {
        if self.components.is_empty() {
            return Err(ChemistryError::InvalidMixtureState(
                "metallurgical composition must contain at least one component".to_string(),
            ));
        }
        let total = self.components.values().sum::<f64>();
        if !total.is_finite() || total <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "metallurgical composition total must be positive and finite".to_string(),
            ));
        }
        for (component, fraction) in &self.components {
            if component.as_str().trim().is_empty() {
                return Err(ChemistryError::InvalidMixtureState(
                    "metallurgical component id must not be empty".to_string(),
                ));
            }
            if !fraction.is_finite() || *fraction < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(format!(
                    "metallurgical component '{}' fraction must be non-negative and finite",
                    component.as_str()
                )));
            }
        }
        if (total - 1.0).abs() > 1.0e-6 {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "metallurgical composition fractions must sum to 1.0, got {total}"
            )));
        }
        Ok(())
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum MetallurgicalPhaseKind {
    Liquid,
    SolidSolution,
    Intermetallic,
    Ferrite,
    Austenite,
    Cementite,
    Graphite,
    Martensite,
    Pearlite,
    Bainite,
    TemperedMartensite,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ComponentLimit {
    pub component: MetallurgicalComponentId,
    pub min_fraction: f64,
    pub max_fraction: f64,
}

impl ComponentLimit {
    pub fn new(
        component: impl Into<MetallurgicalComponentId>,
        min_fraction: f64,
        max_fraction: f64,
    ) -> Self {
        Self {
            component: component.into(),
            min_fraction,
            max_fraction,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CompositionEnergyTerm {
    pub component: MetallurgicalComponentId,
    pub center_fraction: f64,
    pub width_fraction: f64,
    pub penalty_j_per_mol: f64,
}

impl CompositionEnergyTerm {
    pub fn new(
        component: impl Into<MetallurgicalComponentId>,
        center_fraction: f64,
        width_fraction: f64,
        penalty_j_per_mol: f64,
    ) -> Self {
        Self {
            component: component.into(),
            center_fraction,
            width_fraction,
            penalty_j_per_mol,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PhaseFreeEnergyModel {
    pub reference_gibbs_j_per_mol: f64,
    pub entropy_j_per_mol_kelvin: f64,
    pub low_temperature_kelvin: f64,
    pub high_temperature_kelvin: f64,
    pub composition_terms: Vec<CompositionEnergyTerm>,
    pub cooling_rate_stabilization_threshold_kelvin_per_second: Option<f64>,
    pub cooling_rate_stabilization_j_per_mol: f64,
}

impl PhaseFreeEnergyModel {
    pub fn new(reference_gibbs_j_per_mol: f64, entropy_j_per_mol_kelvin: f64) -> Self {
        Self {
            low_temperature_kelvin: 0.0,
            high_temperature_kelvin: f64::MAX,
            reference_gibbs_j_per_mol,
            entropy_j_per_mol_kelvin,
            composition_terms: Vec::new(),
            cooling_rate_stabilization_threshold_kelvin_per_second: None,
            cooling_rate_stabilization_j_per_mol: 0.0,
        }
    }

    pub fn temperature_window(mut self, low: f64, high: f64) -> Self {
        self.low_temperature_kelvin = low;
        self.high_temperature_kelvin = high;
        self
    }

    pub fn composition_term(mut self, term: CompositionEnergyTerm) -> Self {
        self.composition_terms.push(term);
        self
    }

    pub fn cooling_rate_stabilization(
        mut self,
        threshold: f64,
        stabilization_j_per_mol: f64,
    ) -> Self {
        self.cooling_rate_stabilization_threshold_kelvin_per_second = Some(threshold);
        self.cooling_rate_stabilization_j_per_mol = stabilization_j_per_mol;
        self
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalPhasePropertyModel {
    pub hardness_hv: f64,
    pub yield_strength_mpa: f64,
    pub ductility_fraction: f64,
    pub electrical_resistivity_micro_ohm_meter: f64,
    pub thermal_conductivity_w_per_meter_kelvin: f64,
    pub corrosion_resistance_score: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PhaseKineticModel {
    pub diffusion_prefactor_square_meters_per_second: f64,
    pub diffusion_activation_energy_j_per_mol: f64,
    pub transformation_rate_per_second: f64,
    pub grain_growth_rate_micrometers_per_second: f64,
    pub recovery_rate_per_second: f64,
    pub precipitation_rate_per_second: f64,
}

impl PhaseKineticModel {
    pub fn new(
        diffusion_prefactor_square_meters_per_second: f64,
        diffusion_activation_energy_j_per_mol: f64,
        transformation_rate_per_second: f64,
        grain_growth_rate_micrometers_per_second: f64,
        recovery_rate_per_second: f64,
        precipitation_rate_per_second: f64,
    ) -> Self {
        Self {
            diffusion_prefactor_square_meters_per_second,
            diffusion_activation_energy_j_per_mol,
            transformation_rate_per_second,
            grain_growth_rate_micrometers_per_second,
            recovery_rate_per_second,
            precipitation_rate_per_second,
        }
    }

    pub fn for_phase_kind(kind: MetallurgicalPhaseKind) -> Self {
        match kind {
            MetallurgicalPhaseKind::Liquid => Self::new(1.0e-8, 35_000.0, 10.0, 0.0, 2.0, 0.0),
            MetallurgicalPhaseKind::Austenite => {
                Self::new(2.0e-10, 120_000.0, 0.03, 0.012, 2.0e-4, 1.0e-4)
            }
            MetallurgicalPhaseKind::Ferrite | MetallurgicalPhaseKind::SolidSolution => {
                Self::new(8.0e-11, 135_000.0, 0.02, 0.008, 1.5e-4, 8.0e-5)
            }
            MetallurgicalPhaseKind::Pearlite => {
                Self::new(4.0e-11, 150_000.0, 0.015, 0.006, 1.0e-4, 1.2e-4)
            }
            MetallurgicalPhaseKind::Bainite => {
                Self::new(2.0e-11, 155_000.0, 0.025, 0.004, 8.0e-5, 1.5e-4)
            }
            MetallurgicalPhaseKind::Martensite => {
                Self::new(1.0e-13, 180_000.0, 2.0, 0.001, 5.0e-5, 2.0e-5)
            }
            MetallurgicalPhaseKind::TemperedMartensite => {
                Self::new(8.0e-13, 165_000.0, 0.04, 0.003, 1.2e-4, 2.5e-4)
            }
            MetallurgicalPhaseKind::Cementite
            | MetallurgicalPhaseKind::Graphite
            | MetallurgicalPhaseKind::Intermetallic => {
                Self::new(5.0e-13, 180_000.0, 0.005, 0.002, 4.0e-5, 2.0e-4)
            }
        }
    }

    pub(super) fn diffusivity_at_kelvin(&self, temperature_kelvin: f64) -> ChemistryResult<f64> {
        validate_positive_finite(temperature_kelvin, "diffusion temperature")?;
        validate_kinetic_model(self)?;
        let exponent = -self.diffusion_activation_energy_j_per_mol
            / (GAS_CONSTANT_J_PER_MOL_KELVIN * temperature_kelvin);
        Ok(self.diffusion_prefactor_square_meters_per_second * exponent.exp())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalPhaseModel {
    pub id: String,
    pub kind: MetallurgicalPhaseKind,
    pub component_limits: Vec<ComponentLimit>,
    pub free_energy_model: PhaseFreeEnergyModel,
    pub property_model: MetallurgicalPhasePropertyModel,
    pub kinetic_model: PhaseKineticModel,
}

impl MetallurgicalPhaseModel {
    pub fn new(
        id: impl Into<String>,
        kind: MetallurgicalPhaseKind,
        property_model: MetallurgicalPhasePropertyModel,
    ) -> Self {
        Self {
            id: id.into(),
            kind,
            component_limits: Vec::new(),
            free_energy_model: PhaseFreeEnergyModel::new(0.0, 0.0),
            property_model,
            kinetic_model: PhaseKineticModel::for_phase_kind(kind),
        }
    }

    pub fn limit(mut self, limit: ComponentLimit) -> Self {
        self.component_limits.push(limit);
        self
    }

    pub fn free_energy_model(mut self, model: PhaseFreeEnergyModel) -> Self {
        self.free_energy_model = model;
        self
    }

    pub fn kinetic_model(mut self, model: PhaseKineticModel) -> Self {
        self.kinetic_model = model;
        self
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PhaseBoundaryPoint {
    pub composition: BTreeMap<MetallurgicalComponentId, f64>,
    pub solidus_kelvin: f64,
    pub liquidus_kelvin: f64,
}

impl PhaseBoundaryPoint {
    pub fn new(
        composition: impl IntoIterator<Item = (impl Into<MetallurgicalComponentId>, f64)>,
        solidus_kelvin: f64,
        liquidus_kelvin: f64,
    ) -> Self {
        Self {
            composition: composition
                .into_iter()
                .map(|(component, fraction)| (component.into(), fraction))
                .collect(),
            solidus_kelvin,
            liquidus_kelvin,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct PhaseBoundarySnapshot {
    pub solidus_kelvin: f64,
    pub liquidus_kelvin: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalSystem {
    pub id: String,
    pub components: BTreeSet<MetallurgicalComponentId>,
    pub phase_models: Vec<MetallurgicalPhaseModel>,
    pub phase_boundaries: Vec<PhaseBoundaryPoint>,
}

impl MetallurgicalSystem {
    pub fn new(
        id: impl Into<String>,
        components: impl IntoIterator<Item = impl Into<MetallurgicalComponentId>>,
    ) -> Self {
        Self {
            id: id.into(),
            components: components.into_iter().map(Into::into).collect(),
            phase_models: Vec::new(),
            phase_boundaries: Vec::new(),
        }
    }

    pub fn phase_model(mut self, model: MetallurgicalPhaseModel) -> Self {
        self.phase_models.push(model);
        self
    }

    pub fn phase_boundary(mut self, point: PhaseBoundaryPoint) -> Self {
        self.phase_boundaries.push(point);
        self
    }

    pub fn validate(&self) -> ChemistryResult<()> {
        if self.id.trim().is_empty() {
            return Err(ChemistryError::InvalidMixtureState(
                "metallurgical system id must not be empty".to_string(),
            ));
        }
        if self.components.is_empty() {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "metallurgical system '{}' has no components",
                self.id
            )));
        }
        if self.phase_models.is_empty() {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "metallurgical system '{}' has no phase models",
                self.id
            )));
        }
        if self.phase_boundaries.is_empty() {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "metallurgical system '{}' has no phase-boundary data",
                self.id
            )));
        }
        for phase in &self.phase_models {
            validate_phase_model(self, phase)?;
        }
        for point in &self.phase_boundaries {
            validate_phase_boundary_point(self, point)?;
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalPhaseAmount {
    pub phase_id: String,
    pub kind: MetallurgicalPhaseKind,
    pub fraction: f64,
    pub property_model: MetallurgicalPhasePropertyModel,
    pub kinetic_model: PhaseKineticModel,
}

#[derive(Debug, Clone, PartialEq)]
pub struct GrainStructure {
    pub average_grain_size_micrometers: f64,
    pub distribution_width: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DefectState {
    pub vacancy_fraction: f64,
    pub dislocation_density_per_square_meter: f64,
    pub cold_work_fraction: f64,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum MechanicalWorkingMode {
    Forging,
    Rolling,
    Drawing,
    Extrusion,
    Machining,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MechanicalWorkingProcess {
    pub mode: MechanicalWorkingMode,
    pub true_strain: f64,
    pub strain_rate_per_second: f64,
    pub temperature_kelvin: f64,
    pub duration_seconds: f64,
}

impl MechanicalWorkingProcess {
    pub fn new(
        mode: MechanicalWorkingMode,
        true_strain: f64,
        strain_rate_per_second: f64,
        temperature_kelvin: f64,
        duration_seconds: f64,
    ) -> Self {
        Self {
            mode,
            true_strain,
            strain_rate_per_second,
            temperature_kelvin,
            duration_seconds,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MechanicalHistoryState {
    pub accumulated_true_strain: f64,
    pub recent_true_strain: f64,
    pub strain_rate_per_second: f64,
    pub recrystallized_fraction: f64,
    pub deformation_temperature_kelvin: f64,
    pub elapsed_work_seconds: f64,
}

impl MechanicalHistoryState {
    pub(super) fn initial(temperature_kelvin: f64) -> Self {
        Self {
            accumulated_true_strain: 0.0,
            recent_true_strain: 0.0,
            strain_rate_per_second: 0.0,
            recrystallized_fraction: 0.0,
            deformation_temperature_kelvin: temperature_kelvin,
            elapsed_work_seconds: 0.0,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ThermalTreatmentState {
    pub previous_temperature_kelvin: f64,
    pub peak_temperature_kelvin: f64,
    pub cooling_rate_kelvin_per_second: f64,
    pub hold_time_seconds: f64,
    pub elapsed_time_seconds: f64,
}

impl ThermalTreatmentState {
    pub fn initial(temperature_kelvin: f64) -> ChemistryResult<Self> {
        validate_non_negative_finite(temperature_kelvin, "initial metallurgical temperature")?;
        Ok(Self {
            previous_temperature_kelvin: temperature_kelvin,
            peak_temperature_kelvin: temperature_kelvin,
            cooling_rate_kelvin_per_second: 0.0,
            hold_time_seconds: 0.0,
            elapsed_time_seconds: 0.0,
        })
    }

    pub fn advance(&self, temperature_kelvin: f64, delta_seconds: f64) -> ChemistryResult<Self> {
        validate_non_negative_finite(temperature_kelvin, "metallurgical temperature")?;
        validate_positive_finite(delta_seconds, "metallurgical tick duration")?;
        let cooling_rate =
            ((self.previous_temperature_kelvin - temperature_kelvin) / delta_seconds).max(0.0);
        let hold_time_seconds =
            if (self.previous_temperature_kelvin - temperature_kelvin).abs() < 1.0e-6 {
                self.hold_time_seconds + delta_seconds
            } else {
                0.0
            };
        Ok(Self {
            previous_temperature_kelvin: temperature_kelvin,
            peak_temperature_kelvin: self.peak_temperature_kelvin.max(temperature_kelvin),
            cooling_rate_kelvin_per_second: cooling_rate,
            hold_time_seconds,
            elapsed_time_seconds: self.elapsed_time_seconds + delta_seconds,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DiffusionState {
    pub effective_diffusivity_square_meters_per_second: f64,
    pub diffusion_length_micrometers: f64,
    pub homogenization_fraction: f64,
    pub precipitation_fraction: f64,
    pub aging_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlloyPropertySnapshot {
    pub hardness_hv: f64,
    pub yield_strength_mpa: f64,
    pub ductility_fraction: f64,
    pub electrical_resistivity_micro_ohm_meter: f64,
    pub thermal_conductivity_w_per_meter_kelvin: f64,
    pub corrosion_resistance_score: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlloyServicePropertySnapshot {
    pub fracture_toughness_mpa_sqrt_meter: f64,
    pub brittleness_score: f64,
    pub wear_resistance_score: f64,
    pub electrical_conductivity_percent_iacs: f64,
    pub high_temperature_stability_score: f64,
    pub softening_temperature_kelvin: f64,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MetallurgicalUseKind {
    Structural,
    CuttingTool,
    Spring,
    ElectricalConductor,
    ThermalConductor,
    CorrosionResistant,
    HighTemperature,
    WearResistant,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalUseSuitability {
    pub kind: MetallurgicalUseKind,
    pub score: f64,
    pub limiting_factor: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlloyUseProfile {
    pub suitability: Vec<MetallurgicalUseSuitability>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalSystemSelectionDiagnostic {
    pub system_id: String,
    pub covers_composition: bool,
    pub missing_components: Vec<MetallurgicalComponentId>,
    pub composition_distance: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalPhaseDiagnostic {
    pub phase_id: String,
    pub kind: MetallurgicalPhaseKind,
    pub selected: bool,
    pub fraction: f64,
    pub gibbs_j_per_mol: Option<f64>,
    pub energy_above_minimum_j_per_mol: Option<f64>,
    pub reason: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalThermalDiagnostic {
    pub previous_temperature_kelvin: f64,
    pub current_temperature_kelvin: f64,
    pub cooling_rate_kelvin_per_second: f64,
    pub hold_time_seconds: f64,
    pub delta_seconds: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalDiagnosticReport {
    pub selected_system_id: Option<String>,
    pub considered_systems: Vec<MetallurgicalSystemSelectionDiagnostic>,
    pub phase_boundaries: Option<PhaseBoundarySnapshot>,
    pub phase_reasons: Vec<MetallurgicalPhaseDiagnostic>,
    pub thermal_reason: MetallurgicalThermalDiagnostic,
    pub unmodeled_reason: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum MetallurgicalStateKind {
    Modeled { system_id: String },
    Unmodeled { reason: String },
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalState {
    pub kind: MetallurgicalStateKind,
    pub composition: MetallurgicalComposition,
    pub temperature_kelvin: f64,
    pub phase_boundaries: Option<PhaseBoundarySnapshot>,
    pub phases: Vec<MetallurgicalPhaseAmount>,
    pub grain_structure: GrainStructure,
    pub defect_state: DefectState,
    pub mechanical_history: MechanicalHistoryState,
    pub diffusion_state: DiffusionState,
    pub thermal_treatment: ThermalTreatmentState,
    pub properties: AlloyPropertySnapshot,
    pub service_properties: AlloyServicePropertySnapshot,
    pub use_profile: AlloyUseProfile,
    pub diagnostics: MetallurgicalDiagnosticReport,
}
