use std::collections::{BTreeMap, BTreeSet};

use super::alloy::AlloyPhaseSnapshot;
use super::error::{ChemistryError, ChemistryResult};

const DEFAULT_GRAIN_SIZE_MICROMETERS: f64 = 50.0;
const TRACE_COMPONENT_FRACTION: f64 = 1.0e-9;

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
pub struct PhaseWeightRule {
    pub component: Option<MetallurgicalComponentId>,
    pub center_fraction: f64,
    pub width_fraction: f64,
    pub low_temperature_kelvin: f64,
    pub high_temperature_kelvin: f64,
    pub base_weight: f64,
    pub cooling_rate_bonus_threshold_kelvin_per_second: Option<f64>,
    pub cooling_rate_bonus: f64,
}

impl PhaseWeightRule {
    pub fn constant(base_weight: f64) -> Self {
        Self {
            component: None,
            center_fraction: 0.0,
            width_fraction: 1.0,
            low_temperature_kelvin: 0.0,
            high_temperature_kelvin: f64::MAX,
            base_weight,
            cooling_rate_bonus_threshold_kelvin_per_second: None,
            cooling_rate_bonus: 0.0,
        }
    }

    pub fn composition_centered(
        component: impl Into<MetallurgicalComponentId>,
        center_fraction: f64,
        width_fraction: f64,
        base_weight: f64,
    ) -> Self {
        Self {
            component: Some(component.into()),
            center_fraction,
            width_fraction,
            low_temperature_kelvin: 0.0,
            high_temperature_kelvin: f64::MAX,
            base_weight,
            cooling_rate_bonus_threshold_kelvin_per_second: None,
            cooling_rate_bonus: 0.0,
        }
    }

    pub fn temperature_window(mut self, low: f64, high: f64) -> Self {
        self.low_temperature_kelvin = low;
        self.high_temperature_kelvin = high;
        self
    }

    pub fn cooling_rate_bonus(mut self, threshold: f64, bonus: f64) -> Self {
        self.cooling_rate_bonus_threshold_kelvin_per_second = Some(threshold);
        self.cooling_rate_bonus = bonus;
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
pub struct MetallurgicalPhaseModel {
    pub id: String,
    pub kind: MetallurgicalPhaseKind,
    pub component_limits: Vec<ComponentLimit>,
    pub weight_rule: PhaseWeightRule,
    pub property_model: MetallurgicalPhasePropertyModel,
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
            weight_rule: PhaseWeightRule::constant(1.0),
            property_model,
        }
    }

    pub fn limit(mut self, limit: ComponentLimit) -> Self {
        self.component_limits.push(limit);
        self
    }

    pub fn weight_rule(mut self, rule: PhaseWeightRule) -> Self {
        self.weight_rule = rule;
        self
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalSystem {
    pub id: String,
    pub components: BTreeSet<MetallurgicalComponentId>,
    pub phase_models: Vec<MetallurgicalPhaseModel>,
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
        }
    }

    pub fn phase_model(mut self, model: MetallurgicalPhaseModel) -> Self {
        self.phase_models.push(model);
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
        for phase in &self.phase_models {
            validate_phase_model(self, phase)?;
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

#[derive(Debug, Clone, PartialEq)]
pub struct ThermalTreatmentState {
    pub previous_temperature_kelvin: f64,
    pub peak_temperature_kelvin: f64,
    pub cooling_rate_kelvin_per_second: f64,
    pub hold_time_seconds: f64,
}

impl ThermalTreatmentState {
    pub fn initial(temperature_kelvin: f64) -> ChemistryResult<Self> {
        validate_non_negative_finite(temperature_kelvin, "initial metallurgical temperature")?;
        Ok(Self {
            previous_temperature_kelvin: temperature_kelvin,
            peak_temperature_kelvin: temperature_kelvin,
            cooling_rate_kelvin_per_second: 0.0,
            hold_time_seconds: 0.0,
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
        })
    }
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
pub enum MetallurgicalStateKind {
    Modeled { system_id: String },
    Unmodeled { reason: String },
}

#[derive(Debug, Clone, PartialEq)]
pub struct MetallurgicalState {
    pub kind: MetallurgicalStateKind,
    pub composition: MetallurgicalComposition,
    pub temperature_kelvin: f64,
    pub phases: Vec<MetallurgicalPhaseAmount>,
    pub grain_structure: GrainStructure,
    pub defect_state: DefectState,
    pub thermal_treatment: ThermalTreatmentState,
    pub properties: AlloyPropertySnapshot,
}

pub fn metallurgical_state_from_alloy_phase(
    alloy: &AlloyPhaseSnapshot,
    systems: &[MetallurgicalSystem],
    previous: Option<&MetallurgicalState>,
    delta_seconds: f64,
) -> ChemistryResult<MetallurgicalState> {
    validate_positive_finite(delta_seconds, "metallurgical tick duration")?;
    let composition = MetallurgicalComposition::from_alloy_phase(alloy)?;
    let thermal_treatment = match previous {
        Some(previous) => previous
            .thermal_treatment
            .advance(alloy.temperature_kelvin, delta_seconds)?,
        None => ThermalTreatmentState::initial(alloy.temperature_kelvin)?,
    };
    for system in systems {
        system.validate()?;
        if system_matches_composition(system, &composition) {
            return modeled_state(alloy, composition, system, thermal_treatment, previous);
        }
    }
    Ok(unmodeled_state(
        alloy,
        composition,
        thermal_treatment,
        "no registered metallurgical system covers all components",
    ))
}

fn modeled_state(
    alloy: &AlloyPhaseSnapshot,
    composition: MetallurgicalComposition,
    system: &MetallurgicalSystem,
    thermal_treatment: ThermalTreatmentState,
    previous: Option<&MetallurgicalState>,
) -> ChemistryResult<MetallurgicalState> {
    let mut weighted_phases = Vec::new();
    for model in &system.phase_models {
        if !phase_limits_match(model, &composition) {
            continue;
        }
        let weight = phase_weight(
            model,
            &composition,
            alloy.temperature_kelvin,
            &thermal_treatment,
        )?;
        if weight > TRACE_COMPONENT_FRACTION {
            weighted_phases.push((model, weight));
        }
    }
    if weighted_phases.is_empty() {
        return Ok(unmodeled_state(
            alloy,
            composition,
            thermal_treatment,
            "registered system has no stable phase for this composition and temperature",
        ));
    }
    let total_weight = weighted_phases
        .iter()
        .map(|(_, weight)| *weight)
        .sum::<f64>();
    validate_positive_finite(total_weight, "metallurgical phase total weight")?;
    let phases = weighted_phases
        .into_iter()
        .map(|(model, weight)| MetallurgicalPhaseAmount {
            phase_id: model.id.clone(),
            kind: model.kind,
            fraction: weight / total_weight,
            property_model: model.property_model.clone(),
        })
        .collect::<Vec<_>>();
    let grain_structure = estimate_grain_structure(&thermal_treatment, previous);
    let defect_state = estimate_defect_state(&thermal_treatment, previous);
    let properties = estimate_properties(&phases, &grain_structure, &defect_state)?;
    Ok(MetallurgicalState {
        kind: MetallurgicalStateKind::Modeled {
            system_id: system.id.clone(),
        },
        composition,
        temperature_kelvin: alloy.temperature_kelvin,
        phases,
        grain_structure,
        defect_state,
        thermal_treatment,
        properties,
    })
}

fn unmodeled_state(
    alloy: &AlloyPhaseSnapshot,
    composition: MetallurgicalComposition,
    thermal_treatment: ThermalTreatmentState,
    reason: impl Into<String>,
) -> MetallurgicalState {
    MetallurgicalState {
        kind: MetallurgicalStateKind::Unmodeled {
            reason: reason.into(),
        },
        composition,
        temperature_kelvin: alloy.temperature_kelvin,
        phases: Vec::new(),
        grain_structure: GrainStructure {
            average_grain_size_micrometers: DEFAULT_GRAIN_SIZE_MICROMETERS,
            distribution_width: 0.3,
        },
        defect_state: DefectState {
            vacancy_fraction: 0.0,
            dislocation_density_per_square_meter: 1.0e10,
            cold_work_fraction: 0.0,
        },
        thermal_treatment,
        properties: AlloyPropertySnapshot {
            hardness_hv: 0.0,
            yield_strength_mpa: 0.0,
            ductility_fraction: 0.0,
            electrical_resistivity_micro_ohm_meter: 0.0,
            thermal_conductivity_w_per_meter_kelvin: 0.0,
            corrosion_resistance_score: 0.0,
        },
    }
}

fn system_matches_composition(
    system: &MetallurgicalSystem,
    composition: &MetallurgicalComposition,
) -> bool {
    composition
        .components
        .keys()
        .all(|component| system.components.contains(component))
}

fn phase_limits_match(
    phase: &MetallurgicalPhaseModel,
    composition: &MetallurgicalComposition,
) -> bool {
    phase.component_limits.iter().all(|limit| {
        let fraction = composition.fraction_of(&limit.component);
        fraction >= limit.min_fraction && fraction <= limit.max_fraction
    })
}

fn phase_weight(
    phase: &MetallurgicalPhaseModel,
    composition: &MetallurgicalComposition,
    temperature_kelvin: f64,
    thermal_treatment: &ThermalTreatmentState,
) -> ChemistryResult<f64> {
    validate_non_negative_finite(temperature_kelvin, "metallurgical temperature")?;
    let rule = &phase.weight_rule;
    validate_phase_weight_rule(rule)?;
    if temperature_kelvin < rule.low_temperature_kelvin
        || temperature_kelvin > rule.high_temperature_kelvin
    {
        return Ok(0.0);
    }
    let temperature_span = (rule.high_temperature_kelvin - rule.low_temperature_kelvin).max(1.0);
    let middle = (rule.high_temperature_kelvin + rule.low_temperature_kelvin) * 0.5;
    let temperature_factor =
        (1.0 - ((temperature_kelvin - middle).abs() / (temperature_span * 0.5))).max(0.0);
    let composition_factor = match &rule.component {
        Some(component) => {
            let fraction = composition.fraction_of(component);
            (1.0 - ((fraction - rule.center_fraction).abs() / rule.width_fraction)).max(0.0)
        }
        None => 1.0,
    };
    let mut weight = rule.base_weight * temperature_factor.max(0.05) * composition_factor;
    if let Some(threshold) = rule.cooling_rate_bonus_threshold_kelvin_per_second {
        if thermal_treatment.cooling_rate_kelvin_per_second >= threshold {
            weight += rule.cooling_rate_bonus;
        }
    }
    Ok(weight.max(0.0))
}

fn estimate_grain_structure(
    thermal_treatment: &ThermalTreatmentState,
    previous: Option<&MetallurgicalState>,
) -> GrainStructure {
    let previous_size = previous
        .map(|state| state.grain_structure.average_grain_size_micrometers)
        .unwrap_or(DEFAULT_GRAIN_SIZE_MICROMETERS);
    let cooling_refinement = 1.0 / (1.0 + thermal_treatment.cooling_rate_kelvin_per_second / 50.0);
    let hold_growth = 1.0 + thermal_treatment.hold_time_seconds.min(3600.0) / 3600.0;
    let temperature_growth = if thermal_treatment.previous_temperature_kelvin > 1000.0 {
        thermal_treatment.previous_temperature_kelvin / 1000.0
    } else {
        1.0
    };
    GrainStructure {
        average_grain_size_micrometers: (previous_size
            * hold_growth
            * temperature_growth
            * cooling_refinement.max(0.15))
        .clamp(1.0, 5000.0),
        distribution_width: (0.2 + thermal_treatment.cooling_rate_kelvin_per_second / 500.0)
            .clamp(0.1, 1.0),
    }
}

fn estimate_defect_state(
    thermal_treatment: &ThermalTreatmentState,
    previous: Option<&MetallurgicalState>,
) -> DefectState {
    let previous_dislocation = previous
        .map(|state| state.defect_state.dislocation_density_per_square_meter)
        .unwrap_or(1.0e10);
    let quench_factor = (thermal_treatment.cooling_rate_kelvin_per_second / 200.0).clamp(0.0, 1.0);
    let recovery_factor = if thermal_treatment.previous_temperature_kelvin > 0.0 {
        (1.0 - thermal_treatment.hold_time_seconds / 7200.0).clamp(0.2, 1.0)
    } else {
        1.0
    };
    DefectState {
        vacancy_fraction: (quench_factor * 1.0e-4).clamp(0.0, 1.0e-3),
        dislocation_density_per_square_meter: (previous_dislocation * recovery_factor
            + quench_factor * 5.0e13)
            .clamp(1.0e8, 1.0e16),
        cold_work_fraction: previous
            .map(|state| state.defect_state.cold_work_fraction)
            .unwrap_or(0.0),
    }
}

fn estimate_properties(
    phases: &[MetallurgicalPhaseAmount],
    grain_structure: &GrainStructure,
    defect_state: &DefectState,
) -> ChemistryResult<AlloyPropertySnapshot> {
    let total = phases.iter().map(|phase| phase.fraction).sum::<f64>();
    if !total.is_finite() || (total - 1.0).abs() > 1.0e-6 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "metallurgical phase fractions must sum to 1.0, got {total}"
        )));
    }
    let mut properties = AlloyPropertySnapshot {
        hardness_hv: 0.0,
        yield_strength_mpa: 0.0,
        ductility_fraction: 0.0,
        electrical_resistivity_micro_ohm_meter: 0.0,
        thermal_conductivity_w_per_meter_kelvin: 0.0,
        corrosion_resistance_score: 0.0,
    };
    for phase in phases {
        let f = phase.fraction;
        let model = &phase.property_model;
        properties.hardness_hv += f * model.hardness_hv;
        properties.yield_strength_mpa += f * model.yield_strength_mpa;
        properties.ductility_fraction += f * model.ductility_fraction;
        properties.electrical_resistivity_micro_ohm_meter +=
            f * model.electrical_resistivity_micro_ohm_meter;
        properties.thermal_conductivity_w_per_meter_kelvin +=
            f * model.thermal_conductivity_w_per_meter_kelvin;
        properties.corrosion_resistance_score += f * model.corrosion_resistance_score;
    }
    let grain_strengthening = 80.0 / grain_structure.average_grain_size_micrometers.sqrt();
    let dislocation_strengthening =
        (defect_state.dislocation_density_per_square_meter / 1.0e12).sqrt() * 15.0;
    properties.yield_strength_mpa += grain_strengthening + dislocation_strengthening;
    properties.hardness_hv += (grain_strengthening + dislocation_strengthening) * 0.3;
    properties.ductility_fraction = (properties.ductility_fraction
        * (1.0 - defect_state.vacancy_fraction * 1000.0))
        .clamp(0.0, 1.0);
    Ok(properties)
}

fn validate_phase_model(
    system: &MetallurgicalSystem,
    phase: &MetallurgicalPhaseModel,
) -> ChemistryResult<()> {
    if phase.id.trim().is_empty() {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "metallurgical system '{}' has a phase with empty id",
            system.id
        )));
    }
    for limit in &phase.component_limits {
        if !system.components.contains(&limit.component) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "phase '{}' references component '{}' outside metallurgical system '{}'",
                phase.id,
                limit.component.as_str(),
                system.id
            )));
        }
        validate_fraction(limit.min_fraction, "phase component minimum")?;
        validate_fraction(limit.max_fraction, "phase component maximum")?;
        if limit.min_fraction > limit.max_fraction {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "phase '{}' has component limit minimum above maximum",
                phase.id
            )));
        }
    }
    validate_phase_weight_rule(&phase.weight_rule)?;
    validate_property_model(&phase.property_model)?;
    Ok(())
}

fn validate_phase_weight_rule(rule: &PhaseWeightRule) -> ChemistryResult<()> {
    validate_non_negative_finite(rule.low_temperature_kelvin, "phase low temperature")?;
    validate_positive_finite(rule.high_temperature_kelvin, "phase high temperature")?;
    if rule.low_temperature_kelvin >= rule.high_temperature_kelvin {
        return Err(ChemistryError::InvalidMixtureState(
            "phase temperature window low bound must be below high bound".to_string(),
        ));
    }
    validate_non_negative_finite(rule.base_weight, "phase base weight")?;
    validate_fraction(rule.center_fraction, "phase composition center")?;
    validate_positive_finite(rule.width_fraction, "phase composition width")?;
    if let Some(threshold) = rule.cooling_rate_bonus_threshold_kelvin_per_second {
        validate_non_negative_finite(threshold, "phase cooling-rate threshold")?;
        validate_non_negative_finite(rule.cooling_rate_bonus, "phase cooling-rate bonus")?;
    }
    Ok(())
}

fn validate_property_model(model: &MetallurgicalPhasePropertyModel) -> ChemistryResult<()> {
    validate_non_negative_finite(model.hardness_hv, "phase hardness")?;
    validate_non_negative_finite(model.yield_strength_mpa, "phase yield strength")?;
    validate_fraction(model.ductility_fraction, "phase ductility")?;
    validate_non_negative_finite(
        model.electrical_resistivity_micro_ohm_meter,
        "phase electrical resistivity",
    )?;
    validate_non_negative_finite(
        model.thermal_conductivity_w_per_meter_kelvin,
        "phase thermal conductivity",
    )?;
    validate_fraction(
        model.corrosion_resistance_score,
        "phase corrosion resistance",
    )?;
    Ok(())
}

fn validate_fraction(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be finite and within 0.0..=1.0"
        )));
    }
    Ok(())
}

fn validate_positive_finite(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be positive and finite"
        )));
    }
    Ok(())
}

fn validate_non_negative_finite(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be non-negative and finite"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::alloy::alloy_phase_snapshots;
    use crate::chemistry::metallurgy_data::default_metallurgical_systems;
    use crate::chemistry::mixture::Mixture;
    use crate::chemistry::registry::ChemistryRegistryBuilder;
    use crate::chemistry::substance::{
        LiquidPhasePreference, SolventRole, Substance, SubstancePhaseProperties,
        SubstanceRepresentation,
    };

    #[test]
    fn default_metallurgical_systems_are_valid() {
        for system in default_metallurgical_systems() {
            system.validate().unwrap();
        }
    }

    #[test]
    fn iron_carbon_melt_gets_modeled_liquid_state() {
        let registry = test_registry().build().unwrap();
        let mut mixture = Mixture::new(1900.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:test_iron", 0.98)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:carbon", 0.02)
            .unwrap();

        let alloy = alloy_phase_snapshots(&registry, &mixture)
            .unwrap()
            .remove(0);
        let state = metallurgical_state_from_alloy_phase(
            &alloy,
            &default_metallurgical_systems(),
            None,
            1.0,
        )
        .unwrap();

        assert!(matches!(
            state.kind,
            MetallurgicalStateKind::Modeled { ref system_id } if system_id == "metallurgy:fe_c"
        ));
        assert!(state
            .phases
            .iter()
            .any(|phase| phase.kind == MetallurgicalPhaseKind::Liquid && phase.fraction > 0.5));
    }

    #[test]
    fn fast_cooling_promotes_martensite_and_harder_state() {
        let registry = test_registry().build().unwrap();
        let systems = default_metallurgical_systems();
        let mut hot = Mixture::new(1900.0).unwrap();
        hot.add_substance(&registry, "destroy:test_iron", 0.97)
            .unwrap();
        hot.add_substance(&registry, "destroy:carbon", 0.03)
            .unwrap();
        let hot_alloy = alloy_phase_snapshots(&registry, &hot).unwrap().remove(0);
        let hot_state =
            metallurgical_state_from_alloy_phase(&hot_alloy, &systems, None, 1.0).unwrap();

        let mut slow_alloy = hot_alloy.clone();
        slow_alloy.temperature_kelvin = 500.0;
        let slow_state =
            metallurgical_state_from_alloy_phase(&slow_alloy, &systems, Some(&hot_state), 50.0)
                .unwrap();

        let fast_state =
            metallurgical_state_from_alloy_phase(&slow_alloy, &systems, Some(&hot_state), 1.0)
                .unwrap();

        let fast_martensite = phase_fraction(&fast_state, MetallurgicalPhaseKind::Martensite);
        let slow_martensite = phase_fraction(&slow_state, MetallurgicalPhaseKind::Martensite);
        assert!(
            fast_martensite > slow_martensite,
            "fast martensite {fast_martensite}, slow martensite {slow_martensite}"
        );
        assert!(
            fast_state.properties.hardness_hv > slow_state.properties.hardness_hv,
            "fast hardness {}, slow hardness {}",
            fast_state.properties.hardness_hv,
            slow_state.properties.hardness_hv
        );
    }

    #[test]
    fn unknown_metal_system_is_explicitly_unmodeled() {
        let registry = test_registry().build().unwrap();
        let mut mixture = Mixture::new(2000.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:test_lead", 1.0)
            .unwrap();

        let alloy = alloy_phase_snapshots(&registry, &mixture)
            .unwrap()
            .remove(0);
        let state = metallurgical_state_from_alloy_phase(
            &alloy,
            &default_metallurgical_systems(),
            None,
            1.0,
        )
        .unwrap();

        assert!(matches!(
            state.kind,
            MetallurgicalStateKind::Unmodeled { .. }
        ));
        assert!(state.phases.is_empty());
    }

    fn phase_fraction(state: &MetallurgicalState, kind: MetallurgicalPhaseKind) -> f64 {
        state
            .phases
            .iter()
            .find(|phase| phase.kind == kind)
            .map(|phase| phase.fraction)
            .unwrap_or(0.0)
    }

    fn test_registry() -> ChemistryRegistryBuilder {
        ChemistryRegistryBuilder::new()
            .substance(test_metal("destroy:test_iron", "Fe", 55.845, 1811.0))
            .substance(test_metal("destroy:test_lead", "Pb", 207.2, 600.61))
            .substance(
                Substance::new("destroy:carbon", 0, 12.011, 2_200.0, 4300.0, 8.5, 0.0)
                    .with_melting_point_kelvin(1000.0)
                    .with_phase_properties(molten_metal_phase_properties())
                    .with_representation(SubstanceRepresentation::UnspecifiedMaterial {
                        reason: "test metallurgical carbon component".to_string(),
                    }),
            )
    }

    fn test_metal(
        id: &'static str,
        element: &'static str,
        molar_mass: f64,
        melting_point: f64,
    ) -> Substance {
        Substance::new(id, 0, molar_mass, 7_800.0, 3300.0, 25.0, 0.0)
            .with_solid_density_grams_per_bucket(7_800.0)
            .with_melting_point_kelvin(melting_point)
            .with_phase_properties(molten_metal_phase_properties())
            .with_representation(SubstanceRepresentation::Metal {
                element_symbol: element.to_string(),
            })
    }

    fn molten_metal_phase_properties() -> SubstancePhaseProperties {
        SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::MoltenMetal,
            aqueous_solubility_mol_per_bucket: Some(0.0),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: true,
            solvent_role: SolventRole::NotSolvent,
        }
    }
}
