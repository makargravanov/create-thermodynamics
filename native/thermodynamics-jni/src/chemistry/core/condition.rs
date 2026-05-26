use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, MixturePhase, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::registry::ChemistryRegistry;
use super::substance::SubstanceId;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum AcidityCondition {
    Acidic,
    Neutral,
    Basic,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum RedoxCondition {
    Oxidizing,
    Reducing,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum AtmosphereCondition {
    Air,
    Inert,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReactionCondition {
    pub phases: Vec<MixturePhase>,
    pub acidity: Option<AcidityCondition>,
    pub solvent: Option<SubstanceId>,
    pub min_temperature_kelvin: Option<f64>,
    pub max_temperature_kelvin: Option<f64>,
    pub min_water_activity: Option<f64>,
    pub max_water_activity: Option<f64>,
    pub max_oxygen_activity: Option<f64>,
    pub gas_pressure_atm: Option<f64>,
    pub atmosphere: Option<AtmosphereCondition>,
    pub redox: Option<RedoxCondition>,
    pub redox_strength: f64,
    pub rate_multiplier: f64,
    pub reason: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReactionConditionEvaluation {
    pub allowed: bool,
    pub multiplier: f64,
    pub blocked_reasons: Vec<String>,
}

impl ReactionCondition {
    pub fn new(reason: impl Into<String>) -> Self {
        Self {
            phases: Vec::new(),
            acidity: None,
            solvent: None,
            min_temperature_kelvin: None,
            max_temperature_kelvin: None,
            min_water_activity: None,
            max_water_activity: None,
            max_oxygen_activity: None,
            gas_pressure_atm: None,
            atmosphere: None,
            redox: None,
            redox_strength: 0.0,
            rate_multiplier: 1.0,
            reason: reason.into(),
        }
    }

    pub fn in_phase(mut self, phase: MixturePhase) -> Self {
        if !self.phases.contains(&phase) {
            self.phases.push(phase);
        }
        self
    }

    pub fn acidity(mut self, acidity: AcidityCondition) -> Self {
        self.acidity = Some(acidity);
        self
    }

    pub fn solvent(mut self, solvent: impl Into<SubstanceId>) -> Self {
        self.solvent = Some(solvent.into());
        self
    }

    pub fn min_temperature_kelvin(mut self, value: f64) -> Self {
        self.min_temperature_kelvin = Some(value);
        self
    }

    pub fn max_temperature_kelvin(mut self, value: f64) -> Self {
        self.max_temperature_kelvin = Some(value);
        self
    }

    pub fn min_water_activity(mut self, value: f64) -> Self {
        self.min_water_activity = Some(value);
        self
    }

    pub fn max_water_activity(mut self, value: f64) -> Self {
        self.max_water_activity = Some(value);
        self
    }

    pub fn max_oxygen_activity(mut self, value: f64) -> Self {
        self.max_oxygen_activity = Some(value);
        self
    }

    pub fn gas_pressure_atm(mut self, value: f64) -> Self {
        self.gas_pressure_atm = Some(value);
        self
    }

    pub fn atmosphere(mut self, atmosphere: AtmosphereCondition) -> Self {
        self.atmosphere = Some(atmosphere);
        self
    }

    pub fn redox(mut self, condition: RedoxCondition, strength: f64) -> Self {
        self.redox = Some(condition);
        self.redox_strength = strength;
        self
    }

    pub fn rate_multiplier(mut self, value: f64) -> Self {
        self.rate_multiplier = value;
        self
    }

    pub fn validate(&self, reaction_id: &str) -> ChemistryResult<()> {
        if self.reason.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "reaction condition must explain what it represents".to_string(),
            });
        }
        validate_optional_positive(
            self.min_temperature_kelvin,
            reaction_id,
            "minimum temperature",
        )?;
        validate_optional_positive(
            self.max_temperature_kelvin,
            reaction_id,
            "maximum temperature",
        )?;
        if let (Some(min), Some(max)) = (self.min_temperature_kelvin, self.max_temperature_kelvin) {
            if min > max {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.to_string(),
                    reason: "minimum temperature cannot exceed maximum temperature".to_string(),
                });
            }
        }
        validate_optional_unit(
            self.min_water_activity,
            reaction_id,
            "minimum water activity",
        )?;
        validate_optional_unit(
            self.max_water_activity,
            reaction_id,
            "maximum water activity",
        )?;
        validate_optional_unit(
            self.max_oxygen_activity,
            reaction_id,
            "maximum oxygen activity",
        )?;
        validate_optional_positive(self.gas_pressure_atm, reaction_id, "gas pressure")?;
        if let (Some(min), Some(max)) = (self.min_water_activity, self.max_water_activity) {
            if min > max {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.to_string(),
                    reason: "minimum water activity cannot exceed maximum water activity"
                        .to_string(),
                });
            }
        }
        if self.redox.is_some() && (!self.redox_strength.is_finite() || self.redox_strength <= 0.0)
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "redox condition strength must be positive and finite".to_string(),
            });
        }
        if !self.rate_multiplier.is_finite() || self.rate_multiplier < 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "condition rate multiplier must be non-negative and finite".to_string(),
            });
        }
        Ok(())
    }
}

pub fn evaluate_reaction_conditions(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    conditions: &[ReactionCondition],
) -> ChemistryResult<ReactionConditionEvaluation> {
    let mut multiplier = 1.0;
    let mut blocked_reasons = Vec::new();
    for condition in conditions {
        let mut blocked = false;
        if !condition.phases.is_empty()
            && !condition
                .phases
                .iter()
                .any(|phase| mixture.total_in_phase(*phase) > TRACE_CONCENTRATION_MOL_PER_BUCKET)
        {
            blocked = true;
        }
        if let Some(acidity) = condition.acidity {
            let ph = mixture.ph(registry)?;
            let acidity_matches = match (acidity, ph) {
                (AcidityCondition::Acidic, Some(value)) => value < 6.0,
                (AcidityCondition::Neutral, Some(value)) => (6.0..=8.0).contains(&value),
                (AcidityCondition::Basic, Some(value)) => value > 8.0,
                (_, None) => false,
            };
            if !acidity_matches {
                blocked = true;
            }
        }
        if let Some(solvent) = &condition.solvent {
            registry.substance(solvent)?;
            if mixture.concentration_of(solvent) <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                blocked = true;
            }
        }
        let temperature = mixture.temperature_kelvin();
        if condition
            .min_temperature_kelvin
            .is_some_and(|min| temperature < min)
            || condition
                .max_temperature_kelvin
                .is_some_and(|max| temperature > max)
        {
            blocked = true;
        }
        let water_activity = mixture.activity_of(
            registry,
            &SubstanceId::from("destroy:water"),
            MixturePhase::Aqueous,
        )?;
        if condition
            .min_water_activity
            .is_some_and(|min| water_activity < min)
            || condition
                .max_water_activity
                .is_some_and(|max| water_activity > max)
        {
            blocked = true;
        }
        if condition.max_oxygen_activity.is_some()
            || condition.atmosphere == Some(AtmosphereCondition::Inert)
        {
            let oxygen = SubstanceId::from("destroy:oxygen");
            registry.substance(&oxygen)?;
            let oxygen_activity = mixture.activity_of(registry, &oxygen, MixturePhase::Gas)?;
            if condition
                .max_oxygen_activity
                .is_some_and(|max| oxygen_activity > max)
            {
                blocked = true;
            }
            if condition.atmosphere == Some(AtmosphereCondition::Inert)
                && oxygen_activity > TRACE_CONCENTRATION_MOL_PER_BUCKET
            {
                blocked = true;
            }
        }
        if condition.gas_pressure_atm.is_some() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<condition-evaluation>".to_string(),
                reason: "gas pressure conditions require a pressure field in ReactionContext"
                    .to_string(),
            });
        }
        if condition.redox.is_some() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<condition-evaluation>".to_string(),
                reason: "redox reaction conditions must be represented by RedoxAnnotation so electron balance remains explicit".to_string(),
            });
        }
        if blocked {
            blocked_reasons.push(condition.reason.clone());
        } else {
            multiplier *= condition.rate_multiplier;
        }
    }
    Ok(ReactionConditionEvaluation {
        allowed: blocked_reasons.is_empty(),
        multiplier,
        blocked_reasons,
    })
}

fn validate_optional_positive(
    value: Option<f64>,
    reaction_id: &str,
    label: &str,
) -> ChemistryResult<()> {
    if value.is_some_and(|value| !value.is_finite() || value <= 0.0) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: format!("{label} must be positive and finite"),
        });
    }
    Ok(())
}

fn validate_optional_unit(
    value: Option<f64>,
    reaction_id: &str,
    label: &str,
) -> ChemistryResult<()> {
    if value.is_some_and(|value| !value.is_finite() || !(0.0..=1.0).contains(&value)) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: format!("{label} must be finite and within 0.0..=1.0"),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::substance::{Substance, SubstancePhaseProperties};
    use crate::chemistry::ChemistryRegistryBuilder;

    #[test]
    fn acidic_condition_requires_actual_acidic_solution() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1000.0, 373.15, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(Substance::new(
                "destroy:proton",
                1,
                1.0,
                1000.0,
                f64::MAX,
                100.0,
                20_000.0,
            ))
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        let neutral = evaluate_reaction_conditions(
            &registry,
            &mixture,
            &[ReactionCondition::new("acid required").acidity(AcidityCondition::Acidic)],
        )
        .unwrap();
        assert!(!neutral.allowed);

        mixture
            .add_substance(&registry, "destroy:proton", 0.1)
            .unwrap();
        let acidic = evaluate_reaction_conditions(
            &registry,
            &mixture,
            &[ReactionCondition::new("acid required").acidity(AcidityCondition::Acidic)],
        )
        .unwrap();
        assert!(acidic.allowed);
    }
}
