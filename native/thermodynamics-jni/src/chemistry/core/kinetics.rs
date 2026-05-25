use std::collections::BTreeMap;
use std::fmt::{Display, Formatter};

use super::catalysis::CatalystSurfaceId;
use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, MixturePhase, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::reaction::{
    ProductDistribution, ProductDistributionVariant, ReactionId, StoichiometricTerm,
    GAS_CONSTANT_J_PER_MOL_KELVIN,
};
use super::simulation::ReactionContext;
use super::substance::SubstanceId;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReactionChannelId(String);

impl ReactionChannelId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<channel>".to_string(),
                reason: "reaction channel id must not be empty".to_string(),
            });
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for ReactionChannelId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for ReactionChannelId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

impl From<String> for ReactionChannelId {
    fn from(value: String) -> Self {
        Self(value)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionChannelMode {
    Kinetic,
    Thermodynamic,
    Mixed,
    Photochemical,
    Surface,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum LightBand {
    Ultraviolet,
    Visible,
    Infrared,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ChannelConditionEffect {
    Phase {
        phase: MixturePhase,
        multiplier: f64,
    },
    Light {
        band: LightBand,
        minimum_power: f64,
        multiplier: f64,
    },
    Surface {
        surface_id: CatalystSurfaceId,
        multiplier: f64,
    },
}

#[derive(Debug, Clone)]
pub struct ReactionChannel {
    pub id: ReactionChannelId,
    pub products: Vec<StoichiometricTerm>,
    pub activation_gibbs_kj_per_mol: f64,
    pub pre_exponential_factor: f64,
    pub mode: ReactionChannelMode,
    pub condition_effects: Vec<ChannelConditionEffect>,
}

#[derive(Debug, Clone)]
pub struct IsomerEnergy {
    pub substance_id: SubstanceId,
    pub relative_gibbs_kj_per_mol: f64,
    pub phase: Option<MixturePhase>,
    pub surface_id: Option<CatalystSurfaceId>,
}

#[derive(Debug, Clone)]
pub struct TransitionStateEnergy {
    pub channel_id: ReactionChannelId,
    pub activation_gibbs_kj_per_mol: f64,
    pub phase: Option<MixturePhase>,
    pub surface_id: Option<CatalystSurfaceId>,
}

#[derive(Debug, Clone, Default)]
pub struct EnergyModel {
    isomer_energies: BTreeMap<(SubstanceId, Option<MixturePhase>, Option<CatalystSurfaceId>), f64>,
    transition_state_energies: BTreeMap<
        (
            ReactionChannelId,
            Option<MixturePhase>,
            Option<CatalystSurfaceId>,
        ),
        f64,
    >,
}

impl ReactionChannel {
    pub fn new(
        id: impl Into<ReactionChannelId>,
        products: impl IntoIterator<Item = StoichiometricTerm>,
        activation_gibbs_kj_per_mol: f64,
    ) -> Self {
        Self {
            id: id.into(),
            products: products.into_iter().collect(),
            activation_gibbs_kj_per_mol,
            pre_exponential_factor: 10_000.0,
            mode: ReactionChannelMode::Kinetic,
            condition_effects: Vec::new(),
        }
    }

    pub fn with_pre_exponential_factor(mut self, value: f64) -> Self {
        self.pre_exponential_factor = value;
        self
    }

    pub fn with_mode(mut self, mode: ReactionChannelMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn with_condition_effect(mut self, effect: ChannelConditionEffect) -> Self {
        self.condition_effects.push(effect);
        self
    }

    pub fn validate(&self, reaction_id: &ReactionId) -> ChemistryResult<()> {
        if self.id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "reaction channel id must not be empty".to_string(),
            });
        }
        if self.products.is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: format!("reaction channel '{}' must contain products", self.id),
            });
        }
        for term in &self.products {
            if term.coefficient == 0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.to_string(),
                    reason: format!(
                        "reaction channel '{}' product coefficients must be greater than zero",
                        self.id
                    ),
                });
            }
        }
        if !self.activation_gibbs_kj_per_mol.is_finite() || self.activation_gibbs_kj_per_mol < 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: format!(
                    "reaction channel '{}' activation Gibbs energy must be non-negative and finite",
                    self.id
                ),
            });
        }
        if !self.pre_exponential_factor.is_finite() || self.pre_exponential_factor <= 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: format!(
                    "reaction channel '{}' pre-exponential factor must be positive and finite",
                    self.id
                ),
            });
        }
        for effect in &self.condition_effects {
            validate_condition_effect(reaction_id, &self.id, effect)?;
        }
        Ok(())
    }

    pub fn rate_weight(
        &self,
        mixture: &Mixture,
        context: &ReactionContext,
    ) -> ChemistryResult<f64> {
        channel_rate_weight(
            self.pre_exponential_factor,
            self.activation_gibbs_kj_per_mol,
            &self.condition_effects,
            mixture,
            context,
        )
    }
}

impl EnergyModel {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_isomer_energy(mut self, energy: IsomerEnergy) -> ChemistryResult<Self> {
        validate_energy(energy.relative_gibbs_kj_per_mol, "isomer Gibbs energy")?;
        let key = (energy.substance_id, energy.phase, energy.surface_id);
        if self
            .isomer_energies
            .insert(key, energy.relative_gibbs_kj_per_mol)
            .is_some()
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<energy-model>".to_string(),
                reason: "duplicate isomer energy entry".to_string(),
            });
        }
        Ok(self)
    }

    pub fn with_transition_state_energy(
        mut self,
        energy: TransitionStateEnergy,
    ) -> ChemistryResult<Self> {
        validate_energy(
            energy.activation_gibbs_kj_per_mol,
            "transition-state Gibbs energy",
        )?;
        let key = (energy.channel_id, energy.phase, energy.surface_id);
        if self
            .transition_state_energies
            .insert(key, energy.activation_gibbs_kj_per_mol)
            .is_some()
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<energy-model>".to_string(),
                reason: "duplicate transition-state energy entry".to_string(),
            });
        }
        Ok(self)
    }

    pub fn isomer_energy(
        &self,
        substance_id: &SubstanceId,
        phase: Option<MixturePhase>,
        surface_id: Option<&CatalystSurfaceId>,
    ) -> Option<f64> {
        let surface_id = surface_id.cloned();
        self.isomer_energies
            .get(&(substance_id.clone(), phase, surface_id.clone()))
            .copied()
            .or_else(|| {
                self.isomer_energies
                    .get(&(substance_id.clone(), None, surface_id))
                    .copied()
            })
            .or_else(|| {
                self.isomer_energies
                    .get(&(substance_id.clone(), phase, None))
                    .copied()
            })
            .or_else(|| {
                self.isomer_energies
                    .get(&(substance_id.clone(), None, None))
                    .copied()
            })
    }

    pub fn transition_state_energy(
        &self,
        channel_id: &ReactionChannelId,
        phase: Option<MixturePhase>,
        surface_id: Option<&CatalystSurfaceId>,
    ) -> Option<f64> {
        let surface_id = surface_id.cloned();
        self.transition_state_energies
            .get(&(channel_id.clone(), phase, surface_id.clone()))
            .copied()
            .or_else(|| {
                self.transition_state_energies
                    .get(&(channel_id.clone(), None, surface_id))
                    .copied()
            })
            .or_else(|| {
                self.transition_state_energies
                    .get(&(channel_id.clone(), phase, None))
                    .copied()
            })
            .or_else(|| {
                self.transition_state_energies
                    .get(&(channel_id.clone(), None, None))
                    .copied()
            })
    }

    pub fn equilibrium_distribution<I>(
        &self,
        substances: I,
        temperature_kelvin: f64,
        phase: Option<MixturePhase>,
    ) -> ChemistryResult<BTreeMap<SubstanceId, f64>>
    where
        I: IntoIterator<Item = SubstanceId>,
    {
        validate_temperature(temperature_kelvin)?;
        let mut weights = Vec::new();
        for substance_id in substances {
            let Some(gibbs) = self.isomer_energy(&substance_id, phase, None) else {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: "<energy-model>".to_string(),
                    reason: format!("missing isomer energy for '{substance_id}'"),
                });
            };
            let weight =
                (-(gibbs * 1000.0) / (GAS_CONSTANT_J_PER_MOL_KELVIN * temperature_kelvin)).exp();
            weights.push((substance_id, weight));
        }
        normalize_weights("<energy-model>", weights)
    }
}

pub fn channel_product_distribution(
    reaction_id: &ReactionId,
    channels: &[ReactionChannel],
    mixture: &Mixture,
    context: &ReactionContext,
) -> ChemistryResult<ProductDistribution> {
    let weights = channels
        .iter()
        .map(|channel| {
            Ok((
                channel.id.to_string(),
                channel.rate_weight(mixture, context)?,
            ))
        })
        .collect::<ChemistryResult<Vec<_>>>()?;
    let fractions = normalize_weights(reaction_id.as_str(), weights)?;
    Ok(ProductDistribution {
        variants: channels
            .iter()
            .map(|channel| ProductDistributionVariant {
                fraction: fractions
                    .get(channel.id.as_str())
                    .copied()
                    .unwrap_or_default(),
                products: channel.products.clone(),
            })
            .filter(|variant| variant.fraction > 0.0)
            .collect(),
    })
}

pub fn channel_rate_sum_per_second(
    reaction_id: &ReactionId,
    channels: &[ReactionChannel],
    mixture: &Mixture,
    context: &ReactionContext,
) -> ChemistryResult<f64> {
    let mut total = 0.0;
    for channel in channels {
        total += channel.rate_weight(mixture, context)?;
    }
    if !total.is_finite() || total < 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: "calculated channel rate must be non-negative and finite".to_string(),
        });
    }
    Ok(total)
}

fn channel_rate_weight(
    pre_exponential_factor: f64,
    activation_gibbs_kj_per_mol: f64,
    condition_effects: &[ChannelConditionEffect],
    mixture: &Mixture,
    context: &ReactionContext,
) -> ChemistryResult<f64> {
    validate_temperature(mixture.temperature_kelvin())?;
    let mut weight = pre_exponential_factor
        * (-(activation_gibbs_kj_per_mol * 1000.0)
            / (GAS_CONSTANT_J_PER_MOL_KELVIN * mixture.temperature_kelvin()))
        .exp();
    for effect in condition_effects {
        weight *= condition_multiplier(effect, mixture, context)?;
        if weight == 0.0 {
            return Ok(0.0);
        }
    }
    if !weight.is_finite() || weight < 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<channel>".to_string(),
            reason: "channel weight must be non-negative and finite".to_string(),
        });
    }
    Ok(weight)
}

fn condition_multiplier(
    effect: &ChannelConditionEffect,
    mixture: &Mixture,
    context: &ReactionContext,
) -> ChemistryResult<f64> {
    match effect {
        ChannelConditionEffect::Phase { phase, multiplier } => {
            if mixture.total_concentration_in_phase(*phase) <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                Ok(0.0)
            } else {
                Ok(*multiplier)
            }
        }
        ChannelConditionEffect::Light {
            band,
            minimum_power,
            multiplier,
        } => {
            let power = context.light_power(*band);
            if power < *minimum_power {
                Ok(0.0)
            } else {
                Ok(power * *multiplier)
            }
        }
        ChannelConditionEffect::Surface {
            surface_id,
            multiplier,
        } => {
            let free = context.free_sites(surface_id);
            if free <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                Ok(0.0)
            } else {
                Ok(*multiplier)
            }
        }
    }
}

fn validate_condition_effect(
    reaction_id: &ReactionId,
    channel_id: &ReactionChannelId,
    effect: &ChannelConditionEffect,
) -> ChemistryResult<()> {
    match effect {
        ChannelConditionEffect::Phase { multiplier, .. }
        | ChannelConditionEffect::Surface { multiplier, .. } => {
            if !multiplier.is_finite() || *multiplier < 0.0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.to_string(),
                    reason: format!(
                        "reaction channel '{channel_id}' condition multiplier must be non-negative and finite"
                    ),
                });
            }
        }
        ChannelConditionEffect::Light {
            minimum_power,
            multiplier,
            ..
        } => {
            if !minimum_power.is_finite() || *minimum_power < 0.0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.to_string(),
                    reason: format!(
                        "reaction channel '{channel_id}' minimum light power must be non-negative and finite"
                    ),
                });
            }
            if !multiplier.is_finite() || *multiplier < 0.0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.to_string(),
                    reason: format!(
                        "reaction channel '{channel_id}' light multiplier must be non-negative and finite"
                    ),
                });
            }
        }
    }
    Ok(())
}

fn normalize_weights<K>(
    reaction_id: &str,
    weights: Vec<(K, f64)>,
) -> ChemistryResult<BTreeMap<K, f64>>
where
    K: Ord + Display,
{
    let total = weights.iter().map(|(_, weight)| *weight).sum::<f64>();
    if !total.is_finite() || total <= 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: "cannot calculate distribution because all channel weights are zero or invalid"
                .to_string(),
        });
    }
    let mut result = BTreeMap::new();
    for (id, weight) in weights {
        if !weight.is_finite() || weight < 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: format!("weight for '{id}' must be non-negative and finite"),
            });
        }
        if weight > 0.0 {
            result.insert(id, weight / total);
        }
    }
    Ok(result)
}

fn validate_temperature(temperature_kelvin: f64) -> ChemistryResult<()> {
    if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<kinetics>".to_string(),
            reason: "temperature must be positive and finite".to_string(),
        });
    }
    Ok(())
}

fn validate_energy(value: f64, name: &str) -> ChemistryResult<()> {
    if !value.is_finite() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<energy-model>".to_string(),
            reason: format!("{name} must be finite"),
        });
    }
    Ok(())
}
