use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, MixturePhase};
use super::reaction::{Reaction, StoichiometricTerm, GAS_CONSTANT_J_PER_MOL_KELVIN};
use super::registry::ChemistryRegistry;

pub const STANDARD_TEMPERATURE_KELVIN: f64 = 298.15;
const QUOTIENT_ACTIVITY_FLOOR: f64 = 1.0e-30;

#[derive(Debug, Clone, PartialEq)]
pub struct ReactionThermodynamics {
    pub reference_temperature_kelvin: f64,
    pub gibbs_free_energy_change_kj_per_mol: f64,
}

impl ReactionThermodynamics {
    pub fn from_gibbs_free_energy_change_kj_per_mol(value: f64) -> Self {
        Self {
            reference_temperature_kelvin: STANDARD_TEMPERATURE_KELVIN,
            gibbs_free_energy_change_kj_per_mol: value,
        }
    }

    pub fn from_gibbs_free_energy_change_at_kelvin(value: f64, temperature_kelvin: f64) -> Self {
        Self {
            reference_temperature_kelvin: temperature_kelvin,
            gibbs_free_energy_change_kj_per_mol: value,
        }
    }

    pub fn from_equilibrium_constant(value: f64) -> ChemistryResult<Self> {
        Self::from_equilibrium_constant_at_kelvin(value, STANDARD_TEMPERATURE_KELVIN)
    }

    pub fn from_equilibrium_constant_at_kelvin(
        value: f64,
        temperature_kelvin: f64,
    ) -> ChemistryResult<Self> {
        validate_temperature(temperature_kelvin, "thermodynamic reference temperature")?;
        validate_equilibrium_constant(value, "reaction equilibrium constant")?;
        Ok(Self {
            reference_temperature_kelvin: temperature_kelvin,
            gibbs_free_energy_change_kj_per_mol: delta_g_from_equilibrium_constant(
                value,
                temperature_kelvin,
            )?,
        })
    }

    pub fn validate(&self, reaction_id: &str) -> ChemistryResult<()> {
        validate_temperature(
            self.reference_temperature_kelvin,
            "thermodynamic reference temperature",
        )
        .map_err(|error| remap_invalid_mixture_to_reaction(error, reaction_id))?;
        if !self.gibbs_free_energy_change_kj_per_mol.is_finite() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "Gibbs free energy change must be finite".to_string(),
            });
        }
        Ok(())
    }

    pub fn entropy_change_j_per_mol_kelvin(&self, enthalpy_change_kj_per_mol: f64) -> f64 {
        (enthalpy_change_kj_per_mol - self.gibbs_free_energy_change_kj_per_mol) * 1000.0
            / self.reference_temperature_kelvin
    }

    pub fn gibbs_free_energy_change_at_kelvin(
        &self,
        enthalpy_change_kj_per_mol: f64,
        temperature_kelvin: f64,
    ) -> ChemistryResult<f64> {
        validate_temperature(temperature_kelvin, "reaction temperature")?;
        if !enthalpy_change_kj_per_mol.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(
                "reaction enthalpy change must be finite".to_string(),
            ));
        }
        let entropy_change = self.entropy_change_j_per_mol_kelvin(enthalpy_change_kj_per_mol);
        let value = enthalpy_change_kj_per_mol - temperature_kelvin * entropy_change / 1000.0;
        if !value.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(
                "Gibbs free energy change became non-finite".to_string(),
            ));
        }
        Ok(value)
    }

    pub fn equilibrium_constant_at_kelvin(
        &self,
        enthalpy_change_kj_per_mol: f64,
        temperature_kelvin: f64,
    ) -> ChemistryResult<f64> {
        let delta_g =
            self.gibbs_free_energy_change_at_kelvin(enthalpy_change_kj_per_mol, temperature_kelvin)?;
        equilibrium_constant_from_delta_g(delta_g, temperature_kelvin)
    }
}

pub fn equilibrium_constant_from_delta_g(
    delta_g_kj_per_mol: f64,
    temperature_kelvin: f64,
) -> ChemistryResult<f64> {
    validate_temperature(temperature_kelvin, "thermodynamic temperature")?;
    if !delta_g_kj_per_mol.is_finite() {
        return Err(ChemistryError::InvalidMixtureState(
            "Gibbs free energy change must be finite".to_string(),
        ));
    }
    let value =
        (-(delta_g_kj_per_mol * 1000.0) / (GAS_CONSTANT_J_PER_MOL_KELVIN * temperature_kelvin))
            .exp();
    validate_equilibrium_constant(value, "derived equilibrium constant")?;
    Ok(value)
}

pub fn delta_g_from_equilibrium_constant(
    equilibrium_constant: f64,
    temperature_kelvin: f64,
) -> ChemistryResult<f64> {
    validate_temperature(temperature_kelvin, "thermodynamic temperature")?;
    validate_equilibrium_constant(equilibrium_constant, "equilibrium constant")?;
    Ok(-GAS_CONSTANT_J_PER_MOL_KELVIN
        * temperature_kelvin
        * equilibrium_constant.ln()
        / 1000.0)
}

pub fn reaction_thermodynamic_rate_factor(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    reaction: &Reaction,
) -> ChemistryResult<f64> {
    let Some(thermodynamics) = &reaction.thermodynamics else {
        return Ok(1.0);
    };
    let equilibrium_constant = thermodynamics.equilibrium_constant_at_kelvin(
        reaction.enthalpy_change_kj_per_mol,
        mixture.temperature_kelvin(),
    )?;
    let quotient = reaction_quotient(registry, mixture, reaction)?;
    if quotient <= 0.0 {
        return Ok(1.0);
    }
    let ratio = quotient / equilibrium_constant;
    if !ratio.is_finite() || ratio < 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "reaction thermodynamic quotient must be non-negative and finite".to_string(),
        });
    }
    Ok((1.0 - ratio).clamp(0.0, 1.0))
}

pub fn reaction_quotient(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    reaction: &Reaction,
) -> ChemistryResult<f64> {
    if !reaction.channels.is_empty() || reaction.product_distribution.is_some() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "thermodynamic reaction quotient requires a single product set".to_string(),
        });
    }
    let products = activity_product(registry, mixture, reaction, &reaction.products, true)?;
    let reactants = activity_product(registry, mixture, reaction, &reaction.reactants, false)?;
    if reactants <= QUOTIENT_ACTIVITY_FLOOR {
        return Ok(f64::INFINITY);
    }
    let quotient = products / reactants;
    if !quotient.is_finite() && quotient != f64::INFINITY {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "reaction quotient must be finite or positive infinity".to_string(),
        });
    }
    Ok(quotient)
}

fn activity_product(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    reaction: &Reaction,
    terms: &[StoichiometricTerm],
    products: bool,
) -> ChemistryResult<f64> {
    let mut product = 1.0;
    for term in terms {
        let phases = if products {
            vec![reaction
                .product_phases
                .get(&term.substance_id)
                .copied()
                .unwrap_or(MixturePhase::Aqueous)]
        } else {
            reaction
                .phase_access
                .get(&term.substance_id)
                .cloned()
                .unwrap_or_else(super::reaction::ReactionPhaseAccess::liquid)
                .phases
        };
        let activity = phases
            .iter()
            .map(|phase| mixture.activity_of(registry, &term.substance_id, *phase))
            .sum::<ChemistryResult<f64>>()?
            .max(0.0);
        if activity <= QUOTIENT_ACTIVITY_FLOOR {
            return Ok(0.0);
        }
        product *= activity.powi(term.coefficient as i32);
        if !product.is_finite() || product < 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "reaction activity product must be non-negative and finite".to_string(),
            });
        }
    }
    Ok(product)
}

fn validate_temperature(value: f64, name: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{name} must be positive and finite"
        )));
    }
    Ok(())
}

fn validate_equilibrium_constant(value: f64, name: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{name} must be positive and finite"
        )));
    }
    Ok(())
}

fn remap_invalid_mixture_to_reaction(error: ChemistryError, reaction_id: &str) -> ChemistryError {
    match error {
        ChemistryError::InvalidMixtureState(reason) => ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason,
        },
        other => other,
    }
}
