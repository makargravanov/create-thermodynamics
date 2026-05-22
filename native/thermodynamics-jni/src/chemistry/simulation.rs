use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::reaction::Reaction;
use super::registry::ChemistryRegistry;

pub const TICKS_PER_SECOND: f64 = 20.0;
pub const EQUILIBRIUM_EPSILON_MOL_PER_BUCKET: f64 = TRACE_CONCENTRATION_MOL_PER_BUCKET;

#[derive(Debug, Clone)]
pub struct SimulationReport {
    pub ticks: u32,
    pub reached_equilibrium: bool,
}

pub fn react_for_tick(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    cycles: u32,
) -> ChemistryResult<bool> {
    if cycles == 0 {
        return Err(ChemistryError::InvalidMixtureState(
            "cycles must be greater than zero".to_string(),
        ));
    }

    let mut any_changed = false;
    for _ in 0..cycles {
        let before = snapshot(mixture);
        let mut reactions_with_rates = Vec::new();
        for reaction in registry.reactions() {
            if reaction.has_external_context() {
                continue;
            }
            let rate =
                reaction_rate_mol_per_bucket_per_tick(registry, mixture, reaction)? / cycles as f64;
            if rate > 0.0 {
                reactions_with_rates.push((reaction, rate));
            }
        }
        reactions_with_rates.sort_by(|(_, left), (_, right)| right.total_cmp(left));

        for (reaction, rate) in reactions_with_rates {
            let limited = limit_by_reactants(mixture, reaction, rate);
            if limited > 0.0 {
                apply_reaction(registry, mixture, reaction, limited)?;
            }
        }

        if changed_since(mixture, &before) {
            any_changed = true;
        }
    }
    Ok(any_changed)
}

pub fn react_until_equilibrium(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    max_ticks: u32,
    cycles_per_tick: u32,
) -> ChemistryResult<SimulationReport> {
    for tick in 0..max_ticks {
        let changed = react_for_tick(registry, mixture, cycles_per_tick)?;
        if !changed {
            return Ok(SimulationReport {
                ticks: tick + 1,
                reached_equilibrium: true,
            });
        }
    }
    Ok(SimulationReport {
        ticks: max_ticks,
        reached_equilibrium: false,
    })
}

pub fn reaction_rate_mol_per_bucket_per_tick(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    reaction: &Reaction,
) -> ChemistryResult<f64> {
    let mut rate =
        reaction.rate_constant_per_second(mixture.temperature_kelvin())? / TICKS_PER_SECOND;
    for (substance_id, order) in &reaction.orders {
        registry.substance(substance_id)?;
        let concentration = mixture.concentration_of(substance_id);
        if concentration <= 0.0 {
            return Ok(0.0);
        }
        rate *= concentration.powi(*order as i32);
    }
    if !rate.is_finite() || rate < 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "calculated reaction rate must be non-negative and finite".to_string(),
        });
    }
    Ok(rate)
}

fn limit_by_reactants(mixture: &Mixture, reaction: &Reaction, requested_moles: f64) -> f64 {
    reaction
        .reactants
        .iter()
        .fold(requested_moles, |current, reactant| {
            let available =
                mixture.concentration_of(&reactant.substance_id) / reactant.coefficient as f64;
            current.min(available)
        })
}

fn apply_reaction(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    reaction: &Reaction,
    moles_per_bucket: f64,
) -> ChemistryResult<()> {
    if !moles_per_bucket.is_finite() || moles_per_bucket < 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "moles to apply must be non-negative and finite".to_string(),
        });
    }
    for reactant in &reaction.reactants {
        mixture.change_concentration(
            registry,
            &reactant.substance_id,
            -(reactant.coefficient as f64) * moles_per_bucket,
        )?;
    }
    for product in &reaction.products {
        mixture.change_concentration(
            registry,
            &product.substance_id,
            (product.coefficient as f64) * moles_per_bucket,
        )?;
    }
    mixture.heat(
        registry,
        -reaction.enthalpy_change_kj_per_mol * 1000.0 * moles_per_bucket,
    )?;
    Ok(())
}

fn snapshot(mixture: &Mixture) -> Vec<(String, f64)> {
    mixture
        .substances()
        .map(|id| (id.to_string(), mixture.concentration_of(id)))
        .collect()
}

fn changed_since(mixture: &Mixture, before: &[(String, f64)]) -> bool {
    for (id, previous) in before {
        let current = mixture.concentration_of(&id.as_str().into());
        if (current - previous).abs() > EQUILIBRIUM_EPSILON_MOL_PER_BUCKET {
            return true;
        }
    }
    for id in mixture.substances() {
        if !before.iter().any(|(before_id, _)| before_id == id.as_str())
            && mixture.concentration_of(id) > EQUILIBRIUM_EPSILON_MOL_PER_BUCKET
        {
            return true;
        }
    }
    false
}
