use std::collections::BTreeMap;

use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::reaction::Reaction;
use super::registry::{ChemistryRegistry, IndexedReaction, SubstanceIndex};

pub const TICKS_PER_SECOND: f64 = 20.0;
pub const EQUILIBRIUM_EPSILON_MOL_PER_BUCKET: f64 = TRACE_CONCENTRATION_MOL_PER_BUCKET;

#[derive(Debug, Clone)]
pub struct SimulationReport {
    pub ticks: u32,
    pub reached_equilibrium: bool,
    pub reaction_results: BTreeMap<String, f64>,
}

#[derive(Debug, Clone)]
pub struct ReactionContext {
    pub uv_power: f64,
    pub external_reactants: BTreeMap<String, f64>,
    pub external_catalysts: BTreeMap<String, f64>,
    pub reaction_results: BTreeMap<String, f64>,
}

impl Default for ReactionContext {
    fn default() -> Self {
        Self {
            uv_power: 0.0,
            external_reactants: BTreeMap::new(),
            external_catalysts: BTreeMap::new(),
            reaction_results: BTreeMap::new(),
        }
    }
}

impl ReactionContext {
    pub fn with_uv_power(mut self, uv_power: f64) -> ChemistryResult<Self> {
        if !uv_power.is_finite() || uv_power < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "UV power must be non-negative and finite".to_string(),
            ));
        }
        self.uv_power = uv_power;
        Ok(self)
    }

    pub fn add_external_reactant(
        &mut self,
        description: impl Into<String>,
        moles_per_bucket: f64,
    ) -> ChemistryResult<()> {
        add_external(&mut self.external_reactants, description, moles_per_bucket)
    }

    pub fn add_external_catalyst(
        &mut self,
        description: impl Into<String>,
        moles_per_bucket: f64,
    ) -> ChemistryResult<()> {
        add_external(&mut self.external_catalysts, description, moles_per_bucket)
    }
}

pub fn react_for_tick(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    cycles: u32,
) -> ChemistryResult<bool> {
    let mut context = ReactionContext::default();
    react_for_tick_with_context(registry, mixture, &mut context, cycles)
}

pub fn react_for_tick_with_context(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    context: &mut ReactionContext,
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
        for reaction_index in
            registry.reaction_candidate_indices_for_substance_indices(mixture.component_indices())
        {
            let reaction = registry.reaction_by_index(reaction_index)?;
            let indexed_reaction = registry.indexed_reaction(reaction_index)?;
            if !context_allows_reaction(reaction, context) {
                continue;
            }
            let rate = reaction_rate_mol_per_bucket_per_tick_for_indexed_reaction(
                mixture,
                reaction,
                indexed_reaction,
                context,
            )? / cycles as f64;
            if rate > 0.0 {
                reactions_with_rates.push((reaction_index, rate));
            }
        }
        reactions_with_rates.sort_by(|(left_index, left), (right_index, right)| {
            let left_reaction = registry
                .reaction_by_index(*left_index)
                .expect("reaction candidate index must be valid");
            let right_reaction = registry
                .reaction_by_index(*right_index)
                .expect("reaction candidate index must be valid");
            right
                .total_cmp(left)
                .then_with(|| left_reaction.id.as_str().cmp(right_reaction.id.as_str()))
        });

        for (reaction_index, rate) in reactions_with_rates {
            let reaction = registry.reaction_by_index(reaction_index)?;
            let indexed_reaction = registry.indexed_reaction(reaction_index)?;
            let limited =
                limit_by_reactants_and_context(mixture, context, reaction, indexed_reaction, rate);
            if limited > 0.0 {
                apply_reaction(
                    registry,
                    mixture,
                    context,
                    reaction,
                    indexed_reaction,
                    limited,
                )?;
            }
        }

        mixture.validate(registry)?;
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
    let mut context = ReactionContext::default();
    react_until_equilibrium_with_context(
        registry,
        mixture,
        &mut context,
        max_ticks,
        cycles_per_tick,
    )
}

pub fn react_until_equilibrium_with_context(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    context: &mut ReactionContext,
    max_ticks: u32,
    cycles_per_tick: u32,
) -> ChemistryResult<SimulationReport> {
    for tick in 0..max_ticks {
        let changed = react_for_tick_with_context(registry, mixture, context, cycles_per_tick)?;
        if !changed {
            return Ok(SimulationReport {
                ticks: tick + 1,
                reached_equilibrium: true,
                reaction_results: context.reaction_results.clone(),
            });
        }
    }
    Ok(SimulationReport {
        ticks: max_ticks,
        reached_equilibrium: false,
        reaction_results: context.reaction_results.clone(),
    })
}

pub fn reaction_rate_mol_per_bucket_per_tick(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    reaction: &Reaction,
) -> ChemistryResult<f64> {
    reaction_rate_mol_per_bucket_per_tick_with_context(
        registry,
        mixture,
        reaction,
        &ReactionContext::default(),
    )
}

pub fn reaction_rate_mol_per_bucket_per_tick_with_context(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    reaction: &Reaction,
    context: &ReactionContext,
) -> ChemistryResult<f64> {
    let mut rate =
        reaction.rate_constant_per_second(mixture.temperature_kelvin())? / TICKS_PER_SECOND;
    if reaction.requires_uv {
        rate *= context.uv_power;
    }
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

fn reaction_rate_mol_per_bucket_per_tick_for_indexed_reaction(
    mixture: &Mixture,
    reaction: &Reaction,
    indexed_reaction: &IndexedReaction,
    context: &ReactionContext,
) -> ChemistryResult<f64> {
    let mut rate =
        reaction.rate_constant_per_second(mixture.temperature_kelvin())? / TICKS_PER_SECOND;
    if reaction.requires_uv {
        rate *= context.uv_power;
    }
    for (substance, order) in &indexed_reaction.orders {
        let concentration = mixture.concentration_of_index(*substance);
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

fn context_allows_reaction(reaction: &Reaction, context: &ReactionContext) -> bool {
    if reaction.requires_uv && context.uv_power <= 0.0 {
        return false;
    }
    for catalyst in &reaction.external_catalysts {
        if context
            .external_catalysts
            .get(&catalyst.description)
            .copied()
            .unwrap_or(0.0)
            < catalyst.moles_per_reaction
        {
            return false;
        }
    }
    for reactant in &reaction.external_reactants {
        if context
            .external_reactants
            .get(&reactant.description)
            .copied()
            .unwrap_or(0.0)
            <= 0.0
        {
            return false;
        }
    }
    true
}

fn limit_by_reactants_and_context(
    mixture: &Mixture,
    context: &ReactionContext,
    reaction: &Reaction,
    indexed_reaction: &IndexedReaction,
    requested_moles: f64,
) -> f64 {
    let molecule_limited =
        indexed_reaction
            .reactants
            .iter()
            .fold(requested_moles, |current, reactant| {
                let available = mixture.concentration_of_index(reactant.substance)
                    / reactant.coefficient as f64;
                current.min(available)
            });
    reaction
        .external_reactants
        .iter()
        .fold(molecule_limited, |current, reactant| {
            let available = context
                .external_reactants
                .get(&reactant.description)
                .copied()
                .unwrap_or(0.0)
                / reactant.moles_per_reaction;
            current.min(available)
        })
}

fn apply_reaction(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    context: &mut ReactionContext,
    reaction: &Reaction,
    indexed_reaction: &IndexedReaction,
    moles_per_bucket: f64,
) -> ChemistryResult<()> {
    if !moles_per_bucket.is_finite() || moles_per_bucket < 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "moles to apply must be non-negative and finite".to_string(),
        });
    }
    let previous_mixture = mixture.clone();
    let previous_context = context.clone();
    let result = apply_reaction_inner(
        registry,
        mixture,
        context,
        reaction,
        indexed_reaction,
        moles_per_bucket,
    );
    if result.is_err() {
        *mixture = previous_mixture;
        *context = previous_context;
    }
    result
}

fn apply_reaction_inner(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    context: &mut ReactionContext,
    reaction: &Reaction,
    indexed_reaction: &IndexedReaction,
    moles_per_bucket: f64,
) -> ChemistryResult<()> {
    let mut deltas = Vec::new();
    for reactant in &indexed_reaction.reactants {
        add_delta(
            &mut deltas,
            reactant.substance,
            -(reactant.coefficient as f64) * moles_per_bucket,
        );
    }
    for product in &indexed_reaction.products {
        add_delta(
            &mut deltas,
            product.substance,
            (product.coefficient as f64) * moles_per_bucket,
        );
    }
    mixture.apply_concentration_deltas_by_index(registry, &deltas)?;
    apply_external_reactants(context, reaction, moles_per_bucket)?;
    for result in &reaction.reaction_results {
        *context
            .reaction_results
            .entry(result.description.clone())
            .or_insert(0.0) += result.moles_per_reaction * moles_per_bucket;
    }
    mixture.heat(
        registry,
        -reaction.enthalpy_change_kj_per_mol * 1000.0 * moles_per_bucket,
    )?;
    Ok(())
}

fn add_delta(deltas: &mut Vec<(SubstanceIndex, f64)>, substance: SubstanceIndex, delta: f64) {
    if let Some((_, existing)) = deltas
        .iter_mut()
        .find(|(existing_substance, _)| *existing_substance == substance)
    {
        *existing += delta;
    } else {
        deltas.push((substance, delta));
    }
}

fn apply_external_reactants(
    context: &mut ReactionContext,
    reaction: &Reaction,
    moles_per_bucket: f64,
) -> ChemistryResult<()> {
    for external in &reaction.external_reactants {
        let current = context
            .external_reactants
            .get(&external.description)
            .copied()
            .unwrap_or(0.0);
        let next = current - external.moles_per_reaction * moles_per_bucket;
        if next < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "external reactant '{}' would become negative: {next}",
                external.description
            )));
        }
        if next <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            context.external_reactants.remove(&external.description);
        } else {
            context
                .external_reactants
                .insert(external.description.clone(), next);
        }
    }
    Ok(())
}

fn add_external(
    target: &mut BTreeMap<String, f64>,
    description: impl Into<String>,
    moles_per_bucket: f64,
) -> ChemistryResult<()> {
    if !moles_per_bucket.is_finite() || moles_per_bucket < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "external amount must be non-negative and finite".to_string(),
        ));
    }
    let description = description.into();
    if description.trim().is_empty() {
        return Err(ChemistryError::InvalidMixtureState(
            "external description must not be empty".to_string(),
        ));
    }
    *target.entry(description).or_insert(0.0) += moles_per_bucket;
    Ok(())
}

fn snapshot(mixture: &Mixture) -> Vec<(SubstanceIndex, f64)> {
    mixture
        .component_indices()
        .map(|substance| (substance, mixture.concentration_of_index(substance)))
        .collect()
}

fn changed_since(mixture: &Mixture, before: &[(SubstanceIndex, f64)]) -> bool {
    for (substance, previous) in before {
        let current = mixture.concentration_of_index(*substance);
        if (current - previous).abs() > EQUILIBRIUM_EPSILON_MOL_PER_BUCKET {
            return true;
        }
    }
    for substance in mixture.component_indices() {
        if !before
            .iter()
            .any(|(before_substance, _)| *before_substance == substance)
            && mixture.concentration_of_index(substance) > EQUILIBRIUM_EPSILON_MOL_PER_BUCKET
        {
            return true;
        }
    }
    false
}
