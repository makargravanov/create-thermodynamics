use std::collections::BTreeMap;

use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, MixturePhase, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::reaction::Reaction;
use super::registry::{
    ChemistryRegistry, IndexedReaction, ReactionCandidateScratch, SubstanceIndex,
};

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

#[derive(Debug, Copy, Clone, PartialEq)]
struct ReactionApplication {
    max_concentration_delta: f64,
    thermal_changed: bool,
}

impl ReactionApplication {
    fn changed(self) -> bool {
        self.max_concentration_delta > EQUILIBRIUM_EPSILON_MOL_PER_BUCKET || self.thermal_changed
    }
}

#[derive(Debug, Clone)]
struct ContextCheckpoint {
    external_reactants: Vec<(String, Option<f64>)>,
    reaction_results: Vec<(String, Option<f64>)>,
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
    let mut candidate_scratch = ReactionCandidateScratch::new();
    let mut reactions_with_rates = Vec::new();
    for _ in 0..cycles {
        registry.collect_reaction_candidate_indices_for_substance_indices(
            mixture.component_indices(),
            &mut candidate_scratch,
        );
        reactions_with_rates.clear();
        for &reaction_index in candidate_scratch.candidates() {
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

        for (reaction_index, rate) in reactions_with_rates.drain(..) {
            let reaction = registry.reaction_by_index(reaction_index)?;
            let indexed_reaction = registry.indexed_reaction(reaction_index)?;
            let limited =
                limit_by_reactants_and_context(mixture, context, reaction, indexed_reaction, rate);
            if limited > 0.0 {
                let application = apply_reaction(
                    registry,
                    mixture,
                    context,
                    reaction,
                    indexed_reaction,
                    limited,
                )?;
                if application.changed() {
                    any_changed = true;
                }
            }
        }

        mixture.validate(registry)?;
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
        let phases = reaction
            .phase_access
            .get(substance_id)
            .cloned()
            .unwrap_or_else(super::reaction::ReactionPhaseAccess::liquid)
            .phases;
        let concentration = phases
            .iter()
            .map(|phase| mixture.concentration_in_phase(substance_id, *phase))
            .sum::<f64>();
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
    for (substance, order, phases) in &indexed_reaction.orders {
        let concentration = mixture.concentration_of_index_in_phases(*substance, phases);
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
                let available = mixture
                    .concentration_of_index_in_phases(reactant.substance, &reactant.phases)
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
) -> ChemistryResult<ReactionApplication> {
    if !moles_per_bucket.is_finite() || moles_per_bucket < 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "moles to apply must be non-negative and finite".to_string(),
        });
    }
    let deltas = concentration_deltas(indexed_reaction, moles_per_bucket);
    validate_external_reactants(context, reaction, moles_per_bucket)?;
    let mixture_checkpoint = mixture.checkpoint_for_reaction(&deltas);
    let context_checkpoint = checkpoint_context(context, reaction);
    let result = apply_reaction_inner(
        registry,
        mixture,
        context,
        reaction,
        indexed_reaction,
        moles_per_bucket,
        &deltas,
        &mixture_checkpoint,
    );
    if result.is_err() {
        mixture.restore_checkpoint(mixture_checkpoint);
        restore_context(context, context_checkpoint);
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
    _deltas: &[(SubstanceIndex, f64)],
    _mixture_checkpoint: &super::mixture::MixtureCheckpoint,
) -> ChemistryResult<ReactionApplication> {
    let reactants = indexed_reaction
        .reactants
        .iter()
        .map(|term| (term.substance, term.coefficient, term.phases.clone()))
        .collect::<Vec<_>>();
    let products = indexed_reaction
        .products
        .iter()
        .map(|term| {
            (
                term.substance,
                term.coefficient,
                term.phases
                    .first()
                    .copied()
                    .unwrap_or(MixturePhase::Aqueous),
            )
        })
        .collect::<Vec<_>>();
    let max_concentration_delta = mixture.apply_reaction_phase_deltas_by_index(
        registry,
        &reactants,
        &products,
        moles_per_bucket,
    )?;
    apply_external_reactants(context, reaction, moles_per_bucket)?;
    apply_reaction_results(context, reaction, moles_per_bucket);
    mixture.heat(
        registry,
        -reaction.enthalpy_change_kj_per_mol * 1000.0 * moles_per_bucket,
    )?;
    Ok(ReactionApplication {
        max_concentration_delta,
        thermal_changed: reaction.enthalpy_change_kj_per_mol != 0.0 && moles_per_bucket > 0.0,
    })
}

fn concentration_deltas(
    indexed_reaction: &IndexedReaction,
    moles_per_bucket: f64,
) -> Vec<(SubstanceIndex, f64)> {
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
    deltas
}

fn apply_reaction_results(
    context: &mut ReactionContext,
    reaction: &Reaction,
    moles_per_bucket: f64,
) {
    for result in &reaction.reaction_results {
        *context
            .reaction_results
            .entry(result.description.clone())
            .or_insert(0.0) += result.moles_per_reaction * moles_per_bucket;
    }
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

fn validate_external_reactants(
    context: &ReactionContext,
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
        if !next.is_finite() {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "external reactant '{}' would become non-finite: {next}",
                external.description
            )));
        }
        if next < -TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "external reactant '{}' would become negative: {next}",
                external.description
            )));
        }
    }
    Ok(())
}

fn checkpoint_context(context: &ReactionContext, reaction: &Reaction) -> ContextCheckpoint {
    ContextCheckpoint {
        external_reactants: reaction
            .external_reactants
            .iter()
            .map(|external| {
                (
                    external.description.clone(),
                    context
                        .external_reactants
                        .get(&external.description)
                        .copied(),
                )
            })
            .collect(),
        reaction_results: reaction
            .reaction_results
            .iter()
            .map(|result| {
                (
                    result.description.clone(),
                    context.reaction_results.get(&result.description).copied(),
                )
            })
            .collect(),
    }
}

fn restore_context(context: &mut ReactionContext, checkpoint: ContextCheckpoint) {
    for (description, previous) in checkpoint.external_reactants {
        match previous {
            Some(value) => {
                context.external_reactants.insert(description, value);
            }
            None => {
                context.external_reactants.remove(&description);
            }
        }
    }
    for (description, previous) in checkpoint.reaction_results {
        match previous {
            Some(value) => {
                context.reaction_results.insert(description, value);
            }
            None => {
                context.reaction_results.remove(&description);
            }
        }
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::registry::ChemistryRegistryBuilder;
    use crate::chemistry::substance::{LiquidPhasePreference, Substance, SubstancePhaseProperties};

    fn simple_registry(reaction: Reaction) -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:a",
                0,
                10.0,
                10_000.0,
                500.0,
                100.0,
                20_000.0,
            ))
            .substance(Substance::new(
                "destroy:b",
                0,
                10.0,
                10_000.0,
                500.0,
                100.0,
                20_000.0,
            ))
            .substance(Substance::new(
                "destroy:solvent",
                0,
                20.0,
                20_000.0,
                500.0,
                100.0,
                20_000.0,
            ))
            .reaction(reaction)
            .build()
            .unwrap()
    }

    fn precipitate_registry(reaction: Reaction) -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(
                Substance::new(
                    "destroy:solid_reactant",
                    0,
                    10.0,
                    10_000.0,
                    500.0,
                    100.0,
                    20_000.0,
                )
                .with_phase_properties(SubstancePhaseProperties {
                    preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                    aqueous_solubility_mol_per_bucket: Some(0.1),
                    organic_solubility_mol_per_bucket: Some(0.0),
                    can_precipitate: true,
                }),
            )
            .substance(Substance::new(
                "destroy:product",
                0,
                10.0,
                10_000.0,
                500.0,
                100.0,
                20_000.0,
            ))
            .reaction(reaction)
            .build()
            .unwrap()
    }

    #[test]
    fn failed_reaction_application_restores_mixture_and_context() {
        let registry = simple_registry(
            Reaction::builder("destroy:a_to_b")
                .reactant("destroy:a", 1, 1)
                .product("destroy:b", 1)
                .reaction_result("result:test", 1.0)
                .pre_exponential_factor(1.0)
                .activation_energy_kj_per_mol(0.0)
                .build(),
        );
        let reaction_index = registry.reaction_index(&"destroy:a_to_b".into()).unwrap();
        let indexed_reaction = registry.indexed_reaction(reaction_index).unwrap();
        let mut invalid_reaction = registry.reaction_by_index(reaction_index).unwrap().clone();
        invalid_reaction.enthalpy_change_kj_per_mol = f64::INFINITY;
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture.add_substance(&registry, "destroy:a", 1.0).unwrap();
        let mut context = ReactionContext::default();

        let error = apply_reaction(
            &registry,
            &mut mixture,
            &mut context,
            &invalid_reaction,
            indexed_reaction,
            0.25,
        )
        .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidMixtureState(_)));
        assert_eq!(mixture.concentration_of(&"destroy:a".into()), 1.0);
        assert_eq!(mixture.concentration_of(&"destroy:b".into()), 0.0);
        assert_eq!(mixture.temperature_kelvin(), 298.0);
        assert!(context.reaction_results.is_empty());
    }

    #[test]
    fn sub_threshold_concentration_change_does_not_keep_equilibrium_loop_running() {
        let registry = simple_registry(
            Reaction::builder("destroy:a_to_b")
                .reactant("destroy:a", 1, 1)
                .product("destroy:b", 1)
                .pre_exponential_factor(1.0e12)
                .activation_energy_kj_per_mol(0.0)
                .build(),
        );
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(
                &registry,
                "destroy:a",
                EQUILIBRIUM_EPSILON_MOL_PER_BUCKET / 10.0,
            )
            .unwrap();

        let changed = react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(!changed);
    }

    #[test]
    fn sub_threshold_reaction_with_heat_counts_as_changed() {
        let registry = simple_registry(
            Reaction::builder("destroy:a_to_b")
                .reactant("destroy:a", 1, 1)
                .product("destroy:b", 1)
                .pre_exponential_factor(1.0e12)
                .activation_energy_kj_per_mol(0.0)
                .enthalpy_change_kj_per_mol(-10.0)
                .build(),
        );
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(
                &registry,
                "destroy:a",
                EQUILIBRIUM_EPSILON_MOL_PER_BUCKET / 10.0,
            )
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:solvent", 1.0)
            .unwrap();

        let changed = react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(changed);
        assert!(mixture.temperature_kelvin() > 298.0);
    }

    #[test]
    fn default_reaction_access_does_not_consume_solid_precipitate() {
        let registry = precipitate_registry(
            Reaction::builder("destroy:solid_to_product")
                .reactant("destroy:solid_reactant", 1, 1)
                .product("destroy:product", 1)
                .pre_exponential_factor(1.0e12)
                .activation_energy_kj_per_mol(0.0)
                .build(),
        );
        let reactant: crate::chemistry::substance::SubstanceId = "destroy:solid_reactant".into();
        let product: crate::chemistry::substance::SubstanceId = "destroy:product".into();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, reactant.clone(), 1.0)
            .unwrap();

        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!((mixture.concentration_of(&reactant) - 0.9).abs() < 1.0e-9);
        assert_eq!(
            mixture.concentration_in_phase(&reactant, MixturePhase::Aqueous),
            0.1
        );
        assert_eq!(
            mixture.concentration_in_phase(&reactant, MixturePhase::Solid),
            0.8
        );
        assert!((mixture.concentration_of(&product) - 0.1).abs() < 1.0e-9);
    }
}
