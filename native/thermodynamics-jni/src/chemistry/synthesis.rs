use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::dynamic::DynamicChemistryRegistry;
use super::error::{ChemistryError, ChemistryResult};
use super::frowns::write_frowns;
use super::molecule::MolecularStructure;
use super::reaction::{Reaction, ReactionId, StoichiometricTerm};
use super::substance::SubstanceId;

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisPlanner {
    pub max_steps: usize,
    pub max_routes: usize,
    pub allowed_reaction_prefixes: BTreeSet<String>,
    pub safety_policy: SynthesisSafetyPolicy,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SynthesisSafetyPolicy {
    denied_substance_ids: BTreeMap<SubstanceId, String>,
    denied_canonical_structures: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisRoute {
    pub target: SubstanceId,
    pub steps: Vec<SynthesisStep>,
    pub estimated_yield: f64,
    pub score: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisStep {
    pub reaction_id: ReactionId,
    pub reactants: Vec<SubstanceId>,
    pub products: Vec<SubstanceId>,
    pub product_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
struct FrontierState {
    known: BTreeSet<SubstanceId>,
    steps: Vec<SynthesisStep>,
    estimated_yield: f64,
}

impl Default for SynthesisPlanner {
    fn default() -> Self {
        Self {
            max_steps: 6,
            max_routes: 8,
            allowed_reaction_prefixes: BTreeSet::new(),
            safety_policy: SynthesisSafetyPolicy::default(),
        }
    }
}

impl SynthesisSafetyPolicy {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn deny_substance(
        mut self,
        substance_id: impl Into<SubstanceId>,
        reason: impl Into<String>,
    ) -> ChemistryResult<Self> {
        let reason = reason.into();
        if reason.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-policy>".to_string(),
                reason: "denied substance reason must not be empty".to_string(),
            });
        }
        self.denied_substance_ids
            .insert(substance_id.into(), reason);
        Ok(self)
    }

    pub fn deny_structure(
        mut self,
        structure: MolecularStructure,
        reason: impl Into<String>,
    ) -> ChemistryResult<Self> {
        let reason = reason.into();
        if reason.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-policy>".to_string(),
                reason: "denied structure reason must not be empty".to_string(),
            });
        }
        self.denied_canonical_structures
            .insert(write_frowns(&structure)?, reason);
        Ok(self)
    }

    fn check_substance(
        &self,
        registry: &DynamicChemistryRegistry,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<()> {
        if let Some(reason) = self.denied_substance_ids.get(substance_id) {
            return Err(denied_synthesis(substance_id, reason));
        }
        let substance = registry.substance(substance_id)?;
        if let Some(structure) = substance.molecular_structure.as_ref() {
            let canonical = write_frowns(structure)?;
            if let Some(reason) = self.denied_canonical_structures.get(&canonical) {
                return Err(denied_synthesis(substance_id, reason));
            }
        }
        Ok(())
    }

    fn allows_substance(
        &self,
        registry: &DynamicChemistryRegistry,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<bool> {
        match self.check_substance(registry, substance_id) {
            Ok(()) => Ok(true),
            Err(ChemistryError::InvalidReaction { reaction_id, .. })
                if reaction_id == "<synthesis-policy>" =>
            {
                Ok(false)
            }
            Err(error) => Err(error),
        }
    }
}

impl SynthesisPlanner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_max_steps(mut self, max_steps: usize) -> Self {
        self.max_steps = max_steps;
        self
    }

    pub fn with_max_routes(mut self, max_routes: usize) -> Self {
        self.max_routes = max_routes;
        self
    }

    pub fn allow_reaction_prefix(mut self, prefix: impl Into<String>) -> Self {
        self.allowed_reaction_prefixes.insert(prefix.into());
        self
    }

    pub fn with_safety_policy(mut self, safety_policy: SynthesisSafetyPolicy) -> Self {
        self.safety_policy = safety_policy;
        self
    }

    pub fn find_routes(
        &self,
        registry: &mut DynamicChemistryRegistry,
        starting_substances: impl IntoIterator<Item = SubstanceId>,
        target_structure: MolecularStructure,
    ) -> ChemistryResult<Vec<SynthesisRoute>> {
        if self.max_steps == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-planner>".to_string(),
                reason: "max_steps must be greater than zero".to_string(),
            });
        }
        if self.max_routes == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-planner>".to_string(),
                reason: "max_routes must be greater than zero".to_string(),
            });
        }
        let target_canonical = write_frowns(&target_structure)?;
        if let Some(reason) = self
            .safety_policy
            .denied_canonical_structures
            .get(&target_canonical)
        {
            return Err(denied_synthesis_id("<target-structure>", reason));
        }
        let target = registry.resolve_structure(target_structure)?;
        self.safety_policy.check_substance(registry, &target)?;
        let mut initial_known = BTreeSet::new();
        for substance_id in starting_substances {
            registry.substance(&substance_id)?;
            self.safety_policy
                .check_substance(registry, &substance_id)?;
            initial_known.insert(substance_id);
        }
        if initial_known.contains(&target) {
            return Ok(vec![SynthesisRoute {
                target,
                steps: Vec::new(),
                estimated_yield: 1.0,
                score: 0.0,
            }]);
        }

        let mut routes = Vec::new();
        let mut queue = VecDeque::from([FrontierState {
            known: initial_known.clone(),
            steps: Vec::new(),
            estimated_yield: 1.0,
        }]);
        let mut seen_known_sets = BTreeSet::from([known_set_key(&initial_known)]);

        while let Some(state) = queue.pop_front() {
            if state.steps.len() >= self.max_steps {
                continue;
            }
            registry.generate_reactions_for_substances(state.known.iter().cloned(), 1)?;
            let candidates = registry
                .reactions()
                .filter(|reaction| self.reaction_allowed(reaction))
                .filter_map(|reaction| {
                    synthesis_step_from_available_reaction(reaction, &state.known)
                })
                .collect::<Vec<_>>();
            for step in candidates {
                let mut next_known = state.known.clone();
                let mut added_new_product = false;
                if !step_products_allowed(registry, &self.safety_policy, &step)? {
                    continue;
                }
                for product in &step.products {
                    if next_known.insert(product.clone()) {
                        added_new_product = true;
                    }
                }
                if !added_new_product {
                    continue;
                }
                let mut next_steps = state.steps.clone();
                next_steps.push(step.clone());
                let estimated_yield = state.estimated_yield * step.product_fraction;
                if next_known.contains(&target) {
                    routes.push(SynthesisRoute {
                        target: target.clone(),
                        score: route_score(&next_steps, estimated_yield),
                        steps: next_steps.clone(),
                        estimated_yield,
                    });
                    routes.sort_by(|left, right| {
                        left.score
                            .total_cmp(&right.score)
                            .then_with(|| right.estimated_yield.total_cmp(&left.estimated_yield))
                    });
                    routes.truncate(self.max_routes);
                    continue;
                }
                let key = known_set_key(&next_known);
                if seen_known_sets.insert(key) {
                    queue.push_back(FrontierState {
                        known: next_known,
                        steps: next_steps,
                        estimated_yield,
                    });
                }
            }
        }
        Ok(routes)
    }

    fn reaction_allowed(&self, reaction: &Reaction) -> bool {
        self.allowed_reaction_prefixes.is_empty()
            || self
                .allowed_reaction_prefixes
                .iter()
                .any(|prefix| reaction.id.as_str().starts_with(prefix))
    }
}

fn step_products_allowed(
    registry: &DynamicChemistryRegistry,
    safety_policy: &SynthesisSafetyPolicy,
    step: &SynthesisStep,
) -> ChemistryResult<bool> {
    for product in &step.products {
        if !safety_policy.allows_substance(registry, product)? {
            return Ok(false);
        }
    }
    Ok(true)
}

fn synthesis_step_from_available_reaction(
    reaction: &Reaction,
    known: &BTreeSet<SubstanceId>,
) -> Option<SynthesisStep> {
    let reactants = reaction
        .reactants
        .iter()
        .map(|term| term.substance_id.clone())
        .collect::<Vec<_>>();
    if !reactants.iter().all(|substance| known.contains(substance)) {
        return None;
    }
    if !reaction.products.is_empty() {
        return Some(SynthesisStep {
            reaction_id: reaction.id.clone(),
            reactants,
            products: term_ids(&reaction.products),
            product_fraction: 1.0,
        });
    }
    if let Some(distribution) = &reaction.product_distribution {
        let mut products = Vec::new();
        let mut fraction_by_product = BTreeMap::new();
        for variant in &distribution.variants {
            for product in term_ids(&variant.products) {
                products.push(product.clone());
                *fraction_by_product.entry(product).or_insert(0.0) += variant.fraction;
            }
        }
        products.sort();
        products.dedup();
        let product_fraction = products
            .iter()
            .filter_map(|product| fraction_by_product.get(product))
            .copied()
            .fold(0.0_f64, f64::max);
        return Some(SynthesisStep {
            reaction_id: reaction.id.clone(),
            reactants,
            products,
            product_fraction,
        });
    }
    if !reaction.channels.is_empty() {
        let mut products = Vec::new();
        for channel in &reaction.channels {
            products.extend(term_ids(&channel.products));
        }
        products.sort();
        products.dedup();
        let product_fraction = if products.is_empty() {
            0.0
        } else {
            1.0 / products.len() as f64
        };
        return Some(SynthesisStep {
            reaction_id: reaction.id.clone(),
            reactants,
            products,
            product_fraction,
        });
    }
    None
}

fn term_ids(terms: &[StoichiometricTerm]) -> Vec<SubstanceId> {
    terms.iter().map(|term| term.substance_id.clone()).collect()
}

fn route_score(steps: &[SynthesisStep], estimated_yield: f64) -> f64 {
    steps.len() as f64 + (1.0 - estimated_yield.clamp(0.0, 1.0))
}

fn known_set_key(known: &BTreeSet<SubstanceId>) -> Vec<SubstanceId> {
    known.iter().cloned().collect()
}

fn denied_synthesis(substance_id: &SubstanceId, reason: &str) -> ChemistryError {
    denied_synthesis_id(substance_id.as_str(), reason)
}

fn denied_synthesis_id(substance_id: &str, reason: &str) -> ChemistryError {
    ChemistryError::InvalidReaction {
        reaction_id: "<synthesis-policy>".to_string(),
        reason: format!("synthesis route for '{substance_id}' is denied: {reason}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::dynamic::DynamicChemistryRegistry;
    use crate::chemistry::frowns::parse_frowns;

    #[test]
    fn planner_finds_local_route_without_generating_the_whole_space() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let routes = SynthesisPlanner::new()
            .with_max_steps(2)
            .find_routes(
                &mut registry,
                [
                    SubstanceId::from("destroy:chloroethane"),
                    SubstanceId::from("destroy:hydroxide"),
                ],
                parse_frowns("CCO").unwrap(),
            )
            .unwrap();
        assert!(!routes.is_empty());
        assert!(routes[0].steps.iter().any(|step| step
            .reaction_id
            .as_str()
            .starts_with("halide_hydroxide_substitution")));
    }

    #[test]
    fn planner_finds_five_step_chain_across_generated_products() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let routes = SynthesisPlanner::new()
            .with_max_steps(5)
            .with_max_routes(4)
            .allow_reaction_prefix("halide_cyanide_substitution/")
            .allow_reaction_prefix("nitrile_hydrolysis/")
            .allow_reaction_prefix("amide_hydrolysis/")
            .allow_reaction_prefix("acyl_chloride_formation/")
            .allow_reaction_prefix("acyl_chloride_esterification/")
            .find_routes(
                &mut registry,
                [
                    SubstanceId::from("destroy:chloroethane"),
                    SubstanceId::from("destroy:cyanide"),
                    SubstanceId::from("destroy:water"),
                    SubstanceId::from("destroy:phosgene"),
                    SubstanceId::from("destroy:ethanol"),
                ],
                parse_frowns("CCC(=O)OCC").unwrap(),
            )
            .unwrap();

        assert!(!routes.is_empty());
        let route = &routes[0];
        assert_eq!(route.steps.len(), 5);
        let ids = route
            .steps
            .iter()
            .map(|step| step.reaction_id.as_str())
            .collect::<Vec<_>>();
        assert!(ids
            .iter()
            .any(|id| id.starts_with("halide_cyanide_substitution/")));
        assert!(ids.iter().any(|id| id.starts_with("nitrile_hydrolysis/")));
        assert!(ids.iter().any(|id| id.starts_with("amide_hydrolysis/")));
        assert!(ids
            .iter()
            .any(|id| id.starts_with("acyl_chloride_formation/")));
        assert!(ids
            .iter()
            .any(|id| id.starts_with("acyl_chloride_esterification/")));
    }

    #[test]
    fn safety_policy_denies_target_routes_before_search() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let policy = SynthesisSafetyPolicy::new()
            .deny_structure(
                parse_frowns("CCO").unwrap(),
                "test policy denies this target",
            )
            .unwrap();
        let error = SynthesisPlanner::new()
            .with_safety_policy(policy)
            .find_routes(
                &mut registry,
                [
                    SubstanceId::from("destroy:chloroethane"),
                    SubstanceId::from("destroy:hydroxide"),
                ],
                parse_frowns("CCO").unwrap(),
            )
            .unwrap_err();
        assert!(matches!(
            error,
            ChemistryError::InvalidReaction { reaction_id, reason }
                if reaction_id == "<synthesis-policy>" && reason.contains("denied")
        ));
    }

    #[test]
    fn safety_policy_filters_denied_intermediates() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let policy = SynthesisSafetyPolicy::new()
            .deny_structure(
                parse_frowns("CCC(=O)N").unwrap(),
                "test policy denies this intermediate",
            )
            .unwrap();
        let routes = SynthesisPlanner::new()
            .with_safety_policy(policy)
            .with_max_steps(5)
            .with_max_routes(4)
            .allow_reaction_prefix("halide_cyanide_substitution/")
            .allow_reaction_prefix("nitrile_hydrolysis/")
            .allow_reaction_prefix("amide_hydrolysis/")
            .allow_reaction_prefix("acyl_chloride_formation/")
            .allow_reaction_prefix("acyl_chloride_esterification/")
            .find_routes(
                &mut registry,
                [
                    SubstanceId::from("destroy:chloroethane"),
                    SubstanceId::from("destroy:cyanide"),
                    SubstanceId::from("destroy:water"),
                    SubstanceId::from("destroy:phosgene"),
                    SubstanceId::from("destroy:ethanol"),
                ],
                parse_frowns("CCC(=O)OCC").unwrap(),
            )
            .unwrap();
        assert!(routes.is_empty());
    }
}
