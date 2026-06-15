use std::collections::{BTreeMap, BTreeSet, VecDeque};

use super::dynamic::DynamicChemistryRegistry;
use super::error::{ChemistryError, ChemistryResult};
use super::frowns::write_frowns;
use super::molecule::MolecularStructure;
use super::reaction::{ExternalRequirement, Reaction, ReactionId, StoichiometricTerm};
use super::substance::SubstanceId;

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisPlanner {
    pub max_steps: usize,
    pub max_routes: usize,
    pub allowed_reaction_prefixes: BTreeSet<String>,
    pub safety_policy: SynthesisSafetyPolicy,
    pub include_routes_with_missing_inputs: bool,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SynthesisSafetyPolicy {
    denied_substance_ids: BTreeMap<SubstanceId, String>,
    denied_canonical_structures: BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisRequest {
    pub target: SynthesisTarget,
    pub available_substances: BTreeSet<SubstanceId>,
    pub max_steps: Option<usize>,
    pub max_routes: Option<usize>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SynthesisTarget {
    Substance(SubstanceId),
    Structure(MolecularStructure),
    Frowns(String),
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisRoute {
    pub target: SubstanceId,
    pub steps: Vec<SynthesisStep>,
    pub estimated_yield: f64,
    pub score: f64,
    pub required_additions: Vec<SynthesisRequirement>,
    pub condition_hints: Vec<SynthesisConditionHint>,
    pub explanation: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisPlanReport {
    pub target: SubstanceId,
    pub status: SynthesisPlanStatus,
    pub routes: Vec<SynthesisRoute>,
    pub required_additions: Vec<SynthesisRequirement>,
    pub condition_hints: Vec<SynthesisConditionHint>,
    pub unsupported_reason: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisReachabilityRequest {
    pub targets: Vec<SubstanceId>,
    pub available_substances: BTreeSet<SubstanceId>,
    pub generation_iterations: usize,
    pub max_steps: usize,
    pub max_routes_per_target: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisReachabilityReport {
    pub starting_substance_count: usize,
    pub target_count: usize,
    pub reachable_targets: Vec<SubstanceId>,
    pub unreachable_targets: Vec<SubstanceId>,
    pub target_reports: Vec<SynthesisTargetReachability>,
    pub generation_report: super::dynamic::DynamicGenerationReport,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisTargetReachability {
    pub target: SubstanceId,
    pub status: SynthesisPlanStatus,
    pub has_known_producer: bool,
    pub producer_reactions: Vec<ReactionId>,
    pub missing_reactants: Vec<SynthesisRequirement>,
    pub routes: Vec<SynthesisRoute>,
    pub unsupported_reason: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SynthesisPlanStatus {
    TargetAlreadyAvailable,
    RoutesFound,
    RequiresAdditionalInputs,
    SearchLimitReached,
    UnsupportedByCurrentModel,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisStep {
    pub reaction_id: ReactionId,
    pub reactants: Vec<SubstanceId>,
    pub products: Vec<SubstanceId>,
    pub missing_reactants: Vec<SubstanceId>,
    pub external_reactants: Vec<SynthesisExternalHint>,
    pub external_catalysts: Vec<SynthesisExternalHint>,
    pub condition_hints: Vec<SynthesisConditionHint>,
    pub requires_uv: bool,
    pub phase_hints: Vec<String>,
    pub product_fraction: f64,
    pub step_cost: f64,
    pub safety_penalty: f64,
    pub condition_penalty: f64,
    pub selectivity_penalty: f64,
    pub purification_penalty: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisRequirement {
    pub substance_id: SubstanceId,
    pub reason: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisExternalHint {
    pub description: String,
    pub moles_per_reaction: f64,
    pub unchecked_mass_reason: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SynthesisConditionHint {
    pub description: String,
    pub penalty: f64,
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
            include_routes_with_missing_inputs: true,
        }
    }
}

impl SynthesisRequest {
    pub fn for_substance(target: impl Into<SubstanceId>) -> Self {
        Self {
            target: SynthesisTarget::Substance(target.into()),
            available_substances: BTreeSet::new(),
            max_steps: None,
            max_routes: None,
        }
    }

    pub fn for_structure(target: MolecularStructure) -> Self {
        Self {
            target: SynthesisTarget::Structure(target),
            available_substances: BTreeSet::new(),
            max_steps: None,
            max_routes: None,
        }
    }

    pub fn for_frowns(target: impl Into<String>) -> Self {
        Self {
            target: SynthesisTarget::Frowns(target.into()),
            available_substances: BTreeSet::new(),
            max_steps: None,
            max_routes: None,
        }
    }

    pub fn with_available_substance(mut self, substance_id: impl Into<SubstanceId>) -> Self {
        self.available_substances.insert(substance_id.into());
        self
    }

    pub fn with_available_substances<I>(mut self, substance_ids: I) -> Self
    where
        I: IntoIterator<Item = SubstanceId>,
    {
        self.available_substances.extend(substance_ids);
        self
    }

    pub fn with_max_steps(mut self, max_steps: usize) -> Self {
        self.max_steps = Some(max_steps);
        self
    }

    pub fn with_max_routes(mut self, max_routes: usize) -> Self {
        self.max_routes = Some(max_routes);
        self
    }
}

impl SynthesisReachabilityRequest {
    pub fn new(targets: impl IntoIterator<Item = SubstanceId>) -> Self {
        Self {
            targets: targets.into_iter().collect(),
            available_substances: BTreeSet::new(),
            generation_iterations: 3,
            max_steps: 6,
            max_routes_per_target: 3,
        }
    }

    pub fn with_available_substance(mut self, substance_id: impl Into<SubstanceId>) -> Self {
        self.available_substances.insert(substance_id.into());
        self
    }

    pub fn with_available_substances<I>(mut self, substance_ids: I) -> Self
    where
        I: IntoIterator<Item = SubstanceId>,
    {
        self.available_substances.extend(substance_ids);
        self
    }

    pub fn with_generation_iterations(mut self, generation_iterations: usize) -> Self {
        self.generation_iterations = generation_iterations;
        self
    }

    pub fn with_max_steps(mut self, max_steps: usize) -> Self {
        self.max_steps = max_steps;
        self
    }

    pub fn with_max_routes_per_target(mut self, max_routes_per_target: usize) -> Self {
        self.max_routes_per_target = max_routes_per_target;
        self
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

    pub fn include_routes_with_missing_inputs(mut self, value: bool) -> Self {
        self.include_routes_with_missing_inputs = value;
        self
    }

    pub fn plan_routes(
        &self,
        registry: &DynamicChemistryRegistry,
        request: SynthesisRequest,
    ) -> ChemistryResult<Vec<SynthesisRoute>> {
        let mut working_registry = registry.clone();
        let target = resolve_target(&mut working_registry, request.target)?;
        self.find_routes_in_working_registry(
            &mut working_registry,
            request.available_substances,
            target,
            request.max_steps.unwrap_or(self.max_steps),
            request.max_routes.unwrap_or(self.max_routes),
        )
    }

    pub fn plan_report(
        &self,
        registry: &DynamicChemistryRegistry,
        request: SynthesisRequest,
    ) -> ChemistryResult<SynthesisPlanReport> {
        let mut working_registry = registry.clone();
        let target = resolve_target(&mut working_registry, request.target.clone())?;
        let starting_substances = request.available_substances.clone();
        let requested_max_steps = request.max_steps;
        let max_steps = request.max_steps.unwrap_or(self.max_steps);
        let max_routes = request.max_routes.unwrap_or(self.max_routes);
        let routes = self.find_routes_in_working_registry(
            &mut working_registry,
            request.available_substances,
            target.clone(),
            max_steps,
            max_routes,
        )?;
        let mut required_additions = route_required_additions(&routes);
        let condition_hints = routes_condition_hints(&routes);
        let (status, unsupported_reason) = if starting_substances.contains(&target) {
            (SynthesisPlanStatus::TargetAlreadyAvailable, None)
        } else if routes.is_empty() {
            if starting_substances.is_empty() {
                working_registry.generate_reactions(1)?;
            }
            let direct_requirements =
                direct_target_requirements(&working_registry, self, &starting_substances, &target)?;
            if direct_requirements.is_empty() {
                if requested_max_steps.is_some_and(|requested| requested < self.max_steps) {
                    (
                        SynthesisPlanStatus::SearchLimitReached,
                        Some(format!(
                            "no route to '{}' was found within {max_steps} step(s)",
                            target.as_str()
                        )),
                    )
                } else {
                    (
                        SynthesisPlanStatus::UnsupportedByCurrentModel,
                        Some(format!(
                            "no known reaction in the current model produces '{}'",
                            target.as_str()
                        )),
                    )
                }
            } else {
                required_additions = direct_requirements;
                (SynthesisPlanStatus::RequiresAdditionalInputs, None)
            }
        } else if required_additions.is_empty() {
            (SynthesisPlanStatus::RoutesFound, None)
        } else {
            (SynthesisPlanStatus::RequiresAdditionalInputs, None)
        };

        Ok(SynthesisPlanReport {
            target,
            status,
            routes,
            required_additions,
            condition_hints,
            unsupported_reason,
        })
    }

    pub fn analyze_reachability(
        &self,
        registry: &DynamicChemistryRegistry,
        request: SynthesisReachabilityRequest,
    ) -> ChemistryResult<SynthesisReachabilityReport> {
        if request.targets.is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-reachability>".to_string(),
                reason: "reachability target list must not be empty".to_string(),
            });
        }
        if request.generation_iterations == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-reachability>".to_string(),
                reason: "generation_iterations must be greater than zero".to_string(),
            });
        }
        if request.max_steps == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-reachability>".to_string(),
                reason: "max_steps must be greater than zero".to_string(),
            });
        }
        if request.max_routes_per_target == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-reachability>".to_string(),
                reason: "max_routes_per_target must be greater than zero".to_string(),
            });
        }

        let mut working_registry = registry.clone();
        let mut targets = request.targets;
        targets.sort();
        targets.dedup();
        for target in &targets {
            working_registry.substance(target)?;
            self.safety_policy
                .check_substance(&working_registry, target)?;
        }
        for substance in &request.available_substances {
            working_registry.substance(substance)?;
            self.safety_policy
                .check_substance(&working_registry, substance)?;
        }

        let generation_report = if request.available_substances.is_empty() {
            working_registry.generate_reactions(request.generation_iterations)?
        } else {
            working_registry.generate_reactions_for_substances(
                request.available_substances.iter().cloned(),
                request.generation_iterations,
            )?
        };

        let mut target_reports = Vec::new();
        let mut reachable_targets = Vec::new();
        let mut unreachable_targets = Vec::new();
        for target in targets {
            let producer_reactions = producer_reactions_for_target(&working_registry, &target);
            let plan = self
                .clone()
                .with_max_steps(request.max_steps)
                .with_max_routes(request.max_routes_per_target)
                .plan_report(
                    &working_registry,
                    SynthesisRequest::for_substance(target.clone())
                        .with_available_substances(
                            request
                                .available_substances
                                .iter()
                                .cloned()
                                .collect::<Vec<_>>(),
                        )
                        .with_max_steps(request.max_steps)
                        .with_max_routes(request.max_routes_per_target),
                )?;
            let reachable = matches!(
                plan.status,
                SynthesisPlanStatus::TargetAlreadyAvailable | SynthesisPlanStatus::RoutesFound
            );
            if reachable {
                reachable_targets.push(target.clone());
            } else {
                unreachable_targets.push(target.clone());
            }
            target_reports.push(SynthesisTargetReachability {
                target,
                status: plan.status,
                has_known_producer: !producer_reactions.is_empty(),
                producer_reactions,
                missing_reactants: plan.required_additions,
                routes: plan.routes,
                unsupported_reason: plan.unsupported_reason,
            });
        }

        Ok(SynthesisReachabilityReport {
            starting_substance_count: request.available_substances.len(),
            target_count: target_reports.len(),
            reachable_targets,
            unreachable_targets,
            target_reports,
            generation_report,
        })
    }

    pub fn find_routes(
        &self,
        registry: &mut DynamicChemistryRegistry,
        starting_substances: impl IntoIterator<Item = SubstanceId>,
        target_structure: MolecularStructure,
    ) -> ChemistryResult<Vec<SynthesisRoute>> {
        let request = SynthesisRequest::for_structure(target_structure)
            .with_available_substances(starting_substances.into_iter().collect::<Vec<_>>());
        self.plan_routes(registry, request)
    }

    fn find_routes_in_working_registry(
        &self,
        registry: &mut DynamicChemistryRegistry,
        starting_substances: impl IntoIterator<Item = SubstanceId>,
        target: SubstanceId,
        max_steps: usize,
        max_routes: usize,
    ) -> ChemistryResult<Vec<SynthesisRoute>> {
        if max_steps == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-planner>".to_string(),
                reason: "max_steps must be greater than zero".to_string(),
            });
        }
        if max_routes == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<synthesis-planner>".to_string(),
                reason: "max_routes must be greater than zero".to_string(),
            });
        }
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
                required_additions: Vec::new(),
                condition_hints: Vec::new(),
                explanation: "target is already available".to_string(),
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
            if state.steps.len() >= max_steps || state.known.is_empty() {
                continue;
            }
            registry.generate_reactions_for_substances(state.known.iter().cloned(), 1)?;
            let candidates = registry
                .reaction_candidates_for_substances(state.known.iter())
                .into_iter()
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
                    routes.push(build_route(
                        target.clone(),
                        estimated_yield,
                        next_steps.clone(),
                        &initial_known,
                        "route can be run from the available substances".to_string(),
                    ));
                    routes.sort_by(compare_routes);
                    routes.truncate(max_routes);
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

        if self.include_routes_with_missing_inputs && routes.len() < max_routes {
            let mut visiting = BTreeSet::new();
            let mut backward = backward_routes_for_target(
                registry,
                &self.safety_policy,
                &initial_known,
                &target,
                max_steps,
                max_routes - routes.len(),
                &mut visiting,
                self,
            )?;
            routes.append(&mut backward);
            routes.sort_by(compare_routes);
            routes.dedup_by(|left, right| route_key(left) == route_key(right));
            routes.truncate(max_routes);
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

fn resolve_target(
    registry: &mut DynamicChemistryRegistry,
    target: SynthesisTarget,
) -> ChemistryResult<SubstanceId> {
    match target {
        SynthesisTarget::Substance(id) => {
            registry.substance(&id)?;
            Ok(id)
        }
        SynthesisTarget::Structure(structure) => registry.resolve_structure(structure),
        SynthesisTarget::Frowns(code) => registry.resolve_frowns(&code),
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
    let reactants = term_ids(&reaction.reactants);
    if !reactants.iter().all(|substance| known.contains(substance)) {
        return None;
    }
    synthesis_step_from_reaction(reaction, known)
}

fn synthesis_step_from_reaction(
    reaction: &Reaction,
    known: &BTreeSet<SubstanceId>,
) -> Option<SynthesisStep> {
    let reactants = term_ids(&reaction.reactants);
    let missing_reactants = reactants
        .iter()
        .filter(|substance| !known.contains(*substance))
        .cloned()
        .collect::<Vec<_>>();
    let products = reaction_products(reaction);
    if products.is_empty() {
        return None;
    }
    let product_fraction = reaction_product_fraction(reaction, &products);
    let mut step = SynthesisStep {
        reaction_id: reaction.id.clone(),
        reactants,
        products,
        missing_reactants,
        external_reactants: external_hints(&reaction.external_reactants),
        external_catalysts: external_hints(&reaction.external_catalysts),
        condition_hints: condition_hints(reaction),
        requires_uv: reaction.requires_uv,
        phase_hints: phase_hints(reaction),
        product_fraction,
        step_cost: 0.0,
        safety_penalty: safety_penalty(reaction),
        condition_penalty: condition_penalty(reaction),
        selectivity_penalty: selectivity_penalty(reaction, product_fraction),
        purification_penalty: purification_penalty(reaction, product_fraction),
    };
    step.step_cost = synthesis_step_cost(&step);
    Some(step)
}

fn backward_routes_for_target(
    registry: &DynamicChemistryRegistry,
    safety_policy: &SynthesisSafetyPolicy,
    available: &BTreeSet<SubstanceId>,
    target: &SubstanceId,
    remaining_steps: usize,
    max_routes: usize,
    visiting: &mut BTreeSet<SubstanceId>,
    planner: &SynthesisPlanner,
) -> ChemistryResult<Vec<SynthesisRoute>> {
    if remaining_steps == 0 || max_routes == 0 || !visiting.insert(target.clone()) {
        return Ok(Vec::new());
    }
    let mut routes = Vec::new();
    for reaction in registry
        .reactions()
        .filter(|reaction| planner.reaction_allowed(reaction))
        .filter(|reaction| reaction_can_produce(reaction, target))
    {
        let Some(step) = synthesis_step_from_reaction(reaction, available) else {
            continue;
        };
        if !step_products_allowed(registry, safety_policy, &step)? {
            continue;
        }
        let mut sub_steps = Vec::new();
        let mut requirements = Vec::new();
        let mut denied_required_reactant = false;
        for reactant in step.missing_reactants.clone() {
            if !safety_policy.allows_substance(registry, &reactant)? {
                denied_required_reactant = true;
                break;
            }
            let nested = backward_routes_for_target(
                registry,
                safety_policy,
                available,
                &reactant,
                remaining_steps.saturating_sub(1),
                1,
                visiting,
                planner,
            )?;
            if let Some(route) = nested.into_iter().next() {
                sub_steps.extend(route.steps);
                requirements.extend(route.required_additions);
            } else {
                requirements.push(SynthesisRequirement {
                    substance_id: reactant.clone(),
                    reason: format!("required reactant for planned step '{}'", step.reaction_id),
                });
            }
        }
        if denied_required_reactant {
            continue;
        }
        sub_steps.push(step);
        if sub_steps.len() > remaining_steps {
            continue;
        }
        let estimated_yield = sub_steps
            .iter()
            .map(|step| step.product_fraction)
            .product::<f64>();
        let mut route = build_route(
            target.clone(),
            estimated_yield,
            sub_steps,
            available,
            if requirements.is_empty() {
                "route can be completed from available substances and generated intermediates"
                    .to_string()
            } else {
                "route is chemically known, but requires additional input substances".to_string()
            },
        );
        route.required_additions.extend(requirements);
        route.required_additions.sort_by(|left, right| {
            left.substance_id
                .cmp(&right.substance_id)
                .then_with(|| left.reason.cmp(&right.reason))
        });
        route
            .required_additions
            .dedup_by(|left, right| left.substance_id == right.substance_id);
        route.score = route_score(
            &route.steps,
            route.estimated_yield,
            &route.required_additions,
        );
        routes.push(route);
        routes.sort_by(compare_routes);
        routes.truncate(max_routes);
    }
    visiting.remove(target);
    Ok(routes)
}

fn reaction_can_produce(reaction: &Reaction, target: &SubstanceId) -> bool {
    reaction_products(reaction)
        .iter()
        .any(|product| product == target)
}

fn producer_reactions_for_target(
    registry: &DynamicChemistryRegistry,
    target: &SubstanceId,
) -> Vec<ReactionId> {
    registry
        .reactions()
        .filter(|reaction| reaction_can_produce(reaction, target))
        .map(|reaction| reaction.id.clone())
        .collect()
}

fn reaction_products(reaction: &Reaction) -> Vec<SubstanceId> {
    let mut products = term_ids(&reaction.products);
    if let Some(distribution) = &reaction.product_distribution {
        for variant in &distribution.variants {
            products.extend(term_ids(&variant.products));
        }
    }
    for channel in &reaction.channels {
        products.extend(term_ids(&channel.products));
    }
    products.sort();
    products.dedup();
    products
}

fn reaction_product_fraction(reaction: &Reaction, products: &[SubstanceId]) -> f64 {
    if !reaction.products.is_empty() {
        return 1.0;
    }
    if let Some(distribution) = &reaction.product_distribution {
        let mut fraction_by_product = BTreeMap::new();
        for variant in &distribution.variants {
            for product in term_ids(&variant.products) {
                *fraction_by_product.entry(product).or_insert(0.0) += variant.fraction;
            }
        }
        return products
            .iter()
            .filter_map(|product| fraction_by_product.get(product))
            .copied()
            .fold(0.0_f64, f64::max);
    }
    if !reaction.channels.is_empty() && !products.is_empty() {
        return 1.0 / products.len() as f64;
    }
    0.0
}

fn term_ids(terms: &[StoichiometricTerm]) -> Vec<SubstanceId> {
    terms.iter().map(|term| term.substance_id.clone()).collect()
}

fn build_route(
    target: SubstanceId,
    estimated_yield: f64,
    steps: Vec<SynthesisStep>,
    known_at_start: &BTreeSet<SubstanceId>,
    explanation: String,
) -> SynthesisRoute {
    let mut required_additions = Vec::new();
    for step in &steps {
        for reactant in &step.missing_reactants {
            if !known_at_start.contains(reactant) {
                required_additions.push(SynthesisRequirement {
                    substance_id: reactant.clone(),
                    reason: format!("required by reaction '{}'", step.reaction_id),
                });
            }
        }
    }
    required_additions.sort_by(|left, right| left.substance_id.cmp(&right.substance_id));
    required_additions.dedup_by(|left, right| left.substance_id == right.substance_id);
    let condition_hints = route_condition_hints(&steps);
    let score = route_score(&steps, estimated_yield, &required_additions);
    SynthesisRoute {
        target,
        steps,
        estimated_yield,
        score,
        required_additions,
        condition_hints,
        explanation,
    }
}

fn route_required_additions(routes: &[SynthesisRoute]) -> Vec<SynthesisRequirement> {
    let mut required = routes
        .iter()
        .flat_map(|route| route.required_additions.iter().cloned())
        .collect::<Vec<_>>();
    deduplicate_requirements(&mut required);
    required
}

fn routes_condition_hints(routes: &[SynthesisRoute]) -> Vec<SynthesisConditionHint> {
    let mut by_description = BTreeMap::<String, f64>::new();
    for route in routes {
        for hint in &route.condition_hints {
            *by_description
                .entry(hint.description.clone())
                .or_insert(0.0) += hint.penalty;
        }
    }
    by_description
        .into_iter()
        .map(|(description, penalty)| SynthesisConditionHint {
            description,
            penalty,
        })
        .collect()
}

fn direct_target_requirements(
    registry: &DynamicChemistryRegistry,
    planner: &SynthesisPlanner,
    starting_substances: &BTreeSet<SubstanceId>,
    target: &SubstanceId,
) -> ChemistryResult<Vec<SynthesisRequirement>> {
    let mut requirements = Vec::new();
    for reaction in registry
        .reactions()
        .filter(|reaction| planner.reaction_allowed(reaction))
        .filter(|reaction| reaction_can_produce(reaction, target))
    {
        for reactant in term_ids(&reaction.reactants) {
            if starting_substances.contains(&reactant) {
                continue;
            }
            if !planner
                .safety_policy
                .allows_substance(registry, &reactant)?
            {
                continue;
            }
            requirements.push(SynthesisRequirement {
                substance_id: reactant,
                reason: format!("required by reaction '{}'", reaction.id),
            });
        }
    }
    deduplicate_requirements(&mut requirements);
    Ok(requirements)
}

fn deduplicate_requirements(requirements: &mut Vec<SynthesisRequirement>) {
    requirements.sort_by(|left, right| {
        left.substance_id
            .cmp(&right.substance_id)
            .then_with(|| left.reason.cmp(&right.reason))
    });
    requirements.dedup_by(|left, right| left.substance_id == right.substance_id);
}

fn route_score(
    steps: &[SynthesisStep],
    estimated_yield: f64,
    required_additions: &[SynthesisRequirement],
) -> f64 {
    steps.iter().map(|step| step.step_cost).sum::<f64>()
        + required_additions.len() as f64 * 1.5
        + (1.0 - estimated_yield.clamp(0.0, 1.0)) * 4.0
}

fn synthesis_step_cost(step: &SynthesisStep) -> f64 {
    1.0 + step.safety_penalty
        + step.condition_penalty
        + step.selectivity_penalty
        + step.purification_penalty
        + step.external_reactants.len() as f64 * 0.4
        + step.external_catalysts.len() as f64 * 0.2
        + if step.requires_uv { 0.5 } else { 0.0 }
}

fn safety_penalty(reaction: &Reaction) -> f64 {
    let mut penalty = 0.0;
    if reaction.requires_uv {
        penalty += 0.5;
    }
    penalty += reaction
        .external_reactants
        .iter()
        .filter(|external| external.unchecked_mass_reason.is_some())
        .count() as f64
        * 0.75;
    penalty += reaction.surface_requirements.len() as f64 * 0.25;
    if reaction.activation_energy_kj_per_mol > 80.0 {
        penalty += 0.5;
    }
    penalty
}

fn condition_penalty(reaction: &Reaction) -> f64 {
    let mut penalty = reaction.conditions.len() as f64 * 0.25;
    for condition in &reaction.conditions {
        if condition
            .min_temperature_kelvin
            .is_some_and(|value| value > 373.15)
        {
            penalty += 0.5;
        }
        if condition.max_temperature_kelvin.is_some() {
            penalty += 0.25;
        }
        if condition.min_water_activity.is_some() || condition.max_water_activity.is_some() {
            penalty += 0.2;
        }
        if condition.atmosphere.is_some() || condition.max_oxygen_activity.is_some() {
            penalty += 0.3;
        }
        if condition.gas_pressure_atm.is_some_and(|value| value > 1.5) {
            penalty += 0.4;
        }
    }
    penalty
}

fn selectivity_penalty(reaction: &Reaction, product_fraction: f64) -> f64 {
    let mut penalty = if reaction.selectivity_profile.is_some() {
        0.15
    } else {
        0.0
    };
    if reaction.channels.len() > 1 || reaction.product_distribution.is_some() {
        penalty += (1.0 - product_fraction.clamp(0.0, 1.0)) * 2.0;
    }
    penalty
}

fn purification_penalty(reaction: &Reaction, product_fraction: f64) -> f64 {
    let mut penalty = 0.0;
    if reaction.channels.len() > 1 {
        penalty += 0.5;
    }
    if reaction.product_distribution.is_some() {
        penalty += 0.5;
    }
    if product_fraction < 0.75 {
        penalty += 0.5;
    }
    if !reaction.product_phases.is_empty() {
        penalty += 0.15;
    }
    penalty
}

fn external_hints(requirements: &[ExternalRequirement]) -> Vec<SynthesisExternalHint> {
    requirements
        .iter()
        .map(|requirement| SynthesisExternalHint {
            description: requirement.description.clone(),
            moles_per_reaction: requirement.moles_per_reaction,
            unchecked_mass_reason: requirement.unchecked_mass_reason.clone(),
        })
        .collect()
}

fn condition_hints(reaction: &Reaction) -> Vec<SynthesisConditionHint> {
    let mut hints = reaction
        .conditions
        .iter()
        .map(|condition| SynthesisConditionHint {
            description: condition.reason.clone(),
            penalty: 0.25,
        })
        .collect::<Vec<_>>();
    if reaction.requires_uv {
        hints.push(SynthesisConditionHint {
            description: "requires ultraviolet light".to_string(),
            penalty: 0.5,
        });
    }
    for catalyst in &reaction.external_catalysts {
        hints.push(SynthesisConditionHint {
            description: format!("requires catalyst '{}'", catalyst.description),
            penalty: 0.2,
        });
    }
    for requirement in &reaction.surface_requirements {
        hints.push(SynthesisConditionHint {
            description: format!("requires surface '{}'", requirement.surface_id),
            penalty: 0.25,
        });
    }
    hints
}

fn phase_hints(reaction: &Reaction) -> Vec<String> {
    reaction
        .phase_access
        .iter()
        .map(|(substance, access)| {
            let phases = access
                .phases
                .iter()
                .map(|phase| format!("{phase:?}"))
                .collect::<Vec<_>>()
                .join(", ");
            format!("{substance}: {phases}")
        })
        .collect()
}

fn route_condition_hints(steps: &[SynthesisStep]) -> Vec<SynthesisConditionHint> {
    let mut by_description = BTreeMap::<String, f64>::new();
    for step in steps {
        for hint in &step.condition_hints {
            *by_description
                .entry(hint.description.clone())
                .or_insert(0.0) += hint.penalty;
        }
    }
    by_description
        .into_iter()
        .map(|(description, penalty)| SynthesisConditionHint {
            description,
            penalty,
        })
        .collect()
}

fn known_set_key(known: &BTreeSet<SubstanceId>) -> Vec<SubstanceId> {
    known.iter().cloned().collect()
}

fn route_key(route: &SynthesisRoute) -> Vec<String> {
    route
        .steps
        .iter()
        .map(|step| step.reaction_id.to_string())
        .collect()
}

fn compare_routes(left: &SynthesisRoute, right: &SynthesisRoute) -> std::cmp::Ordering {
    left.score
        .total_cmp(&right.score)
        .then_with(|| right.estimated_yield.total_cmp(&left.estimated_yield))
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
    fn planner_uses_private_working_registry() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let before = registry.reactions().count();
        let routes = SynthesisPlanner::new()
            .with_max_steps(2)
            .plan_routes(
                &registry,
                SynthesisRequest::for_frowns("CCO")
                    .with_available_substance("destroy:chloroethane")
                    .with_available_substance("destroy:hydroxide"),
            )
            .unwrap();
        assert!(!routes.is_empty());
        assert_eq!(registry.reactions().count(), before);
    }

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
    fn planner_reports_missing_inputs_for_known_route() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let routes = SynthesisPlanner::new()
            .with_max_steps(1)
            .plan_routes(&registry, SynthesisRequest::for_substance("destroy:water"))
            .unwrap();
        assert!(!routes.is_empty());
        assert!(!routes[0].required_additions.is_empty());
    }

    #[test]
    fn planner_report_marks_found_route_without_mutating_registry() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let before = registry.reactions().count();
        let report = SynthesisPlanner::new()
            .with_max_steps(2)
            .plan_report(
                &registry,
                SynthesisRequest::for_frowns("CCO")
                    .with_available_substance("destroy:chloroethane")
                    .with_available_substance("destroy:hydroxide"),
            )
            .unwrap();

        assert_eq!(report.status, SynthesisPlanStatus::RoutesFound);
        assert!(!report.routes.is_empty());
        assert!(report.required_additions.is_empty());
        assert!(report.unsupported_reason.is_none());
        assert_eq!(registry.reactions().count(), before);
    }

    #[test]
    fn planner_report_marks_known_route_with_missing_inputs() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let report = SynthesisPlanner::new()
            .with_max_steps(1)
            .plan_report(&registry, SynthesisRequest::for_substance("destroy:water"))
            .unwrap();

        assert_eq!(report.status, SynthesisPlanStatus::RequiresAdditionalInputs);
        assert!(!report.routes.is_empty());
        assert!(!report.required_additions.is_empty());
        assert!(report.unsupported_reason.is_none());
    }

    #[test]
    fn planner_report_with_empty_inputs_suggests_dynamic_reactants() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let before = registry.reactions().count();
        let report = SynthesisPlanner::new()
            .with_max_steps(2)
            .allow_reaction_prefix("halide_hydroxide_substitution/")
            .plan_report(&registry, SynthesisRequest::for_frowns("CCO"))
            .unwrap();

        assert_eq!(report.status, SynthesisPlanStatus::RequiresAdditionalInputs);
        assert!(report.routes.is_empty());
        assert!(report
            .required_additions
            .iter()
            .any(|requirement| requirement.substance_id.as_str() == "destroy:chloroethane"));
        assert!(report
            .required_additions
            .iter()
            .any(|requirement| requirement.substance_id.as_str() == "destroy:hydroxide"));
        assert!(report.unsupported_reason.is_none());
        assert_eq!(registry.reactions().count(), before);
    }

    #[test]
    fn planner_report_marks_target_without_known_producer() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let report = SynthesisPlanner::new()
            .with_max_steps(2)
            .plan_report(&registry, SynthesisRequest::for_substance("destroy:argon"))
            .unwrap();

        assert_eq!(
            report.status,
            SynthesisPlanStatus::UnsupportedByCurrentModel
        );
        assert!(report.routes.is_empty());
        assert!(report.required_additions.is_empty());
        assert!(report
            .unsupported_reason
            .as_deref()
            .is_some_and(|reason| reason.contains("destroy:argon")));
    }

    #[test]
    fn reachability_report_distinguishes_reachable_targets_from_model_gaps() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let before_reactions = registry.reactions().count();
        let report = SynthesisPlanner::new()
            .analyze_reachability(
                &registry,
                SynthesisReachabilityRequest::new([
                    SubstanceId::from("destroy:ethanol"),
                    SubstanceId::from("destroy:argon"),
                ])
                .with_available_substance("destroy:chloroethane")
                .with_available_substance("destroy:hydroxide")
                .with_generation_iterations(1)
                .with_max_steps(2)
                .with_max_routes_per_target(2),
            )
            .unwrap();

        assert_eq!(report.starting_substance_count, 2);
        assert_eq!(report.target_count, 2);
        assert!(report
            .reachable_targets
            .contains(&SubstanceId::from("destroy:ethanol")));
        assert!(report
            .unreachable_targets
            .contains(&SubstanceId::from("destroy:argon")));
        assert!(report
            .target_reports
            .iter()
            .find(|target| target.target == SubstanceId::from("destroy:ethanol"))
            .is_some_and(|target| {
                target.status == SynthesisPlanStatus::RoutesFound
                    && target.has_known_producer
                    && target.producer_reactions.iter().any(|reaction| {
                        reaction
                            .as_str()
                            .starts_with("halide_hydroxide_substitution/")
                    })
            }));
        assert!(report
            .target_reports
            .iter()
            .find(|target| target.target == SubstanceId::from("destroy:argon"))
            .is_some_and(|target| {
                target.status == SynthesisPlanStatus::UnsupportedByCurrentModel
                    && !target.has_known_producer
                    && target.missing_reactants.is_empty()
            }));
        assert_eq!(registry.reactions().count(), before_reactions);
    }

    #[test]
    fn planner_report_distinguishes_search_limit_from_unsupported_target() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let report = SynthesisPlanner::new()
            .with_max_steps(5)
            .with_max_routes(4)
            .allow_reaction_prefix("halide_cyanide_substitution/")
            .allow_reaction_prefix("nitrile_hydrolysis/")
            .allow_reaction_prefix("amide_hydrolysis/")
            .allow_reaction_prefix("acyl_chloride_formation/")
            .allow_reaction_prefix("acyl_chloride_esterification/")
            .plan_report(
                &registry,
                SynthesisRequest::for_frowns("CCC(=O)OCC")
                    .with_available_substance("destroy:chloroethane")
                    .with_available_substance("destroy:cyanide")
                    .with_available_substance("destroy:water")
                    .with_available_substance("destroy:phosgene")
                    .with_available_substance("destroy:ethanol")
                    .with_max_steps(2),
            )
            .unwrap();

        assert_eq!(report.status, SynthesisPlanStatus::SearchLimitReached);
        assert!(report.routes.is_empty());
        assert!(report
            .unsupported_reason
            .as_deref()
            .is_some_and(|reason| reason.contains("2 step")));
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

    #[test]
    fn safety_policy_filters_denied_missing_reactants_in_backward_search() {
        let registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let policy = SynthesisSafetyPolicy::new()
            .deny_substance("destroy:hydroxide", "test policy denies hydroxide input")
            .unwrap();
        let routes = SynthesisPlanner::new()
            .with_safety_policy(policy)
            .with_max_steps(1)
            .allow_reaction_prefix("destroy:neutralization")
            .plan_routes(&registry, SynthesisRequest::for_substance("destroy:water"))
            .unwrap();

        assert!(routes.is_empty());
    }

    #[test]
    fn planner_uses_tms_protection_as_synthesis_step() {
        let mut setup_registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        setup_registry
            .generate_reactions_for_substances(
                [
                    SubstanceId::from("destroy:ethanol"),
                    SubstanceId::from("destroy:trimethylsilyl_chloride"),
                ],
                1,
            )
            .unwrap();
        let protected_id = setup_registry
            .reactions()
            .find(|reaction| {
                reaction
                    .id
                    .as_str()
                    .starts_with("alcohol_silyl_protection/")
                    && reaction
                        .reactants
                        .iter()
                        .any(|term| term.substance_id == SubstanceId::from("destroy:ethanol"))
            })
            .and_then(|reaction| reaction.products.first())
            .map(|term| term.substance_id.clone())
            .expect("ethanol must generate a TMS-protected product");
        let target_structure = setup_registry
            .substance(&protected_id)
            .unwrap()
            .molecular_structure
            .as_ref()
            .unwrap()
            .clone();

        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let routes = SynthesisPlanner::new()
            .with_max_steps(1)
            .allow_reaction_prefix("alcohol_silyl_protection/")
            .find_routes(
                &mut registry,
                [
                    SubstanceId::from("destroy:ethanol"),
                    SubstanceId::from("destroy:trimethylsilyl_chloride"),
                ],
                target_structure,
            )
            .unwrap();

        assert!(!routes.is_empty());
        assert!(routes[0].steps.iter().any(|step| step
            .reaction_id
            .as_str()
            .starts_with("alcohol_silyl_protection/")));
    }

    #[test]
    fn planner_uses_tms_deprotection_to_restore_alcohol() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        registry
            .generate_reactions_for_substances(
                [
                    SubstanceId::from("destroy:ethanol"),
                    SubstanceId::from("destroy:trimethylsilyl_chloride"),
                ],
                1,
            )
            .unwrap();
        let protected_id = registry
            .reactions()
            .find(|reaction| {
                reaction
                    .id
                    .as_str()
                    .starts_with("alcohol_silyl_protection/")
                    && reaction
                        .reactants
                        .iter()
                        .any(|term| term.substance_id == SubstanceId::from("destroy:ethanol"))
            })
            .and_then(|reaction| reaction.products.first())
            .map(|term| term.substance_id.clone())
            .expect("ethanol must generate a TMS-protected product");

        let routes = SynthesisPlanner::new()
            .with_max_steps(1)
            .allow_reaction_prefix("silyl_ether_deprotection/")
            .find_routes(
                &mut registry,
                [
                    protected_id,
                    SubstanceId::from("destroy:fluoride"),
                    SubstanceId::from("destroy:proton"),
                ],
                parse_frowns("CCO").unwrap(),
            )
            .unwrap();

        assert!(!routes.is_empty());
        assert!(routes[0].steps.iter().any(|step| step
            .reaction_id
            .as_str()
            .starts_with("silyl_ether_deprotection/")));
    }

    #[test]
    fn planner_finds_methane_to_acetylene_route_without_prefix_filters() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let routes = SynthesisPlanner::new()
            .with_max_steps(4)
            .with_max_routes(5)
            .find_routes(
                &mut registry,
                [SubstanceId::from("destroy:methane")],
                parse_frowns("C#C").unwrap(),
            )
            .unwrap();

        assert!(!routes.is_empty());
        let best = &routes[0];
        assert_eq!(best.steps.len(), 3);
        assert!(best.required_additions.is_empty());
        assert!(best
            .explanation
            .contains("route can be run from the available substances"));

        let reaction_ids = best
            .steps
            .iter()
            .map(|step| step.reaction_id.as_str())
            .collect::<Vec<_>>();
        assert!(
            reaction_ids[0].starts_with("dehydrogenative_coupling/destroy_methane/destroy_methane")
        );
        assert!(reaction_ids[1].starts_with("pyrolysis/destroy_linear_C_C_"));
        assert!(reaction_ids[2].starts_with("pyrolysis/destroy_ethene/"));
        assert!(best.steps.last().is_some_and(|step| step
            .products
            .contains(&SubstanceId::from("destroy:acetylene"))));
    }

    #[test]
    fn planner_reaches_iodomethane_through_general_alcohol_hydrohalogenation() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let routes = SynthesisPlanner::new()
            .with_max_steps(1)
            .allow_reaction_prefix("alcohol_hydrohalogenation/")
            .find_routes(
                &mut registry,
                [
                    SubstanceId::from("destroy:methanol"),
                    SubstanceId::from("destroy:hydroiodic_acid"),
                ],
                parse_frowns("CI").unwrap(),
            )
            .unwrap();

        assert!(!routes.is_empty());
        assert!(routes[0].steps.iter().any(|step| step
            .reaction_id
            .as_str()
            .starts_with("alcohol_hydrohalogenation/")));
        assert!(routes[0].steps.iter().any(|step| step
            .products
            .contains(&SubstanceId::from("destroy:iodomethane"))));
    }

    #[test]
    fn planner_reaches_trimethyl_borate_by_repeating_borate_esterification() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let routes = SynthesisPlanner::new()
            .with_max_steps(3)
            .allow_reaction_prefix("borate_esterification/")
            .find_routes(
                &mut registry,
                [
                    SubstanceId::from("destroy:boric_acid"),
                    SubstanceId::from("destroy:methanol"),
                ],
                parse_frowns("COB(OC)OC").unwrap(),
            )
            .unwrap();

        assert!(!routes.is_empty());
        assert_eq!(routes[0].steps.len(), 3);
        assert!(routes[0].steps.iter().all(|step| step
            .reaction_id
            .as_str()
            .starts_with("borate_esterification/")));
        assert!(routes[0].steps.last().is_some_and(|step| step
            .products
            .contains(&SubstanceId::from("destroy:trimethyl_borate"))));
    }
}
