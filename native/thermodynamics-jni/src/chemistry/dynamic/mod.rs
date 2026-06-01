use std::collections::BTreeMap;
use std::collections::VecDeque;

use super::error::{ChemistryError, ChemistryResult};
use super::frowns::{parse_frowns, write_frowns};
use super::functional_group::{FunctionalGroup, FunctionalGroupType};
use super::kinetics::EnergyModel;
use super::molecule::{MolecularEditor, MolecularStructure, Stereochemistry};
use super::organic::{self, OrganicGenerationSpace};
use super::reaction::{Reaction, ReactionId, StoichiometricTerm};
use super::reactive_site::{
    try_find_reactive_sites, ReactiveRole, ReactiveSiteKey, ReactiveSiteKind,
};
use super::registry::{ChemistryRegistry, ChemistryRegistryBuilder, SubstanceIndex};
use super::solution::AcidBaseSpec;
use super::substance::{
    LiquidPhasePreference, SolventRole, Substance, SubstanceId, SubstancePhaseProperties,
    SubstanceTagId,
};

const DEFAULT_DYNAMIC_DENSITY: f64 = 1000.0;
const DEFAULT_DYNAMIC_HEAT_CAPACITY: f64 = 100.0;
const DEFAULT_DYNAMIC_LATENT_HEAT: f64 = 20_000.0;
const DEFAULT_DYNAMIC_COLOR: u32 = 0x20FF_FFFF;
const MAX_DYNAMIC_ATOMS: usize = 100;
const MAX_DYNAMIC_WORK_ITEMS: usize = 1_000_000;
const MAX_DYNAMIC_QUEUE_ITEMS: usize = 100_000;
const MASS_TOLERANCE_GRAMS_PER_MOL: f64 = 1.0e-6;

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct DynamicReactionIndex(usize);

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum KnownSubstanceIndex {
    Static(SubstanceIndex),
    Dynamic(usize),
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
enum OrganicGeneratorKind {
    HalideHydroxideSubstitution,
    HalideAmmoniaSubstitution,
    HalideCyanideSubstitution,
    HalideAmineSubstitution,
    AlcoholOxidation,
    AlcoholDehydration,
    ThionylChlorideSubstitution,
    CarboxylicAcidEsterification,
    AcylChlorideEsterification,
    AlkoxideProtonation,
    NitrileHydrolysis,
    NitrileHydrogenation,
    NitroHydrogenation,
    AcylChlorideHydrolysis,
    AcylChlorideFormation,
    AldehydeOxidation,
    CyanideNucleophilicAddition,
    WolffKishnerReduction,
    AmideHydrolysis,
    AminePhosgenation,
    CyanamideAddition,
    IsocyanateHydrolysis,
    BoraneOxidation,
    BorateEsterHydrolysis,
    AlkeneChlorination,
    AlkeneChlorohydrination,
    AlkeneHydrolysis,
    AlkeneHydroborationWithBorane,
    AlkeneHydrochlorination,
    AlkeneHydrogenation,
    AlkeneHydroiodination,
    AlkeneIodination,
    AlkyneChlorination,
    AlkyneChlorohydrination,
    AlkyneHydrolysis,
    AlkyneHydroborationWithBorane,
    AlkyneHydrochlorination,
    AlkyneHydrogenation,
    AlkyneHydroiodination,
    AlkyneIodination,
    AromaticNitration,
    EpoxideHydrolysis,
    OrganometallicCarbonylAddition,
    AldolAddition,
    AlphaHalogenation,
    AldolDehydration,
    EnamineFormation,
    EnolateAlkylation,
    MichaelAddition,
    ClaisenCondensation,
    PhosphoniumSaltFormation,
    PhosphoniumYlideFormation,
    WittigOlefination,
    HornerWadsworthEmmonsOlefination,
    JuliaOlefination,
    AlcoholSilylProtection,
    SilylEtherDeprotection,
    AcetalDeprotection,
    AmineBocProtection,
    BocDeprotection,
    AmineCbzProtection,
    CbzDeprotection,
}

impl OrganicGeneratorKind {
    fn bit(self) -> u64 {
        1_u64 << self.ordinal()
    }

    fn ordinal(self) -> u32 {
        match self {
            OrganicGeneratorKind::HalideHydroxideSubstitution => 0,
            OrganicGeneratorKind::HalideAmmoniaSubstitution => 1,
            OrganicGeneratorKind::HalideCyanideSubstitution => 2,
            OrganicGeneratorKind::HalideAmineSubstitution => 3,
            OrganicGeneratorKind::AlcoholOxidation => 4,
            OrganicGeneratorKind::AlcoholDehydration => 5,
            OrganicGeneratorKind::ThionylChlorideSubstitution => 6,
            OrganicGeneratorKind::CarboxylicAcidEsterification => 7,
            OrganicGeneratorKind::AcylChlorideEsterification => 8,
            OrganicGeneratorKind::AlkoxideProtonation => 9,
            OrganicGeneratorKind::NitrileHydrolysis => 10,
            OrganicGeneratorKind::NitrileHydrogenation => 11,
            OrganicGeneratorKind::NitroHydrogenation => 12,
            OrganicGeneratorKind::AcylChlorideHydrolysis => 13,
            OrganicGeneratorKind::AcylChlorideFormation => 14,
            OrganicGeneratorKind::AldehydeOxidation => 15,
            OrganicGeneratorKind::CyanideNucleophilicAddition => 16,
            OrganicGeneratorKind::WolffKishnerReduction => 17,
            OrganicGeneratorKind::AmideHydrolysis => 18,
            OrganicGeneratorKind::AminePhosgenation => 19,
            OrganicGeneratorKind::CyanamideAddition => 20,
            OrganicGeneratorKind::IsocyanateHydrolysis => 21,
            OrganicGeneratorKind::BoraneOxidation => 22,
            OrganicGeneratorKind::BorateEsterHydrolysis => 23,
            OrganicGeneratorKind::AlkeneChlorination => 24,
            OrganicGeneratorKind::AlkeneChlorohydrination => 25,
            OrganicGeneratorKind::AlkeneHydrolysis => 26,
            OrganicGeneratorKind::AlkeneHydroborationWithBorane => 27,
            OrganicGeneratorKind::AlkeneHydrochlorination => 28,
            OrganicGeneratorKind::AlkeneHydrogenation => 29,
            OrganicGeneratorKind::AlkeneHydroiodination => 30,
            OrganicGeneratorKind::AlkeneIodination => 31,
            OrganicGeneratorKind::AlkyneChlorination => 32,
            OrganicGeneratorKind::AlkyneChlorohydrination => 33,
            OrganicGeneratorKind::AlkyneHydrolysis => 34,
            OrganicGeneratorKind::AlkyneHydroborationWithBorane => 35,
            OrganicGeneratorKind::AlkyneHydrochlorination => 36,
            OrganicGeneratorKind::AlkyneHydrogenation => 37,
            OrganicGeneratorKind::AlkyneHydroiodination => 38,
            OrganicGeneratorKind::AlkyneIodination => 39,
            OrganicGeneratorKind::AromaticNitration => 40,
            OrganicGeneratorKind::EpoxideHydrolysis => 41,
            OrganicGeneratorKind::OrganometallicCarbonylAddition => 42,
            OrganicGeneratorKind::AldolAddition => 43,
            OrganicGeneratorKind::AlphaHalogenation => 44,
            OrganicGeneratorKind::AldolDehydration => 45,
            OrganicGeneratorKind::EnamineFormation => 46,
            OrganicGeneratorKind::EnolateAlkylation => 47,
            OrganicGeneratorKind::MichaelAddition => 48,
            OrganicGeneratorKind::ClaisenCondensation => 49,
            OrganicGeneratorKind::PhosphoniumSaltFormation => 50,
            OrganicGeneratorKind::PhosphoniumYlideFormation => 51,
            OrganicGeneratorKind::WittigOlefination => 52,
            OrganicGeneratorKind::HornerWadsworthEmmonsOlefination => 53,
            OrganicGeneratorKind::JuliaOlefination => 54,
            OrganicGeneratorKind::AlcoholSilylProtection => 55,
            OrganicGeneratorKind::SilylEtherDeprotection => 56,
            OrganicGeneratorKind::AcetalDeprotection => 57,
            OrganicGeneratorKind::AmineBocProtection => 58,
            OrganicGeneratorKind::BocDeprotection => 59,
            OrganicGeneratorKind::AmineCbzProtection => 60,
            OrganicGeneratorKind::CbzDeprotection => 61,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
struct SiteHandle {
    substance: KnownSubstanceIndex,
    site_index: usize,
}

#[derive(Debug, Clone)]
struct SiteBucket {
    site_kind: ReactiveSiteKind,
    handles: Vec<SiteHandle>,
}

#[derive(Debug, Clone)]
pub struct DynamicChemistryRegistry {
    static_registry: ChemistryRegistry,
    dynamic_substances: Vec<Substance>,
    dynamic_substance_id_to_index: BTreeMap<SubstanceId, usize>,
    dynamic_reactions: Vec<Reaction>,
    dynamic_reaction_id_to_index: BTreeMap<ReactionId, DynamicReactionIndex>,
    dynamic_reaction_index_by_substance: Vec<Vec<DynamicReactionIndex>>,
    dynamic_unindexed_reaction_indices: Vec<DynamicReactionIndex>,
    canonical_to_id: BTreeMap<String, SubstanceId>,
    canonical_by_id: BTreeMap<SubstanceId, String>,
    site_index: Vec<SiteBucket>,
    site_handles_by_substance: Vec<Vec<SiteHandle>>,
    processed_generation_masks: BTreeMap<(usize, ReactiveSiteKey), u64>,
    dynamic_acid_base_specs: Vec<AcidBaseSpec>,
    energy_model: EnergyModel,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DynamicGenerationReport {
    pub iterations: usize,
    pub added_substances: usize,
    pub added_reactions: usize,
    pub processed_work_items: usize,
    pub skipped_duplicates: usize,
    pub remaining_queue: usize,
    pub reached_fixed_point: bool,
    pub generator_errors: Vec<String>,
}

impl DynamicChemistryRegistry {
    pub fn from_destroy_catalog() -> ChemistryResult<Self> {
        Self::from_registry(super::destroy_registry_builder()?.build()?)
    }

    pub fn from_registry(registry: ChemistryRegistry) -> ChemistryResult<Self> {
        let static_substance_count = registry.substance_count();
        let mut result = Self {
            static_registry: registry,
            dynamic_substances: Vec::new(),
            dynamic_substance_id_to_index: BTreeMap::new(),
            dynamic_reactions: Vec::new(),
            dynamic_reaction_id_to_index: BTreeMap::new(),
            dynamic_reaction_index_by_substance: vec![Vec::new(); static_substance_count],
            dynamic_unindexed_reaction_indices: Vec::new(),
            canonical_to_id: BTreeMap::new(),
            canonical_by_id: BTreeMap::new(),
            site_index: Vec::new(),
            site_handles_by_substance: vec![Vec::new(); static_substance_count],
            processed_generation_masks: BTreeMap::new(),
            dynamic_acid_base_specs: Vec::new(),
            energy_model: EnergyModel::new(),
        };
        result.rebuild_canonical_index()?;
        result.rebuild_site_index()?;
        Ok(result)
    }

    pub fn static_registry(&self) -> &ChemistryRegistry {
        &self.static_registry
    }

    pub fn energy_model(&self) -> &EnergyModel {
        &self.energy_model
    }

    pub fn energy_model_mut(&mut self) -> &mut EnergyModel {
        &mut self.energy_model
    }

    pub fn with_energy_model(mut self, energy_model: EnergyModel) -> Self {
        self.energy_model = energy_model;
        self
    }

    pub fn to_registry(&self) -> ChemistryResult<ChemistryRegistry> {
        let mut builder = ChemistryRegistryBuilder::from_registry(&self.static_registry);
        for substance in &self.dynamic_substances {
            builder = builder.substance(substance.clone());
        }
        for reaction in &self.dynamic_reactions {
            builder = builder.reaction(reaction.clone());
        }
        for spec in &self.dynamic_acid_base_specs {
            builder = builder.acid_base_pair(spec.clone());
        }
        builder.build()
    }

    pub fn substance(&self, id: &SubstanceId) -> ChemistryResult<&Substance> {
        self.dynamic_substance_id_to_index
            .get(id)
            .and_then(|index| self.dynamic_substances.get(*index))
            .or_else(|| self.static_registry.substance(id).ok())
            .ok_or_else(|| ChemistryError::InvalidMixtureState(format!("unknown substance '{id}'")))
    }

    pub fn reaction(&self, id: &ReactionId) -> ChemistryResult<&Reaction> {
        self.dynamic_reaction_id_to_index
            .get(id)
            .and_then(|index| self.dynamic_reactions.get(index.0))
            .or_else(|| self.static_registry.reaction(id).ok())
            .ok_or_else(|| ChemistryError::UnknownReaction(id.to_string()))
    }

    pub fn substances(&self) -> impl Iterator<Item = &Substance> {
        self.static_registry
            .substances()
            .chain(self.dynamic_substances.iter())
    }

    pub fn reactions(&self) -> impl Iterator<Item = &Reaction> {
        self.static_registry
            .reactions()
            .chain(self.dynamic_reactions.iter())
    }

    pub fn validate_substance_can_enter_mixture(
        &self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<()> {
        let substance = self.substance(substance_id)?;
        if substance
            .tags
            .iter()
            .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
        {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "hypothetical substance '{substance_id}' cannot be added to a mixture"
            )));
        }
        Ok(())
    }

    pub fn reaction_candidates_for_substances<'registry, 'substances, I>(
        &'registry self,
        substances: I,
    ) -> Vec<&'registry Reaction>
    where
        I: IntoIterator<Item = &'substances SubstanceId>,
    {
        let substance_ids = substances.into_iter().collect::<Vec<_>>();
        let mut seen = vec![false; self.dynamic_reactions.len()];
        let mut dynamic_reaction_indices = Vec::new();
        for reaction_index in &self.dynamic_unindexed_reaction_indices {
            mark_dynamic_reaction_candidate(
                &mut seen,
                &mut dynamic_reaction_indices,
                *reaction_index,
            );
        }
        for substance_id in &substance_ids {
            let Some(substance_index) = self.known_substance_index(substance_id) else {
                continue;
            };
            let slot =
                known_substance_slot(self.static_registry.substance_count(), substance_index);
            if let Some(indexed_reactions) = self.dynamic_reaction_index_by_substance.get(slot) {
                for reaction_index in indexed_reactions {
                    mark_dynamic_reaction_candidate(
                        &mut seen,
                        &mut dynamic_reaction_indices,
                        *reaction_index,
                    );
                }
            }
        }
        let mut result = self
            .static_registry
            .reaction_candidates_for_substances(substance_ids)
            .into_iter()
            .collect::<Vec<_>>();
        result.extend(
            dynamic_reaction_indices
                .into_iter()
                .filter_map(|reaction_index| self.dynamic_reaction_by_index(reaction_index)),
        );
        result
    }

    pub fn resolve_frowns(&mut self, code: &str) -> ChemistryResult<SubstanceId> {
        self.resolve_structure(parse_frowns(code)?)
    }

    pub fn resolve_structure(
        &mut self,
        structure: MolecularStructure,
    ) -> ChemistryResult<SubstanceId> {
        let canonical = write_frowns(&structure)?;
        if let Some(id) = self.canonical_to_id.get(&canonical) {
            return Ok(id.clone());
        }
        let substance = build_dynamic_substance(canonical.clone(), structure)?;
        let id = substance.id.clone();
        self.add_dynamic_substance_with_canonical(substance, Some(canonical))?;
        Ok(id)
    }

    pub fn generate_reactions_for(
        &mut self,
        substance_id: &SubstanceId,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        self.substance(substance_id)?;
        let seeds = vec![self.known_substance_index_or_error(substance_id)?];
        self.generate_reactions_from_scope(seeds, Some(max_iterations))
    }

    pub fn generate_reactions_for_to_fixed_point(
        &mut self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<DynamicGenerationReport> {
        self.substance(substance_id)?;
        let seeds = vec![self.known_substance_index_or_error(substance_id)?];
        self.generate_reactions_from_scope(seeds, None)
    }

    pub fn generate_reactions_for_substances(
        &mut self,
        substance_ids: impl IntoIterator<Item = SubstanceId>,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        let seeds = self.validated_substance_indices(substance_ids)?;
        self.generate_reactions_from_scope(seeds, Some(max_iterations))
    }

    pub fn generate_reactions_for_substances_to_fixed_point(
        &mut self,
        substance_ids: impl IntoIterator<Item = SubstanceId>,
    ) -> ChemistryResult<DynamicGenerationReport> {
        let seeds = self.validated_substance_indices(substance_ids)?;
        self.generate_reactions_from_scope(seeds, None)
    }

    pub fn generate_reactions(
        &mut self,
        max_iterations: usize,
    ) -> ChemistryResult<DynamicGenerationReport> {
        self.generate_reactions_from_scope(self.all_known_substance_indices(), Some(max_iterations))
    }

    pub fn generate_to_fixed_point(&mut self) -> ChemistryResult<DynamicGenerationReport> {
        self.generate_reactions_from_scope(self.all_known_substance_indices(), None)
    }

    fn generate_reactions_from_scope(
        &mut self,
        seeds: Vec<KnownSubstanceIndex>,
        max_iterations: Option<usize>,
    ) -> ChemistryResult<DynamicGenerationReport> {
        if max_iterations == Some(0) {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<dynamic-generation>".to_string(),
                reason: "max_iterations must be greater than zero".to_string(),
            });
        }

        let mut added_substances = 0usize;
        let mut added_reactions = 0usize;
        let mut processed_work_items = 0usize;
        let mut skipped_duplicates = 0usize;
        let mut queue = VecDeque::new();
        let mut queued = vec![false; self.known_substance_count()];
        let mut scope = vec![false; self.known_substance_count()];
        for seed in seeds {
            self.mark_known_substance(&mut scope, seed, true);
            enqueue_known_substance(
                self.static_registry.substance_count(),
                &mut queue,
                &mut queued,
                seed,
            );
        }
        let mut iteration = 0usize;
        loop {
            if max_iterations.is_some_and(|max| iteration >= max) {
                return Ok(DynamicGenerationReport {
                    iterations: iteration,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: queue.len(),
                    reached_fixed_point: false,
                    generator_errors: Vec::new(),
                });
            }
            if queue.is_empty() {
                return Ok(DynamicGenerationReport {
                    iterations: iteration,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: 0,
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }
            let batch_len = queue.len();
            let mut unprocessed_seeds = Vec::new();
            let mut pending_generation_mask_updates = Vec::new();
            for _ in 0..batch_len {
                let seed = queue.pop_front().expect("batch length was measured");
                self.mark_known_substance(&mut queued, seed, false);
                let mask_updates = self.generation_mask_updates_for_substance(seed)?;
                if mask_updates.is_empty() {
                    skipped_duplicates += 1;
                    continue;
                }
                pending_generation_mask_updates.extend(mask_updates);
                unprocessed_seeds.push(seed);
                processed_work_items += 1;
                if processed_work_items > MAX_DYNAMIC_WORK_ITEMS {
                    let seed_id = self.known_substance_id(seed)?.to_string();
                    return Err(ChemistryError::GenerationInvariantViolation {
                        generator: "<dynamic-generation>".to_string(),
                        substance_id: seed_id,
                        reason: format!(
                            "processed more than {MAX_DYNAMIC_WORK_ITEMS} work items without reaching a fixed point"
                        ),
                    });
                }
            }
            if unprocessed_seeds.is_empty() {
                return Ok(DynamicGenerationReport {
                    iterations: iteration + 1,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: 0,
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }

            let generated =
                self.generate_organic_for_known_substances(&unprocessed_seeds, &scope)?;
            let mut changed = false;

            let mut staged = self.clone();
            staged.apply_generation_mask_updates(&pending_generation_mask_updates);
            let mut generated_id_remap = BTreeMap::new();
            let mut new_substance_ids = Vec::new();
            for substance in generated.substances {
                if staged.substance(&substance.id).is_ok() {
                    skipped_duplicates += 1;
                    continue;
                }
                let generated_id = substance.id.clone();
                let canonical = substance
                    .molecular_structure
                    .as_ref()
                    .map(write_frowns)
                    .transpose()?
                    .ok_or_else(|| ChemistryError::InvalidSubstance {
                        substance_id: substance.id.to_string(),
                        reason: "generated dynamic substance has no structure".to_string(),
                    })?;
                if let Some(existing) = staged.canonical_to_id.get(&canonical) {
                    generated_id_remap.insert(generated_id, existing.clone());
                    skipped_duplicates += 1;
                    continue;
                }
                let substance_id = substance.id.clone();
                staged.add_dynamic_substance_with_canonical(substance, Some(canonical))?;
                staged.known_substance_index(&substance_id).ok_or_else(|| {
                    ChemistryError::InvalidMixtureState(format!(
                        "generated substance '{substance_id}' was not indexed",
                    ))
                })?;
                new_substance_ids.push(substance_id);
                added_substances += 1;
                changed = true;
            }

            for reaction in generated.reactions {
                let reaction = remap_reaction_substances(reaction, &generated_id_remap);
                if staged.reaction(&reaction.id).is_ok() {
                    skipped_duplicates += 1;
                    continue;
                }
                staged.add_dynamic_reaction(reaction)?;
                added_reactions += 1;
                changed = true;
            }

            let mut staged_scope = scope.clone();
            let mut staged_queue = queue.clone();
            let mut staged_queued = queued.clone();
            for substance_id in &new_substance_ids {
                let new_index = staged.known_substance_index(substance_id).ok_or_else(|| {
                    ChemistryError::InvalidMixtureState(format!(
                        "generated substance '{substance_id}' was not indexed",
                    ))
                })?;
                staged.mark_known_substance(&mut staged_scope, new_index, true);
                enqueue_known_substance(
                    staged.static_registry.substance_count(),
                    &mut staged_queue,
                    &mut staged_queued,
                    new_index,
                );
                if staged_queue.len() > MAX_DYNAMIC_QUEUE_ITEMS {
                    return Err(ChemistryError::GenerationInvariantViolation {
                        generator: "<dynamic-generation>".to_string(),
                        substance_id: substance_id.to_string(),
                        reason: format!(
                            "dynamic generation queue exceeded {MAX_DYNAMIC_QUEUE_ITEMS} substances before reaching a fixed point"
                        ),
                    });
                }
            }

            *self = staged;
            scope = staged_scope;
            queue = staged_queue;
            queued = staged_queued;

            if !changed {
                return Ok(DynamicGenerationReport {
                    iterations: iteration + 1,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: queue.len(),
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }
            if queue.is_empty() {
                return Ok(DynamicGenerationReport {
                    iterations: iteration + 1,
                    added_substances,
                    added_reactions,
                    processed_work_items,
                    skipped_duplicates,
                    remaining_queue: 0,
                    reached_fixed_point: true,
                    generator_errors: Vec::new(),
                });
            }
            iteration += 1;
        }
    }

    fn rebuild_canonical_index(&mut self) -> ChemistryResult<()> {
        self.canonical_to_id.clear();
        self.canonical_by_id.clear();
        let canonical_entries = self
            .substances()
            .filter_map(|substance| {
                substance
                    .molecular_structure
                    .as_ref()
                    .map(|structure| (substance.id.clone(), structure))
            })
            .map(|(id, structure)| write_frowns(structure).map(|canonical| (canonical, id)))
            .collect::<ChemistryResult<Vec<_>>>()?;
        for (canonical, id) in canonical_entries {
            self.canonical_to_id
                .entry(canonical.clone())
                .or_insert_with(|| id.clone());
            self.canonical_by_id.insert(id, canonical);
        }
        Ok(())
    }

    fn rebuild_site_index(&mut self) -> ChemistryResult<()> {
        self.site_index.clear();
        self.site_handles_by_substance = vec![
            Vec::new();
            self.static_registry.substance_count()
                + self.dynamic_substances.len()
        ];
        let substance_indices = self.static_registry.substance_indices().collect::<Vec<_>>();
        for substance_index in substance_indices {
            self.add_site_handles_for_substance(KnownSubstanceIndex::Static(substance_index))?;
        }
        Ok(())
    }

    fn add_dynamic_substance_with_canonical(
        &mut self,
        substance: Substance,
        canonical: Option<String>,
    ) -> ChemistryResult<()> {
        substance.validate()?;
        if self.static_registry.substance(&substance.id).is_ok()
            || self
                .dynamic_substance_id_to_index
                .contains_key(&substance.id)
        {
            return Err(ChemistryError::DuplicateSubstance(substance.id.to_string()));
        }
        for tag in &substance.tags {
            if !self.static_registry.has_substance_tag(tag) {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: substance.id.to_string(),
                    reason: format!("unknown substance tag '{tag}'"),
                });
            }
        }
        let canonical = match (canonical, substance.molecular_structure.as_ref()) {
            (Some(canonical), _) => Some(canonical),
            (None, Some(structure)) => Some(write_frowns(structure)?),
            (None, None) => None,
        };
        if let Some(canonical) = canonical {
            self.canonical_to_id
                .entry(canonical.clone())
                .or_insert_with(|| substance.id.clone());
            self.canonical_by_id.insert(substance.id.clone(), canonical);
        }
        let substance_index = self.dynamic_substances.len();
        self.dynamic_substance_id_to_index
            .insert(substance.id.clone(), substance_index);
        self.dynamic_substances.push(substance);
        self.dynamic_reaction_index_by_substance.push(Vec::new());
        self.site_handles_by_substance.push(Vec::new());
        self.add_site_handles_for_substance(KnownSubstanceIndex::Dynamic(substance_index))?;
        self.register_dynamic_acidity_for_substance(substance_index)?;
        Ok(())
    }

    fn register_dynamic_acidity_for_substance(
        &mut self,
        substance_index: usize,
    ) -> ChemistryResult<()> {
        let acid = self
            .dynamic_substances
            .get(substance_index)
            .ok_or_else(|| {
                ChemistryError::InvalidMixtureState(format!(
                    "invalid dynamic substance index {substance_index}"
                ))
            })?;
        if acid
            .tags
            .iter()
            .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
        {
            return Ok(());
        }
        let Some(structure) = acid.molecular_structure.as_ref() else {
            return Ok(());
        };
        let acid_id = acid.id.clone();
        let structure = structure.clone();
        let groups = acid.functional_groups.clone();
        for (group_index, group) in groups.iter().enumerate() {
            let Some((base_structure, pka)) =
                conjugate_base_for_group(&structure, group, &acid_id)?
            else {
                continue;
            };
            let canonical = write_frowns(&base_structure)?;
            let base_id = if let Some(existing) = self.canonical_to_id.get(&canonical) {
                existing.clone()
            } else {
                let base = build_dynamic_substance(canonical.clone(), base_structure)?;
                let base_id = base.id.clone();
                self.add_dynamic_substance_with_canonical(base, Some(canonical))?;
                base_id
            };
            if self
                .static_registry
                .acid_base_specs()
                .any(|spec| spec.acid == acid_id && spec.conjugate_base == base_id)
                || self
                    .dynamic_acid_base_specs
                    .iter()
                    .any(|spec| spec.acid == acid_id && spec.conjugate_base == base_id)
            {
                continue;
            }
            self.dynamic_acid_base_specs.push(AcidBaseSpec::new(
                format!(
                    "{}/acid_center_{}",
                    acid_id.as_str().replace(':', "_"),
                    group_index
                ),
                acid_id.clone(),
                base_id,
                pka,
            ));
        }
        Ok(())
    }

    fn add_dynamic_reaction(&mut self, reaction: Reaction) -> ChemistryResult<()> {
        if self.static_registry.reaction(&reaction.id).is_ok()
            || self.dynamic_reaction_id_to_index.contains_key(&reaction.id)
        {
            return Err(ChemistryError::DuplicateReaction(reaction.id.to_string()));
        }
        self.validate_dynamic_reaction(&reaction)?;
        let reaction_index = DynamicReactionIndex(self.dynamic_reactions.len());
        let indexed_substances = self.indexed_substances_for_reaction(&reaction);
        index_dynamic_reaction(
            &mut self.dynamic_reaction_index_by_substance,
            &mut self.dynamic_unindexed_reaction_indices,
            self.static_registry.substance_count(),
            reaction_index,
            indexed_substances,
        );
        self.dynamic_reaction_id_to_index
            .insert(reaction.id.clone(), reaction_index);
        self.dynamic_reactions.push(reaction);
        Ok(())
    }

    fn dynamic_reaction_by_index(&self, index: DynamicReactionIndex) -> Option<&Reaction> {
        self.dynamic_reactions.get(index.0)
    }

    fn known_substance_index(&self, substance_id: &SubstanceId) -> Option<KnownSubstanceIndex> {
        self.dynamic_substance_id_to_index
            .get(substance_id)
            .copied()
            .map(KnownSubstanceIndex::Dynamic)
            .or_else(|| {
                self.static_registry
                    .substance_index(substance_id)
                    .map(KnownSubstanceIndex::Static)
            })
    }

    fn indexed_substances_for_reaction(&self, reaction: &Reaction) -> Vec<KnownSubstanceIndex> {
        let mut substances = Vec::new();
        for reactant in &reaction.reactants {
            if let Some(substance_index) = self.known_substance_index(&reactant.substance_id) {
                insert_sorted_unique(&mut substances, substance_index);
            }
        }
        for ordered_substance in reaction.orders.keys() {
            if let Some(substance_index) = self.known_substance_index(ordered_substance) {
                insert_sorted_unique(&mut substances, substance_index);
            }
        }
        substances
    }

    fn validate_dynamic_reaction(&self, reaction: &Reaction) -> ChemistryResult<()> {
        reaction.validate_shape()?;
        for term in reaction
            .reactants
            .iter()
            .chain(reaction.products.iter())
            .chain(
                reaction
                    .product_distribution
                    .iter()
                    .flat_map(|distribution| distribution.variants.iter())
                    .flat_map(|variant| variant.products.iter()),
            )
            .chain(
                reaction
                    .channels
                    .iter()
                    .flat_map(|channel| channel.products.iter()),
            )
        {
            self.substance(&term.substance_id)
                .map_err(|_| ChemistryError::UnknownSubstance {
                    reaction_id: reaction.id.to_string(),
                    substance_id: term.substance_id.to_string(),
                })?;
        }
        for ordered_substance in reaction.orders.keys() {
            self.substance(ordered_substance)
                .map_err(|_| ChemistryError::UnknownSubstance {
                    reaction_id: reaction.id.to_string(),
                    substance_id: ordered_substance.to_string(),
                })?;
        }
        let external_reactant_charge = reaction
            .external_reactants
            .iter()
            .filter_map(|requirement| {
                requirement
                    .charge
                    .map(|charge| charge as f64 * requirement.moles_per_reaction)
            })
            .sum::<f64>();
        let reactant_charge = reaction
            .reactants
            .iter()
            .map(|term| {
                self.substance(&term.substance_id)
                    .map(|substance| substance.charge as f64 * term.coefficient as f64)
            })
            .sum::<ChemistryResult<f64>>()?
            + external_reactant_charge;
        let product_charge = self.dynamic_product_charge(reaction)?;
        if (reactant_charge - product_charge).abs() > 1.0e-9
            && !reaction.allow_charge_imbalance
        {
            return Err(ChemistryError::ChargeNotConserved {
                reaction_id: reaction.id.to_string(),
                reactants: reactant_charge.round() as i32,
                products: product_charge.round() as i32,
            });
        }

        let external_reactant_mass = reaction
            .external_reactants
            .iter()
            .filter_map(|requirement| {
                requirement
                    .molar_mass_grams
                    .map(|mass| mass * requirement.moles_per_reaction)
            })
            .sum::<f64>();
        let reactant_mass = reaction
            .reactants
            .iter()
            .map(|term| {
                self.substance(&term.substance_id)
                    .map(|substance| substance.molar_mass_grams * term.coefficient as f64)
            })
            .sum::<ChemistryResult<f64>>()?
            + external_reactant_mass;
        let product_mass = self.dynamic_product_mass(reaction)?;
        if (reactant_mass - product_mass).abs() > MASS_TOLERANCE_GRAMS_PER_MOL
            && !reaction.allow_mass_imbalance
        {
            return Err(ChemistryError::MassNotConserved {
                reaction_id: reaction.id.to_string(),
                reactants: reactant_mass,
                products: product_mass,
            });
        }
        Ok(())
    }

    fn dynamic_product_charge(&self, reaction: &Reaction) -> ChemistryResult<f64> {
        if !reaction.channels.is_empty() {
            let external_reactant_charge = reaction
                .external_reactants
                .iter()
                .filter_map(|requirement| {
                    requirement
                        .charge
                        .map(|charge| charge as f64 * requirement.moles_per_reaction)
                })
                .sum::<f64>();
            let external_product_charge = reaction
                .external_products
                .iter()
                .filter_map(|requirement| {
                    requirement
                        .charge
                        .map(|charge| charge as f64 * requirement.moles_per_reaction)
                })
                .sum::<f64>();
            let reactant_charge = reaction
                .reactants
                .iter()
                .map(|term| {
                    self.substance(&term.substance_id)
                        .map(|substance| substance.charge as f64 * term.coefficient as f64)
                })
                .sum::<ChemistryResult<f64>>()?
                + external_reactant_charge;
            let mut first_channel_charge = None;
            for channel in &reaction.channels {
                let channel_charge = channel
                    .products
                    .iter()
                    .map(|term| {
                        self.substance(&term.substance_id)
                            .map(|substance| substance.charge as f64 * term.coefficient as f64)
                    })
                    .sum::<ChemistryResult<f64>>()?
                    + external_product_charge;
                if first_channel_charge.is_none() {
                    first_channel_charge = Some(channel_charge);
                }
                if (reactant_charge - channel_charge).abs() > 1.0e-9
                    && !reaction.allow_charge_imbalance
                {
                    return Err(ChemistryError::ChargeNotConserved {
                        reaction_id: format!("{}:{}", reaction.id, channel.id),
                        reactants: reactant_charge.round() as i32,
                        products: channel_charge.round() as i32,
                    });
                }
            }
            return Ok(first_channel_charge.unwrap_or(0.0));
        }
        if let Some(distribution) = &reaction.product_distribution {
            let external_product_charge = reaction
                .external_products
                .iter()
                .filter_map(|requirement| {
                    requirement
                        .charge
                        .map(|charge| charge as f64 * requirement.moles_per_reaction)
                })
                .sum::<f64>();
            return distribution
                .variants
                .iter()
                .map(|variant| {
                    let charge = variant
                        .products
                        .iter()
                        .map(|term| {
                            self.substance(&term.substance_id)
                                .map(|substance| substance.charge as f64 * term.coefficient as f64)
                        })
                        .sum::<ChemistryResult<f64>>()?;
                    Ok((charge + external_product_charge) * variant.fraction)
                })
                .sum();
        }
        let external_product_charge = reaction
            .external_products
            .iter()
            .filter_map(|requirement| {
                requirement
                    .charge
                    .map(|charge| charge as f64 * requirement.moles_per_reaction)
            })
            .sum::<f64>();
        let product_charge = reaction
            .products
            .iter()
            .map(|term| {
                self.substance(&term.substance_id)
                    .map(|substance| substance.charge as f64 * term.coefficient as f64)
            })
            .sum::<ChemistryResult<f64>>()?;
        Ok(product_charge + external_product_charge)
    }

    fn dynamic_product_mass(&self, reaction: &Reaction) -> ChemistryResult<f64> {
        if !reaction.channels.is_empty() {
            let external_reactant_mass = reaction
                .external_reactants
                .iter()
                .filter_map(|requirement| {
                    requirement
                        .molar_mass_grams
                        .map(|mass| mass * requirement.moles_per_reaction)
                })
                .sum::<f64>();
            let external_product_mass = reaction
                .external_products
                .iter()
                .filter_map(|requirement| {
                    requirement
                        .molar_mass_grams
                        .map(|mass| mass * requirement.moles_per_reaction)
                })
                .sum::<f64>();
            let reactant_mass = reaction
                .reactants
                .iter()
                .map(|term| {
                    self.substance(&term.substance_id)
                        .map(|substance| substance.molar_mass_grams * term.coefficient as f64)
                })
                .sum::<ChemistryResult<f64>>()?
                + external_reactant_mass;
            let mut first_channel_mass = None;
            for channel in &reaction.channels {
                let channel_mass = channel
                    .products
                    .iter()
                    .map(|term| {
                        self.substance(&term.substance_id)
                            .map(|substance| substance.molar_mass_grams * term.coefficient as f64)
                    })
                    .sum::<ChemistryResult<f64>>()?
                    + external_product_mass;
                if first_channel_mass.is_none() {
                    first_channel_mass = Some(channel_mass);
                }
                if (reactant_mass - channel_mass).abs() > MASS_TOLERANCE_GRAMS_PER_MOL
                    && !reaction.allow_mass_imbalance
                {
                    return Err(ChemistryError::MassNotConserved {
                        reaction_id: format!("{}:{}", reaction.id, channel.id),
                        reactants: reactant_mass,
                        products: channel_mass,
                    });
                }
            }
            return Ok(first_channel_mass.unwrap_or(0.0));
        }
        if let Some(distribution) = &reaction.product_distribution {
            let external_product_mass = reaction
                .external_products
                .iter()
                .filter_map(|requirement| {
                    requirement
                        .molar_mass_grams
                        .map(|mass| mass * requirement.moles_per_reaction)
                })
                .sum::<f64>();
            return distribution
                .variants
                .iter()
                .map(|variant| {
                    let mass = variant
                        .products
                        .iter()
                        .map(|term| {
                            self.substance(&term.substance_id).map(|substance| {
                                substance.molar_mass_grams * term.coefficient as f64
                            })
                        })
                        .sum::<ChemistryResult<f64>>()?;
                    Ok((mass + external_product_mass) * variant.fraction)
                })
                .sum();
        }
        let external_product_mass = reaction
            .external_products
            .iter()
            .filter_map(|requirement| {
                requirement
                    .molar_mass_grams
                    .map(|mass| mass * requirement.moles_per_reaction)
            })
            .sum::<f64>();
        let product_mass = reaction
            .products
            .iter()
            .map(|term| {
                self.substance(&term.substance_id)
                    .map(|substance| substance.molar_mass_grams * term.coefficient as f64)
            })
            .sum::<ChemistryResult<f64>>()?;
        Ok(product_mass + external_product_mass)
    }

    fn validated_substance_indices(
        &self,
        substance_ids: impl IntoIterator<Item = SubstanceId>,
    ) -> ChemistryResult<Vec<KnownSubstanceIndex>> {
        let mut result = Vec::new();
        for substance_id in substance_ids {
            let index = self.known_substance_index_or_error(&substance_id)?;
            insert_sorted_unique(&mut result, index);
        }
        Ok(result)
    }

    fn generation_mask_updates_for_substance(
        &self,
        substance_index: KnownSubstanceIndex,
    ) -> ChemistryResult<Vec<(usize, ReactiveSiteKey, u64)>> {
        let site_keys = self.generation_site_keys_for_substance(substance_index)?;
        let slot = known_substance_slot(self.static_registry.substance_count(), substance_index);
        let mut updates = Vec::new();
        for site_key in site_keys {
            let current_mask = self
                .processed_generation_masks
                .get(&(slot, site_key.clone()))
                .copied()
                .unwrap_or(0);
            let mut next_mask = current_mask;
            for generator in generators_for_site(&site_key.kind, &site_key.roles) {
                let bit = generator.bit();
                if next_mask & bit == 0 {
                    next_mask |= bit;
                }
            }
            if next_mask != current_mask {
                updates.push((slot, site_key, next_mask));
            }
        }
        Ok(updates)
    }

    fn apply_generation_mask_updates(&mut self, updates: &[(usize, ReactiveSiteKey, u64)]) {
        for (slot, site_key, mask) in updates {
            self.processed_generation_masks
                .insert((*slot, site_key.clone()), *mask);
        }
    }

    fn generation_site_keys_for_substance(
        &self,
        substance_index: KnownSubstanceIndex,
    ) -> ChemistryResult<Vec<ReactiveSiteKey>> {
        let substance = self.substance_by_known_index(substance_index)?;
        let Some(structure) = substance.molecular_structure.as_ref() else {
            return Ok(Vec::new());
        };
        Ok(try_find_reactive_sites(structure)?
            .into_iter()
            .map(|site| site.key())
            .collect())
    }

    fn generate_organic_for_known_substances(
        &self,
        seeds: &[KnownSubstanceIndex],
        scope: &[bool],
    ) -> ChemistryResult<organic::GeneratedOrganicCatalog> {
        let seed_ids = seeds
            .iter()
            .map(|seed| self.known_substance_id(*seed).cloned())
            .collect::<ChemistryResult<std::collections::BTreeSet<_>>>()?;
        let mut scope_ids = std::collections::BTreeSet::new();
        let mut scoped_substances = Vec::new();
        for substance in self.all_known_substance_indices() {
            let slot = known_substance_slot(self.static_registry.substance_count(), substance);
            if scope.get(slot).copied().unwrap_or(false) {
                let substance = self.substance_by_known_index(substance)?;
                scope_ids.insert(substance.id.clone());
                scoped_substances.push(substance);
            }
        }
        let space =
            OrganicGenerationSpace::from_substances_for_scope(scoped_substances, &scope_ids)?;
        organic::generate_organic_reactions_for_seed_substances(
            &space,
            &seed_ids,
            self.canonical_to_id.clone(),
            &crate::chemistry::selectivity::types::SelectivityContext::default(),
        )
    }

    fn add_site_handles_for_substance(
        &mut self,
        substance: KnownSubstanceIndex,
    ) -> ChemistryResult<()> {
        let Some(substance_data) = self.substance_by_known_index(substance).ok() else {
            return Ok(());
        };
        let Some(structure) = substance_data.molecular_structure.as_ref() else {
            return Ok(());
        };
        let site_kinds = try_find_reactive_sites(structure)?
            .into_iter()
            .map(|site| site.kind)
            .collect::<Vec<_>>();
        for (site_index, site_kind) in site_kinds.into_iter().enumerate() {
            let handle = SiteHandle {
                substance,
                site_index,
            };
            let slot = known_substance_slot(self.static_registry.substance_count(), substance);
            if self.site_handles_by_substance.len() <= slot {
                self.site_handles_by_substance.resize(slot + 1, Vec::new());
            }
            self.site_handles_by_substance[slot].push(handle);
            push_site_handle(&mut self.site_index, site_kind, handle);
        }
        Ok(())
    }

    fn all_known_substance_indices(&self) -> Vec<KnownSubstanceIndex> {
        self.static_registry
            .substance_indices()
            .map(KnownSubstanceIndex::Static)
            .chain((0..self.dynamic_substances.len()).map(KnownSubstanceIndex::Dynamic))
            .collect()
    }

    fn known_substance_count(&self) -> usize {
        self.static_registry.substance_count() + self.dynamic_substances.len()
    }

    fn known_substance_index_or_error(
        &self,
        substance_id: &SubstanceId,
    ) -> ChemistryResult<KnownSubstanceIndex> {
        self.known_substance_index(substance_id).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!("unknown substance '{substance_id}'"))
        })
    }

    fn substance_by_known_index(&self, index: KnownSubstanceIndex) -> ChemistryResult<&Substance> {
        match index {
            KnownSubstanceIndex::Static(index) => self.static_registry.substance_by_index(index),
            KnownSubstanceIndex::Dynamic(index) => {
                self.dynamic_substances.get(index).ok_or_else(|| {
                    ChemistryError::InvalidMixtureState(format!(
                        "invalid dynamic substance index {index}"
                    ))
                })
            }
        }
    }

    fn known_substance_id(&self, index: KnownSubstanceIndex) -> ChemistryResult<&SubstanceId> {
        Ok(&self.substance_by_known_index(index)?.id)
    }

    fn mark_known_substance(
        &self,
        marks: &mut Vec<bool>,
        substance: KnownSubstanceIndex,
        value: bool,
    ) {
        let slot = known_substance_slot(self.static_registry.substance_count(), substance);
        if marks.len() <= slot {
            marks.resize(slot + 1, false);
        }
        marks[slot] = value;
    }
}

fn remap_reaction_substances(
    mut reaction: Reaction,
    remap: &BTreeMap<SubstanceId, SubstanceId>,
) -> Reaction {
    if remap.is_empty() {
        return reaction;
    }
    for term in &mut reaction.reactants {
        remap_term(term, remap);
    }
    for term in &mut reaction.products {
        remap_term(term, remap);
    }
    for substance_id in reaction.orders.keys().cloned().collect::<Vec<_>>() {
        if let Some(replacement) = remap.get(&substance_id) {
            if let Some(order) = reaction.orders.remove(&substance_id) {
                reaction.orders.insert(replacement.clone(), order);
            }
        }
    }
    if let Some(distribution) = &mut reaction.product_distribution {
        for variant in &mut distribution.variants {
            for term in &mut variant.products {
                remap_term(term, remap);
            }
        }
    }
    for channel in &mut reaction.channels {
        for term in &mut channel.products {
            remap_term(term, remap);
        }
    }
    let phase_entries = reaction.phase_access.clone();
    for (substance_id, access) in phase_entries {
        if let Some(replacement) = remap.get(&substance_id) {
            reaction.phase_access.remove(&substance_id);
            reaction.phase_access.insert(replacement.clone(), access);
        }
    }
    let product_phase_entries = reaction.product_phases.clone();
    for (substance_id, phase) in product_phase_entries {
        if let Some(replacement) = remap.get(&substance_id) {
            reaction.product_phases.remove(&substance_id);
            reaction.product_phases.insert(replacement.clone(), phase);
        }
    }
    reaction
}

fn remap_term(term: &mut StoichiometricTerm, remap: &BTreeMap<SubstanceId, SubstanceId>) {
    if let Some(replacement) = remap.get(&term.substance_id) {
        term.substance_id = replacement.clone();
    }
}

fn generators_for_site(
    site_kind: &ReactiveSiteKind,
    roles: &[ReactiveRole],
) -> &'static [OrganicGeneratorKind] {
    match site_kind {
        ReactiveSiteKind::Halide => &[
            OrganicGeneratorKind::HalideHydroxideSubstitution,
            OrganicGeneratorKind::HalideAmmoniaSubstitution,
            OrganicGeneratorKind::HalideCyanideSubstitution,
            OrganicGeneratorKind::HalideAmineSubstitution,
            OrganicGeneratorKind::EnolateAlkylation,
        ],
        ReactiveSiteKind::Alcohol => &[
            OrganicGeneratorKind::AlcoholOxidation,
            OrganicGeneratorKind::AlcoholDehydration,
            OrganicGeneratorKind::ThionylChlorideSubstitution,
            OrganicGeneratorKind::AlcoholSilylProtection,
            OrganicGeneratorKind::CarboxylicAcidEsterification,
            OrganicGeneratorKind::AcylChlorideEsterification,
        ],
        ReactiveSiteKind::SilylEther => &[OrganicGeneratorKind::SilylEtherDeprotection],
        ReactiveSiteKind::Acetal | ReactiveSiteKind::Ketal => {
            &[OrganicGeneratorKind::AcetalDeprotection]
        }
        ReactiveSiteKind::Alkoxide => &[OrganicGeneratorKind::AlkoxideProtonation],
        ReactiveSiteKind::Ester => &[OrganicGeneratorKind::ClaisenCondensation],
        ReactiveSiteKind::Nitrile => &[
            OrganicGeneratorKind::NitrileHydrolysis,
            OrganicGeneratorKind::NitrileHydrogenation,
        ],
        ReactiveSiteKind::Nitro => &[OrganicGeneratorKind::NitroHydrogenation],
        ReactiveSiteKind::AcylChloride => &[
            OrganicGeneratorKind::AcylChlorideHydrolysis,
            OrganicGeneratorKind::AcylChlorideEsterification,
        ],
        ReactiveSiteKind::CarboxylicAcid => &[
            OrganicGeneratorKind::AcylChlorideFormation,
            OrganicGeneratorKind::CarboxylicAcidEsterification,
        ],
        ReactiveSiteKind::Aldehyde => &[
            OrganicGeneratorKind::AldehydeOxidation,
            OrganicGeneratorKind::CyanideNucleophilicAddition,
            OrganicGeneratorKind::WolffKishnerReduction,
            OrganicGeneratorKind::OrganometallicCarbonylAddition,
            OrganicGeneratorKind::AldolAddition,
            OrganicGeneratorKind::EnamineFormation,
        ],
        ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => &[
            OrganicGeneratorKind::CyanideNucleophilicAddition,
            OrganicGeneratorKind::WolffKishnerReduction,
            OrganicGeneratorKind::OrganometallicCarbonylAddition,
            OrganicGeneratorKind::AldolAddition,
            OrganicGeneratorKind::EnamineFormation,
        ],
        ReactiveSiteKind::Amide => &[OrganicGeneratorKind::AmideHydrolysis],
        ReactiveSiteKind::PrimaryAmine => &[
            OrganicGeneratorKind::AminePhosgenation,
            OrganicGeneratorKind::AmineBocProtection,
            OrganicGeneratorKind::AmineCbzProtection,
        ],
        ReactiveSiteKind::Phosphine => &[OrganicGeneratorKind::PhosphoniumSaltFormation],
        ReactiveSiteKind::PhosphoniumSalt => &[OrganicGeneratorKind::PhosphoniumYlideFormation],
        ReactiveSiteKind::PhosphorusYlide => &[OrganicGeneratorKind::WittigOlefination],
        ReactiveSiteKind::PhosphonateCarbanion => {
            &[OrganicGeneratorKind::HornerWadsworthEmmonsOlefination]
        }
        ReactiveSiteKind::SulfoneCarbanion => &[OrganicGeneratorKind::JuliaOlefination],
        ReactiveSiteKind::NonTertiaryAmine => &[
            OrganicGeneratorKind::CyanamideAddition,
            OrganicGeneratorKind::HalideAmineSubstitution,
            OrganicGeneratorKind::EnamineFormation,
            OrganicGeneratorKind::AmineBocProtection,
            OrganicGeneratorKind::AmineCbzProtection,
        ],
        ReactiveSiteKind::BocCarbamate => &[OrganicGeneratorKind::BocDeprotection],
        ReactiveSiteKind::CbzCarbamate => &[OrganicGeneratorKind::CbzDeprotection],
        ReactiveSiteKind::Isocyanate => &[OrganicGeneratorKind::IsocyanateHydrolysis],
        ReactiveSiteKind::Borane => &[OrganicGeneratorKind::BoraneOxidation],
        ReactiveSiteKind::BorateEster => &[OrganicGeneratorKind::BorateEsterHydrolysis],
        ReactiveSiteKind::Alkene => &[
            OrganicGeneratorKind::AlkeneChlorination,
            OrganicGeneratorKind::AlkeneChlorohydrination,
            OrganicGeneratorKind::AlkeneHydrolysis,
            OrganicGeneratorKind::AlkeneHydroborationWithBorane,
            OrganicGeneratorKind::AlkeneHydrochlorination,
            OrganicGeneratorKind::AlkeneHydrogenation,
            OrganicGeneratorKind::AlkeneHydroiodination,
            OrganicGeneratorKind::AlkeneIodination,
            OrganicGeneratorKind::MichaelAddition,
        ],
        ReactiveSiteKind::Alkyne => &[
            OrganicGeneratorKind::AlkyneChlorination,
            OrganicGeneratorKind::AlkyneChlorohydrination,
            OrganicGeneratorKind::AlkyneHydrolysis,
            OrganicGeneratorKind::AlkyneHydroborationWithBorane,
            OrganicGeneratorKind::AlkyneHydrochlorination,
            OrganicGeneratorKind::AlkyneHydrogenation,
            OrganicGeneratorKind::AlkyneHydroiodination,
            OrganicGeneratorKind::AlkyneIodination,
        ],
        ReactiveSiteKind::AromaticRing => &[OrganicGeneratorKind::AromaticNitration],
        ReactiveSiteKind::Epoxide => &[OrganicGeneratorKind::EpoxideHydrolysis],
        ReactiveSiteKind::Organomagnesium
        | ReactiveSiteKind::Organolithium
        | ReactiveSiteKind::Organocopper => {
            if roles.contains(&ReactiveRole::Nucleophile) {
                &[OrganicGeneratorKind::OrganometallicCarbonylAddition]
            } else {
                &[]
            }
        }
        ReactiveSiteKind::Enol | ReactiveSiteKind::Enolate => &[
            OrganicGeneratorKind::AldolAddition,
            OrganicGeneratorKind::AlphaHalogenation,
            OrganicGeneratorKind::AldolDehydration,
            OrganicGeneratorKind::EnolateAlkylation,
            OrganicGeneratorKind::MichaelAddition,
            OrganicGeneratorKind::ClaisenCondensation,
        ],
        _ => &[],
    }
}

fn index_dynamic_reaction(
    by_substance: &mut Vec<Vec<DynamicReactionIndex>>,
    unindexed: &mut Vec<DynamicReactionIndex>,
    static_substance_count: usize,
    reaction_index: DynamicReactionIndex,
    substances: Vec<KnownSubstanceIndex>,
) {
    if substances.is_empty() {
        unindexed.push(reaction_index);
        return;
    }

    for substance in substances {
        let slot = known_substance_slot(static_substance_count, substance);
        if by_substance.len() <= slot {
            by_substance.resize(slot + 1, Vec::new());
        }
        by_substance[slot].push(reaction_index);
    }
}

fn known_substance_slot(static_substance_count: usize, substance: KnownSubstanceIndex) -> usize {
    match substance {
        KnownSubstanceIndex::Static(index) => index.as_usize(),
        KnownSubstanceIndex::Dynamic(index) => static_substance_count + index,
    }
}

fn insert_sorted_unique<T: Ord + Copy>(values: &mut Vec<T>, value: T) {
    match values.binary_search(&value) {
        Ok(_) => {}
        Err(index) => values.insert(index, value),
    }
}

fn enqueue_known_substance(
    static_substance_count: usize,
    queue: &mut VecDeque<KnownSubstanceIndex>,
    queued: &mut Vec<bool>,
    substance: KnownSubstanceIndex,
) {
    let slot = known_substance_slot(static_substance_count, substance);
    if queued.len() <= slot {
        queued.resize(slot + 1, false);
    }
    if !queued[slot] {
        queued[slot] = true;
        queue.push_back(substance);
    }
}

fn push_site_handle(
    buckets: &mut Vec<SiteBucket>,
    site_kind: ReactiveSiteKind,
    handle: SiteHandle,
) {
    if let Some(bucket) = buckets
        .iter_mut()
        .find(|bucket| bucket.site_kind == site_kind)
    {
        bucket.handles.push(handle);
    } else {
        buckets.push(SiteBucket {
            site_kind,
            handles: vec![handle],
        });
    }
}

fn mark_dynamic_reaction_candidate(
    seen: &mut [bool],
    result: &mut Vec<DynamicReactionIndex>,
    reaction_index: DynamicReactionIndex,
) {
    if let Some(slot) = seen.get_mut(reaction_index.0) {
        if !*slot {
            *slot = true;
            result.push(reaction_index);
        }
    }
}

fn conjugate_base_for_group(
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    substance_id: &SubstanceId,
) -> ChemistryResult<Option<(MolecularStructure, f64)>> {
    match group.group_type {
        FunctionalGroupType::CarboxylicAcid => {
            let oxygen = group_atom(group, 2, substance_id, "carboxylic acid oxygen")?;
            let proton = group_atom(group, 3, substance_id, "carboxylic acid proton")?;
            deprotonated_structure(structure, oxygen, proton, -1.0).map(|structure| {
                Some((
                    structure,
                    estimated_acid_pka(FunctionalGroupType::CarboxylicAcid),
                ))
            })
        }
        FunctionalGroupType::BoricAcid => {
            let oxygen = group_atom(group, 1, substance_id, "boric acid oxygen")?;
            let proton = group_atom(group, 2, substance_id, "boric acid proton")?;
            deprotonated_structure(structure, oxygen, proton, -1.0).map(|structure| {
                Some((
                    structure,
                    estimated_acid_pka(FunctionalGroupType::BoricAcid),
                ))
            })
        }
        _ => Ok(None),
    }
}

fn group_atom(
    group: &FunctionalGroup,
    index: usize,
    substance_id: &SubstanceId,
    label: &str,
) -> ChemistryResult<usize> {
    group
        .atoms
        .get(index)
        .copied()
        .ok_or_else(|| ChemistryError::InvalidSubstance {
            substance_id: substance_id.to_string(),
            reason: format!("acidic functional group is missing {label}"),
        })
}

fn deprotonated_structure(
    structure: &MolecularStructure,
    charged_atom: usize,
    proton: usize,
    charged_atom_charge: f64,
) -> ChemistryResult<MolecularStructure> {
    let mut editor = MolecularEditor::new(structure);
    editor.replace_atom(
        charged_atom,
        &structure.atoms[charged_atom].element,
        charged_atom_charge,
    )?;
    editor.remove_atom(proton)?;
    editor.finish()
}

fn estimated_acid_pka(group_type: FunctionalGroupType) -> f64 {
    match group_type {
        FunctionalGroupType::CarboxylicAcid => 4.8,
        FunctionalGroupType::BoricAcid => 9.2,
        _ => unreachable!("only acid functional groups have estimated pKa values"),
    }
}

fn build_dynamic_substance(
    canonical_frowns: String,
    structure: MolecularStructure,
) -> ChemistryResult<Substance> {
    if structure.atom_count() >= MAX_DYNAMIC_ATOMS {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: canonical_frowns,
            reason: "dynamic substance has too many atoms".to_string(),
        });
    }
    let summary = structure.summary()?;
    let mut tags = Vec::new();
    if structure.atoms.iter().any(|atom| atom.element == "R") {
        tags.push(SubstanceTagId::from("destroy:hypothetical"));
    }
    if structure
        .stereochemistry
        .iter()
        .any(|stereo| matches!(stereo, Stereochemistry::Mixture { .. }))
    {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: canonical_frowns,
            reason: "stereo mixtures are not substances; generators must distribute products into concrete stereoisomers".to_string(),
        });
    }
    Ok(Substance::new(
        SubstanceId::new(canonical_frowns.clone())?,
        summary.charge,
        summary.molar_mass_grams,
        DEFAULT_DYNAMIC_DENSITY,
        if summary.charge == 0 {
            estimate_dynamic_boiling_point(summary.molar_mass_grams)
        } else {
            f64::MAX
        },
        DEFAULT_DYNAMIC_HEAT_CAPACITY,
        DEFAULT_DYNAMIC_LATENT_HEAT,
    )
    .with_phase_properties(estimate_dynamic_phase_properties(
        summary.charge,
        summary.molar_mass_grams,
        &structure,
    ))
    .with_catalog_metadata(Some(canonical_frowns), None, DEFAULT_DYNAMIC_COLOR, tags)
    .with_molecular_structure(structure))
}

fn estimate_dynamic_boiling_point(molar_mass_grams: f64) -> f64 {
    2.042_598_921_281_41 * molar_mass_grams + 178.176_866_128_713
}

fn estimate_dynamic_phase_properties(
    charge: i32,
    molar_mass_grams: f64,
    structure: &MolecularStructure,
) -> SubstancePhaseProperties {
    let estimate = estimate_dynamic_phase_profile(charge, molar_mass_grams, structure);
    if estimate.ionic {
        return SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: Some(10.0),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: false,
            solvent_role: SolventRole::NotSolvent,
        };
    }

    let can_precipitate = estimate.solid_forming_tendency >= 0.7;
    let can_be_solvent = conservative_dynamic_solvent_candidate(&estimate, molar_mass_grams);
    let solvent_role = if can_be_solvent {
        SolventRole::ConservativePredictedSolvent
    } else {
        SolventRole::NotSolvent
    };
    if estimate.estimated_log_p <= -0.5 || estimate.polarity_score >= 4.0 {
        let organic_solubility = (0.35
            - estimate.polarity_score * 0.04
            - estimate.hydrogen_bond_donor_count as f64 * 0.03)
            .clamp(0.02, 0.35);
        SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: None,
            organic_solubility_mol_per_bucket: Some(organic_solubility),
            can_precipitate,
            can_form_liquid_phase: can_be_solvent,
            solvent_role,
        }
    } else if estimate.estimated_log_p >= 1.0 {
        let aqueous_solubility = (0.08
            + estimate.polarity_score * 0.04
            + estimate.hydrogen_bond_acceptor_count as f64 * 0.02
            - estimate.estimated_log_p * 0.02)
            .clamp(0.005, 0.5);
        SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Organic,
            aqueous_solubility_mol_per_bucket: Some(aqueous_solubility),
            organic_solubility_mol_per_bucket: None,
            can_precipitate,
            can_form_liquid_phase: can_be_solvent,
            solvent_role,
        }
    } else if estimate.carbon_count >= estimate.hetero_atom_count {
        SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Organic,
            aqueous_solubility_mol_per_bucket: Some(
                (0.2 + estimate.polarity_score * 0.05).clamp(0.05, 0.7),
            ),
            organic_solubility_mol_per_bucket: None,
            can_precipitate,
            can_form_liquid_phase: can_be_solvent,
            solvent_role,
        }
    } else {
        SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: None,
            organic_solubility_mol_per_bucket: Some(0.2),
            can_precipitate,
            can_form_liquid_phase: can_be_solvent,
            solvent_role,
        }
    }
}

fn conservative_dynamic_solvent_candidate(
    estimate: &DynamicPhaseEstimate,
    molar_mass_grams: f64,
) -> bool {
    !estimate.ionic
        && molar_mass_grams <= 150.0
        && estimate.solid_forming_tendency < 0.35
        && estimate.carbon_count > 0
        && estimate.hetero_atom_count <= 4
        && estimate.polarity_score <= 4.0
}

#[derive(Debug, Clone, PartialEq)]
struct DynamicPhaseEstimate {
    polarity_score: f64,
    hydrogen_bond_donor_count: usize,
    hydrogen_bond_acceptor_count: usize,
    carbon_count: usize,
    hetero_atom_count: usize,
    ionic: bool,
    estimated_log_p: f64,
    solid_forming_tendency: f64,
}

fn estimate_dynamic_phase_profile(
    charge: i32,
    molar_mass_grams: f64,
    structure: &MolecularStructure,
) -> DynamicPhaseEstimate {
    let carbon_count = count_atoms(structure, "C");
    let halogen_count = structure
        .atoms
        .iter()
        .filter(|atom| matches!(atom.element.as_str(), "F" | "Cl" | "Br" | "I"))
        .count();
    let hetero_atom_count = structure
        .atoms
        .iter()
        .filter(|atom| !matches!(atom.element.as_str(), "C" | "H" | "R"))
        .count();
    let hydrogen_bond_donor_count = structure
        .atoms
        .iter()
        .enumerate()
        .filter(|(index, atom)| {
            matches!(atom.element.as_str(), "O" | "N" | "S")
                && structure.explicit_hydrogen_count(*index) > 0
        })
        .count();
    let hydrogen_bond_acceptor_count = structure
        .atoms
        .iter()
        .filter(|atom| matches!(atom.element.as_str(), "O" | "N" | "S" | "P"))
        .filter(|atom| atom.charge <= 0.0)
        .count();
    let ionic = charge != 0 || structure.atoms.iter().any(|atom| atom.charge.abs() >= 0.5);
    let polarity_score = (charge.abs() as f64 * 8.0)
        + hydrogen_bond_donor_count as f64 * 1.4
        + hydrogen_bond_acceptor_count as f64 * 0.9
        + hetero_atom_count as f64 * 0.45
        + halogen_count as f64 * 0.15
        - carbon_count as f64 * 0.18;
    let estimated_log_p = carbon_count as f64 * 0.52 + halogen_count as f64 * 0.28
        - hetero_atom_count as f64 * 0.42
        - hydrogen_bond_donor_count as f64 * 0.9
        - hydrogen_bond_acceptor_count as f64 * 0.45
        - charge.abs() as f64 * 6.0;
    let heavy_atom_count = structure
        .atoms
        .iter()
        .filter(|atom| !matches!(atom.element.as_str(), "H" | "R"))
        .count();
    let solid_forming_tendency = ((molar_mass_grams - 160.0) / 220.0
        + heavy_atom_count.saturating_sub(12) as f64 * 0.04
        + hydrogen_bond_donor_count as f64 * 0.08
        + hydrogen_bond_acceptor_count as f64 * 0.04)
        .clamp(0.0, 1.0);

    DynamicPhaseEstimate {
        polarity_score,
        hydrogen_bond_donor_count,
        hydrogen_bond_acceptor_count,
        carbon_count,
        hetero_atom_count,
        ionic,
        estimated_log_p,
        solid_forming_tendency,
    }
}

fn count_atoms(structure: &MolecularStructure, element: &str) -> usize {
    structure
        .atoms
        .iter()
        .filter(|atom| atom.element == element)
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::mixture::Mixture;
    use crate::chemistry::simulation::react_for_tick;

    #[test]
    fn generation_tracking_uses_dense_generator_masks() {
        assert!(std::mem::size_of::<OrganicGeneratorKind>() <= 1);
        assert!(OrganicGeneratorKind::AlkyneIodination.ordinal() < u64::BITS);
    }

    #[test]
    fn resolving_dynamic_substance_updates_site_index() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let id = registry.resolve_frowns("CCCC=C").unwrap();
        let index = registry.known_substance_index(&id).unwrap();

        assert!(registry.site_index.iter().any(|bucket| bucket.site_kind
            == ReactiveSiteKind::Alkene
            && bucket
                .handles
                .iter()
                .any(|handle| handle.substance == index)));
    }

    #[test]
    fn resolves_known_substance_by_frowns() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acetone = registry.resolve_frowns("CC(=O)C").unwrap();
        assert_eq!(acetone.as_str(), "destroy:acetone");
    }

    #[test]
    fn creates_stable_dynamic_substance_without_duplicates() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let first = registry.resolve_frowns("CCCCCCCC").unwrap();
        let second = registry.resolve_frowns("CCCCCCCC").unwrap();
        assert_eq!(first, second);
        assert!(first.as_str().starts_with("destroy:linear:"));
    }

    #[test]
    fn dynamic_substances_distinguish_stereoisomers() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let clockwise = registry
            .resolve_frowns(
                "destroy:graph:atoms=C.H.Cl.F.I;bonds=0-s-1,0-s-2,0-s-3,0-s-4;stereo=t:0:1.2.3.4:cw",
            )
            .unwrap();
        let repeated = registry
            .resolve_frowns(
                "destroy:graph:atoms=C.H.Cl.F.I;bonds=0-s-1,0-s-2,0-s-3,0-s-4;stereo=t:0:1.2.3.4:cw",
            )
            .unwrap();
        let counter_clockwise = registry
            .resolve_frowns(
                "destroy:graph:atoms=C.H.Cl.F.I;bonds=0-s-1,0-s-2,0-s-3,0-s-4;stereo=t:0:1.2.3.4:ccw",
            )
            .unwrap();

        assert_eq!(clockwise, repeated);
        assert_ne!(clockwise, counter_clockwise);
    }

    #[test]
    fn dynamic_stereo_mixture_is_not_a_substance() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let error = registry
            .resolve_frowns(
                "destroy:graph:atoms=C.H.Cl.F.I;bonds=0-s-1,0-s-2,0-s-3,0-s-4;stereo=mix:tetra:0.1.2.3.4",
            )
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }

    #[test]
    fn canonical_codes_are_cached_for_static_and_dynamic_substances() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        assert_eq!(
            registry.canonical_to_id["destroy:linear:C(=O)(C)(C)"].as_str(),
            "destroy:acetone"
        );
        assert_eq!(
            registry.canonical_by_id[&SubstanceId::from("destroy:acetone")],
            "destroy:linear:C(=O)(C)(C)"
        );

        let dynamic = registry.resolve_frowns("CCCCCCCC").unwrap();
        let cached = registry.canonical_by_id[&dynamic].clone();
        let repeated = registry.resolve_frowns("CCCCCCCC").unwrap();

        assert_eq!(dynamic, repeated);
        assert_eq!(registry.canonical_to_id[&cached], dynamic);
    }

    #[test]
    fn hypothetical_dynamic_substance_cannot_enter_mixture() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let id = registry.resolve_frowns("R1C").unwrap();
        let substance = registry.substance(&id).unwrap();
        assert!(substance
            .tags
            .iter()
            .any(|tag| tag.as_str() == "destroy:hypothetical"));
        assert!(registry.validate_substance_can_enter_mixture(&id).is_err());
    }

    #[test]
    fn dynamic_substance_gets_phase_properties_from_structure() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let hydrocarbon = registry.resolve_frowns("CCCCCCCC").unwrap();
        let substance = registry.substance(&hydrocarbon).unwrap();

        assert_eq!(
            substance.phase_properties.preferred_liquid_phase,
            LiquidPhasePreference::Organic
        );
        assert_eq!(
            substance.phase_properties.aqueous_solubility_mol_per_bucket,
            Some(0.005)
        );
        assert_eq!(
            substance.phase_properties.organic_solubility_mol_per_bucket,
            None
        );
        assert!(substance.phase_properties.can_form_liquid_phase);
        assert_eq!(
            substance.phase_properties.solvent_role,
            SolventRole::ConservativePredictedSolvent
        );
    }

    #[test]
    fn dynamic_phase_profile_distinguishes_hydrocarbon_alcohol_and_sugar_like_molecule() {
        let octane = parse_frowns("CCCCCCCC").unwrap();
        let octane_summary = octane.summary().unwrap();
        let octane_estimate =
            estimate_dynamic_phase_profile(0, octane_summary.molar_mass_grams, &octane);

        let ethanol = parse_frowns("CCO").unwrap();
        let ethanol_summary = ethanol.summary().unwrap();
        let ethanol_estimate =
            estimate_dynamic_phase_profile(0, ethanol_summary.molar_mass_grams, &ethanol);

        let glycerol_like = parse_frowns("C(CO)(O)O").unwrap();
        let glycerol_summary = glycerol_like.summary().unwrap();
        let glycerol_estimate =
            estimate_dynamic_phase_profile(0, glycerol_summary.molar_mass_grams, &glycerol_like);

        assert!(octane_estimate.estimated_log_p > ethanol_estimate.estimated_log_p);
        assert!(ethanol_estimate.polarity_score > octane_estimate.polarity_score);
        assert!(glycerol_estimate.polarity_score > ethanol_estimate.polarity_score);
        assert_eq!(glycerol_estimate.hydrogen_bond_donor_count, 3);
        assert_eq!(glycerol_estimate.hydrogen_bond_acceptor_count, 3);
    }

    #[test]
    fn dynamic_phase_properties_make_polar_neutral_molecule_aqueous() {
        let structure = parse_frowns("C(CO)(O)O").unwrap();
        let summary = structure.summary().unwrap();
        let properties =
            estimate_dynamic_phase_properties(summary.charge, summary.molar_mass_grams, &structure);

        assert_eq!(
            properties.preferred_liquid_phase,
            LiquidPhasePreference::Aqueous
        );
        assert_eq!(properties.aqueous_solubility_mol_per_bucket, None);
        assert!(properties
            .organic_solubility_mol_per_bucket
            .is_some_and(|value| value < 0.2));
        assert_eq!(properties.solvent_role, SolventRole::NotSolvent);
    }

    #[test]
    fn dynamic_ethanol_like_molecule_is_conservative_predicted_solvent() {
        let structure = parse_frowns("CCO").unwrap();
        let summary = structure.summary().unwrap();
        let properties =
            estimate_dynamic_phase_properties(summary.charge, summary.molar_mass_grams, &structure);

        assert!(properties.can_form_liquid_phase);
        assert_eq!(
            properties.solvent_role,
            SolventRole::ConservativePredictedSolvent
        );
    }

    #[test]
    fn dynamic_phase_properties_allow_large_neutral_molecules_to_precipitate() {
        let structure = parse_frowns("CCCCCCCCCCCCCCCCCCCCO").unwrap();
        let summary = structure.summary().unwrap();
        let properties =
            estimate_dynamic_phase_properties(summary.charge, summary.molar_mass_grams, &structure);

        assert_eq!(
            properties.preferred_liquid_phase,
            LiquidPhasePreference::Organic
        );
        assert!(properties.can_precipitate);
        assert_eq!(properties.solvent_role, SolventRole::NotSolvent);
        assert!(properties
            .aqueous_solubility_mol_per_bucket
            .is_some_and(|value| value < 0.1));
    }

    #[test]
    fn dynamic_carboxylic_acid_creates_conjugate_base_and_acid_equilibrium() {
        let mut dynamic = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acid = dynamic.resolve_frowns("CCC(=O)O").unwrap();
        let spec = dynamic
            .dynamic_acid_base_specs
            .iter()
            .find(|spec| spec.acid == acid)
            .expect("dynamic carboxylic acid must register acidity")
            .clone();
        let base = dynamic.substance(&spec.conjugate_base).unwrap();

        assert_eq!(base.charge, -1);
        assert_ne!(spec.conjugate_base, acid);
        assert!((spec.pka - 4.8).abs() < 1.0e-9);

        let registry = dynamic.to_registry().unwrap();
        assert!(registry
            .acid_base_specs()
            .any(|registered| registered.acid == acid
                && registered.conjugate_base == spec.conjugate_base));

        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        mixture.add_substance(&registry, acid.clone(), 0.1).unwrap();
        react_for_tick(&registry, &mut mixture, 1).unwrap();

        assert!(mixture.ph(&registry).unwrap().is_some_and(|ph| ph < 7.0));
        assert!(mixture.concentration_of(&spec.conjugate_base) > 0.0);
    }

    #[test]
    fn dynamic_products_are_available_for_later_generation() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let alkene = registry.resolve_frowns("CCCC=C").unwrap();
        let report = registry.generate_reactions_for(&alkene, 2).unwrap();
        assert!(report.added_substances > 0);
        assert!(report.added_reactions > 0);
        assert!(report.processed_work_items > 0);
        assert!(registry.reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("alkene_hydrogenation/destroy_linear_CCCC_C/")
                || reaction.id.as_str().starts_with("alkene_hydrogenation/")
        }));
    }

    #[test]
    fn repeated_dynamic_generation_reaches_fixed_point_without_duplicates() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let alkene = registry.resolve_frowns("CCCC=C").unwrap();
        let first = registry.generate_reactions_for(&alkene, 1).unwrap();
        let second = registry.generate_reactions_for(&alkene, 1).unwrap();

        assert!(first.added_substances > 0);
        assert!(second.reached_fixed_point);
        assert_eq!(second.added_substances, 0);
        assert_eq!(second.added_reactions, 0);
        assert!(second.skipped_duplicates > 0 || second.processed_work_items > 0);
    }

    #[test]
    fn dynamic_generation_can_run_to_fixed_point_without_iteration_limit() {
        let mut registry = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let methane = registry.resolve_frowns("C").unwrap();
        let first = registry
            .generate_reactions_for_to_fixed_point(&methane)
            .unwrap();
        let second = registry
            .generate_reactions_for_to_fixed_point(&methane)
            .unwrap();

        assert!(first.reached_fixed_point);
        assert_eq!(first.remaining_queue, 0);
        assert!(first.generator_errors.is_empty());
        assert!(second.reached_fixed_point);
        assert_eq!(second.added_substances, 0);
        assert_eq!(second.added_reactions, 0);
        assert_eq!(second.processed_work_items, 0);
    }

    #[test]
    fn pair_generators_use_generation_scope_not_whole_registry() {
        let mut acid_only = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acetic_acid = acid_only.resolve_frowns("CC(=O)O").unwrap();
        acid_only.generate_reactions_for(&acetic_acid, 1).unwrap();
        assert!(!acid_only.reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("carboxylic_acid_esterification/destroy_acetic_acid/destroy_ethanol/")
        }));

        let mut acid_and_alcohol = DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let acetic_acid = acid_and_alcohol.resolve_frowns("CC(=O)O").unwrap();
        let ethanol = acid_and_alcohol.resolve_frowns("CCO").unwrap();
        acid_and_alcohol
            .generate_reactions_for_substances([acetic_acid, ethanol], 1)
            .unwrap();
        assert!(acid_and_alcohol.reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("carboxylic_acid_esterification/destroy_acetic_acid/destroy_ethanol/")
        }));
    }
}
