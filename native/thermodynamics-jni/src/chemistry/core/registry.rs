use std::collections::{BTreeMap, BTreeSet};

use super::alloy::AlloyPhaseSnapshot;
use super::catalysis::{CatalystSurfaceId, CatalystSurfaceSpec};
use super::complex::ComplexSpec;
use super::error::{ChemistryError, ChemistryResult};
use super::kinetics::{ChannelConditionEffect, ReactionChannel};
use super::metallurgy::generation::{
    validate_compound_phases, validate_element_data, validate_pair_interactions,
};
use super::metallurgy::{
    metallurgical_state_from_alloy_phase, MetallurgicalCompoundPhaseData, MetallurgicalElementData,
    MetallurgicalPairInteractionData, MetallurgicalState, MetallurgicalSystem,
};
use super::mixture::MixturePhase;
use super::reaction::{ProductDistribution, Reaction, ReactionId, StoichiometricTerm};
use super::redox::{
    build_redox_pair_reaction, validate_half_reaction_conservation, validate_half_reaction_shape,
    validate_redox_annotation, validate_redox_pair, RedoxHalfReaction, RedoxPair, RedoxPairSpec,
};
use super::solution::{
    AcidBaseSpec, EquilibriumSpec, IndexedEquilibrium, IndexedEquilibriumTerm,
    IndexedPrecipitation, PrecipitationSpec,
};
use super::substance::{
    LiquidPhasePreference, MaterialFormulaUnit, SolventRole, Substance, SubstanceAggregateState,
    SubstanceId, SubstanceRepresentation, SubstanceTagId,
};

const MASS_TOLERANCE_GRAMS_PER_MOL: f64 = 1.0e-6;
const THERMO_TOLERANCE: f64 = 1.0e-6;

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct SubstanceIndex(usize);

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct ReactionIndex(usize);

impl SubstanceIndex {
    pub(crate) fn as_usize(self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone, Default)]
pub(crate) struct ReactionCandidateScratch {
    marks: Vec<u32>,
    generation: u32,
    candidates: Vec<ReactionIndex>,
}

impl ReactionCandidateScratch {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn candidates(&self) -> &[ReactionIndex] {
        &self.candidates
    }

    fn prepare(&mut self, reaction_count: usize) {
        if self.marks.len() < reaction_count {
            self.marks.resize(reaction_count, 0);
        }
        self.candidates.clear();
        if self.generation == u32::MAX {
            self.marks.fill(0);
            self.generation = 1;
        } else {
            self.generation += 1;
        }
    }
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedStoichiometricTerm {
    pub substance: SubstanceIndex,
    pub coefficient: u32,
    pub phases: Vec<MixturePhase>,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedProductDistributionVariant {
    pub fraction: f64,
    pub products: Vec<IndexedStoichiometricTerm>,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedReactionChannel {
    pub channel: ReactionChannel,
    pub products: Vec<IndexedStoichiometricTerm>,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedReaction {
    pub reaction: ReactionIndex,
    pub reactants: Vec<IndexedStoichiometricTerm>,
    pub products: Vec<IndexedStoichiometricTerm>,
    pub product_distribution: Option<Vec<IndexedProductDistributionVariant>>,
    pub channels: Vec<IndexedReactionChannel>,
    pub orders: Vec<(SubstanceIndex, u32, Vec<MixturePhase>)>,
}

#[derive(Debug, Clone)]
pub(crate) struct SubstancePropertiesTable {
    pub charge: Vec<i32>,
    pub molar_mass_grams: Vec<f64>,
    pub liquid_density_grams_per_bucket: Vec<f64>,
    pub solid_density_grams_per_bucket: Vec<f64>,
    pub melting_point_kelvin: Vec<f64>,
    pub boiling_point_kelvin: Vec<f64>,
    pub molar_heat_capacity_j_per_mol_kelvin: Vec<f64>,
    pub fusion_heat_j_per_mol: Vec<f64>,
    pub latent_heat_j_per_mol: Vec<f64>,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum SolventMiscibility {
    FullyMiscible,
    PartiallyMiscible { limit_mol_per_bucket: f64 },
    Immiscible,
}

#[derive(Debug, Clone, PartialEq)]
pub enum GasSolubilityModel {
    Henry {
        henry_mol_per_bucket_pascal: f64,
        temperature_kelvin: f64,
        salting_out_coefficient: f64,
        transfer_coefficient_per_tick: f64,
        estimated: bool,
    },
}

#[derive(Debug, Clone)]
pub struct ChemistryRegistry {
    substances_by_index: Vec<Substance>,
    substance_id_to_index: BTreeMap<SubstanceId, SubstanceIndex>,
    reactions_by_index: Vec<Reaction>,
    reaction_id_to_index: BTreeMap<ReactionId, ReactionIndex>,
    indexed_reactions: Vec<IndexedReaction>,
    reaction_index_by_substance: Vec<Vec<ReactionIndex>>,
    unindexed_reaction_indices: Vec<ReactionIndex>,
    substance_properties: SubstancePropertiesTable,
    substance_tags: BTreeSet<SubstanceTagId>,
    solvent_miscibility: BTreeMap<(SubstanceIndex, SubstanceIndex), SolventMiscibility>,
    gas_solubility: BTreeMap<SubstanceIndex, GasSolubilityModel>,
    acid_base_specs: Vec<AcidBaseSpec>,
    indexed_equilibria: Vec<IndexedEquilibrium>,
    indexed_precipitations: Vec<IndexedPrecipitation>,
    redox_half_reactions: BTreeMap<String, RedoxHalfReaction>,
    catalyst_surface_specs: BTreeMap<CatalystSurfaceId, CatalystSurfaceSpec>,
    complex_specs: Vec<ComplexSpec>,
    metallurgical_systems: Vec<MetallurgicalSystem>,
    metallurgical_element_data: Vec<MetallurgicalElementData>,
    metallurgical_pair_interactions: Vec<MetallurgicalPairInteractionData>,
    metallurgical_compound_phases: Vec<MetallurgicalCompoundPhaseData>,
}

impl ChemistryRegistry {
    pub fn substance(&self, id: &SubstanceId) -> ChemistryResult<&Substance> {
        self.substance_index(id)
            .and_then(|index| self.substance_by_index(index).ok())
            .ok_or_else(|| ChemistryError::InvalidMixtureState(format!("unknown substance '{id}'")))
    }

    pub fn reaction(&self, id: &ReactionId) -> ChemistryResult<&Reaction> {
        self.reaction_index(id)
            .and_then(|index| self.reaction_by_index(index).ok())
            .ok_or_else(|| ChemistryError::UnknownReaction(id.to_string()))
    }

    pub fn reactions(&self) -> impl Iterator<Item = &Reaction> {
        self.reactions_by_index.iter()
    }

    pub fn reaction_candidates_for_substances<'registry, 'substances, I>(
        &'registry self,
        substances: I,
    ) -> Vec<&'registry Reaction>
    where
        I: IntoIterator<Item = &'substances SubstanceId>,
    {
        let substance_indices = substances
            .into_iter()
            .filter_map(|substance_id| self.substance_index(substance_id));
        self.reaction_candidate_indices_for_substance_indices(substance_indices)
            .into_iter()
            .filter_map(|reaction_index| self.reaction_by_index(reaction_index).ok())
            .collect()
    }

    pub(crate) fn reaction_candidate_indices_for_substance_indices<I>(
        &self,
        substances: I,
    ) -> Vec<ReactionIndex>
    where
        I: IntoIterator<Item = SubstanceIndex>,
    {
        let mut scratch = ReactionCandidateScratch::new();
        self.collect_reaction_candidate_indices_for_substance_indices(substances, &mut scratch);
        scratch.candidates
    }

    pub(crate) fn collect_reaction_candidate_indices_for_substance_indices<I>(
        &self,
        substances: I,
        scratch: &mut ReactionCandidateScratch,
    ) where
        I: IntoIterator<Item = SubstanceIndex>,
    {
        scratch.prepare(self.reactions_by_index.len());
        for reaction_index in &self.unindexed_reaction_indices {
            mark_reaction_candidate(scratch, *reaction_index);
        }
        for substance_index in substances {
            if let Some(indexed_reactions) = self.reaction_index_by_substance.get(substance_index.0)
            {
                for reaction_index in indexed_reactions {
                    mark_reaction_candidate(scratch, *reaction_index);
                }
            }
        }
    }

    pub fn substances(&self) -> impl Iterator<Item = &Substance> {
        self.substances_by_index.iter()
    }

    pub fn substance_tags(&self) -> impl Iterator<Item = &SubstanceTagId> {
        self.substance_tags.iter()
    }

    pub fn substance_count(&self) -> usize {
        self.substances_by_index.len()
    }

    pub(crate) fn substance_indices(&self) -> impl Iterator<Item = SubstanceIndex> + '_ {
        (0..self.substances_by_index.len()).map(SubstanceIndex)
    }

    pub fn has_substance_tag(&self, id: &SubstanceTagId) -> bool {
        self.substance_tags.contains(id)
    }

    pub(crate) fn substance_index(&self, id: &SubstanceId) -> Option<SubstanceIndex> {
        self.substance_id_to_index.get(id).copied()
    }

    pub(crate) fn reaction_index(&self, id: &ReactionId) -> Option<ReactionIndex> {
        self.reaction_id_to_index.get(id).copied()
    }

    pub(crate) fn substance_by_index(&self, index: SubstanceIndex) -> ChemistryResult<&Substance> {
        self.substances_by_index.get(index.0).ok_or_else(|| {
            ChemistryError::InvalidMixtureState(format!("invalid substance index {}", index.0))
        })
    }

    pub(crate) fn reaction_by_index(&self, index: ReactionIndex) -> ChemistryResult<&Reaction> {
        self.reactions_by_index
            .get(index.0)
            .ok_or_else(|| ChemistryError::UnknownReaction(format!("<reaction-index:{}>", index.0)))
    }

    pub(crate) fn indexed_reaction(
        &self,
        index: ReactionIndex,
    ) -> ChemistryResult<&IndexedReaction> {
        self.indexed_reactions
            .get(index.0)
            .ok_or_else(|| ChemistryError::UnknownReaction(format!("<reaction-index:{}>", index.0)))
    }

    pub(crate) fn substance_properties(&self) -> &SubstancePropertiesTable {
        &self.substance_properties
    }

    pub(crate) fn solvent_miscibility(
        &self,
        left: SubstanceIndex,
        right: SubstanceIndex,
    ) -> SolventMiscibility {
        if left == right {
            return SolventMiscibility::FullyMiscible;
        }
        if self.same_default_miscible_material_melt(left, right) {
            return SolventMiscibility::FullyMiscible;
        }
        ordered_pair(left, right)
            .and_then(|key| self.solvent_miscibility.get(&key).copied())
            .unwrap_or(SolventMiscibility::Immiscible)
    }

    fn same_default_miscible_material_melt(
        &self,
        left: SubstanceIndex,
        right: SubstanceIndex,
    ) -> bool {
        let (Ok(left), Ok(right)) = (
            self.substance_by_index(left),
            self.substance_by_index(right),
        ) else {
            return false;
        };
        matches!(
            (
                left.phase_properties.preferred_liquid_phase,
                right.phase_properties.preferred_liquid_phase,
            ),
            (
                LiquidPhasePreference::MoltenMetal,
                LiquidPhasePreference::MoltenMetal
            ) | (
                LiquidPhasePreference::MoltenSlag,
                LiquidPhasePreference::MoltenSlag
            )
        ) && left.phase_properties.can_form_liquid_phase
            && right.phase_properties.can_form_liquid_phase
            && left.phase_properties.solvent_role == SolventRole::NotSolvent
            && right.phase_properties.solvent_role == SolventRole::NotSolvent
    }

    pub(crate) fn gas_solubility(&self, substance: SubstanceIndex) -> Option<&GasSolubilityModel> {
        self.gas_solubility.get(&substance)
    }

    pub fn acid_base_specs(&self) -> impl Iterator<Item = &AcidBaseSpec> {
        self.acid_base_specs.iter()
    }

    pub(crate) fn indexed_equilibria(&self) -> &[IndexedEquilibrium] {
        &self.indexed_equilibria
    }

    pub(crate) fn indexed_precipitations(&self) -> &[IndexedPrecipitation] {
        &self.indexed_precipitations
    }

    pub fn redox_half_reactions(&self) -> impl Iterator<Item = &RedoxHalfReaction> {
        self.redox_half_reactions.values()
    }

    pub fn redox_half_reaction(&self, id: &str) -> Option<&RedoxHalfReaction> {
        self.redox_half_reactions.get(id)
    }

    pub fn catalyst_surface_spec(&self, id: &CatalystSurfaceId) -> Option<&CatalystSurfaceSpec> {
        self.catalyst_surface_specs.get(id)
    }

    pub fn catalyst_surface_specs(&self) -> impl Iterator<Item = &CatalystSurfaceSpec> {
        self.catalyst_surface_specs.values()
    }

    pub fn complex_specs(&self) -> impl Iterator<Item = &ComplexSpec> {
        self.complex_specs.iter()
    }

    pub fn metallurgical_systems(&self) -> &[MetallurgicalSystem] {
        &self.metallurgical_systems
    }

    pub fn metallurgical_element_data(&self) -> &[MetallurgicalElementData] {
        &self.metallurgical_element_data
    }

    pub fn metallurgical_pair_interactions(&self) -> &[MetallurgicalPairInteractionData] {
        &self.metallurgical_pair_interactions
    }

    pub fn metallurgical_compound_phases(&self) -> &[MetallurgicalCompoundPhaseData] {
        &self.metallurgical_compound_phases
    }

    pub fn metallurgical_state_from_alloy_phase(
        &self,
        alloy: &AlloyPhaseSnapshot,
        previous: Option<&MetallurgicalState>,
        delta_seconds: f64,
    ) -> ChemistryResult<MetallurgicalState> {
        metallurgical_state_from_alloy_phase(
            alloy,
            &self.metallurgical_systems,
            &self.metallurgical_element_data,
            &self.metallurgical_pair_interactions,
            &self.metallurgical_compound_phases,
            previous,
            delta_seconds,
        )
    }
}

#[derive(Default)]
pub struct ChemistryRegistryBuilder {
    substances: Vec<Substance>,
    reactions: Vec<Reaction>,
    substance_tags: BTreeSet<SubstanceTagId>,
    solvent_miscibility: Vec<(SubstanceId, SubstanceId, SolventMiscibility)>,
    gas_solubility: Vec<(SubstanceId, GasSolubilityModel)>,
    acid_base_specs: Vec<AcidBaseSpec>,
    equilibria: Vec<EquilibriumSpec>,
    precipitations: Vec<PrecipitationSpec>,
    redox_half_reactions: Vec<RedoxHalfReaction>,
    redox_pairs: Vec<RedoxPair>,
    redox_pair_specs: Vec<RedoxPairSpec>,
    catalyst_surface_specs: Vec<CatalystSurfaceSpec>,
    complex_specs: Vec<ComplexSpec>,
    metallurgical_systems: Vec<MetallurgicalSystem>,
    metallurgical_element_data: Vec<MetallurgicalElementData>,
    metallurgical_pair_interactions: Vec<MetallurgicalPairInteractionData>,
    metallurgical_compound_phases: Vec<MetallurgicalCompoundPhaseData>,
}

impl ChemistryRegistryBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_registry(registry: &ChemistryRegistry) -> Self {
        let complex_equilibrium_ids = registry
            .complex_specs
            .iter()
            .map(|spec| format!("{}.formation", spec.id))
            .collect::<BTreeSet<_>>();
        Self {
            substances: registry.substances_by_index.clone(),
            reactions: registry.reactions_by_index.clone(),
            substance_tags: registry.substance_tags.clone(),
            solvent_miscibility: registry
                .solvent_miscibility
                .iter()
                .map(|((left, right), miscibility)| {
                    (
                        registry.substances_by_index[left.0].id.clone(),
                        registry.substances_by_index[right.0].id.clone(),
                        *miscibility,
                    )
                })
                .collect(),
            gas_solubility: registry
                .gas_solubility
                .iter()
                .map(|(substance, model)| {
                    (
                        registry.substances_by_index[substance.0].id.clone(),
                        model.clone(),
                    )
                })
                .collect(),
            acid_base_specs: registry.acid_base_specs.clone(),
            equilibria: registry
                .indexed_equilibria
                .iter()
                .filter(|equilibrium| {
                    !equilibrium.spec.id.ends_with(".acid_base_equilibrium")
                        && !equilibrium.spec.id.ends_with(".neutralization_equilibrium")
                        && !equilibrium
                            .spec
                            .id
                            .ends_with(".base_hydrolysis_equilibrium")
                        && !complex_equilibrium_ids.contains(&equilibrium.spec.id)
                })
                .map(|equilibrium| equilibrium.spec.clone())
                .collect(),
            precipitations: registry
                .indexed_precipitations
                .iter()
                .map(|precipitation| precipitation.spec.clone())
                .collect(),
            redox_half_reactions: registry.redox_half_reactions.values().cloned().collect(),
            redox_pairs: Vec::new(),
            redox_pair_specs: Vec::new(),
            catalyst_surface_specs: registry.catalyst_surface_specs.values().cloned().collect(),
            complex_specs: registry.complex_specs.clone(),
            metallurgical_systems: registry.metallurgical_systems.clone(),
            metallurgical_element_data: registry.metallurgical_element_data.clone(),
            metallurgical_pair_interactions: registry.metallurgical_pair_interactions.clone(),
            metallurgical_compound_phases: registry.metallurgical_compound_phases.clone(),
        }
    }

    pub fn substance(mut self, substance: Substance) -> Self {
        self.substances.push(substance);
        self
    }

    pub fn reaction(mut self, reaction: Reaction) -> Self {
        self.reactions.push(reaction);
        self
    }

    pub fn substance_tag(mut self, tag_id: impl Into<SubstanceTagId>) -> Self {
        self.substance_tags.insert(tag_id.into());
        self
    }

    pub fn solvent_miscibility(
        mut self,
        left: impl Into<SubstanceId>,
        right: impl Into<SubstanceId>,
        miscibility: SolventMiscibility,
    ) -> Self {
        self.solvent_miscibility
            .push((left.into(), right.into(), miscibility));
        self
    }

    pub fn gas_solubility(
        mut self,
        substance: impl Into<SubstanceId>,
        model: GasSolubilityModel,
    ) -> Self {
        self.gas_solubility.push((substance.into(), model));
        self
    }

    pub fn acid_base_pair(mut self, spec: AcidBaseSpec) -> Self {
        self.acid_base_specs.push(spec);
        self
    }

    pub fn equilibrium(mut self, spec: EquilibriumSpec) -> Self {
        self.equilibria.push(spec);
        self
    }

    pub fn precipitation(mut self, spec: PrecipitationSpec) -> Self {
        self.precipitations.push(spec);
        self
    }

    pub fn redox_half_reaction(mut self, half_reaction: RedoxHalfReaction) -> Self {
        self.redox_half_reactions.push(half_reaction);
        self
    }

    pub fn redox_pair(mut self, pair: RedoxPair) -> Self {
        self.redox_pairs.push(pair);
        self
    }

    pub fn redox_pair_from_halves(
        mut self,
        reaction_id: impl Into<ReactionId>,
        oxidation_half_id: impl Into<String>,
        reduction_half_id: impl Into<String>,
    ) -> Self {
        self.redox_pair_specs.push(RedoxPairSpec::new(
            reaction_id,
            oxidation_half_id,
            reduction_half_id,
        ));
        self
    }

    pub fn catalyst_surface_spec(mut self, spec: CatalystSurfaceSpec) -> Self {
        self.catalyst_surface_specs.push(spec);
        self
    }

    pub fn complex_spec(mut self, spec: ComplexSpec) -> Self {
        self.complex_specs.push(spec);
        self
    }

    pub fn metallurgical_system(mut self, system: MetallurgicalSystem) -> Self {
        self.metallurgical_systems.push(system);
        self
    }

    pub fn metallurgical_systems(
        mut self,
        systems: impl IntoIterator<Item = MetallurgicalSystem>,
    ) -> Self {
        self.metallurgical_systems.extend(systems);
        self
    }

    pub fn metallurgical_element_data(mut self, data: MetallurgicalElementData) -> Self {
        self.metallurgical_element_data.push(data);
        self
    }

    pub fn metallurgical_elements(
        mut self,
        data: impl IntoIterator<Item = MetallurgicalElementData>,
    ) -> Self {
        self.metallurgical_element_data.extend(data);
        self
    }

    pub fn metallurgical_pair_interaction(
        mut self,
        interaction: MetallurgicalPairInteractionData,
    ) -> Self {
        self.metallurgical_pair_interactions.push(interaction);
        self
    }

    pub fn metallurgical_pair_interactions(
        mut self,
        interactions: impl IntoIterator<Item = MetallurgicalPairInteractionData>,
    ) -> Self {
        self.metallurgical_pair_interactions.extend(interactions);
        self
    }

    pub fn metallurgical_compound_phase(mut self, phase: MetallurgicalCompoundPhaseData) -> Self {
        self.metallurgical_compound_phases.push(phase);
        self
    }

    pub fn metallurgical_compound_phases(
        mut self,
        phases: impl IntoIterator<Item = MetallurgicalCompoundPhaseData>,
    ) -> Self {
        self.metallurgical_compound_phases.extend(phases);
        self
    }

    pub fn build(self) -> ChemistryResult<ChemistryRegistry> {
        let mut redox_half_reactions = BTreeMap::new();
        for half in self.redox_half_reactions {
            validate_half_reaction_shape(&half)?;
            if redox_half_reactions.insert(half.id.clone(), half).is_some() {
                return Err(ChemistryError::DuplicateReaction(
                    "<redox-half>".to_string(),
                ));
            }
        }
        let mut reactions = self.reactions;
        for pair in self.redox_pairs {
            validate_redox_pair(&pair, &redox_half_reactions)?;
            reactions.push(pair.reaction);
        }
        for spec in self.redox_pair_specs {
            reactions.push(build_redox_pair_reaction(&spec, &redox_half_reactions)?);
        }
        let mut substances = self.substances;
        let mut equilibria = self.equilibria;
        let complex_specs = build_complex_specs(&self.complex_specs, &substances)?;
        for (spec, substance) in &complex_specs {
            if let Some(substance) = substance {
                substances.push(substance.clone());
            }
            equilibria.push(spec.to_equilibrium());
        }

        let catalyst_surface_specs = build_catalyst_surface_specs(&self.catalyst_surface_specs)?;
        let metallurgical_systems = build_metallurgical_systems(self.metallurgical_systems)?;
        let metallurgical_element_data = validate_element_data(self.metallurgical_element_data)?;
        let metallurgical_pair_interactions =
            validate_pair_interactions(self.metallurgical_pair_interactions)?;
        let metallurgical_compound_phases =
            validate_compound_phases(self.metallurgical_compound_phases)?;

        let mut substance_map = BTreeMap::new();
        for substance in substances {
            substance.validate()?;
            let id = substance.id.clone();
            if substance_map.insert(id.clone(), substance).is_some() {
                return Err(ChemistryError::DuplicateSubstance(id.to_string()));
            }
        }

        let mut substances_by_index = Vec::new();
        let mut substance_id_to_index = BTreeMap::new();
        for (id, substance) in substance_map {
            let substance_index = SubstanceIndex(substances_by_index.len());
            substance_id_to_index.insert(id, substance_index);
            substances_by_index.push(substance);
        }
        let substance_properties = build_substance_properties_table(&substances_by_index);

        let mut reaction_map = BTreeMap::new();
        for reaction in reactions {
            reaction.validate_shape()?;
            let id = reaction.id.clone();
            if reaction_map.insert(id.clone(), reaction).is_some() {
                return Err(ChemistryError::DuplicateReaction(id.to_string()));
            }
        }

        let mut reactions_by_index = Vec::new();
        let mut reaction_id_to_index = BTreeMap::new();
        for (id, reaction) in reaction_map {
            let reaction_index = ReactionIndex(reactions_by_index.len());
            reaction_id_to_index.insert(id, reaction_index);
            reactions_by_index.push(reaction);
        }

        let indexed_reactions =
            build_indexed_reactions(&reactions_by_index, &substance_id_to_index)?;
        let (reaction_index_by_substance, unindexed_reaction_indices) =
            build_reaction_index(&indexed_reactions, substances_by_index.len());
        let solvent_miscibility = build_solvent_miscibility(
            &self.solvent_miscibility,
            &substances_by_index,
            &substance_id_to_index,
        )?;
        let gas_solubility = build_gas_solubility(
            &self.gas_solubility,
            &substances_by_index,
            &substance_id_to_index,
        )?;
        let acid_base_specs = validate_acid_base_specs(
            &self.acid_base_specs,
            &substances_by_index,
            &substance_id_to_index,
        )?;
        for spec in &acid_base_specs {
            equilibria.extend(spec.to_equilibria());
        }
        let indexed_equilibria =
            build_indexed_equilibria(&equilibria, &substances_by_index, &substance_id_to_index)?;
        let indexed_precipitations = build_indexed_precipitations(
            &self.precipitations,
            &substances_by_index,
            &substance_id_to_index,
        )?;
        let registry = ChemistryRegistry {
            substances_by_index,
            substance_id_to_index,
            reactions_by_index,
            reaction_id_to_index,
            indexed_reactions,
            reaction_index_by_substance,
            unindexed_reaction_indices,
            substance_properties,
            substance_tags: self.substance_tags,
            solvent_miscibility,
            gas_solubility,
            acid_base_specs,
            indexed_equilibria,
            indexed_precipitations,
            redox_half_reactions,
            catalyst_surface_specs,
            complex_specs: complex_specs.into_iter().map(|(spec, _)| spec).collect(),
            metallurgical_systems,
            metallurgical_element_data,
            metallurgical_pair_interactions,
            metallurgical_compound_phases,
        };
        registry.validate_redox_half_reactions()?;
        registry.validate_substance_tags()?;
        registry.validate_material_representations()?;
        registry.validate_reactions()?;
        Ok(registry)
    }
}

impl ChemistryRegistry {
    fn validate_substance_tags(&self) -> ChemistryResult<()> {
        for substance in &self.substances_by_index {
            for tag in &substance.tags {
                if !self.substance_tags.contains(tag) {
                    return Err(ChemistryError::InvalidSubstance {
                        substance_id: substance.id.to_string(),
                        reason: format!("unknown substance tag '{tag}'"),
                    });
                }
            }
        }
        Ok(())
    }

    fn validate_material_representations(&self) -> ChemistryResult<()> {
        for substance in &self.substances_by_index {
            match &substance.representation {
                SubstanceRepresentation::IonicSolid { formula_units }
                | SubstanceRepresentation::Oxide { formula_units } => {
                    self.validate_material_formula(substance, formula_units, &[])?;
                }
                SubstanceRepresentation::Hydrate {
                    formula_units,
                    water_count,
                } => {
                    let water = SubstanceId::from("destroy:water");
                    let hydrate_water = [MaterialFormulaUnit::new(water, *water_count)];
                    self.validate_material_formula(substance, formula_units, &hydrate_water)?;
                }
                SubstanceRepresentation::Molecular
                | SubstanceRepresentation::Ion { .. }
                | SubstanceRepresentation::Metal { .. }
                | SubstanceRepresentation::MetallurgicalSolute { .. }
                | SubstanceRepresentation::SurfaceMaterial { .. }
                | SubstanceRepresentation::UnspecifiedMaterial { .. } => {}
            }
        }
        Ok(())
    }

    fn validate_material_formula(
        &self,
        material: &Substance,
        formula_units: &[MaterialFormulaUnit],
        extra_units: &[MaterialFormulaUnit],
    ) -> ChemistryResult<()> {
        let mut mass = 0.0;
        let mut charge = 0;
        for unit in formula_units.iter().chain(extra_units.iter()) {
            let component = self.substance(&unit.substance_id)?;
            mass += component.molar_mass_grams * unit.coefficient as f64;
            charge += component.charge * unit.coefficient as i32;
        }
        if charge != material.charge {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: material.id.to_string(),
                reason: format!(
                    "material formula charge {charge} does not match substance charge {}",
                    material.charge
                ),
            });
        }
        if (mass - material.molar_mass_grams).abs() > MASS_TOLERANCE_GRAMS_PER_MOL {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: material.id.to_string(),
                reason: format!(
                    "material formula mass {mass} does not match substance mass {}",
                    material.molar_mass_grams
                ),
            });
        }
        Ok(())
    }

    fn validate_reactions(&self) -> ChemistryResult<()> {
        let known_reaction_ids = self
            .reaction_id_to_index
            .iter()
            .map(|(id, index)| (id.clone(), index.0))
            .collect::<BTreeMap<_, _>>();
        for reaction in &self.reactions_by_index {
            validate_redox_annotation(reaction, &known_reaction_ids, &self.redox_half_reactions)?;
            if let Some(redox) = &reaction.redox {
                for participant in [redox.oxidant.as_ref(), redox.reductant.as_ref()]
                    .into_iter()
                    .flatten()
                {
                    if !participant.as_str().starts_with("external:")
                        && !self.substance_id_to_index.contains_key(participant)
                    {
                        return Err(ChemistryError::UnknownSubstance {
                            reaction_id: reaction.id.to_string(),
                            substance_id: participant.to_string(),
                        });
                    }
                }
            }
            for term in reaction
                .reactants
                .iter()
                .chain(reaction.products.iter())
                .chain(distributed_product_terms(reaction))
                .chain(channel_product_terms(reaction))
            {
                if !self.substance_id_to_index.contains_key(&term.substance_id) {
                    return Err(ChemistryError::UnknownSubstance {
                        reaction_id: reaction.id.to_string(),
                        substance_id: term.substance_id.to_string(),
                    });
                }
            }
            for ordered_substance in reaction.orders.keys() {
                if !self.substance_id_to_index.contains_key(ordered_substance) {
                    return Err(ChemistryError::UnknownSubstance {
                        reaction_id: reaction.id.to_string(),
                        substance_id: ordered_substance.to_string(),
                    });
                }
            }
            for substance_id in reaction
                .phase_access
                .keys()
                .chain(reaction.product_phases.keys())
            {
                if !self.substance_id_to_index.contains_key(substance_id) {
                    return Err(ChemistryError::UnknownSubstance {
                        reaction_id: reaction.id.to_string(),
                        substance_id: substance_id.to_string(),
                    });
                }
            }
            for requirement in &reaction.surface_requirements {
                let Some(surface_spec) = self.catalyst_surface_specs.get(&requirement.surface_id)
                else {
                    return Err(ChemistryError::InvalidReaction {
                        reaction_id: reaction.id.to_string(),
                        reason: format!("unknown catalyst surface '{}'", requirement.surface_id),
                    });
                };
                if !requirement
                    .phases
                    .iter()
                    .any(|phase| surface_spec.accessible_phases.contains(phase))
                {
                    return Err(ChemistryError::InvalidReaction {
                        reaction_id: reaction.id.to_string(),
                        reason: format!(
                            "surface '{}' is not accessible from any required phase",
                            requirement.surface_id
                        ),
                    });
                }
            }
            for step in &reaction.surface_steps {
                if !self.catalyst_surface_specs.contains_key(step.surface_id()) {
                    return Err(ChemistryError::InvalidReaction {
                        reaction_id: reaction.id.to_string(),
                        reason: format!("unknown catalyst surface '{}'", step.surface_id()),
                    });
                }
            }
            for channel in &reaction.channels {
                for effect in &channel.condition_effects {
                    if let ChannelConditionEffect::Surface { surface_id, .. } = effect {
                        if !self.catalyst_surface_specs.contains_key(surface_id) {
                            return Err(ChemistryError::InvalidReaction {
                                reaction_id: reaction.id.to_string(),
                                reason: format!(
                                    "reaction channel '{}' references unknown catalyst surface '{}'",
                                    channel.id, surface_id
                                ),
                            });
                        }
                    }
                }
            }

            for requirement in reaction
                .external_reactants
                .iter()
                .chain(reaction.external_products.iter())
                .chain(reaction.external_catalysts.iter())
            {
                if requirement.description.trim().is_empty() {
                    return Err(ChemistryError::InvalidReaction {
                        reaction_id: reaction.id.to_string(),
                        reason: "external requirements must have a description".to_string(),
                    });
                }
            }
            for result in &reaction.reaction_results {
                if result.description.trim().is_empty() {
                    return Err(ChemistryError::InvalidReaction {
                        reaction_id: reaction.id.to_string(),
                        reason: "reaction results must have a description".to_string(),
                    });
                }
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
                    self.substance_index(&term.substance_id)
                        .map(|index| {
                            self.substance_properties.charge[index.0] as f64
                                * term.coefficient as f64
                        })
                        .ok_or_else(|| ChemistryError::UnknownSubstance {
                            reaction_id: reaction.id.to_string(),
                            substance_id: term.substance_id.to_string(),
                        })
                })
                .sum::<ChemistryResult<f64>>()?
                + external_reactant_charge;
            let product_charge = product_charge(reaction, self)?;
            if (reactant_charge - product_charge).abs() > 1.0e-9 && !reaction.allow_charge_imbalance
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
            let product_mass = product_mass(reaction, self)?;
            if (reactant_mass - product_mass).abs() > MASS_TOLERANCE_GRAMS_PER_MOL
                && !reaction.allow_mass_imbalance
            {
                return Err(ChemistryError::MassNotConserved {
                    reaction_id: reaction.id.to_string(),
                    reactants: reactant_mass,
                    products: product_mass,
                });
            }

            if let Some(reverse_id) = &reaction.reverse_reaction_id {
                let reverse = self.reaction(reverse_id)?;
                if reverse.reverse_reaction_id.as_ref() != Some(&reaction.id) {
                    return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                        reaction_id: reaction.id.to_string(),
                        reverse_id: reverse_id.to_string(),
                        reason: "reverse reaction must point back to the forward reaction"
                            .to_string(),
                    });
                }
                if stoichiometric_map(&reaction.reactants) != stoichiometric_map(&reverse.products)
                    || stoichiometric_map(&reaction.products)
                        != stoichiometric_map(&reverse.reactants)
                {
                    return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                        reaction_id: reaction.id.to_string(),
                        reverse_id: reverse_id.to_string(),
                        reason: "reverse reaction must mirror closed reactants and products"
                            .to_string(),
                    });
                }
                if reaction.requires_uv != reverse.requires_uv {
                    return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                        reaction_id: reaction.id.to_string(),
                        reverse_id: reverse_id.to_string(),
                        reason: "reverse reaction must carry the same UV requirement".to_string(),
                    });
                }
                if (reaction.enthalpy_change_kj_per_mol + reverse.enthalpy_change_kj_per_mol).abs()
                    > THERMO_TOLERANCE
                {
                    return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                        reaction_id: reaction.id.to_string(),
                        reverse_id: reverse_id.to_string(),
                        reason: "enthalpy changes must sum to zero".to_string(),
                    });
                }
                match (&reaction.thermodynamics, &reverse.thermodynamics) {
                    (Some(forward), Some(backward)) => {
                        if (forward.reference_temperature_kelvin
                            - backward.reference_temperature_kelvin)
                            .abs()
                            > THERMO_TOLERANCE
                        {
                            return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                                reaction_id: reaction.id.to_string(),
                                reverse_id: reverse_id.to_string(),
                                reason: "thermodynamic reference temperatures must match"
                                    .to_string(),
                            });
                        }
                        if (forward.gibbs_free_energy_change_kj_per_mol
                            + backward.gibbs_free_energy_change_kj_per_mol)
                            .abs()
                            > THERMO_TOLERANCE
                        {
                            return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                                reaction_id: reaction.id.to_string(),
                                reverse_id: reverse_id.to_string(),
                                reason: "Gibbs free energy changes must sum to zero".to_string(),
                            });
                        }
                    }
                    (None, None) => {}
                    _ => {
                        return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                            reaction_id: reaction.id.to_string(),
                            reverse_id: reverse_id.to_string(),
                            reason: "paired reverse reactions must either both define thermodynamics or both omit it".to_string(),
                        });
                    }
                }
                let expected_reverse_activation =
                    reaction.activation_energy_kj_per_mol - reaction.enthalpy_change_kj_per_mol;
                if (expected_reverse_activation - reverse.activation_energy_kj_per_mol).abs()
                    > THERMO_TOLERANCE
                {
                    return Err(ChemistryError::ReversibleThermodynamicsMismatch {
                        reaction_id: reaction.id.to_string(),
                        reverse_id: reverse_id.to_string(),
                        reason: "activation energies must match Hess relation".to_string(),
                    });
                }
            }
        }
        Ok(())
    }

    fn validate_redox_half_reactions(&self) -> ChemistryResult<()> {
        for half in self.redox_half_reactions.values() {
            validate_half_reaction_conservation(half, |substance_id| {
                let substance = self.substance(substance_id)?;
                Ok((substance.molar_mass_grams, substance.charge))
            })?;
        }
        Ok(())
    }
}

fn stoichiometric_map(terms: &[StoichiometricTerm]) -> BTreeMap<SubstanceId, u32> {
    let mut result = BTreeMap::new();
    for term in terms {
        *result.entry(term.substance_id.clone()).or_insert(0) += term.coefficient;
    }
    result
}

fn distributed_product_terms(reaction: &Reaction) -> impl Iterator<Item = &StoichiometricTerm> {
    reaction
        .product_distribution
        .iter()
        .flat_map(|distribution| distribution.variants.iter())
        .flat_map(|variant| variant.products.iter())
}

fn channel_product_terms(reaction: &Reaction) -> impl Iterator<Item = &StoichiometricTerm> {
    reaction
        .channels
        .iter()
        .flat_map(|channel| channel.products.iter())
}

fn build_metallurgical_systems(
    systems: Vec<MetallurgicalSystem>,
) -> ChemistryResult<Vec<MetallurgicalSystem>> {
    let mut ids = BTreeSet::new();
    let mut result = Vec::new();
    for system in systems {
        system.validate()?;
        if !ids.insert(system.id.clone()) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "metallurgical system '{}' is registered more than once",
                system.id
            )));
        }
        result.push(system);
    }
    Ok(result)
}

fn product_charge(reaction: &Reaction, registry: &ChemistryRegistry) -> ChemistryResult<f64> {
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
                registry
                    .substance_index(&term.substance_id)
                    .map(|index| {
                        registry.substance_properties.charge[index.0] as f64
                            * term.coefficient as f64
                    })
                    .ok_or_else(|| ChemistryError::UnknownSubstance {
                        reaction_id: reaction.id.to_string(),
                        substance_id: term.substance_id.to_string(),
                    })
            })
            .sum::<ChemistryResult<f64>>()?
            + external_reactant_charge;
        let mut first_channel_charge = None;
        for channel in &reaction.channels {
            let channel_charge = channel
                .products
                .iter()
                .map(|term| {
                    registry
                        .substance_index(&term.substance_id)
                        .map(|index| {
                            registry.substance_properties.charge[index.0] as f64
                                * term.coefficient as f64
                        })
                        .ok_or_else(|| ChemistryError::UnknownSubstance {
                            reaction_id: reaction.id.to_string(),
                            substance_id: term.substance_id.to_string(),
                        })
                })
                .sum::<ChemistryResult<f64>>()?
                + external_product_charge;
            if first_channel_charge.is_none() {
                first_channel_charge = Some(channel_charge);
            }
            if (reactant_charge - channel_charge).abs() > 1.0e-9 && !reaction.allow_charge_imbalance
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
                        registry
                            .substance_index(&term.substance_id)
                            .map(|index| {
                                registry.substance_properties.charge[index.0] as f64
                                    * term.coefficient as f64
                            })
                            .ok_or_else(|| ChemistryError::UnknownSubstance {
                                reaction_id: reaction.id.to_string(),
                                substance_id: term.substance_id.to_string(),
                            })
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
            registry
                .substance_index(&term.substance_id)
                .map(|index| {
                    registry.substance_properties.charge[index.0] as f64 * term.coefficient as f64
                })
                .ok_or_else(|| ChemistryError::UnknownSubstance {
                    reaction_id: reaction.id.to_string(),
                    substance_id: term.substance_id.to_string(),
                })
        })
        .sum::<ChemistryResult<f64>>()?;
    Ok(product_charge + external_product_charge)
}

fn product_mass(reaction: &Reaction, registry: &ChemistryRegistry) -> ChemistryResult<f64> {
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
                registry
                    .substance(&term.substance_id)
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
                    registry
                        .substance(&term.substance_id)
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
                        registry
                            .substance(&term.substance_id)
                            .map(|substance| substance.molar_mass_grams * term.coefficient as f64)
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
            registry
                .substance(&term.substance_id)
                .map(|substance| substance.molar_mass_grams * term.coefficient as f64)
        })
        .sum::<ChemistryResult<f64>>()?;
    Ok(product_mass + external_product_mass)
}

fn build_substance_properties_table(substances: &[Substance]) -> SubstancePropertiesTable {
    SubstancePropertiesTable {
        charge: substances
            .iter()
            .map(|substance| substance.charge)
            .collect(),
        molar_mass_grams: substances
            .iter()
            .map(|substance| substance.molar_mass_grams)
            .collect(),
        liquid_density_grams_per_bucket: substances
            .iter()
            .map(|substance| substance.liquid_density_grams_per_bucket)
            .collect(),
        solid_density_grams_per_bucket: substances
            .iter()
            .map(|substance| substance.solid_density_grams_per_bucket)
            .collect(),
        melting_point_kelvin: substances
            .iter()
            .map(|substance| substance.melting_point_kelvin)
            .collect(),
        boiling_point_kelvin: substances
            .iter()
            .map(|substance| substance.boiling_point_kelvin)
            .collect(),
        molar_heat_capacity_j_per_mol_kelvin: substances
            .iter()
            .map(|substance| substance.molar_heat_capacity_j_per_mol_kelvin)
            .collect(),
        fusion_heat_j_per_mol: substances
            .iter()
            .map(|substance| substance.fusion_heat_j_per_mol)
            .collect(),
        latent_heat_j_per_mol: substances
            .iter()
            .map(|substance| substance.latent_heat_j_per_mol)
            .collect(),
    }
}

fn build_complex_specs(
    specs: &[ComplexSpec],
    base_substances: &[Substance],
) -> ChemistryResult<Vec<(ComplexSpec, Option<Substance>)>> {
    let mut substances_by_id = base_substances
        .iter()
        .map(|substance| (substance.id.clone(), substance.clone()))
        .collect::<BTreeMap<_, _>>();
    let mut seen = BTreeSet::new();
    let mut result = Vec::new();
    for spec in specs {
        spec.validate_shape()?;
        if !seen.insert(spec.id.clone()) {
            return Err(ChemistryError::DuplicateSubstance(spec.id.to_string()));
        }
        let central = substances_by_id.get(&spec.central_ion).ok_or_else(|| {
            ChemistryError::UnknownSubstance {
                reaction_id: format!("{}.formation", spec.id),
                substance_id: spec.central_ion.to_string(),
            }
        })?;
        if central.charge == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: format!("{}.formation", spec.id),
                reason: "complex central substance must be an ion with explicit charge".to_string(),
            });
        }
        let mut charge = central.charge;
        let mut mass = central.molar_mass_grams;
        for ligand in &spec.ligands {
            let ligand_substance = substances_by_id.get(&ligand.substance_id).ok_or_else(|| {
                ChemistryError::UnknownSubstance {
                    reaction_id: format!("{}.formation", spec.id),
                    substance_id: ligand.substance_id.to_string(),
                }
            })?;
            let ligand_solubility = match spec.phase {
                MixturePhase::Aqueous => {
                    ligand_substance
                        .phase_properties
                        .aqueous_solubility_mol_per_bucket
                }
                MixturePhase::Organic => {
                    ligand_substance
                        .phase_properties
                        .organic_solubility_mol_per_bucket
                }
                MixturePhase::MoltenMetal
                | MixturePhase::MoltenSlag
                | MixturePhase::Gas
                | MixturePhase::Solid => None,
            };
            if ligand_substance.phase_properties.can_precipitate && ligand_solubility == Some(0.0) {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: format!("{}.formation", spec.id),
                    reason: format!(
                        "complex ligand '{}' is not available as a dissolved ligand",
                        ligand.substance_id
                    ),
                });
            }
            charge += ligand_substance.charge * ligand.count as i32;
            mass += ligand_substance.molar_mass_grams * ligand.count as f64;
        }
        if let Some(existing) = substances_by_id.get(&spec.id) {
            if existing.charge != spec.charge {
                return Err(ChemistryError::ChargeNotConserved {
                    reaction_id: format!("{}.formation", spec.id),
                    reactants: charge,
                    products: existing.charge,
                });
            }
            if (existing.molar_mass_grams - mass).abs() > MASS_TOLERANCE_GRAMS_PER_MOL {
                return Err(ChemistryError::MassNotConserved {
                    reaction_id: format!("{}.formation", spec.id),
                    reactants: mass,
                    products: existing.molar_mass_grams,
                });
            }
            result.push((spec.clone(), None));
        } else {
            let substance = spec.to_substance(mass, charge)?;
            substances_by_id.insert(spec.id.clone(), substance.clone());
            result.push((spec.clone(), Some(substance)));
        }
    }
    Ok(result)
}

fn build_catalyst_surface_specs(
    specs: &[CatalystSurfaceSpec],
) -> ChemistryResult<BTreeMap<CatalystSurfaceId, CatalystSurfaceSpec>> {
    let mut result = BTreeMap::new();
    for spec in specs {
        spec.validate()?;
        if result.insert(spec.id.clone(), spec.clone()).is_some() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.to_string(),
                reason: "duplicate catalyst surface spec".to_string(),
            });
        }
    }
    Ok(result)
}

fn build_solvent_miscibility(
    entries: &[(SubstanceId, SubstanceId, SolventMiscibility)],
    substances: &[Substance],
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<BTreeMap<(SubstanceIndex, SubstanceIndex), SolventMiscibility>> {
    let mut result = BTreeMap::new();
    for (left_id, right_id, miscibility) in entries {
        let left = *substance_id_to_index.get(left_id).ok_or_else(|| {
            ChemistryError::InvalidSubstance {
                substance_id: left_id.to_string(),
                reason: "unknown solvent in miscibility table".to_string(),
            }
        })?;
        let right = *substance_id_to_index.get(right_id).ok_or_else(|| {
            ChemistryError::InvalidSubstance {
                substance_id: right_id.to_string(),
                reason: "unknown solvent in miscibility table".to_string(),
            }
        })?;
        validate_solvent_for_miscibility(&substances[left.0])?;
        validate_solvent_for_miscibility(&substances[right.0])?;
        validate_miscibility(*miscibility, left_id)?;
        let Some(key) = ordered_pair(left, right) else {
            continue;
        };
        if result.insert(key, *miscibility).is_some() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: format!("{left_id}+{right_id}"),
                reason: "duplicate solvent miscibility entry".to_string(),
            });
        }
    }
    Ok(result)
}

fn validate_solvent_for_miscibility(substance: &Substance) -> ChemistryResult<()> {
    if substance.phase_properties.solvent_role == SolventRole::NotSolvent
        || !substance.phase_properties.can_form_liquid_phase
    {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: substance.id.to_string(),
            reason: "miscibility table may only reference liquid-forming solvents".to_string(),
        });
    }
    Ok(())
}

fn validate_miscibility(
    miscibility: SolventMiscibility,
    substance_id: &SubstanceId,
) -> ChemistryResult<()> {
    if let SolventMiscibility::PartiallyMiscible {
        limit_mol_per_bucket,
    } = miscibility
    {
        if !limit_mol_per_bucket.is_finite() || limit_mol_per_bucket < 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "partial miscibility limit must be non-negative and finite".to_string(),
            });
        }
    }
    Ok(())
}

fn build_gas_solubility(
    entries: &[(SubstanceId, GasSolubilityModel)],
    substances: &[Substance],
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<BTreeMap<SubstanceIndex, GasSolubilityModel>> {
    let mut result = BTreeMap::new();
    for (substance_id, model) in entries {
        let substance = *substance_id_to_index.get(substance_id).ok_or_else(|| {
            ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "unknown substance in gas solubility table".to_string(),
            }
        })?;
        if substances[substance.0].aggregate_state_at(298.0)? != SubstanceAggregateState::Gas {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "gas solubility table may only reference substances that are gases near room temperature".to_string(),
            });
        }
        validate_gas_solubility_model(substance_id, model)?;
        if result.insert(substance, model.clone()).is_some() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "duplicate gas solubility entry".to_string(),
            });
        }
    }
    Ok(result)
}

fn validate_acid_base_specs(
    specs: &[AcidBaseSpec],
    substances: &[Substance],
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<Vec<AcidBaseSpec>> {
    let mut seen = BTreeSet::new();
    let mut result = Vec::new();
    for spec in specs {
        if spec.id.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<acid-base>".to_string(),
                reason: "acid-base spec id must not be empty".to_string(),
            });
        }
        if !spec.pka.is_finite() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "pKa must be finite".to_string(),
            });
        }
        if !seen.insert(spec.id.clone()) {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "duplicate acid-base spec".to_string(),
            });
        }
        let acid = *substance_id_to_index.get(&spec.acid).ok_or_else(|| {
            ChemistryError::UnknownSubstance {
                reaction_id: spec.id.clone(),
                substance_id: spec.acid.to_string(),
            }
        })?;
        let base = *substance_id_to_index
            .get(&spec.conjugate_base)
            .ok_or_else(|| ChemistryError::UnknownSubstance {
                reaction_id: spec.id.clone(),
                substance_id: spec.conjugate_base.to_string(),
            })?;
        let proton = *substance_id_to_index.get(&spec.proton).ok_or_else(|| {
            ChemistryError::UnknownSubstance {
                reaction_id: spec.id.clone(),
                substance_id: spec.proton.to_string(),
            }
        })?;
        let acid_charge = substances[acid.0].charge;
        let base_charge = substances[base.0].charge;
        let proton_charge = substances[proton.0].charge;
        if acid_charge != base_charge + proton_charge {
            return Err(ChemistryError::ChargeNotConserved {
                reaction_id: spec.id.clone(),
                reactants: acid_charge,
                products: base_charge + proton_charge,
            });
        }
        result.push(spec.clone());
    }
    Ok(result)
}

fn build_indexed_equilibria(
    specs: &[EquilibriumSpec],
    substances: &[Substance],
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<Vec<IndexedEquilibrium>> {
    let mut seen = BTreeSet::new();
    let mut result = Vec::new();
    for spec in specs {
        validate_equilibrium_spec_shape(spec)?;
        if !seen.insert(spec.id.clone()) {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "duplicate equilibrium spec".to_string(),
            });
        }
        let reactants = index_equilibrium_terms(&spec.id, &spec.reactants, substance_id_to_index)?;
        let products = index_equilibrium_terms(&spec.id, &spec.products, substance_id_to_index)?;
        validate_equilibrium_conservation(spec, substances, &reactants, &products)?;
        result.push(IndexedEquilibrium {
            spec: spec.clone(),
            reactants,
            products,
        });
    }
    Ok(result)
}

fn build_indexed_precipitations(
    specs: &[PrecipitationSpec],
    substances: &[Substance],
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<Vec<IndexedPrecipitation>> {
    let mut seen = BTreeSet::new();
    let mut result = Vec::new();
    for spec in specs {
        validate_precipitation_spec_shape(spec)?;
        if !seen.insert(spec.id.clone()) {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "duplicate precipitation spec".to_string(),
            });
        }
        let solid = *substance_id_to_index.get(&spec.solid).ok_or_else(|| {
            ChemistryError::UnknownSubstance {
                reaction_id: spec.id.clone(),
                substance_id: spec.solid.to_string(),
            }
        })?;
        if !substances[solid.0].phase_properties.can_precipitate {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: spec.solid.to_string(),
                reason: "precipitation solid must be allowed to precipitate".to_string(),
            });
        }
        let ions = index_equilibrium_terms(&spec.id, &spec.ions, substance_id_to_index)?;
        validate_precipitation_conservation(spec, substances, solid, &ions)?;
        result.push(IndexedPrecipitation {
            spec: spec.clone(),
            solid,
            ions,
        });
    }
    Ok(result)
}

fn validate_precipitation_spec_shape(spec: &PrecipitationSpec) -> ChemistryResult<()> {
    if spec.id.trim().is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<precipitation>".to_string(),
            reason: "precipitation id must not be empty".to_string(),
        });
    }
    if spec.ions.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "precipitation must have dissolved ions".to_string(),
        });
    }
    if !spec.solubility_product.is_finite() || spec.solubility_product <= 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "solubility product must be positive and finite".to_string(),
        });
    }
    if !spec.reference_temperature_kelvin.is_finite() || spec.reference_temperature_kelvin <= 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "precipitation reference temperature must be positive and finite".to_string(),
        });
    }
    if !spec.enthalpy_change_kj_per_mol.is_finite() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "precipitation enthalpy change must be finite".to_string(),
        });
    }
    for term in &spec.ions {
        if term.coefficient == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "precipitation coefficients must be greater than zero".to_string(),
            });
        }
        if term.phase != MixturePhase::Aqueous {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "precipitation ions must be aqueous".to_string(),
            });
        }
    }
    Ok(())
}

fn validate_precipitation_conservation(
    spec: &PrecipitationSpec,
    substances: &[Substance],
    solid: SubstanceIndex,
    ions: &[IndexedEquilibriumTerm],
) -> ChemistryResult<()> {
    let ion_charge = ions
        .iter()
        .map(|term| substances[term.substance.0].charge * term.coefficient as i32)
        .sum::<i32>();
    let solid_charge = substances[solid.0].charge;
    if solid_charge != ion_charge {
        return Err(ChemistryError::ChargeNotConserved {
            reaction_id: spec.id.clone(),
            reactants: solid_charge,
            products: ion_charge,
        });
    }
    let ion_mass = ions
        .iter()
        .map(|term| substances[term.substance.0].molar_mass_grams * term.coefficient as f64)
        .sum::<f64>();
    let solid_mass = substances[solid.0].molar_mass_grams;
    if (solid_mass - ion_mass).abs() > MASS_TOLERANCE_GRAMS_PER_MOL {
        return Err(ChemistryError::MassNotConserved {
            reaction_id: spec.id.clone(),
            reactants: solid_mass,
            products: ion_mass,
        });
    }
    Ok(())
}

fn validate_equilibrium_spec_shape(spec: &EquilibriumSpec) -> ChemistryResult<()> {
    if spec.id.trim().is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<equilibrium>".to_string(),
            reason: "equilibrium id must not be empty".to_string(),
        });
    }
    if spec.reactants.is_empty() || spec.products.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "equilibrium must have reactants and products".to_string(),
        });
    }
    if !spec.equilibrium_constant.is_finite() || spec.equilibrium_constant <= 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "equilibrium constant must be positive and finite".to_string(),
        });
    }
    if !spec.reference_temperature_kelvin.is_finite() || spec.reference_temperature_kelvin <= 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "equilibrium reference temperature must be positive and finite".to_string(),
        });
    }
    if !spec.enthalpy_change_kj_per_mol.is_finite() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.id.clone(),
            reason: "equilibrium enthalpy change must be finite".to_string(),
        });
    }
    for term in spec.reactants.iter().chain(spec.products.iter()) {
        if term.coefficient == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "equilibrium coefficients must be greater than zero".to_string(),
            });
        }
        if term.phase == MixturePhase::Gas || term.phase == MixturePhase::Solid {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: spec.id.clone(),
                reason: "solution equilibria may only use liquid phases".to_string(),
            });
        }
    }
    Ok(())
}

fn index_equilibrium_terms(
    equilibrium_id: &str,
    terms: &[super::solution::EquilibriumTerm],
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<Vec<IndexedEquilibriumTerm>> {
    terms
        .iter()
        .map(|term| {
            let substance = *substance_id_to_index
                .get(&term.substance_id)
                .ok_or_else(|| ChemistryError::UnknownSubstance {
                    reaction_id: equilibrium_id.to_string(),
                    substance_id: term.substance_id.to_string(),
                })?;
            Ok(IndexedEquilibriumTerm {
                substance,
                coefficient: term.coefficient,
                phase: term.phase,
            })
        })
        .collect()
}

fn validate_equilibrium_conservation(
    spec: &EquilibriumSpec,
    substances: &[Substance],
    reactants: &[IndexedEquilibriumTerm],
    products: &[IndexedEquilibriumTerm],
) -> ChemistryResult<()> {
    let reactant_charge = reactants
        .iter()
        .map(|term| substances[term.substance.0].charge * term.coefficient as i32)
        .sum::<i32>();
    let product_charge = products
        .iter()
        .map(|term| substances[term.substance.0].charge * term.coefficient as i32)
        .sum::<i32>();
    if reactant_charge != product_charge {
        return Err(ChemistryError::ChargeNotConserved {
            reaction_id: spec.id.clone(),
            reactants: reactant_charge,
            products: product_charge,
        });
    }
    let reactant_mass = reactants
        .iter()
        .map(|term| substances[term.substance.0].molar_mass_grams * term.coefficient as f64)
        .sum::<f64>();
    let product_mass = products
        .iter()
        .map(|term| substances[term.substance.0].molar_mass_grams * term.coefficient as f64)
        .sum::<f64>();
    if (reactant_mass - product_mass).abs() > MASS_TOLERANCE_GRAMS_PER_MOL {
        return Err(ChemistryError::MassNotConserved {
            reaction_id: spec.id.clone(),
            reactants: reactant_mass,
            products: product_mass,
        });
    }
    Ok(())
}

fn validate_gas_solubility_model(
    substance_id: &SubstanceId,
    model: &GasSolubilityModel,
) -> ChemistryResult<()> {
    match model {
        GasSolubilityModel::Henry {
            henry_mol_per_bucket_pascal,
            temperature_kelvin,
            salting_out_coefficient,
            transfer_coefficient_per_tick,
            estimated: _,
        } => {
            if !henry_mol_per_bucket_pascal.is_finite() || *henry_mol_per_bucket_pascal < 0.0 {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: substance_id.to_string(),
                    reason: "Henry constant must be non-negative and finite".to_string(),
                });
            }
            if !temperature_kelvin.is_finite() || *temperature_kelvin <= 0.0 {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: substance_id.to_string(),
                    reason: "Henry reference temperature must be positive and finite".to_string(),
                });
            }
            if !salting_out_coefficient.is_finite() || *salting_out_coefficient < 0.0 {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: substance_id.to_string(),
                    reason: "salting-out coefficient must be non-negative and finite".to_string(),
                });
            }
            if !transfer_coefficient_per_tick.is_finite() || *transfer_coefficient_per_tick < 0.0 {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: substance_id.to_string(),
                    reason: "gas transfer coefficient must be non-negative and finite".to_string(),
                });
            }
        }
    }
    Ok(())
}

fn build_indexed_reactions(
    reactions: &[Reaction],
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<Vec<IndexedReaction>> {
    reactions
        .iter()
        .enumerate()
        .map(|(reaction_index, reaction)| {
            let reaction_index = ReactionIndex(reaction_index);
            let reactants = reaction
                .reactants
                .iter()
                .map(|term| indexed_stoichiometric_term(reaction, term, substance_id_to_index))
                .collect::<ChemistryResult<Vec<_>>>()?;
            let products = reaction
                .products
                .iter()
                .map(|term| indexed_product_term(reaction, term, substance_id_to_index))
                .collect::<ChemistryResult<Vec<_>>>()?;
            let product_distribution = reaction
                .product_distribution
                .as_ref()
                .map(|distribution| {
                    indexed_product_distribution(reaction, distribution, substance_id_to_index)
                })
                .transpose()?;
            let channels = reaction
                .channels
                .iter()
                .map(|channel| {
                    Ok(IndexedReactionChannel {
                        channel: channel.clone(),
                        products: channel
                            .products
                            .iter()
                            .map(|term| indexed_product_term(reaction, term, substance_id_to_index))
                            .collect::<ChemistryResult<Vec<_>>>()?,
                    })
                })
                .collect::<ChemistryResult<Vec<_>>>()?;
            let orders = reaction
                .orders
                .iter()
                .map(|(substance_id, order)| {
                    let substance = substance_id_to_index
                        .get(substance_id)
                        .copied()
                        .ok_or_else(|| ChemistryError::UnknownSubstance {
                            reaction_id: reaction.id.to_string(),
                            substance_id: substance_id.to_string(),
                        })?;
                    Ok((
                        substance,
                        *order,
                        reaction
                            .phase_access
                            .get(substance_id)
                            .cloned()
                            .unwrap_or_else(super::reaction::ReactionPhaseAccess::liquid)
                            .phases,
                    ))
                })
                .collect::<ChemistryResult<Vec<_>>>()?;
            Ok(IndexedReaction {
                reaction: reaction_index,
                reactants,
                products,
                product_distribution,
                channels,
                orders,
            })
        })
        .collect()
}

fn indexed_product_distribution(
    reaction: &Reaction,
    distribution: &ProductDistribution,
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<Vec<IndexedProductDistributionVariant>> {
    distribution
        .variants
        .iter()
        .map(|variant| {
            Ok(IndexedProductDistributionVariant {
                fraction: variant.fraction,
                products: variant
                    .products
                    .iter()
                    .map(|term| indexed_product_term(reaction, term, substance_id_to_index))
                    .collect::<ChemistryResult<Vec<_>>>()?,
            })
        })
        .collect()
}

fn indexed_stoichiometric_term(
    reaction: &Reaction,
    term: &StoichiometricTerm,
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<IndexedStoichiometricTerm> {
    let substance = substance_id_to_index
        .get(&term.substance_id)
        .copied()
        .ok_or_else(|| ChemistryError::UnknownSubstance {
            reaction_id: reaction.id.to_string(),
            substance_id: term.substance_id.to_string(),
        })?;
    Ok(IndexedStoichiometricTerm {
        substance,
        coefficient: term.coefficient,
        phases: reaction
            .phase_access
            .get(&term.substance_id)
            .cloned()
            .unwrap_or_else(super::reaction::ReactionPhaseAccess::liquid)
            .phases,
    })
}

fn indexed_product_term(
    reaction: &Reaction,
    term: &StoichiometricTerm,
    substance_id_to_index: &BTreeMap<SubstanceId, SubstanceIndex>,
) -> ChemistryResult<IndexedStoichiometricTerm> {
    let substance = substance_id_to_index
        .get(&term.substance_id)
        .copied()
        .ok_or_else(|| ChemistryError::UnknownSubstance {
            reaction_id: reaction.id.to_string(),
            substance_id: term.substance_id.to_string(),
        })?;
    Ok(IndexedStoichiometricTerm {
        substance,
        coefficient: term.coefficient,
        phases: vec![reaction
            .product_phases
            .get(&term.substance_id)
            .copied()
            .unwrap_or(MixturePhase::Aqueous)],
    })
}

fn build_reaction_index(
    indexed_reactions: &[IndexedReaction],
    substance_count: usize,
) -> (Vec<Vec<ReactionIndex>>, Vec<ReactionIndex>) {
    let mut by_substance = vec![Vec::new(); substance_count];
    let mut unindexed = Vec::new();
    for indexed_reaction in indexed_reactions {
        let mut indexed_substances = Vec::new();
        for reactant in &indexed_reaction.reactants {
            insert_sorted_unique(&mut indexed_substances, reactant.substance);
        }
        for (substance, _, _) in &indexed_reaction.orders {
            insert_sorted_unique(&mut indexed_substances, *substance);
        }

        if indexed_substances.is_empty() {
            unindexed.push(indexed_reaction.reaction);
            continue;
        }

        for substance in indexed_substances {
            by_substance[substance.0].push(indexed_reaction.reaction);
        }
    }
    (by_substance, unindexed)
}

fn insert_sorted_unique<T: Ord + Copy>(values: &mut Vec<T>, value: T) {
    match values.binary_search(&value) {
        Ok(_) => {}
        Err(index) => values.insert(index, value),
    }
}

fn ordered_pair(
    left: SubstanceIndex,
    right: SubstanceIndex,
) -> Option<(SubstanceIndex, SubstanceIndex)> {
    if left == right {
        None
    } else if left < right {
        Some((left, right))
    } else {
        Some((right, left))
    }
}

fn mark_reaction_candidate(scratch: &mut ReactionCandidateScratch, reaction_index: ReactionIndex) {
    if let Some(slot) = scratch.marks.get_mut(reaction_index.0) {
        if *slot != scratch.generation {
            *slot = scratch.generation;
            scratch.candidates.push(reaction_index);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::substance::{LiquidPhasePreference, SubstancePhaseProperties};

    fn aqueous_ion(id: &str, charge: i32, mass: f64) -> Substance {
        Substance::new(id, charge, mass, 1000.0, f64::MAX, 50.0, 0.0).with_phase_properties(
            SubstancePhaseProperties {
                preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                aqueous_solubility_mol_per_bucket: Some(10.0),
                organic_solubility_mol_per_bucket: Some(0.0),
                can_precipitate: false,
                can_form_liquid_phase: false,
                solvent_role: SolventRole::NotSolvent,
            },
        )
    }

    #[test]
    fn material_formula_mass_and_charge_are_checked_by_registry() {
        let bad_salt = Substance::new(
            "destroy:bad_sodium_chloride",
            0,
            100.0,
            2_160_000.0,
            f64::MAX,
            50.0,
            0.0,
        )
        .with_phase_properties(SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: Some(0.01),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: false,
            solvent_role: SolventRole::NotSolvent,
        })
        .with_representation(SubstanceRepresentation::IonicSolid {
            formula_units: vec![
                MaterialFormulaUnit::new("destroy:sodium_ion", 1),
                MaterialFormulaUnit::new("destroy:chloride", 1),
            ],
        });

        let error = ChemistryRegistryBuilder::new()
            .substance(aqueous_ion("destroy:sodium_ion", 1, 22.99))
            .substance(aqueous_ion("destroy:chloride", -1, 35.45))
            .substance(bad_salt)
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }
}
