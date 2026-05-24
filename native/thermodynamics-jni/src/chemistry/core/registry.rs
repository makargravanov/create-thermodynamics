use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::mixture::MixturePhase;
use super::reaction::{Reaction, ReactionId, StoichiometricTerm};
use super::substance::{Substance, SubstanceAggregateState, SubstanceId, SubstanceTagId};

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
pub(crate) struct IndexedReaction {
    pub reaction: ReactionIndex,
    pub reactants: Vec<IndexedStoichiometricTerm>,
    pub products: Vec<IndexedStoichiometricTerm>,
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
        ordered_pair(left, right)
            .and_then(|key| self.solvent_miscibility.get(&key).copied())
            .unwrap_or(SolventMiscibility::Immiscible)
    }

    pub(crate) fn gas_solubility(&self, substance: SubstanceIndex) -> Option<&GasSolubilityModel> {
        self.gas_solubility.get(&substance)
    }
}

#[derive(Default)]
pub struct ChemistryRegistryBuilder {
    substances: Vec<Substance>,
    reactions: Vec<Reaction>,
    substance_tags: BTreeSet<SubstanceTagId>,
    solvent_miscibility: Vec<(SubstanceId, SubstanceId, SolventMiscibility)>,
    gas_solubility: Vec<(SubstanceId, GasSolubilityModel)>,
}

impl ChemistryRegistryBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn from_registry(registry: &ChemistryRegistry) -> Self {
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

    pub fn build(self) -> ChemistryResult<ChemistryRegistry> {
        let mut substance_map = BTreeMap::new();
        for substance in self.substances {
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
        for reaction in self.reactions {
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
        };
        registry.validate_substance_tags()?;
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

    fn validate_reactions(&self) -> ChemistryResult<()> {
        for reaction in &self.reactions_by_index {
            for term in reaction.reactants.iter().chain(reaction.products.iter()) {
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

            for requirement in reaction
                .external_reactants
                .iter()
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
                        .map(|charge| charge * requirement.moles_per_reaction.round() as i32)
                })
                .sum::<i32>();

            let reactant_charge = reaction
                .reactants
                .iter()
                .map(|term| {
                    self.substance_index(&term.substance_id)
                        .map(|index| {
                            self.substance_properties.charge[index.0] * term.coefficient as i32
                        })
                        .ok_or_else(|| ChemistryError::UnknownSubstance {
                            reaction_id: reaction.id.to_string(),
                            substance_id: term.substance_id.to_string(),
                        })
                })
                .sum::<ChemistryResult<i32>>()?
                + external_reactant_charge;
            let product_charge = reaction
                .products
                .iter()
                .map(|term| {
                    self.substance_index(&term.substance_id)
                        .map(|index| {
                            self.substance_properties.charge[index.0] * term.coefficient as i32
                        })
                        .ok_or_else(|| ChemistryError::UnknownSubstance {
                            reaction_id: reaction.id.to_string(),
                            substance_id: term.substance_id.to_string(),
                        })
                })
                .sum::<ChemistryResult<i32>>()?;
            if reactant_charge != product_charge && !reaction.allow_charge_imbalance {
                return Err(ChemistryError::ChargeNotConserved {
                    reaction_id: reaction.id.to_string(),
                    reactants: reactant_charge,
                    products: product_charge,
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
            let product_mass = reaction
                .products
                .iter()
                .map(|term| {
                    self.substance(&term.substance_id)
                        .map(|substance| substance.molar_mass_grams * term.coefficient as f64)
                })
                .sum::<ChemistryResult<f64>>()?;
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
}

fn stoichiometric_map(terms: &[StoichiometricTerm]) -> BTreeMap<SubstanceId, u32> {
    let mut result = BTreeMap::new();
    for term in terms {
        *result.entry(term.substance_id.clone()).or_insert(0) += term.coefficient;
    }
    result
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
    if !substance.phase_properties.can_form_liquid_phase {
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
                orders,
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
