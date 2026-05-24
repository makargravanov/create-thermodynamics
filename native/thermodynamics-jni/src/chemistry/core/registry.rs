use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::reaction::{Reaction, ReactionId, StoichiometricTerm};
use super::substance::{Substance, SubstanceId, SubstanceTagId};

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

#[derive(Debug, Clone)]
pub(crate) struct IndexedStoichiometricTerm {
    pub substance: SubstanceIndex,
    pub coefficient: u32,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedReaction {
    pub reaction: ReactionIndex,
    pub reactants: Vec<IndexedStoichiometricTerm>,
    pub products: Vec<IndexedStoichiometricTerm>,
    pub orders: Vec<(SubstanceIndex, u32)>,
}

#[derive(Debug, Clone)]
pub(crate) struct SubstancePropertiesTable {
    pub charge: Vec<i32>,
    pub molar_mass_grams: Vec<f64>,
    pub liquid_density_grams_per_bucket: Vec<f64>,
    pub boiling_point_kelvin: Vec<f64>,
    pub molar_heat_capacity_j_per_mol_kelvin: Vec<f64>,
    pub latent_heat_j_per_mol: Vec<f64>,
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
        let mut seen = vec![false; self.reactions_by_index.len()];
        let mut result = Vec::new();
        for reaction_index in &self.unindexed_reaction_indices {
            mark_reaction_candidate(&mut seen, &mut result, *reaction_index);
        }
        for substance_index in substances {
            if let Some(indexed_reactions) = self.reaction_index_by_substance.get(substance_index.0)
            {
                for reaction_index in indexed_reactions {
                    mark_reaction_candidate(&mut seen, &mut result, *reaction_index);
                }
            }
        }
        result
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
}

#[derive(Default)]
pub struct ChemistryRegistryBuilder {
    substances: Vec<Substance>,
    reactions: Vec<Reaction>,
    substance_tags: BTreeSet<SubstanceTagId>,
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
        boiling_point_kelvin: substances
            .iter()
            .map(|substance| substance.boiling_point_kelvin)
            .collect(),
        molar_heat_capacity_j_per_mol_kelvin: substances
            .iter()
            .map(|substance| substance.molar_heat_capacity_j_per_mol_kelvin)
            .collect(),
        latent_heat_j_per_mol: substances
            .iter()
            .map(|substance| substance.latent_heat_j_per_mol)
            .collect(),
    }
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
                .map(|term| indexed_stoichiometric_term(reaction, term, substance_id_to_index))
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
                    Ok((substance, *order))
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
        for (substance, _) in &indexed_reaction.orders {
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

fn mark_reaction_candidate(
    seen: &mut [bool],
    result: &mut Vec<ReactionIndex>,
    reaction_index: ReactionIndex,
) {
    if let Some(slot) = seen.get_mut(reaction_index.0) {
        if !*slot {
            *slot = true;
            result.push(reaction_index);
        }
    }
}
