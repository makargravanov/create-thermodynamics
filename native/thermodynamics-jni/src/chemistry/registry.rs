use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::reaction::{Reaction, ReactionId, StoichiometricTerm};
use super::substance::{Substance, SubstanceId, SubstanceTagId};

const MASS_TOLERANCE_GRAMS_PER_MOL: f64 = 1.0e-6;
const THERMO_TOLERANCE: f64 = 1.0e-6;

#[derive(Debug, Clone)]
pub struct ChemistryRegistry {
    substances: BTreeMap<SubstanceId, Substance>,
    reactions: BTreeMap<ReactionId, Reaction>,
    substance_tags: BTreeSet<SubstanceTagId>,
}

impl ChemistryRegistry {
    pub fn substance(&self, id: &SubstanceId) -> ChemistryResult<&Substance> {
        self.substances
            .get(id)
            .ok_or_else(|| ChemistryError::InvalidMixtureState(format!("unknown substance '{id}'")))
    }

    pub fn reaction(&self, id: &ReactionId) -> ChemistryResult<&Reaction> {
        self.reactions
            .get(id)
            .ok_or_else(|| ChemistryError::UnknownReaction(id.to_string()))
    }

    pub fn reactions(&self) -> impl Iterator<Item = &Reaction> {
        self.reactions.values()
    }

    pub fn substances(&self) -> impl Iterator<Item = &Substance> {
        self.substances.values()
    }

    pub fn substance_tags(&self) -> impl Iterator<Item = &SubstanceTagId> {
        self.substance_tags.iter()
    }

    pub fn substance_count(&self) -> usize {
        self.substances.len()
    }

    pub fn has_substance_tag(&self, id: &SubstanceTagId) -> bool {
        self.substance_tags.contains(id)
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
            substances: registry.substances.values().cloned().collect(),
            reactions: registry.reactions.values().cloned().collect(),
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
        let mut substances = BTreeMap::new();
        for substance in self.substances {
            substance.validate()?;
            let id = substance.id.clone();
            if substances.insert(id.clone(), substance).is_some() {
                return Err(ChemistryError::DuplicateSubstance(id.to_string()));
            }
        }

        let mut reactions = BTreeMap::new();
        for reaction in self.reactions {
            reaction.validate_shape()?;
            let id = reaction.id.clone();
            if reactions.insert(id.clone(), reaction).is_some() {
                return Err(ChemistryError::DuplicateReaction(id.to_string()));
            }
        }

        let registry = ChemistryRegistry {
            substances,
            reactions,
            substance_tags: self.substance_tags,
        };
        registry.validate_substance_tags()?;
        registry.validate_reactions()?;
        Ok(registry)
    }
}

impl ChemistryRegistry {
    fn validate_substance_tags(&self) -> ChemistryResult<()> {
        for substance in self.substances.values() {
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
        for reaction in self.reactions.values() {
            for term in reaction.reactants.iter().chain(reaction.products.iter()) {
                if !self.substances.contains_key(&term.substance_id) {
                    return Err(ChemistryError::UnknownSubstance {
                        reaction_id: reaction.id.to_string(),
                        substance_id: term.substance_id.to_string(),
                    });
                }
            }
            for ordered_substance in reaction.orders.keys() {
                if !self.substances.contains_key(ordered_substance) {
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
                .map(|term| self.substances[&term.substance_id].charge * term.coefficient as i32)
                .sum::<i32>()
                + external_reactant_charge;
            let product_charge = reaction
                .products
                .iter()
                .map(|term| self.substances[&term.substance_id].charge * term.coefficient as i32)
                .sum::<i32>();
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
                    self.substances[&term.substance_id].molar_mass_grams * term.coefficient as f64
                })
                .sum::<f64>()
                + external_reactant_mass;
            let product_mass = reaction
                .products
                .iter()
                .map(|term| {
                    self.substances[&term.substance_id].molar_mass_grams * term.coefficient as f64
                })
                .sum::<f64>();
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
                let reverse = self
                    .reactions
                    .get(reverse_id)
                    .ok_or_else(|| ChemistryError::UnknownReaction(reverse_id.to_string()))?;
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
