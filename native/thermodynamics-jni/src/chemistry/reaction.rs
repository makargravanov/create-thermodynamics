use std::collections::BTreeMap;
use std::fmt::{Display, Formatter};

use super::error::{ChemistryError, ChemistryResult};
use super::substance::SubstanceId;

pub const GAS_CONSTANT_J_PER_MOL_KELVIN: f64 = 8.314_462_618_153_24;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReactionId(String);

impl ReactionId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: value,
                reason: "id must not be empty".to_string(),
            });
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for ReactionId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for ReactionId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

#[derive(Debug, Clone)]
pub struct StoichiometricTerm {
    pub substance_id: SubstanceId,
    pub coefficient: u32,
}

impl StoichiometricTerm {
    pub fn new(substance_id: impl Into<SubstanceId>, coefficient: u32) -> Self {
        Self {
            substance_id: substance_id.into(),
            coefficient,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Reaction {
    pub id: ReactionId,
    pub reactants: Vec<StoichiometricTerm>,
    pub products: Vec<StoichiometricTerm>,
    pub orders: BTreeMap<SubstanceId, u32>,
    pub pre_exponential_factor: f64,
    pub activation_energy_kj_per_mol: f64,
    pub enthalpy_change_kj_per_mol: f64,
    pub reverse_reaction_id: Option<ReactionId>,
    pub allow_mass_imbalance: bool,
    pub allow_charge_imbalance: bool,
}

impl Reaction {
    pub fn builder(id: impl Into<ReactionId>) -> ReactionBuilder {
        ReactionBuilder {
            reaction: Reaction {
                id: id.into(),
                reactants: Vec::new(),
                products: Vec::new(),
                orders: BTreeMap::new(),
                pre_exponential_factor: 10_000.0,
                activation_energy_kj_per_mol: 25.0,
                enthalpy_change_kj_per_mol: 0.0,
                reverse_reaction_id: None,
                allow_mass_imbalance: false,
                allow_charge_imbalance: false,
            },
        }
    }

    pub fn rate_constant_per_second(&self, temperature_kelvin: f64) -> ChemistryResult<f64> {
        if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: self.id.to_string(),
                reason: "temperature must be positive and finite".to_string(),
            });
        }
        Ok(self.pre_exponential_factor
            * (-(self.activation_energy_kj_per_mol * 1000.0)
                / (GAS_CONSTANT_J_PER_MOL_KELVIN * temperature_kelvin))
                .exp())
    }

    pub fn validate_shape(&self) -> ChemistryResult<()> {
        let id = self.id.to_string();
        if self.reactants.is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: id,
                reason: "reaction must have at least one reactant".to_string(),
            });
        }
        if self.products.is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: id,
                reason: "reaction must have at least one product".to_string(),
            });
        }
        for term in self.reactants.iter().chain(self.products.iter()) {
            if term.coefficient == 0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: id,
                    reason: "stoichiometric coefficients must be greater than zero".to_string(),
                });
            }
        }
        if !self.pre_exponential_factor.is_finite() || self.pre_exponential_factor <= 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: id,
                reason: "pre-exponential factor must be positive and finite".to_string(),
            });
        }
        if !self.activation_energy_kj_per_mol.is_finite() || self.activation_energy_kj_per_mol < 0.0
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: id,
                reason: "activation energy must be non-negative and finite".to_string(),
            });
        }
        if !self.enthalpy_change_kj_per_mol.is_finite() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: id,
                reason: "enthalpy change must be finite".to_string(),
            });
        }
        Ok(())
    }
}

pub struct ReactionBuilder {
    reaction: Reaction,
}

impl ReactionBuilder {
    pub fn reactant(
        mut self,
        substance_id: impl Into<SubstanceId>,
        coefficient: u32,
        order: u32,
    ) -> Self {
        let substance_id = substance_id.into();
        self.reaction
            .reactants
            .push(StoichiometricTerm::new(substance_id.clone(), coefficient));
        self.reaction.orders.insert(substance_id, order);
        self
    }

    pub fn product(mut self, substance_id: impl Into<SubstanceId>, coefficient: u32) -> Self {
        self.reaction
            .products
            .push(StoichiometricTerm::new(substance_id, coefficient));
        self
    }

    pub fn catalyst_order(mut self, substance_id: impl Into<SubstanceId>, order: u32) -> Self {
        self.reaction.orders.insert(substance_id.into(), order);
        self
    }

    pub fn pre_exponential_factor(mut self, value: f64) -> Self {
        self.reaction.pre_exponential_factor = value;
        self
    }

    pub fn activation_energy_kj_per_mol(mut self, value: f64) -> Self {
        self.reaction.activation_energy_kj_per_mol = value;
        self
    }

    pub fn enthalpy_change_kj_per_mol(mut self, value: f64) -> Self {
        self.reaction.enthalpy_change_kj_per_mol = value;
        self
    }

    pub fn reverse_reaction_id(mut self, id: impl Into<ReactionId>) -> Self {
        self.reaction.reverse_reaction_id = Some(id.into());
        self
    }

    pub fn allow_mass_imbalance(mut self) -> Self {
        self.reaction.allow_mass_imbalance = true;
        self
    }

    pub fn allow_charge_imbalance(mut self) -> Self {
        self.reaction.allow_charge_imbalance = true;
        self
    }

    pub fn build(self) -> Reaction {
        self.reaction
    }
}
