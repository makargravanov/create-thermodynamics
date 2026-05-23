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

impl From<String> for ReactionId {
    fn from(value: String) -> Self {
        Self(value)
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
    pub external_reactants: Vec<ExternalRequirement>,
    pub external_catalysts: Vec<ExternalRequirement>,
    pub reaction_results: Vec<ReactionResult>,
    pub pre_exponential_factor: f64,
    pub activation_energy_kj_per_mol: f64,
    pub enthalpy_change_kj_per_mol: f64,
    pub reverse_reaction_id: Option<ReactionId>,
    pub requires_uv: bool,
    pub display_as_reversible: bool,
    pub show_in_jei: bool,
    pub show_in_jei_condition: Option<String>,
    pub allow_mass_imbalance: bool,
    pub allow_charge_imbalance: bool,
}

#[derive(Debug, Clone)]
pub struct ExternalRequirement {
    pub description: String,
    pub moles_per_reaction: f64,
    pub molar_mass_grams: Option<f64>,
    pub charge: Option<i32>,
    pub unchecked_mass_reason: Option<String>,
}

#[derive(Debug, Clone)]
pub struct ReactionResult {
    pub description: String,
    pub moles_per_reaction: f64,
}

impl Reaction {
    pub fn builder(id: impl Into<ReactionId>) -> ReactionBuilder {
        ReactionBuilder {
            reaction: Reaction {
                id: id.into(),
                reactants: Vec::new(),
                products: Vec::new(),
                orders: BTreeMap::new(),
                external_reactants: Vec::new(),
                external_catalysts: Vec::new(),
                reaction_results: Vec::new(),
                pre_exponential_factor: 10_000.0,
                activation_energy_kj_per_mol: 25.0,
                enthalpy_change_kj_per_mol: 0.0,
                reverse_reaction_id: None,
                requires_uv: false,
                display_as_reversible: false,
                show_in_jei: true,
                show_in_jei_condition: None,
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

    pub fn has_external_context(&self) -> bool {
        self.requires_uv
            || !self.external_reactants.is_empty()
            || !self.external_catalysts.is_empty()
            || !self.reaction_results.is_empty()
    }

    pub fn requires_context_to_proceed(&self) -> bool {
        self.requires_uv
            || !self.external_reactants.is_empty()
            || !self.external_catalysts.is_empty()
    }

    pub fn validate_shape(&self) -> ChemistryResult<()> {
        let reaction_id = self.id.to_string();
        if !self.has_external_context() {
            if self.reactants.is_empty() {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.clone(),
                    reason: "reaction must have at least one reactant".to_string(),
                });
            }
            if self.products.is_empty() {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.clone(),
                    reason: "reaction must have at least one product".to_string(),
                });
            }
        }
        for term in self.reactants.iter().chain(self.products.iter()) {
            if term.coefficient == 0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.clone(),
                    reason: "stoichiometric coefficients must be greater than zero".to_string(),
                });
            }
        }
        for requirement in self
            .external_reactants
            .iter()
            .chain(self.external_catalysts.iter())
        {
            if !requirement.moles_per_reaction.is_finite() || requirement.moles_per_reaction <= 0.0
            {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.clone(),
                    reason: "external requirements must be positive and finite".to_string(),
                });
            }
            if requirement.description.trim().is_empty() {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.clone(),
                    reason: "external requirements must have a description".to_string(),
                });
            }
            match (requirement.molar_mass_grams, requirement.charge) {
                (Some(mass), Some(_)) if mass.is_finite() && mass >= 0.0 => {}
                (None, None) if requirement.unchecked_mass_reason.is_some() => {}
                (Some(_), Some(_)) => {
                    return Err(ChemistryError::InvalidReaction {
                        reaction_id: reaction_id.clone(),
                        reason: "external requirement mass must be non-negative and finite"
                            .to_string(),
                    });
                }
                _ => {
                    return Err(ChemistryError::InvalidReaction {
                        reaction_id: reaction_id.clone(),
                        reason: "external requirement must either provide both mass and charge or an unchecked mass reason".to_string(),
                    });
                }
            }
        }
        for result in &self.reaction_results {
            if !result.moles_per_reaction.is_finite() || result.moles_per_reaction < 0.0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.clone(),
                    reason: "reaction results must be non-negative and finite".to_string(),
                });
            }
            if result.description.trim().is_empty() {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: reaction_id.clone(),
                    reason: "reaction results must have a description".to_string(),
                });
            }
        }
        if !self.pre_exponential_factor.is_finite() || self.pre_exponential_factor <= 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.clone(),
                reason: "pre-exponential factor must be positive and finite".to_string(),
            });
        }
        if !self.activation_energy_kj_per_mol.is_finite() || self.activation_energy_kj_per_mol < 0.0
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.clone(),
                reason: "activation energy must be non-negative and finite".to_string(),
            });
        }
        if !self.enthalpy_change_kj_per_mol.is_finite() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.clone(),
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

    pub fn external_reactant(mut self, description: impl Into<String>, moles: f64) -> Self {
        let description = description.into();
        self.reaction.external_reactants.push(ExternalRequirement {
            unchecked_mass_reason: Some(format!(
                "legacy external reactant '{description}' has no chemical formula in the model"
            )),
            description,
            moles_per_reaction: moles,
            molar_mass_grams: None,
            charge: None,
        });
        self
    }

    pub fn external_catalyst(mut self, description: impl Into<String>, moles: f64) -> Self {
        let description = description.into();
        self.reaction.external_catalysts.push(ExternalRequirement {
            unchecked_mass_reason: Some(format!(
                "legacy external catalyst '{description}' has no chemical formula in the model"
            )),
            description,
            moles_per_reaction: moles,
            molar_mass_grams: None,
            charge: None,
        });
        self
    }

    pub fn chemical_external_reactant(
        mut self,
        description: impl Into<String>,
        moles: f64,
        molar_mass_grams: f64,
        charge: i32,
    ) -> Self {
        self.reaction.external_reactants.push(ExternalRequirement {
            description: description.into(),
            moles_per_reaction: moles,
            molar_mass_grams: Some(molar_mass_grams),
            charge: Some(charge),
            unchecked_mass_reason: None,
        });
        self
    }

    pub fn chemical_external_catalyst(
        mut self,
        description: impl Into<String>,
        moles: f64,
        molar_mass_grams: f64,
        charge: i32,
    ) -> Self {
        self.reaction.external_catalysts.push(ExternalRequirement {
            description: description.into(),
            moles_per_reaction: moles,
            molar_mass_grams: Some(molar_mass_grams),
            charge: Some(charge),
            unchecked_mass_reason: None,
        });
        self
    }

    pub fn reaction_result(mut self, description: impl Into<String>, moles: f64) -> Self {
        self.reaction.reaction_results.push(ReactionResult {
            description: description.into(),
            moles_per_reaction: moles,
        });
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

    pub fn requires_uv(mut self) -> Self {
        self.reaction.requires_uv = true;
        self
    }

    pub fn display_as_reversible(mut self) -> Self {
        self.reaction.display_as_reversible = true;
        self
    }

    pub fn show_in_jei(mut self, value: bool) -> Self {
        self.reaction.show_in_jei = value;
        self
    }

    pub fn show_in_jei_condition(mut self, value: impl Into<String>) -> Self {
        self.reaction.show_in_jei_condition = Some(value.into());
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
