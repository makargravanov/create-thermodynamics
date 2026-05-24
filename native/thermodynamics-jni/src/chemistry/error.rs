use std::error::Error;
use std::fmt::{Display, Formatter};

pub type ChemistryResult<T> = Result<T, ChemistryError>;

#[derive(Debug, Clone, PartialEq)]
pub enum ChemistryError {
    DuplicateSubstance(String),
    DuplicateReaction(String),
    UnknownSubstance {
        reaction_id: String,
        substance_id: String,
    },
    UnknownReaction(String),
    InvalidSubstance {
        substance_id: String,
        reason: String,
    },
    InvalidReaction {
        reaction_id: String,
        reason: String,
    },
    ChargeNotConserved {
        reaction_id: String,
        reactants: i32,
        products: i32,
    },
    MassNotConserved {
        reaction_id: String,
        reactants: f64,
        products: f64,
    },
    ReversibleThermodynamicsMismatch {
        reaction_id: String,
        reverse_id: String,
        reason: String,
    },
    GenerationInvariantViolation {
        generator: String,
        substance_id: String,
        reason: String,
    },
    EquilibriumInvariantViolation {
        equilibrium_id: String,
        reason: String,
    },
    InvalidMixtureState(String),
}

impl Display for ChemistryError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            ChemistryError::DuplicateSubstance(id) => {
                write!(f, "substance '{id}' is registered more than once")
            }
            ChemistryError::DuplicateReaction(id) => {
                write!(f, "reaction '{id}' is registered more than once")
            }
            ChemistryError::UnknownSubstance {
                reaction_id,
                substance_id,
            } => {
                write!(
                    f,
                    "reaction '{reaction_id}' refers to unknown substance '{substance_id}'"
                )
            }
            ChemistryError::UnknownReaction(id) => write!(f, "unknown reaction '{id}'"),
            ChemistryError::InvalidSubstance {
                substance_id,
                reason,
            } => {
                write!(f, "substance '{substance_id}' is invalid: {reason}")
            }
            ChemistryError::InvalidReaction {
                reaction_id,
                reason,
            } => {
                write!(f, "reaction '{reaction_id}' is invalid: {reason}")
            }
            ChemistryError::ChargeNotConserved {
                reaction_id,
                reactants,
                products,
            } => {
                write!(f, "reaction '{reaction_id}' does not conserve charge: reactants={reactants}, products={products}")
            }
            ChemistryError::MassNotConserved {
                reaction_id,
                reactants,
                products,
            } => {
                write!(f, "reaction '{reaction_id}' does not conserve mass: reactants={reactants}, products={products}")
            }
            ChemistryError::ReversibleThermodynamicsMismatch {
                reaction_id,
                reverse_id,
                reason,
            } => {
                write!(f, "reversible reactions '{reaction_id}' and '{reverse_id}' are inconsistent: {reason}")
            }
            ChemistryError::GenerationInvariantViolation {
                generator,
                substance_id,
                reason,
            } => {
                write!(
                    f,
                    "generation invariant failed in '{generator}' for '{substance_id}': {reason}"
                )
            }
            ChemistryError::EquilibriumInvariantViolation {
                equilibrium_id,
                reason,
            } => {
                write!(
                    f,
                    "equilibrium invariant failed in '{equilibrium_id}': {reason}"
                )
            }
            ChemistryError::InvalidMixtureState(reason) => {
                write!(f, "invalid mixture state: {reason}")
            }
        }
    }
}

impl Error for ChemistryError {}
