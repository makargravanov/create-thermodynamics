mod constants;
mod logic;
mod types;
mod validation;

pub use logic::{apply_mechanical_working, metallurgical_state_from_alloy_phase};
pub use types::*;

#[cfg(test)]
mod tests;
