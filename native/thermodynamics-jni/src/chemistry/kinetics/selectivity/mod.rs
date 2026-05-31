//! Kinetic selectivity module for organic reactions
//!
//! This module provides centralized evaluation of reaction selectivity based on:
//! - Substitution degree (1°, 2°, 3°, benzylic, allylic)
//! - Electronic environment (electron-donating/withdrawing groups)
//! - Steric hindrance
//! - Reaction conditions (temperature, pH, solvent)

pub mod carbonyl_addition;
pub mod elimination;
pub mod engine;
pub mod esterification;
pub mod nucleophilic_substitution;
pub mod types;

pub use engine::SelectivityEngine;
pub use engine::SiteDescriptorBuilder;
pub use types::*;
