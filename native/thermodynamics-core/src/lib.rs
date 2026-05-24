pub mod activity;
pub mod adiabatic;
pub mod analysis;
pub mod chemistry;
pub mod equilibrium;
pub mod registry;
pub mod thermal;

pub use activity::{davies_log10_gamma, DAVIES_MAX_IONIC_STRENGTH_MOLAL};
pub use adiabatic::{
    solve_adiabatic_equilibrium, AdiabaticEquilibriumError, AdiabaticEquilibriumProblem,
    AdiabaticEquilibriumResult,
};
pub use analysis::{
    analyze_aqueous_equilibrium, AqueousEquilibriumSummary, AqueousSpeciesSummary,
    EquilibriumAnalysisError,
};
pub use chemistry::{
    ActivityModel, ConstantPressureHeatCapacity, DataSource, Element, ElementId, PhaseId,
    PhaseKind, Species, SpeciesAmount, SpeciesId, StandardEnthalpyOfFormation, StandardGibbsEnergy,
    StandardThermo, TemperatureRange,
};
pub use equilibrium::{
    solve_equilibrium, EquilibriumDiagnostic, EquilibriumError, EquilibriumProblem,
    EquilibriumResiduals, EquilibriumResult,
};
pub use registry::{SpeciesRegistry, SpeciesRegistryError};
pub use thermal::{
    mixture_enthalpy_joule, mixture_heat_capacity_joule_per_kelvin, solve_temperature_for_enthalpy,
    thermal_state_for_composition, MixtureThermalState, ThermalError,
};
