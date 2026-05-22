pub mod activity;
pub mod adiabatic;
pub mod analysis;
pub mod candidates;
pub mod chemistry;
pub mod equilibrium;
pub mod gas;
pub mod registry;
pub mod thermal;

pub use activity::{davies_log10_gamma, DAVIES_MAX_IONIC_STRENGTH_MOLAL};
pub use adiabatic::{
    solve_adiabatic_equilibrium, AdiabaticEquilibriumError, AdiabaticEquilibriumProblem,
    AdiabaticEquilibriumResult,
};
pub use analysis::{
    analyze_aqueous_equilibrium, analyze_phase_equilibrium, AqueousEquilibriumSummary,
    AqueousSpeciesSummary, EquilibriumAnalysisError, PhaseAmountSummary, PhaseEquilibriumSummary,
};
pub use candidates::{
    select_candidate_species, CandidateExclusion, CandidateExclusionReason, CandidatePhaseFilter,
    CandidateSelection, CandidateSelectionError, CandidateSelectionRequest,
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
pub use gas::{
    solve_closed_gas_equilibrium, ClosedGasEquilibriumError, ClosedGasEquilibriumProblem,
    ClosedGasEquilibriumResult,
};
pub use registry::{SpeciesRegistry, SpeciesRegistryError};
pub use thermal::{
    mixture_enthalpy_joule, mixture_heat_capacity_joule_per_kelvin, solve_temperature_for_enthalpy,
    thermal_state_for_composition, MixtureThermalState, ThermalError,
};
