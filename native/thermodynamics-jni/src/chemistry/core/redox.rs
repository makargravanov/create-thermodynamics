use std::collections::BTreeMap;

use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, MixturePhase, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::molecule::MolecularStructure;
use super::reaction::{
    Reaction, ReactionBuilder, ReactionId, StoichiometricTerm, GAS_CONSTANT_J_PER_MOL_KELVIN,
};
use super::registry::ChemistryRegistry;
use super::substance::{LiquidPhasePreference, SubstanceAggregateState, SubstanceId};

pub const ELECTRON_EXTERNAL_ID: &str = "redox:electron";
pub const FARADAY_CONSTANT_COULOMBS_PER_MOL: f64 = 96_485.332_123_310_02;
const DEFAULT_ELECTRODE_EFFICIENCY: f64 = 1.0;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum OxidationStateRule {
    Electronegativity,
    Explicit,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AtomOxidationState {
    pub atom_index: usize,
    pub element: String,
    pub state: f64,
    pub rule: OxidationStateRule,
}

#[derive(Debug, Clone, PartialEq)]
pub struct OxidationStateAssignment {
    pub atom_states: Vec<AtomOxidationState>,
    pub total_charge: i32,
}

impl OxidationStateAssignment {
    pub fn state_sum_for_element(&self, element: &str) -> f64 {
        self.atom_states
            .iter()
            .filter(|state| state.element == element)
            .map(|state| state.state)
            .sum()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RedoxRole {
    Oxidant,
    Reductant,
    OxidantAndReductant,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ExplicitOxidationState {
    pub element: String,
    pub state: f64,
    pub atom_count: u32,
}

impl ExplicitOxidationState {
    pub fn new(element: impl Into<String>, state: f64, atom_count: u32) -> Self {
        Self {
            element: element.into(),
            state,
            atom_count,
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum RedoxEnvironment {
    Acidic,
    Basic,
    Neutral,
    Any,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ElectronSide {
    Reactant,
    Product,
}

#[derive(Debug, Clone)]
pub struct RedoxHalfReaction {
    pub id: String,
    pub reactants: Vec<StoichiometricTerm>,
    pub products: Vec<StoichiometricTerm>,
    pub electron_count: u32,
    pub electron_side: ElectronSide,
    pub environment: RedoxEnvironment,
    pub standard_potential_volts: Option<f64>,
}

impl RedoxHalfReaction {
    pub fn oxidation(
        id: impl Into<String>,
        reactants: impl IntoIterator<Item = (SubstanceId, u32)>,
        products: impl IntoIterator<Item = (SubstanceId, u32)>,
        electron_count: u32,
        environment: RedoxEnvironment,
    ) -> Self {
        Self::new(
            id,
            reactants,
            products,
            electron_count,
            ElectronSide::Product,
            environment,
        )
    }

    pub fn reduction(
        id: impl Into<String>,
        reactants: impl IntoIterator<Item = (SubstanceId, u32)>,
        products: impl IntoIterator<Item = (SubstanceId, u32)>,
        electron_count: u32,
        environment: RedoxEnvironment,
    ) -> Self {
        Self::new(
            id,
            reactants,
            products,
            electron_count,
            ElectronSide::Reactant,
            environment,
        )
    }

    pub fn new(
        id: impl Into<String>,
        reactants: impl IntoIterator<Item = (SubstanceId, u32)>,
        products: impl IntoIterator<Item = (SubstanceId, u32)>,
        electron_count: u32,
        electron_side: ElectronSide,
        environment: RedoxEnvironment,
    ) -> Self {
        Self {
            id: id.into(),
            reactants: reactants
                .into_iter()
                .map(|(substance_id, coefficient)| {
                    StoichiometricTerm::new(substance_id, coefficient)
                })
                .collect(),
            products: products
                .into_iter()
                .map(|(substance_id, coefficient)| {
                    StoichiometricTerm::new(substance_id, coefficient)
                })
                .collect(),
            electron_count,
            electron_side,
            environment,
            standard_potential_volts: None,
        }
    }

    pub fn with_standard_potential_volts(mut self, value: f64) -> Self {
        self.standard_potential_volts = Some(value);
        self
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RedoxPotentialEvaluation {
    pub oxidation_half_id: String,
    pub reduction_half_id: String,
    pub oxidation_potential_volts: f64,
    pub reduction_potential_volts: f64,
    pub cell_potential_volts: f64,
    pub equilibrium_constant: f64,
    pub thermodynamic_rate_factor: f64,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum ElectrodePolarity {
    Anode,
    Cathode,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ElectrodeProcess {
    pub half_reaction_id: String,
    pub polarity: ElectrodePolarity,
    pub phase: MixturePhase,
    pub current_efficiency: f64,
    pub overpotential_volts: f64,
}

impl ElectrodeProcess {
    pub fn anode(half_reaction_id: impl Into<String>) -> Self {
        Self::new(half_reaction_id, ElectrodePolarity::Anode)
    }

    pub fn cathode(half_reaction_id: impl Into<String>) -> Self {
        Self::new(half_reaction_id, ElectrodePolarity::Cathode)
    }

    pub fn new(half_reaction_id: impl Into<String>, polarity: ElectrodePolarity) -> Self {
        Self {
            half_reaction_id: half_reaction_id.into(),
            polarity,
            phase: MixturePhase::Aqueous,
            current_efficiency: DEFAULT_ELECTRODE_EFFICIENCY,
            overpotential_volts: 0.0,
        }
    }

    pub fn in_phase(mut self, phase: MixturePhase) -> Self {
        self.phase = phase;
        self
    }

    pub fn with_current_efficiency(mut self, value: f64) -> Self {
        self.current_efficiency = value;
        self
    }

    pub fn with_overpotential_volts(mut self, value: f64) -> Self {
        self.overpotential_volts = value;
        self
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ElectrolysisCell {
    pub anode: ElectrodeProcess,
    pub cathode: ElectrodeProcess,
    pub applied_voltage_volts: f64,
}

impl ElectrolysisCell {
    pub fn new(
        anode: ElectrodeProcess,
        cathode: ElectrodeProcess,
        applied_voltage_volts: f64,
    ) -> Self {
        Self {
            anode,
            cathode,
            applied_voltage_volts,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ElectrolysisReport {
    pub requested_charge_coulombs: f64,
    pub transferred_electrons_mol_per_bucket: f64,
    pub anode_extent_mol_per_bucket: f64,
    pub cathode_extent_mol_per_bucket: f64,
    pub anode_potential_volts: f64,
    pub cathode_potential_volts: f64,
    pub reversible_voltage_volts: f64,
    pub required_voltage_volts: f64,
}

#[derive(Debug, Clone)]
pub struct RedoxPair {
    pub reaction: Reaction,
    pub oxidation_half_id: String,
    pub reduction_half_id: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RedoxPairSpec {
    pub reaction_id: ReactionId,
    pub oxidation_half_id: String,
    pub reduction_half_id: String,
}

impl RedoxPairSpec {
    pub fn new(
        reaction_id: impl Into<ReactionId>,
        oxidation_half_id: impl Into<String>,
        reduction_half_id: impl Into<String>,
    ) -> Self {
        Self {
            reaction_id: reaction_id.into(),
            oxidation_half_id: oxidation_half_id.into(),
            reduction_half_id: reduction_half_id.into(),
        }
    }
}

impl RedoxPair {
    pub fn new(
        mut reaction: Reaction,
        oxidation_half_id: impl Into<String>,
        reduction_half_id: impl Into<String>,
        electron_count: u32,
        environment: RedoxEnvironment,
    ) -> Self {
        let oxidation_half_id = oxidation_half_id.into();
        let reduction_half_id = reduction_half_id.into();
        reaction.redox = Some(RedoxAnnotation {
            oxidant: None,
            reductant: None,
            transferred_electrons: electron_count,
            environment,
            oxidation_half_id: Some(oxidation_half_id.clone()),
            reduction_half_id: Some(reduction_half_id.clone()),
            electron_balance_checked: true,
        });
        Self {
            reaction,
            oxidation_half_id,
            reduction_half_id,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RedoxAnnotation {
    pub oxidant: Option<SubstanceId>,
    pub reductant: Option<SubstanceId>,
    pub transferred_electrons: u32,
    pub environment: RedoxEnvironment,
    pub oxidation_half_id: Option<String>,
    pub reduction_half_id: Option<String>,
    pub electron_balance_checked: bool,
}

impl RedoxAnnotation {
    pub fn checked(
        oxidant: impl Into<SubstanceId>,
        reductant: impl Into<SubstanceId>,
        transferred_electrons: u32,
        environment: RedoxEnvironment,
    ) -> Self {
        Self {
            oxidant: Some(oxidant.into()),
            reductant: Some(reductant.into()),
            transferred_electrons,
            environment,
            oxidation_half_id: None,
            reduction_half_id: None,
            electron_balance_checked: true,
        }
    }

    pub fn from_halves(
        transferred_electrons: u32,
        environment: RedoxEnvironment,
        oxidation_half_id: impl Into<String>,
        reduction_half_id: impl Into<String>,
    ) -> Self {
        Self {
            oxidant: None,
            reductant: None,
            transferred_electrons,
            environment,
            oxidation_half_id: Some(oxidation_half_id.into()),
            reduction_half_id: Some(reduction_half_id.into()),
            electron_balance_checked: true,
        }
    }
}

pub fn assign_oxidation_states(
    structure: &MolecularStructure,
) -> ChemistryResult<OxidationStateAssignment> {
    structure.validate()?;
    let mut states = structure
        .atoms
        .iter()
        .enumerate()
        .map(|(atom_index, atom)| AtomOxidationState {
            atom_index,
            element: atom.element.clone(),
            state: atom.charge,
            rule: OxidationStateRule::Electronegativity,
        })
        .collect::<Vec<_>>();

    for bond in &structure.bonds {
        let left = &structure.atoms[bond.from];
        let right = &structure.atoms[bond.to];
        let left_electronegativity = electronegativity(&left.element)?;
        let right_electronegativity = electronegativity(&right.element)?;
        if (left_electronegativity - right_electronegativity).abs() <= 1.0e-9 {
            continue;
        }
        if left_electronegativity > right_electronegativity {
            states[bond.from].state -= bond.order;
            states[bond.to].state += bond.order;
        } else {
            states[bond.from].state += bond.order;
            states[bond.to].state -= bond.order;
        }
    }

    let total = states.iter().map(|state| state.state).sum::<f64>();
    if !total.is_finite() {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: structure.source_code.clone(),
            reason: "oxidation states became non-finite".to_string(),
        });
    }
    Ok(OxidationStateAssignment {
        atom_states: states,
        total_charge: total.round() as i32,
    })
}

pub fn explicit_oxidation_assignment(
    substance_id: &SubstanceId,
    charge: i32,
    states: &[ExplicitOxidationState],
) -> ChemistryResult<OxidationStateAssignment> {
    if states.is_empty() {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: substance_id.to_string(),
            reason: "explicit oxidation states must not be empty".to_string(),
        });
    }
    let mut atom_states = Vec::new();
    let mut total = 0.0;
    for state in states {
        if state.element.trim().is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "oxidation state element must not be empty".to_string(),
            });
        }
        if !state.state.is_finite() || state.atom_count == 0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "oxidation state must be finite and atom count must be positive"
                    .to_string(),
            });
        }
        for _ in 0..state.atom_count {
            atom_states.push(AtomOxidationState {
                atom_index: atom_states.len(),
                element: state.element.clone(),
                state: state.state,
                rule: OxidationStateRule::Explicit,
            });
            total += state.state;
        }
    }
    if (total - charge as f64).abs() > 1.0e-6 {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: substance_id.to_string(),
            reason: format!("explicit oxidation states sum to {total}, expected charge {charge}"),
        });
    }
    Ok(OxidationStateAssignment {
        atom_states,
        total_charge: charge,
    })
}

pub(crate) fn validate_half_reaction_shape(half: &RedoxHalfReaction) -> ChemistryResult<()> {
    if half.id.trim().is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<redox-half>".to_string(),
            reason: "redox half reaction id must not be empty".to_string(),
        });
    }
    if half.reactants.is_empty() || half.products.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: half.id.clone(),
            reason: "redox half reaction must have reactants and products".to_string(),
        });
    }
    if half.electron_count == 0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: half.id.clone(),
            reason: "redox half reaction must transfer at least one electron".to_string(),
        });
    }
    if half
        .standard_potential_volts
        .is_some_and(|potential| !potential.is_finite())
    {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: half.id.clone(),
            reason: "standard redox potential must be finite when present".to_string(),
        });
    }
    for term in half.reactants.iter().chain(half.products.iter()) {
        if term.coefficient == 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: half.id.clone(),
                reason: "redox half reaction coefficients must be positive".to_string(),
            });
        }
    }
    Ok(())
}

pub(crate) fn validate_redox_annotation(
    reaction: &Reaction,
    known_reaction_ids: &BTreeMap<ReactionId, usize>,
    half_reactions: &BTreeMap<String, RedoxHalfReaction>,
) -> ChemistryResult<()> {
    let Some(redox) = &reaction.redox else {
        return Ok(());
    };
    if reaction.allow_charge_imbalance {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "redox reactions may not allow charge imbalance".to_string(),
        });
    }
    if redox.transferred_electrons == 0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "redox annotation must transfer at least one electron".to_string(),
        });
    }
    if !redox.electron_balance_checked {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "redox annotation must be explicitly electron-balanced".to_string(),
        });
    }
    if reaction
        .external_reactants
        .iter()
        .chain(reaction.external_products.iter())
        .any(|external| external.description == ELECTRON_EXTERNAL_ID)
    {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "closed redox reaction must not contain free electrons".to_string(),
        });
    }
    if let Some(oxidation_half_id) = &redox.oxidation_half_id {
        let oxidation = half_reactions
            .get(oxidation_half_id)
            .ok_or_else(|| ChemistryError::UnknownReaction(oxidation_half_id.clone()))?;
        if oxidation.electron_side != ElectronSide::Product {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "oxidation half reaction must produce electrons".to_string(),
            });
        }
        if redox.transferred_electrons % oxidation.electron_count != 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "oxidation half electron count must divide redox annotation electron count"
                    .to_string(),
            });
        }
        if oxidation.standard_potential_volts.is_some() && redox.reduction_half_id.is_none() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "redox potentials require both oxidation and reduction halves".to_string(),
            });
        }
    }
    if let Some(reduction_half_id) = &redox.reduction_half_id {
        let reduction = half_reactions
            .get(reduction_half_id)
            .ok_or_else(|| ChemistryError::UnknownReaction(reduction_half_id.clone()))?;
        if reduction.electron_side != ElectronSide::Reactant {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "reduction half reaction must consume electrons".to_string(),
            });
        }
        if redox.transferred_electrons % reduction.electron_count != 0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "reduction half electron count must divide redox annotation electron count"
                    .to_string(),
            });
        }
        if reduction.standard_potential_volts.is_some() && redox.oxidation_half_id.is_none() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "redox potentials require both oxidation and reduction halves".to_string(),
            });
        }
    }
    if let (Some(oxidation_half_id), Some(reduction_half_id)) =
        (&redox.oxidation_half_id, &redox.reduction_half_id)
    {
        let oxidation = half_reactions
            .get(oxidation_half_id)
            .ok_or_else(|| ChemistryError::UnknownReaction(oxidation_half_id.clone()))?;
        let reduction = half_reactions
            .get(reduction_half_id)
            .ok_or_else(|| ChemistryError::UnknownReaction(reduction_half_id.clone()))?;
        if oxidation.standard_potential_volts.is_some()
            != reduction.standard_potential_volts.is_some()
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "paired redox halves must either both define potentials or both omit them"
                    .to_string(),
            });
        }
    }
    if let Some(reverse_id) = &reaction.reverse_reaction_id {
        if !known_reaction_ids.contains_key(reverse_id) {
            return Err(ChemistryError::UnknownReaction(reverse_id.to_string()));
        }
    }
    Ok(())
}

pub fn evaluate_redox_potential(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    reaction: &Reaction,
) -> ChemistryResult<Option<RedoxPotentialEvaluation>> {
    let Some(redox) = &reaction.redox else {
        return Ok(None);
    };
    let (Some(oxidation_half_id), Some(reduction_half_id)) =
        (&redox.oxidation_half_id, &redox.reduction_half_id)
    else {
        return Ok(None);
    };
    let oxidation = registry
        .redox_half_reaction(oxidation_half_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(oxidation_half_id.clone()))?;
    let reduction = registry
        .redox_half_reaction(reduction_half_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(reduction_half_id.clone()))?;
    let (Some(oxidation_standard), Some(reduction_standard)) = (
        oxidation.standard_potential_volts,
        reduction.standard_potential_volts,
    ) else {
        return Ok(None);
    };
    let oxidation_potential =
        half_potential_volts(registry, mixture, oxidation, oxidation_standard)?;
    let reduction_potential =
        half_potential_volts(registry, mixture, reduction, reduction_standard)?;
    let cell_potential = oxidation_potential + reduction_potential;
    let electron_count = redox.transferred_electrons as f64;
    let exponent = electron_count * FARADAY_CONSTANT_COULOMBS_PER_MOL * cell_potential
        / (GAS_CONSTANT_J_PER_MOL_KELVIN * mixture.temperature_kelvin());
    let equilibrium_constant = if exponent > 700.0 {
        f64::INFINITY
    } else {
        exponent.exp()
    };
    let thermodynamic_rate_factor = if exponent > 40.0 {
        1.0
    } else if exponent < -40.0 {
        0.0
    } else {
        let value = 1.0 / (1.0 + (-exponent).exp());
        if value <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            0.0
        } else {
            value
        }
    };
    if !cell_potential.is_finite()
        || !thermodynamic_rate_factor.is_finite()
        || !(0.0..=1.0).contains(&thermodynamic_rate_factor)
    {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction.id.to_string(),
            reason: "redox potential evaluation became invalid".to_string(),
        });
    }
    Ok(Some(RedoxPotentialEvaluation {
        oxidation_half_id: oxidation_half_id.clone(),
        reduction_half_id: reduction_half_id.clone(),
        oxidation_potential_volts: oxidation_potential,
        reduction_potential_volts: reduction_potential,
        cell_potential_volts: cell_potential,
        equilibrium_constant,
        thermodynamic_rate_factor,
    }))
}

pub fn apply_electrolysis_cell(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    cell: &ElectrolysisCell,
    current_amperes: f64,
    duration_seconds: f64,
) -> ChemistryResult<ElectrolysisReport> {
    validate_electrolysis_cell(cell, current_amperes, duration_seconds)?;
    let anode_half = registry
        .redox_half_reaction(&cell.anode.half_reaction_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(cell.anode.half_reaction_id.clone()))?;
    let cathode_half = registry
        .redox_half_reaction(&cell.cathode.half_reaction_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(cell.cathode.half_reaction_id.clone()))?;
    validate_electrode_half(cell, &cell.anode, anode_half)?;
    validate_electrode_half(cell, &cell.cathode, cathode_half)?;
    if !environments_compatible(anode_half.environment, cathode_half.environment) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<electrolysis-cell>".to_string(),
            reason: "electrode half reaction environments are incompatible".to_string(),
        });
    }

    let anode_potential_volts = half_potential_volts_in_phase(
        registry,
        mixture,
        anode_half,
        required_standard_potential(anode_half)?,
        cell.anode.phase,
    )?;
    let cathode_potential_volts = half_potential_volts_in_phase(
        registry,
        mixture,
        cathode_half,
        required_standard_potential(cathode_half)?,
        cell.cathode.phase,
    )?;
    let cell_potential_volts = anode_potential_volts + cathode_potential_volts;
    let reversible_voltage_volts = (-cell_potential_volts).max(0.0);
    let required_voltage_volts = reversible_voltage_volts
        + cell.anode.overpotential_volts
        + cell.cathode.overpotential_volts;
    if cell.applied_voltage_volts + 1.0e-12 < required_voltage_volts {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "applied voltage {} V is below required electrolysis voltage {} V",
            cell.applied_voltage_volts, required_voltage_volts
        )));
    }

    let requested_charge_coulombs = current_amperes * duration_seconds;
    let electron_moles = requested_charge_coulombs / FARADAY_CONSTANT_COULOMBS_PER_MOL;
    let anode_extent_mol_per_bucket =
        electron_moles * cell.anode.current_efficiency / anode_half.electron_count as f64;
    let cathode_extent_mol_per_bucket =
        electron_moles * cell.cathode.current_efficiency / cathode_half.electron_count as f64;

    let anode_reactants = electrode_reactants(registry, anode_half, cell.anode.phase)?;
    let anode_products = electrode_products(registry, mixture, anode_half, cell.anode.phase)?;
    let cathode_reactants = electrode_reactants(registry, cathode_half, cell.cathode.phase)?;
    let cathode_products = electrode_products(registry, mixture, cathode_half, cell.cathode.phase)?;
    let produced_credits = electrode_product_credits(
        &anode_products,
        anode_extent_mol_per_bucket,
        &cathode_products,
        cathode_extent_mol_per_bucket,
    );

    ensure_electrode_reactants_available(
        registry,
        mixture,
        anode_half,
        &anode_reactants,
        anode_extent_mol_per_bucket,
        &produced_credits,
    )?;
    ensure_electrode_reactants_available(
        registry,
        mixture,
        cathode_half,
        &cathode_reactants,
        cathode_extent_mol_per_bucket,
        &produced_credits,
    )?;

    let mut staged = mixture.clone();
    let deltas = electrode_phase_amount_deltas(
        &anode_reactants,
        &anode_products,
        anode_extent_mol_per_bucket,
        &cathode_reactants,
        &cathode_products,
        cathode_extent_mol_per_bucket,
    )?;
    staged.apply_phase_amount_deltas_by_index(registry, &deltas)?;
    *mixture = staged;

    Ok(ElectrolysisReport {
        requested_charge_coulombs,
        transferred_electrons_mol_per_bucket: electron_moles,
        anode_extent_mol_per_bucket,
        cathode_extent_mol_per_bucket,
        anode_potential_volts,
        cathode_potential_volts,
        reversible_voltage_volts,
        required_voltage_volts,
    })
}

fn half_potential_volts(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    half: &RedoxHalfReaction,
    standard_potential_volts: f64,
) -> ChemistryResult<f64> {
    half_potential_volts_in_phase(
        registry,
        mixture,
        half,
        standard_potential_volts,
        MixturePhase::Aqueous,
    )
}

fn half_potential_volts_in_phase(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    half: &RedoxHalfReaction,
    standard_potential_volts: f64,
    phase: MixturePhase,
) -> ChemistryResult<f64> {
    let quotient = half_reaction_quotient_in_phase(registry, mixture, half, phase)?;
    if quotient <= 0.0 || !quotient.is_finite() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: half.id.clone(),
            reason: "redox half reaction quotient must be positive and finite".to_string(),
        });
    }
    let value = standard_potential_volts
        - GAS_CONSTANT_J_PER_MOL_KELVIN * mixture.temperature_kelvin()
            / (half.electron_count as f64 * FARADAY_CONSTANT_COULOMBS_PER_MOL)
            * quotient.ln();
    if !value.is_finite() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: half.id.clone(),
            reason: "redox half potential became non-finite".to_string(),
        });
    }
    Ok(value)
}

fn validate_electrolysis_cell(
    cell: &ElectrolysisCell,
    current_amperes: f64,
    duration_seconds: f64,
) -> ChemistryResult<()> {
    if cell.anode.polarity != ElectrodePolarity::Anode {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<electrolysis-cell>".to_string(),
            reason: "anode process must have anode polarity".to_string(),
        });
    }
    if cell.cathode.polarity != ElectrodePolarity::Cathode {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<electrolysis-cell>".to_string(),
            reason: "cathode process must have cathode polarity".to_string(),
        });
    }
    validate_electrode_numbers(&cell.anode)?;
    validate_electrode_numbers(&cell.cathode)?;
    if !cell.applied_voltage_volts.is_finite() || cell.applied_voltage_volts < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "applied voltage must be non-negative and finite".to_string(),
        ));
    }
    if !current_amperes.is_finite() || current_amperes < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "electrolysis current must be non-negative and finite".to_string(),
        ));
    }
    if !duration_seconds.is_finite() || duration_seconds < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "electrolysis duration must be non-negative and finite".to_string(),
        ));
    }
    Ok(())
}

fn validate_electrode_numbers(process: &ElectrodeProcess) -> ChemistryResult<()> {
    if process.half_reaction_id.trim().is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<electrode-process>".to_string(),
            reason: "electrode half reaction id must not be empty".to_string(),
        });
    }
    if !process.current_efficiency.is_finite()
        || process.current_efficiency <= 0.0
        || process.current_efficiency > 1.0
    {
        return Err(ChemistryError::InvalidMixtureState(
            "electrode current efficiency must be within 0.0..=1.0 and greater than zero"
                .to_string(),
        ));
    }
    if !process.overpotential_volts.is_finite() || process.overpotential_volts < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "electrode overpotential must be non-negative and finite".to_string(),
        ));
    }
    Ok(())
}

fn validate_electrode_half(
    _cell: &ElectrolysisCell,
    process: &ElectrodeProcess,
    half: &RedoxHalfReaction,
) -> ChemistryResult<()> {
    let expected_side = match process.polarity {
        ElectrodePolarity::Anode => ElectronSide::Product,
        ElectrodePolarity::Cathode => ElectronSide::Reactant,
    };
    if half.electron_side != expected_side {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: process.half_reaction_id.clone(),
            reason: match process.polarity {
                ElectrodePolarity::Anode => {
                    "anode process must use an oxidation half reaction".to_string()
                }
                ElectrodePolarity::Cathode => {
                    "cathode process must use a reduction half reaction".to_string()
                }
            },
        });
    }
    if half.standard_potential_volts.is_none() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: process.half_reaction_id.clone(),
            reason: "electrode process requires a standard potential".to_string(),
        });
    }
    Ok(())
}

fn required_standard_potential(half: &RedoxHalfReaction) -> ChemistryResult<f64> {
    half.standard_potential_volts
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: half.id.clone(),
            reason: "electrode process requires a standard potential".to_string(),
        })
}

fn electrode_reactants(
    registry: &ChemistryRegistry,
    half: &RedoxHalfReaction,
    phase: MixturePhase,
) -> ChemistryResult<Vec<(super::registry::SubstanceIndex, u32, Vec<MixturePhase>)>> {
    half.reactants
        .iter()
        .map(|term| {
            let substance = registry
                .substance_index(&term.substance_id)
                .ok_or_else(|| {
                    ChemistryError::InvalidMixtureState(format!(
                        "unknown electrode reactant '{}'",
                        term.substance_id
                    ))
                })?;
            Ok((substance, term.coefficient, vec![phase]))
        })
        .collect()
}

fn electrode_products(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    half: &RedoxHalfReaction,
    phase: MixturePhase,
) -> ChemistryResult<Vec<(super::registry::SubstanceIndex, f64, MixturePhase)>> {
    half.products
        .iter()
        .map(|term| {
            let substance = registry
                .substance_index(&term.substance_id)
                .ok_or_else(|| {
                    ChemistryError::InvalidMixtureState(format!(
                        "unknown electrode product '{}'",
                        term.substance_id
                    ))
                })?;
            let product_phase = electrode_term_phase(registry, mixture, &term.substance_id, phase)?;
            Ok((substance, term.coefficient as f64, product_phase))
        })
        .collect()
}

fn ensure_electrode_reactants_available(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    half: &RedoxHalfReaction,
    reactants: &[(super::registry::SubstanceIndex, u32, Vec<MixturePhase>)],
    extent_mol_per_bucket: f64,
    produced_credits: &BTreeMap<(super::registry::SubstanceIndex, MixturePhase), f64>,
) -> ChemistryResult<()> {
    for (substance, coefficient, phases) in reactants {
        let required = *coefficient as f64 * extent_mol_per_bucket;
        let available = mixture.concentration_of_index_in_phases(*substance, phases);
        let internally_produced = phases
            .iter()
            .map(|phase| {
                produced_credits
                    .get(&(*substance, *phase))
                    .copied()
                    .unwrap_or(0.0)
            })
            .sum::<f64>();
        if available + internally_produced + TRACE_CONCENTRATION_MOL_PER_BUCKET < required {
            let substance_id = &registry.substance_by_index(*substance)?.id;
            return Err(ChemistryError::InvalidMixtureState(format!(
                "electrode half '{}' needs {} mol/bucket of '{}' but only {} is available",
                half.id,
                required,
                substance_id,
                available + internally_produced
            )));
        }
    }
    Ok(())
}

fn electrode_product_credits(
    anode_products: &[(super::registry::SubstanceIndex, f64, MixturePhase)],
    anode_extent_mol_per_bucket: f64,
    cathode_products: &[(super::registry::SubstanceIndex, f64, MixturePhase)],
    cathode_extent_mol_per_bucket: f64,
) -> BTreeMap<(super::registry::SubstanceIndex, MixturePhase), f64> {
    let mut credits = BTreeMap::new();
    for (substance, coefficient, phase) in anode_products {
        *credits.entry((*substance, *phase)).or_insert(0.0) +=
            *coefficient * anode_extent_mol_per_bucket;
    }
    for (substance, coefficient, phase) in cathode_products {
        *credits.entry((*substance, *phase)).or_insert(0.0) +=
            *coefficient * cathode_extent_mol_per_bucket;
    }
    credits
}

fn electrode_phase_amount_deltas(
    anode_reactants: &[(super::registry::SubstanceIndex, u32, Vec<MixturePhase>)],
    anode_products: &[(super::registry::SubstanceIndex, f64, MixturePhase)],
    anode_extent_mol_per_bucket: f64,
    cathode_reactants: &[(super::registry::SubstanceIndex, u32, Vec<MixturePhase>)],
    cathode_products: &[(super::registry::SubstanceIndex, f64, MixturePhase)],
    cathode_extent_mol_per_bucket: f64,
) -> ChemistryResult<Vec<(super::registry::SubstanceIndex, MixturePhase, f64)>> {
    let mut deltas = Vec::new();
    push_electrode_reactant_deltas(&mut deltas, anode_reactants, anode_extent_mol_per_bucket)?;
    push_electrode_product_deltas(&mut deltas, anode_products, anode_extent_mol_per_bucket);
    push_electrode_reactant_deltas(
        &mut deltas,
        cathode_reactants,
        cathode_extent_mol_per_bucket,
    )?;
    push_electrode_product_deltas(&mut deltas, cathode_products, cathode_extent_mol_per_bucket);
    Ok(deltas)
}

fn push_electrode_reactant_deltas(
    deltas: &mut Vec<(super::registry::SubstanceIndex, MixturePhase, f64)>,
    reactants: &[(super::registry::SubstanceIndex, u32, Vec<MixturePhase>)],
    extent_mol_per_bucket: f64,
) -> ChemistryResult<()> {
    for (substance, coefficient, phases) in reactants {
        let [phase] = phases.as_slice() else {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<electrolysis-cell>".to_string(),
                reason: "electrode reactant must have exactly one concrete phase".to_string(),
            });
        };
        deltas.push((
            *substance,
            *phase,
            -(*coefficient as f64) * extent_mol_per_bucket,
        ));
    }
    Ok(())
}

fn push_electrode_product_deltas(
    deltas: &mut Vec<(super::registry::SubstanceIndex, MixturePhase, f64)>,
    products: &[(super::registry::SubstanceIndex, f64, MixturePhase)],
    extent_mol_per_bucket: f64,
) {
    for (substance, coefficient, phase) in products {
        deltas.push((*substance, *phase, *coefficient * extent_mol_per_bucket));
    }
}

fn half_reaction_quotient_in_phase(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    half: &RedoxHalfReaction,
    phase: MixturePhase,
) -> ChemistryResult<f64> {
    let products = terms_activity_product_in_phase(registry, mixture, &half.products, phase)?;
    let reactants = terms_activity_product_in_phase(registry, mixture, &half.reactants, phase)?;
    let quotient = products / reactants;
    if !quotient.is_finite() || quotient <= 0.0 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: half.id.clone(),
            reason: "redox half reaction quotient must be positive and finite".to_string(),
        });
    }
    Ok(quotient)
}

fn terms_activity_product_in_phase(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    terms: &[StoichiometricTerm],
    phase: MixturePhase,
) -> ChemistryResult<f64> {
    let mut product = 1.0;
    for term in terms {
        let term_phase = electrode_term_phase(registry, mixture, &term.substance_id, phase)?;
        let activity = mixture
            .activity_of(registry, &term.substance_id, term_phase)?
            .max(TRACE_CONCENTRATION_MOL_PER_BUCKET);
        product *= activity.powi(term.coefficient as i32);
    }
    if !product.is_finite() || product <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(
            "redox activity product must be positive and finite".to_string(),
        ));
    }
    Ok(product)
}

fn electrode_term_phase(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    substance_id: &SubstanceId,
    electrode_phase: MixturePhase,
) -> ChemistryResult<MixturePhase> {
    let substance = registry.substance(substance_id)?;
    if matches!(
        electrode_phase,
        MixturePhase::MoltenMetal | MixturePhase::MoltenSlag
    ) && substance.charge != 0
    {
        return Ok(electrode_phase);
    }
    if substance.aggregate_state_at(mixture.temperature_kelvin())? == SubstanceAggregateState::Gas {
        return Ok(MixturePhase::Gas);
    }
    Ok(match substance.phase_properties.preferred_liquid_phase {
        LiquidPhasePreference::MoltenMetal => MixturePhase::MoltenMetal,
        LiquidPhasePreference::MoltenSlag => MixturePhase::MoltenSlag,
        LiquidPhasePreference::Aqueous => MixturePhase::Aqueous,
        LiquidPhasePreference::Organic => MixturePhase::Organic,
    })
}

pub(crate) fn validate_half_reaction_conservation<F>(
    half: &RedoxHalfReaction,
    mut mass_and_charge: F,
) -> ChemistryResult<()>
where
    F: FnMut(&SubstanceId) -> ChemistryResult<(f64, i32)>,
{
    validate_half_reaction_shape(half)?;
    let (reactant_mass, mut reactant_charge) = sum_terms(&half.reactants, &mut mass_and_charge)?;
    let (product_mass, mut product_charge) = sum_terms(&half.products, &mut mass_and_charge)?;
    match half.electron_side {
        ElectronSide::Reactant => {
            reactant_charge -= half.electron_count as i32;
        }
        ElectronSide::Product => {
            product_charge -= half.electron_count as i32;
        }
    }
    if (reactant_mass - product_mass).abs() > 1.0e-6 {
        return Err(ChemistryError::MassNotConserved {
            reaction_id: half.id.clone(),
            reactants: reactant_mass,
            products: product_mass,
        });
    }
    if reactant_charge != product_charge {
        return Err(ChemistryError::ChargeNotConserved {
            reaction_id: half.id.clone(),
            reactants: reactant_charge,
            products: product_charge,
        });
    }
    Ok(())
}

pub(crate) fn validate_redox_pair(
    pair: &RedoxPair,
    halves: &BTreeMap<String, RedoxHalfReaction>,
) -> ChemistryResult<()> {
    let oxidation = halves
        .get(&pair.oxidation_half_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(pair.oxidation_half_id.clone()))?;
    let reduction = halves
        .get(&pair.reduction_half_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(pair.reduction_half_id.clone()))?;
    if oxidation.electron_side != ElectronSide::Product
        || reduction.electron_side != ElectronSide::Reactant
    {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: pair.reaction.id.to_string(),
            reason: "redox pair must combine oxidation and reduction halves".to_string(),
        });
    }
    if oxidation.electron_count != reduction.electron_count {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: pair.reaction.id.to_string(),
            reason: "redox pair electron counts must match before closed reaction registration"
                .to_string(),
        });
    }
    if !environments_compatible(oxidation.environment, reduction.environment) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: pair.reaction.id.to_string(),
            reason: "redox half reaction environments are incompatible".to_string(),
        });
    }
    Ok(())
}

pub(crate) fn build_redox_pair_reaction(
    spec: &RedoxPairSpec,
    halves: &BTreeMap<String, RedoxHalfReaction>,
) -> ChemistryResult<Reaction> {
    let oxidation = halves
        .get(&spec.oxidation_half_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(spec.oxidation_half_id.clone()))?;
    let reduction = halves
        .get(&spec.reduction_half_id)
        .ok_or_else(|| ChemistryError::UnknownReaction(spec.reduction_half_id.clone()))?;
    if oxidation.electron_side != ElectronSide::Product
        || reduction.electron_side != ElectronSide::Reactant
    {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.reaction_id.to_string(),
            reason: "redox pair must combine oxidation and reduction halves".to_string(),
        });
    }
    if !environments_compatible(oxidation.environment, reduction.environment) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.reaction_id.to_string(),
            reason: "redox half reaction environments are incompatible".to_string(),
        });
    }
    if oxidation.standard_potential_volts.is_some() != reduction.standard_potential_volts.is_some()
    {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.reaction_id.to_string(),
            reason: "paired redox halves must either both define potentials or both omit them"
                .to_string(),
        });
    }

    let transferred_electrons = least_common_multiple(
        oxidation.electron_count,
        reduction.electron_count,
        &spec.reaction_id,
    )?;
    let oxidation_scale = transferred_electrons / oxidation.electron_count;
    let reduction_scale = transferred_electrons / reduction.electron_count;

    let mut left = BTreeMap::<SubstanceId, u64>::new();
    let mut right = BTreeMap::<SubstanceId, u64>::new();
    add_scaled_terms(
        &mut left,
        &oxidation.reactants,
        oxidation_scale,
        &spec.reaction_id,
    )?;
    add_scaled_terms(
        &mut right,
        &oxidation.products,
        oxidation_scale,
        &spec.reaction_id,
    )?;
    add_scaled_terms(
        &mut left,
        &reduction.reactants,
        reduction_scale,
        &spec.reaction_id,
    )?;
    add_scaled_terms(
        &mut right,
        &reduction.products,
        reduction_scale,
        &spec.reaction_id,
    )?;
    cancel_common_terms(&mut left, &mut right);
    if left.is_empty() || right.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: spec.reaction_id.to_string(),
            reason: "closed redox reaction must have reactants and products after cancellation"
                .to_string(),
        });
    }

    let environment = combined_environment(oxidation.environment, reduction.environment);
    let mut builder = Reaction::builder(spec.reaction_id.clone());
    for (substance_id, coefficient) in left {
        builder = builder.reactant(
            substance_id,
            checked_u32(coefficient, &spec.reaction_id)?,
            1,
        );
    }
    for (substance_id, coefficient) in right {
        builder = builder.product(substance_id, checked_u32(coefficient, &spec.reaction_id)?);
    }
    Ok(builder
        .redox_annotation(RedoxAnnotation::from_halves(
            transferred_electrons,
            environment,
            spec.oxidation_half_id.clone(),
            spec.reduction_half_id.clone(),
        ))
        .build())
}

pub fn electron_reactant(builder: ReactionBuilder, count: u32) -> ReactionBuilder {
    builder.chemical_external_reactant(ELECTRON_EXTERNAL_ID, count as f64, 0.0, -1)
}

pub fn electron_product(builder: ReactionBuilder, count: u32) -> ReactionBuilder {
    builder.chemical_external_product(ELECTRON_EXTERNAL_ID, count as f64, 0.0, -1)
}

fn sum_terms<F>(
    terms: &[StoichiometricTerm],
    mass_and_charge: &mut F,
) -> ChemistryResult<(f64, i32)>
where
    F: FnMut(&SubstanceId) -> ChemistryResult<(f64, i32)>,
{
    let mut mass = 0.0;
    let mut charge = 0;
    for term in terms {
        let (substance_mass, substance_charge) = mass_and_charge(&term.substance_id)?;
        mass += substance_mass * term.coefficient as f64;
        charge += substance_charge * term.coefficient as i32;
    }
    Ok((mass, charge))
}

fn environments_compatible(left: RedoxEnvironment, right: RedoxEnvironment) -> bool {
    left == RedoxEnvironment::Any || right == RedoxEnvironment::Any || left == right
}

fn combined_environment(left: RedoxEnvironment, right: RedoxEnvironment) -> RedoxEnvironment {
    match (left, right) {
        (RedoxEnvironment::Any, other) | (other, RedoxEnvironment::Any) => other,
        (same, _) => same,
    }
}

fn add_scaled_terms(
    target: &mut BTreeMap<SubstanceId, u64>,
    terms: &[StoichiometricTerm],
    scale: u32,
    reaction_id: &ReactionId,
) -> ChemistryResult<()> {
    for term in terms {
        let coefficient = u64::from(term.coefficient)
            .checked_mul(u64::from(scale))
            .ok_or_else(|| ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "redox coefficient overflow while scaling half reaction".to_string(),
            })?;
        *target.entry(term.substance_id.clone()).or_insert(0) += coefficient;
    }
    Ok(())
}

fn cancel_common_terms(
    left: &mut BTreeMap<SubstanceId, u64>,
    right: &mut BTreeMap<SubstanceId, u64>,
) {
    let common = left
        .keys()
        .filter(|substance_id| right.contains_key(*substance_id))
        .cloned()
        .collect::<Vec<_>>();
    for substance_id in common {
        let cancelled = left[&substance_id].min(right[&substance_id]);
        if let Some(value) = left.get_mut(&substance_id) {
            *value -= cancelled;
        }
        if let Some(value) = right.get_mut(&substance_id) {
            *value -= cancelled;
        }
        if left.get(&substance_id).copied() == Some(0) {
            left.remove(&substance_id);
        }
        if right.get(&substance_id).copied() == Some(0) {
            right.remove(&substance_id);
        }
    }
}

fn least_common_multiple(left: u32, right: u32, reaction_id: &ReactionId) -> ChemistryResult<u32> {
    let divisor = greatest_common_divisor(left, right);
    left.checked_div(divisor)
        .and_then(|quotient| quotient.checked_mul(right))
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: "redox electron count overflow while combining half reactions".to_string(),
        })
}

fn greatest_common_divisor(mut left: u32, mut right: u32) -> u32 {
    while right != 0 {
        let next = left % right;
        left = right;
        right = next;
    }
    left
}

fn checked_u32(value: u64, reaction_id: &ReactionId) -> ChemistryResult<u32> {
    u32::try_from(value).map_err(|_| ChemistryError::InvalidReaction {
        reaction_id: reaction_id.to_string(),
        reason: "redox coefficient does not fit into u32".to_string(),
    })
}

fn electronegativity(element: &str) -> ChemistryResult<f64> {
    let value = match element {
        "H" => 2.20,
        "B" => 2.04,
        "C" => 2.55,
        "N" => 3.04,
        "O" => 3.44,
        "F" => 3.98,
        "Al" => 1.61,
        "Li" => 0.98,
        "Na" => 0.93,
        "Mg" => 1.31,
        "P" => 2.19,
        "S" => 2.58,
        "Cl" => 3.16,
        "K" => 0.82,
        "Ca" => 1.00,
        "Br" => 2.96,
        "Cr" => 1.66,
        "Fe" => 1.83,
        "Ni" => 1.91,
        "Cu" => 1.90,
        "Zn" => 1.65,
        "Zr" => 1.33,
        "I" => 2.66,
        "Pt" => 2.28,
        "Au" => 2.54,
        "Hg" => 2.00,
        "Pb" => 2.33,
        "Si" => 1.90,
        "Ar" => 0.0,
        "R" => 2.55,
        _ => {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: "<oxidation-state>".to_string(),
                reason: format!("no electronegativity for element '{element}'"),
            });
        }
    };
    Ok(value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;
    use crate::chemistry::mixture::Mixture;
    use crate::chemistry::molecule::{
        MolecularAtom, MolecularBond, MolecularStructure, ValenceSaturation,
    };
    use crate::chemistry::reaction::Reaction;
    use crate::chemistry::registry::ChemistryRegistryBuilder;
    use crate::chemistry::simulation::{
        reaction_rate_mol_per_bucket_per_tick, reaction_rate_mol_per_bucket_per_tick_with_context,
        ReactionContext,
    };
    use crate::chemistry::substance::Substance;

    fn redox_test_registry() -> ChemistryRegistryBuilder {
        ChemistryRegistryBuilder::new()
            .substance(Substance::new(
                "destroy:water",
                0,
                18.03,
                18_000.0,
                373.0,
                75.0,
                40_650.0,
            ))
            .substance(Substance::new(
                "destroy:proton",
                1,
                1.01,
                1_000.0,
                10_000.0,
                0.0,
                0.0,
            ))
            .substance(Substance::new(
                "destroy:hydroxide",
                -1,
                17.01,
                17_000.0,
                10_000.0,
                75.0,
                0.0,
            ))
            .substance(Substance::new(
                "destroy:hydrogen_peroxide",
                0,
                34.04,
                34_000.0,
                423.0,
                90.0,
                20_000.0,
            ))
            .substance(
                Substance::new("destroy:iron_ii", 2, 55.85, 55_850.0, 10_000.0, 80.0, 0.0)
                    .with_explicit_oxidation_states(vec![ExplicitOxidationState::new(
                        "Fe", 2.0, 1,
                    )]),
            )
            .substance(
                Substance::new("destroy:iron_iii", 3, 55.85, 55_850.0, 10_000.0, 80.0, 0.0)
                    .with_explicit_oxidation_states(vec![ExplicitOxidationState::new(
                        "Fe", 3.0, 1,
                    )]),
            )
    }

    fn water_structure() -> MolecularStructure {
        MolecularStructure {
            source_code: "test:water".to_string(),
            atoms: vec![
                MolecularAtom {
                    element: "O".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
            ],
            bonds: vec![
                MolecularBond {
                    from: 0,
                    to: 1,
                    order: 1.0,
                },
                MolecularBond {
                    from: 0,
                    to: 2,
                    order: 1.0,
                },
            ],
            stereochemistry: Vec::new(),
        }
    }

    fn peroxide_structure() -> MolecularStructure {
        MolecularStructure {
            source_code: "test:peroxide".to_string(),
            atoms: vec![
                MolecularAtom {
                    element: "O".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
                MolecularAtom {
                    element: "O".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                    valence_saturation: ValenceSaturation::Saturate,
                },
            ],
            bonds: vec![
                MolecularBond {
                    from: 0,
                    to: 1,
                    order: 1.0,
                },
                MolecularBond {
                    from: 0,
                    to: 2,
                    order: 1.0,
                },
                MolecularBond {
                    from: 1,
                    to: 3,
                    order: 1.0,
                },
            ],
            stereochemistry: Vec::new(),
        }
    }

    #[test]
    fn oxidation_states_cover_water_peroxide_and_organic_examples() {
        let water = assign_oxidation_states(&water_structure()).unwrap();
        assert!((water.state_sum_for_element("O") + 2.0).abs() < 1.0e-9);
        assert!((water.state_sum_for_element("H") - 2.0).abs() < 1.0e-9);

        let peroxide = assign_oxidation_states(&peroxide_structure()).unwrap();
        assert!((peroxide.state_sum_for_element("O") + 2.0).abs() < 1.0e-9);

        let ethanol = assign_oxidation_states(&parse_frowns("CCO").unwrap()).unwrap();
        let carbon_states = ethanol
            .atom_states
            .iter()
            .filter(|state| state.element == "C")
            .map(|state| state.state)
            .collect::<Vec<_>>();
        assert!(carbon_states
            .iter()
            .any(|state| (*state + 3.0).abs() < 1.0e-9));
        assert!(carbon_states
            .iter()
            .any(|state| (*state + 1.0).abs() < 1.0e-9));

        let acetone = assign_oxidation_states(&parse_frowns("CC(=O)C").unwrap()).unwrap();
        assert!(acetone
            .atom_states
            .iter()
            .any(|state| state.element == "C" && (state.state - 2.0).abs() < 1.0e-9));
    }

    #[test]
    fn explicit_inorganic_oxidation_states_are_checked_against_charge() {
        let assignment = explicit_oxidation_assignment(
            &"destroy:chromate".into(),
            -2,
            &[
                ExplicitOxidationState::new("Cr", 6.0, 1),
                ExplicitOxidationState::new("O", -2.0, 4),
            ],
        )
        .unwrap();
        assert_eq!(assignment.total_charge, -2);

        let error = explicit_oxidation_assignment(
            &"destroy:bad_chromate".into(),
            -1,
            &[
                ExplicitOxidationState::new("Cr", 6.0, 1),
                ExplicitOxidationState::new("O", -2.0, 4),
            ],
        )
        .unwrap_err();
        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }

    #[test]
    fn redox_half_reactions_conserve_mass_and_charge_with_electrons() {
        let registry = redox_test_registry()
            .redox_half_reaction(RedoxHalfReaction::oxidation(
                "iron_ii_to_iron_iii",
                [("destroy:iron_ii".into(), 2)],
                [("destroy:iron_iii".into(), 2)],
                2,
                RedoxEnvironment::Acidic,
            ))
            .redox_half_reaction(RedoxHalfReaction::reduction(
                "peroxide_to_water",
                [
                    ("destroy:hydrogen_peroxide".into(), 1),
                    ("destroy:proton".into(), 2),
                ],
                [("destroy:water".into(), 2)],
                2,
                RedoxEnvironment::Acidic,
            ))
            .build()
            .unwrap();
        assert_eq!(registry.redox_half_reactions().count(), 2);
    }

    #[test]
    fn invalid_redox_half_reaction_fails_registry_build() {
        let error = redox_test_registry()
            .redox_half_reaction(RedoxHalfReaction::oxidation(
                "broken_iron_half",
                [("destroy:iron_ii".into(), 1)],
                [("destroy:iron_iii".into(), 1)],
                2,
                RedoxEnvironment::Acidic,
            ))
            .build()
            .unwrap_err();
        assert!(matches!(error, ChemistryError::ChargeNotConserved { .. }));
    }

    #[test]
    fn redox_pair_registers_closed_reaction_without_free_electrons() {
        let oxidation = RedoxHalfReaction::oxidation(
            "iron_ii_to_iron_iii",
            [("destroy:iron_ii".into(), 2)],
            [("destroy:iron_iii".into(), 2)],
            2,
            RedoxEnvironment::Acidic,
        );
        let reduction = RedoxHalfReaction::reduction(
            "peroxide_to_water",
            [
                ("destroy:hydrogen_peroxide".into(), 1),
                ("destroy:proton".into(), 2),
            ],
            [("destroy:water".into(), 2)],
            2,
            RedoxEnvironment::Acidic,
        );
        let reaction = Reaction::builder("destroy:iron_peroxide_redox")
            .reactant("destroy:iron_ii", 2, 1)
            .reactant("destroy:hydrogen_peroxide", 1, 1)
            .reactant("destroy:proton", 2, 1)
            .product("destroy:iron_iii", 2)
            .product("destroy:water", 2)
            .build();
        let registry = redox_test_registry()
            .redox_half_reaction(oxidation)
            .redox_half_reaction(reduction)
            .redox_pair(RedoxPair::new(
                reaction,
                "iron_ii_to_iron_iii",
                "peroxide_to_water",
                2,
                RedoxEnvironment::Acidic,
            ))
            .build()
            .unwrap();
        let reaction = registry
            .reaction(&"destroy:iron_peroxide_redox".into())
            .unwrap();
        assert!(reaction.redox.is_some());
        assert!(reaction
            .external_reactants
            .iter()
            .all(|external| external.description != ELECTRON_EXTERNAL_ID));
    }

    #[test]
    fn redox_pair_from_halves_scales_and_cancels_electrons() {
        let registry = redox_test_registry()
            .redox_half_reaction(
                RedoxHalfReaction::oxidation(
                    "iron_ii_to_iron_iii",
                    [("destroy:iron_ii".into(), 1)],
                    [("destroy:iron_iii".into(), 1)],
                    1,
                    RedoxEnvironment::Acidic,
                )
                .with_standard_potential_volts(-0.771),
            )
            .redox_half_reaction(
                RedoxHalfReaction::reduction(
                    "peroxide_to_water",
                    [
                        ("destroy:hydrogen_peroxide".into(), 1),
                        ("destroy:proton".into(), 2),
                    ],
                    [("destroy:water".into(), 2)],
                    2,
                    RedoxEnvironment::Acidic,
                )
                .with_standard_potential_volts(1.776),
            )
            .redox_pair_from_halves(
                "destroy:generated_iron_peroxide_redox",
                "iron_ii_to_iron_iii",
                "peroxide_to_water",
            )
            .build()
            .unwrap();

        let reaction = registry
            .reaction(&"destroy:generated_iron_peroxide_redox".into())
            .unwrap();
        assert_eq!(
            reaction
                .reactants
                .iter()
                .find(|term| term.substance_id == "destroy:iron_ii".into())
                .unwrap()
                .coefficient,
            2
        );
        assert_eq!(
            reaction
                .reactants
                .iter()
                .find(|term| term.substance_id == "destroy:hydrogen_peroxide".into())
                .unwrap()
                .coefficient,
            1
        );
        assert_eq!(
            reaction
                .products
                .iter()
                .find(|term| term.substance_id == "destroy:iron_iii".into())
                .unwrap()
                .coefficient,
            2
        );
        assert_eq!(reaction.redox.as_ref().unwrap().transferred_electrons, 2);
        assert!(reaction.external_reactants.is_empty());
        assert!(reaction.external_products.is_empty());

        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_ii", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_iii", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:hydrogen_peroxide", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:proton", 0.1)
            .unwrap();
        let evaluation = evaluate_redox_potential(&registry, &mixture, reaction)
            .unwrap()
            .unwrap();
        assert!(evaluation.cell_potential_volts > 0.0);
        assert!(evaluation.thermodynamic_rate_factor > 0.0);
    }

    #[test]
    fn redox_reaction_rejects_charge_imbalance_escape_hatch_and_free_electrons() {
        let bad_charge = redox_test_registry()
            .reaction(
                Reaction::builder("destroy:bad_redox")
                    .reactant("destroy:iron_ii", 1, 1)
                    .product("destroy:iron_iii", 1)
                    .redox_annotation(RedoxAnnotation::checked(
                        "destroy:iron_iii",
                        "destroy:iron_ii",
                        1,
                        RedoxEnvironment::Any,
                    ))
                    .allow_charge_imbalance()
                    .build(),
            )
            .build()
            .unwrap_err();
        assert!(matches!(bad_charge, ChemistryError::InvalidReaction { .. }));

        let free_electron = redox_test_registry()
            .reaction(
                Reaction::builder("destroy:free_electron_redox")
                    .reactant("destroy:iron_iii", 1, 1)
                    .product("destroy:iron_ii", 1)
                    .electron_reactant(1)
                    .redox_annotation(RedoxAnnotation::checked(
                        "destroy:iron_iii",
                        "destroy:iron_ii",
                        1,
                        RedoxEnvironment::Any,
                    ))
                    .build(),
            )
            .build()
            .unwrap_err();
        assert!(matches!(
            free_electron,
            ChemistryError::InvalidReaction { .. }
        ));
    }

    #[test]
    fn redox_environment_controls_reaction_rate() {
        let registry = redox_test_registry()
            .reaction(
                Reaction::builder("destroy:acidic_redox")
                    .reactant("destroy:iron_ii", 2, 1)
                    .reactant("destroy:hydrogen_peroxide", 1, 1)
                    .reactant("destroy:proton", 2, 1)
                    .product("destroy:iron_iii", 2)
                    .product("destroy:water", 2)
                    .redox_annotation(RedoxAnnotation::checked(
                        "destroy:hydrogen_peroxide",
                        "destroy:iron_ii",
                        2,
                        RedoxEnvironment::Acidic,
                    ))
                    .build(),
            )
            .build()
            .unwrap();
        let reaction = registry.reaction(&"destroy:acidic_redox".into()).unwrap();

        let mut without_acid = Mixture::new(298.0).unwrap();
        without_acid
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        without_acid
            .add_substance(&registry, "destroy:iron_ii", 0.2)
            .unwrap();
        without_acid
            .add_substance(&registry, "destroy:hydrogen_peroxide", 0.2)
            .unwrap();
        assert_eq!(
            reaction_rate_mol_per_bucket_per_tick(&registry, &without_acid, reaction).unwrap(),
            0.0
        );

        let mut acidic = without_acid.clone();
        acidic
            .add_substance(&registry, "destroy:proton", 0.2)
            .unwrap();
        assert!(reaction_rate_mol_per_bucket_per_tick(&registry, &acidic, reaction).unwrap() > 0.0);
    }

    #[test]
    fn redox_potential_uses_nernst_activity_and_ph() {
        let oxidation = RedoxHalfReaction::oxidation(
            "iron_ii_to_iron_iii",
            [("destroy:iron_ii".into(), 2)],
            [("destroy:iron_iii".into(), 2)],
            2,
            RedoxEnvironment::Acidic,
        )
        .with_standard_potential_volts(-0.771);
        let reduction = RedoxHalfReaction::reduction(
            "peroxide_to_water",
            [
                ("destroy:hydrogen_peroxide".into(), 1),
                ("destroy:proton".into(), 2),
            ],
            [("destroy:water".into(), 2)],
            2,
            RedoxEnvironment::Acidic,
        )
        .with_standard_potential_volts(1.776);
        let reaction = Reaction::builder("destroy:iron_peroxide_redox")
            .reactant("destroy:iron_ii", 2, 1)
            .reactant("destroy:hydrogen_peroxide", 1, 1)
            .reactant("destroy:proton", 2, 1)
            .product("destroy:iron_iii", 2)
            .product("destroy:water", 2)
            .redox_annotation(RedoxAnnotation::from_halves(
                2,
                RedoxEnvironment::Acidic,
                "iron_ii_to_iron_iii",
                "peroxide_to_water",
            ))
            .build();
        let registry = redox_test_registry()
            .redox_half_reaction(oxidation)
            .redox_half_reaction(reduction)
            .reaction(reaction)
            .build()
            .unwrap();
        let reaction = registry
            .reaction(&"destroy:iron_peroxide_redox".into())
            .unwrap();

        let mut acidic = Mixture::new(298.0).unwrap();
        acidic
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        acidic
            .add_substance(&registry, "destroy:iron_ii", 0.1)
            .unwrap();
        acidic
            .add_substance(&registry, "destroy:iron_iii", 0.1)
            .unwrap();
        acidic
            .add_substance(&registry, "destroy:hydrogen_peroxide", 0.1)
            .unwrap();
        acidic
            .add_substance(&registry, "destroy:proton", 0.1)
            .unwrap();

        let mut weakly_acidic = Mixture::new(298.0).unwrap();
        weakly_acidic
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        weakly_acidic
            .add_substance(&registry, "destroy:iron_ii", 0.1)
            .unwrap();
        weakly_acidic
            .add_substance(&registry, "destroy:iron_iii", 0.1)
            .unwrap();
        weakly_acidic
            .add_substance(&registry, "destroy:hydrogen_peroxide", 0.1)
            .unwrap();
        weakly_acidic
            .add_substance(&registry, "destroy:proton", 1.0e-6)
            .unwrap();

        let acidic_potential = evaluate_redox_potential(&registry, &acidic, reaction)
            .unwrap()
            .unwrap();
        let weak_potential = evaluate_redox_potential(&registry, &weakly_acidic, reaction)
            .unwrap()
            .unwrap();

        assert!(acidic_potential.cell_potential_volts > weak_potential.cell_potential_volts);
        assert!(acidic_potential.cell_potential_volts > 0.0);
    }

    #[test]
    fn redox_potential_suppresses_unfavorable_closed_reaction_without_context_hack() {
        let oxidation = RedoxHalfReaction::oxidation(
            "iron_ii_to_iron_iii",
            [("destroy:iron_ii".into(), 1)],
            [("destroy:iron_iii".into(), 1)],
            1,
            RedoxEnvironment::Any,
        )
        .with_standard_potential_volts(-0.771);
        let reduction = RedoxHalfReaction::reduction(
            "iron_iii_to_iron_ii",
            [("destroy:iron_iii".into(), 1)],
            [("destroy:iron_ii".into(), 1)],
            1,
            RedoxEnvironment::Any,
        )
        .with_standard_potential_volts(0.1);
        let reaction = Reaction::builder("destroy:unfavorable_electron_exchange")
            .reactant("destroy:iron_ii", 1, 1)
            .reactant("destroy:iron_iii", 1, 1)
            .product("destroy:iron_iii", 1)
            .product("destroy:iron_ii", 1)
            .redox_annotation(RedoxAnnotation::from_halves(
                1,
                RedoxEnvironment::Any,
                "iron_ii_to_iron_iii",
                "iron_iii_to_iron_ii",
            ))
            .pre_exponential_factor(1.0e12)
            .activation_energy_kj_per_mol(0.0)
            .build();
        let registry = redox_test_registry()
            .redox_half_reaction(oxidation)
            .redox_half_reaction(reduction)
            .reaction(reaction)
            .build()
            .unwrap();
        let reaction = registry
            .reaction(&"destroy:unfavorable_electron_exchange".into())
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_ii", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_iii", 0.1)
            .unwrap();

        let evaluation = evaluate_redox_potential(&registry, &mixture, reaction)
            .unwrap()
            .unwrap();
        let rate = reaction_rate_mol_per_bucket_per_tick_with_context(
            &registry,
            &mixture,
            reaction,
            &ReactionContext::default(),
        )
        .unwrap();

        assert!(evaluation.cell_potential_volts < 0.0);
        assert_eq!(evaluation.thermodynamic_rate_factor, 0.0);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn paired_redox_potentials_must_be_complete() {
        let error = redox_test_registry()
            .redox_half_reaction(
                RedoxHalfReaction::oxidation(
                    "iron_ii_to_iron_iii",
                    [("destroy:iron_ii".into(), 1)],
                    [("destroy:iron_iii".into(), 1)],
                    1,
                    RedoxEnvironment::Any,
                )
                .with_standard_potential_volts(-0.771),
            )
            .redox_half_reaction(RedoxHalfReaction::reduction(
                "iron_iii_to_iron_ii",
                [("destroy:iron_iii".into(), 1)],
                [("destroy:iron_ii".into(), 1)],
                1,
                RedoxEnvironment::Any,
            ))
            .reaction(
                Reaction::builder("destroy:incomplete_redox_potential")
                    .reactant("destroy:iron_ii", 1, 1)
                    .reactant("destroy:iron_iii", 1, 1)
                    .product("destroy:iron_iii", 1)
                    .product("destroy:iron_ii", 1)
                    .redox_annotation(RedoxAnnotation::from_halves(
                        1,
                        RedoxEnvironment::Any,
                        "iron_ii_to_iron_iii",
                        "iron_iii_to_iron_ii",
                    ))
                    .build(),
            )
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidReaction { .. }));
    }

    #[test]
    fn electrolysis_cell_applies_half_reactions_by_faraday_law() {
        let oxidation = RedoxHalfReaction::oxidation(
            "iron_ii_to_iron_iii",
            [("destroy:iron_ii".into(), 1)],
            [("destroy:iron_iii".into(), 1)],
            1,
            RedoxEnvironment::Acidic,
        )
        .with_standard_potential_volts(-0.771);
        let reduction = RedoxHalfReaction::reduction(
            "iron_iii_to_iron_ii",
            [("destroy:iron_iii".into(), 1)],
            [("destroy:iron_ii".into(), 1)],
            1,
            RedoxEnvironment::Acidic,
        )
        .with_standard_potential_volts(-0.4);
        let registry = redox_test_registry()
            .redox_half_reaction(oxidation)
            .redox_half_reaction(reduction)
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_ii", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_iii", 0.01)
            .unwrap();

        let cell = ElectrolysisCell::new(
            ElectrodeProcess::anode("iron_ii_to_iron_iii"),
            ElectrodeProcess::cathode("iron_iii_to_iron_ii").with_current_efficiency(0.5),
            1.3,
        );
        let report =
            apply_electrolysis_cell(&registry, &mut mixture, &cell, 96.485_332_123_310_02, 1.0)
                .unwrap();

        assert!((report.transferred_electrons_mol_per_bucket - 0.001).abs() < 1.0e-12);
        assert!((mixture.concentration_of(&"destroy:iron_ii".into()) - 0.0995).abs() < 1.0e-9);
        assert!((mixture.concentration_of(&"destroy:iron_iii".into()) - 0.0105).abs() < 1.0e-9);
        assert!(report.reversible_voltage_volts > 1.1);
    }

    #[test]
    fn electrolysis_cell_rejects_missing_voltage_without_mutating_mixture() {
        let oxidation = RedoxHalfReaction::oxidation(
            "iron_ii_to_iron_iii",
            [("destroy:iron_ii".into(), 1)],
            [("destroy:iron_iii".into(), 1)],
            1,
            RedoxEnvironment::Any,
        )
        .with_standard_potential_volts(-0.8);
        let reduction = RedoxHalfReaction::reduction(
            "iron_iii_to_iron_ii",
            [("destroy:iron_iii".into(), 1)],
            [("destroy:iron_ii".into(), 1)],
            1,
            RedoxEnvironment::Any,
        )
        .with_standard_potential_volts(-0.4);
        let registry = redox_test_registry()
            .redox_half_reaction(oxidation)
            .redox_half_reaction(reduction)
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_ii", 0.1)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_iii", 0.1)
            .unwrap();
        let iron_ii_before = mixture.concentration_of(&"destroy:iron_ii".into());
        let iron_iii_before = mixture.concentration_of(&"destroy:iron_iii".into());
        let cell = ElectrolysisCell::new(
            ElectrodeProcess::anode("iron_ii_to_iron_iii"),
            ElectrodeProcess::cathode("iron_iii_to_iron_ii"),
            0.5,
        );

        let error = apply_electrolysis_cell(&registry, &mut mixture, &cell, 1.0, 1.0).unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidMixtureState(_)));
        assert_eq!(
            mixture.concentration_of(&"destroy:iron_ii".into()),
            iron_ii_before
        );
        assert_eq!(
            mixture.concentration_of(&"destroy:iron_iii".into()),
            iron_iii_before
        );
    }

    #[test]
    fn electrolysis_cell_rejects_wrong_electrode_polarity() {
        let oxidation = RedoxHalfReaction::oxidation(
            "iron_ii_to_iron_iii",
            [("destroy:iron_ii".into(), 1)],
            [("destroy:iron_iii".into(), 1)],
            1,
            RedoxEnvironment::Any,
        )
        .with_standard_potential_volts(-0.771);
        let registry = redox_test_registry()
            .redox_half_reaction(oxidation)
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:iron_ii", 0.1)
            .unwrap();
        let cell = ElectrolysisCell::new(
            ElectrodeProcess::cathode("iron_ii_to_iron_iii"),
            ElectrodeProcess::cathode("iron_ii_to_iron_iii"),
            1.0,
        );

        let error = apply_electrolysis_cell(&registry, &mut mixture, &cell, 1.0, 1.0).unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidReaction { .. }));
    }
}
