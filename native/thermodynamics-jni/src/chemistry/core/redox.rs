use std::collections::BTreeMap;

use super::error::{ChemistryError, ChemistryResult};
use super::molecule::MolecularStructure;
use super::reaction::{Reaction, ReactionBuilder, ReactionId, StoichiometricTerm};
use super::substance::SubstanceId;

pub const ELECTRON_EXTERNAL_ID: &str = "redox:electron";

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
        }
    }
}

#[derive(Debug, Clone)]
pub struct RedoxPair {
    pub reaction: Reaction,
    pub oxidation_half_id: String,
    pub reduction_half_id: String,
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
        if oxidation.electron_count != redox.transferred_electrons {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "oxidation half electron count does not match redox annotation".to_string(),
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
        if reduction.electron_count != redox.transferred_electrons {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction.id.to_string(),
                reason: "reduction half electron count does not match redox annotation".to_string(),
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
    use crate::chemistry::molecule::{MolecularAtom, MolecularBond, MolecularStructure};
    use crate::chemistry::reaction::Reaction;
    use crate::chemistry::registry::ChemistryRegistryBuilder;
    use crate::chemistry::simulation::reaction_rate_mol_per_bucket_per_tick;
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
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
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
                },
                MolecularAtom {
                    element: "O".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
                },
                MolecularAtom {
                    element: "H".to_string(),
                    charge: 0.0,
                    r_group_number: 0,
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
}
