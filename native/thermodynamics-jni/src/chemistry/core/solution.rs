use std::collections::BTreeMap;

use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{Mixture, MixturePhase, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::registry::{ChemistryRegistry, SubstanceIndex};
use super::substance::SubstanceId;

const DAVIES_A_AT_298_K: f64 = 0.509;
const PH_PROTON_ACTIVITY_FLOOR: f64 = 1.0e-14;
const EQUILIBRIUM_QUOTIENT_TOLERANCE: f64 = 1.0e-8;
const EQUILIBRIUM_MIN_EXTENT_MOL_PER_BUCKET: f64 = 1.0e-12;
const EQUILIBRIUM_MAX_PASSES: usize = 256;

#[derive(Debug, Clone, PartialEq)]
pub struct SolutionState {
    pub aqueous_ionic_strength_mol_per_bucket: f64,
    pub proton_activity_mol_per_bucket: Option<f64>,
    pub ph: Option<f64>,
    pub activity_coefficients: BTreeMap<SubstanceId, f64>,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum ActivityModel {
    Davies,
}

impl Default for ActivityModel {
    fn default() -> Self {
        Self::Davies
    }
}

impl ActivityModel {
    pub fn coefficient(
        self,
        charge: i32,
        ionic_strength_mol_per_bucket: f64,
        temperature_kelvin: f64,
    ) -> ChemistryResult<f64> {
        if !ionic_strength_mol_per_bucket.is_finite() || ionic_strength_mol_per_bucket < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "ionic strength must be non-negative and finite".to_string(),
            ));
        }
        if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "temperature must be positive and finite for activity".to_string(),
            ));
        }
        if charge == 0 || ionic_strength_mol_per_bucket == 0.0 {
            return Ok(1.0);
        }
        match self {
            ActivityModel::Davies => {
                let sqrt_i = ionic_strength_mol_per_bucket.sqrt();
                let raw_term = sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_strength_mol_per_bucket;
                let term = raw_term.max(0.0);
                let temperature_factor = (298.0 / temperature_kelvin).sqrt();
                let log10_gamma =
                    -DAVIES_A_AT_298_K * temperature_factor * (charge * charge) as f64 * term;
                let coefficient = 10.0_f64.powf(log10_gamma);
                if !coefficient.is_finite() || !(0.0..=1.0).contains(&coefficient) {
                    return Err(ChemistryError::InvalidMixtureState(format!(
                        "activity coefficient must be within 0.0..=1.0: {coefficient}"
                    )));
                }
                Ok(coefficient)
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AcidBaseSpec {
    pub id: String,
    pub acid: SubstanceId,
    pub conjugate_base: SubstanceId,
    pub pka: f64,
    pub proton: SubstanceId,
}

impl AcidBaseSpec {
    pub fn new(
        id: impl Into<String>,
        acid: impl Into<SubstanceId>,
        conjugate_base: impl Into<SubstanceId>,
        pka: f64,
    ) -> Self {
        Self {
            id: id.into(),
            acid: acid.into(),
            conjugate_base: conjugate_base.into(),
            pka,
            proton: "destroy:proton".into(),
        }
    }

    pub fn with_proton(mut self, proton: impl Into<SubstanceId>) -> Self {
        self.proton = proton.into();
        self
    }

    pub(crate) fn to_equilibria(&self) -> [EquilibriumSpec; 2] {
        let acid_constant = 10.0_f64.powf(-self.pka);
        [
            EquilibriumSpec::new(
                format!("{}.acid_base_equilibrium", self.id),
                [(self.acid.clone(), 1, MixturePhase::Aqueous)],
                [
                    (self.proton.clone(), 1, MixturePhase::Aqueous),
                    (self.conjugate_base.clone(), 1, MixturePhase::Aqueous),
                ],
                acid_constant,
            ),
            EquilibriumSpec::new(
                format!("{}.neutralization_equilibrium", self.id),
                [
                    (self.acid.clone(), 1, MixturePhase::Aqueous),
                    (
                        SubstanceId::from("destroy:hydroxide"),
                        1,
                        MixturePhase::Aqueous,
                    ),
                ],
                [
                    (self.conjugate_base.clone(), 1, MixturePhase::Aqueous),
                    (SubstanceId::from("destroy:water"), 1, MixturePhase::Aqueous),
                ],
                acid_constant / 1.0e-14,
            ),
        ]
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumTerm {
    pub substance_id: SubstanceId,
    pub coefficient: u32,
    pub phase: MixturePhase,
}

impl EquilibriumTerm {
    pub fn new(
        substance_id: impl Into<SubstanceId>,
        coefficient: u32,
        phase: MixturePhase,
    ) -> Self {
        Self {
            substance_id: substance_id.into(),
            coefficient,
            phase,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumSpec {
    pub id: String,
    pub reactants: Vec<EquilibriumTerm>,
    pub products: Vec<EquilibriumTerm>,
    pub equilibrium_constant: f64,
    pub reference_temperature_kelvin: f64,
    pub enthalpy_change_kj_per_mol: f64,
}

impl EquilibriumSpec {
    pub fn new(
        id: impl Into<String>,
        reactants: impl IntoIterator<Item = (SubstanceId, u32, MixturePhase)>,
        products: impl IntoIterator<Item = (SubstanceId, u32, MixturePhase)>,
        equilibrium_constant: f64,
    ) -> Self {
        Self {
            id: id.into(),
            reactants: reactants
                .into_iter()
                .map(|(substance_id, coefficient, phase)| {
                    EquilibriumTerm::new(substance_id, coefficient, phase)
                })
                .collect(),
            products: products
                .into_iter()
                .map(|(substance_id, coefficient, phase)| {
                    EquilibriumTerm::new(substance_id, coefficient, phase)
                })
                .collect(),
            equilibrium_constant,
            reference_temperature_kelvin: 298.0,
            enthalpy_change_kj_per_mol: 0.0,
        }
    }

    pub fn with_reference_temperature_kelvin(mut self, value: f64) -> Self {
        self.reference_temperature_kelvin = value;
        self
    }

    pub fn with_enthalpy_change_kj_per_mol(mut self, value: f64) -> Self {
        self.enthalpy_change_kj_per_mol = value;
        self
    }

    pub(crate) fn constant_at(&self, temperature_kelvin: f64) -> ChemistryResult<f64> {
        if !temperature_kelvin.is_finite() || temperature_kelvin <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "temperature must be positive and finite for equilibrium".to_string(),
            ));
        }
        if self.enthalpy_change_kj_per_mol == 0.0 {
            return Ok(self.equilibrium_constant);
        }
        let factor = (-self.enthalpy_change_kj_per_mol * 1000.0
            / super::reaction::GAS_CONSTANT_J_PER_MOL_KELVIN
            * (1.0 / temperature_kelvin - 1.0 / self.reference_temperature_kelvin))
            .exp();
        let value = self.equilibrium_constant * factor;
        if !value.is_finite() || value <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "equilibrium constant '{}' became invalid at current temperature",
                self.id
            )));
        }
        Ok(value)
    }
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedEquilibriumTerm {
    pub substance: SubstanceIndex,
    pub coefficient: u32,
    pub phase: MixturePhase,
}

#[derive(Debug, Clone)]
pub(crate) struct IndexedEquilibrium {
    pub spec: EquilibriumSpec,
    pub reactants: Vec<IndexedEquilibriumTerm>,
    pub products: Vec<IndexedEquilibriumTerm>,
}

pub fn solution_state(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
) -> ChemistryResult<SolutionState> {
    let ionic_strength = mixture.aqueous_ionic_strength(registry)?;
    let mut coefficients = BTreeMap::new();
    for substance in registry.substances() {
        let coefficient = ActivityModel::default().coefficient(
            substance.charge,
            ionic_strength,
            mixture.temperature_kelvin(),
        )?;
        coefficients.insert(substance.id.clone(), coefficient);
    }

    let proton: SubstanceId = "destroy:proton".into();
    let has_aqueous_phase =
        mixture.total_in_phase(MixturePhase::Aqueous) > TRACE_CONCENTRATION_MOL_PER_BUCKET;
    if !has_aqueous_phase {
        return Ok(SolutionState {
            aqueous_ionic_strength_mol_per_bucket: ionic_strength,
            proton_activity_mol_per_bucket: None,
            ph: None,
            activity_coefficients: coefficients,
        });
    }
    let proton_activity = activity_of(registry, mixture, &proton, MixturePhase::Aqueous)?
        .max(PH_PROTON_ACTIVITY_FLOOR);
    Ok(SolutionState {
        aqueous_ionic_strength_mol_per_bucket: ionic_strength,
        proton_activity_mol_per_bucket: Some(proton_activity),
        ph: Some(-proton_activity.log10()),
        activity_coefficients: coefficients,
    })
}

pub fn activity_of(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    substance_id: &SubstanceId,
    phase: MixturePhase,
) -> ChemistryResult<f64> {
    let substance = registry.substance(substance_id)?;
    if phase == MixturePhase::Aqueous && substance_id.as_str() == "destroy:water" {
        return Ok((mixture.concentration_in_phase(substance_id, phase)
            > TRACE_CONCENTRATION_MOL_PER_BUCKET)
            .then_some(1.0)
            .unwrap_or(0.0));
    }
    let concentration = mixture.concentration_in_phase(substance_id, phase);
    if phase != MixturePhase::Aqueous || substance.charge == 0 {
        return Ok(concentration);
    }
    let coefficient = ActivityModel::default().coefficient(
        substance.charge,
        mixture.aqueous_ionic_strength(registry)?,
        mixture.temperature_kelvin(),
    )?;
    Ok(concentration * coefficient)
}

pub(crate) fn equilibrate_solution_equilibria(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
) -> ChemistryResult<f64> {
    let mut max_delta = 0.0_f64;
    for _ in 0..EQUILIBRIUM_MAX_PASSES {
        let mut pass_delta = 0.0_f64;
        for equilibrium in registry.indexed_equilibria() {
            pass_delta = pass_delta.max(apply_equilibrium(registry, mixture, equilibrium)?);
        }
        max_delta = max_delta.max(pass_delta);
        if pass_delta <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Ok(max_delta);
        }
    }
    Err(ChemistryError::EquilibriumInvariantViolation {
        equilibrium_id: "<all>".to_string(),
        reason: "equilibrium solver did not reach a fixed point".to_string(),
    })
}

fn apply_equilibrium(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    equilibrium: &IndexedEquilibrium,
) -> ChemistryResult<f64> {
    let constant = equilibrium.spec.constant_at(mixture.temperature_kelvin())?;
    let quotient = reaction_quotient(registry, mixture, equilibrium)?;
    if quotient.is_none() {
        return Ok(0.0);
    }
    let quotient = quotient.unwrap_or(0.0);
    let relative_error = ((quotient / constant).ln()).abs();
    if relative_error <= EQUILIBRIUM_QUOTIENT_TOLERANCE {
        return Ok(0.0);
    }
    let forward = quotient < constant;
    let max_extent = max_equilibrium_extent(mixture, equilibrium, forward);
    if max_extent <= EQUILIBRIUM_MIN_EXTENT_MOL_PER_BUCKET {
        return Ok(0.0);
    }
    let mut low = 0.0;
    let mut high = max_extent;
    let initial_distance = equilibrium_distance(constant, quotient);
    for _ in 0..80 {
        let mid = (low + high) * 0.5;
        let trial = trial_quotient(registry, mixture, equilibrium, forward, mid)?;
        let distance = equilibrium_distance(constant, trial);
        if distance < initial_distance && quotient_is_on_target_side(trial, constant, forward) {
            low = mid;
        } else {
            high = mid;
        }
    }
    let extent = low;
    if extent <= EQUILIBRIUM_MIN_EXTENT_MOL_PER_BUCKET {
        return Ok(0.0);
    }
    apply_equilibrium_extent(registry, mixture, equilibrium, forward, extent)
}

fn reaction_quotient(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    equilibrium: &IndexedEquilibrium,
) -> ChemistryResult<Option<f64>> {
    let products = term_activity_product(registry, mixture, &equilibrium.products)?;
    let reactants = term_activity_product(registry, mixture, &equilibrium.reactants)?;
    if products == 0.0 && reactants == 0.0 {
        return Ok(None);
    }
    if reactants == 0.0 {
        return Ok(Some(f64::INFINITY));
    }
    Ok(Some(products / reactants))
}

fn term_activity_product(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    terms: &[IndexedEquilibriumTerm],
) -> ChemistryResult<f64> {
    let mut product = 1.0;
    for term in terms {
        let substance = registry.substance_by_index(term.substance)?;
        let activity = activity_of(registry, mixture, &substance.id, term.phase)?;
        product *= activity.powi(term.coefficient as i32);
    }
    Ok(product)
}

fn max_equilibrium_extent(
    mixture: &Mixture,
    equilibrium: &IndexedEquilibrium,
    forward: bool,
) -> f64 {
    let limiting_terms = if forward {
        &equilibrium.reactants
    } else {
        &equilibrium.products
    };
    limiting_terms
        .iter()
        .map(|term| {
            mixture.concentration_of_index_in_phases(term.substance, &[term.phase])
                / term.coefficient as f64
        })
        .fold(f64::INFINITY, f64::min)
}

fn trial_quotient(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    equilibrium: &IndexedEquilibrium,
    forward: bool,
    extent: f64,
) -> ChemistryResult<f64> {
    let mut cloned = mixture.clone();
    apply_equilibrium_extent(registry, &mut cloned, equilibrium, forward, extent)?;
    Ok(reaction_quotient(registry, &cloned, equilibrium)?.unwrap_or(0.0))
}

fn quotient_is_on_target_side(value: f64, constant: f64, forward: bool) -> bool {
    if forward {
        value <= constant
    } else {
        value >= constant
    }
}

fn equilibrium_distance(constant: f64, quotient: f64) -> f64 {
    if quotient == 0.0 {
        return f64::INFINITY;
    }
    if quotient.is_infinite() {
        return f64::INFINITY;
    }
    ((quotient / constant).ln()).abs()
}

fn apply_equilibrium_extent(
    registry: &ChemistryRegistry,
    mixture: &mut Mixture,
    equilibrium: &IndexedEquilibrium,
    forward: bool,
    extent: f64,
) -> ChemistryResult<f64> {
    let (reactants, products) = if forward {
        (&equilibrium.reactants, &equilibrium.products)
    } else {
        (&equilibrium.products, &equilibrium.reactants)
    };
    let reactants = reactants
        .iter()
        .map(|term| (term.substance, term.coefficient, vec![term.phase]))
        .collect::<Vec<_>>();
    let products = products
        .iter()
        .map(|term| (term.substance, term.coefficient, term.phase))
        .collect::<Vec<_>>();
    mixture.apply_reaction_phase_deltas_by_index(registry, &reactants, &products, extent)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::registry::ChemistryRegistryBuilder;
    use crate::chemistry::substance::{Substance, SubstancePhaseProperties};

    fn aqueous_substance(id: &str, charge: i32, mass: f64) -> Substance {
        Substance::new(id, charge, mass, 1_000.0, 10_000.0, 75.0, 20_000.0)
            .with_phase_properties(SubstancePhaseProperties::aqueous_unlimited())
    }

    fn acid_registry() -> ChemistryRegistry {
        ChemistryRegistryBuilder::new()
            .substance(aqueous_substance("destroy:water", 0, 18.0))
            .substance(aqueous_substance("destroy:proton", 1, 1.0))
            .substance(aqueous_substance("destroy:hydroxide", -1, 17.0))
            .substance(aqueous_substance("destroy:strong_acid", 0, 11.0))
            .substance(aqueous_substance("destroy:strong_base", -1, 10.0))
            .substance(aqueous_substance("destroy:weak_acid", 0, 11.0))
            .substance(aqueous_substance("destroy:weak_base", -1, 10.0))
            .substance(aqueous_substance("destroy:salt_cation", 1, 22.0))
            .substance(aqueous_substance("destroy:salt_anion", -1, 35.0))
            .equilibrium(EquilibriumSpec::new(
                "destroy:water.autoionization",
                [(SubstanceId::from("destroy:water"), 1, MixturePhase::Aqueous)],
                [
                    (
                        SubstanceId::from("destroy:proton"),
                        1,
                        MixturePhase::Aqueous,
                    ),
                    (
                        SubstanceId::from("destroy:hydroxide"),
                        1,
                        MixturePhase::Aqueous,
                    ),
                ],
                1.0e-14,
            ))
            .acid_base_pair(AcidBaseSpec::new(
                "destroy:strong_acid",
                "destroy:strong_acid",
                "destroy:strong_base",
                -6.0,
            ))
            .acid_base_pair(AcidBaseSpec::new(
                "destroy:weak_acid",
                "destroy:weak_acid",
                "destroy:weak_base",
                4.76,
            ))
            .build()
            .unwrap()
    }

    #[test]
    fn pure_water_reaches_neutral_ph_after_autoionization() {
        let registry = acid_registry();
        let mut mixture = Mixture::new(298.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();

        equilibrate_solution_equilibria(&registry, &mut mixture).unwrap();

        let ph = mixture.ph(&registry).unwrap().unwrap();
        assert!((ph - 7.0).abs() < 0.1, "pH was {ph}");
    }

    #[test]
    fn strong_acid_and_strong_base_move_ph_in_opposite_directions() {
        let registry = acid_registry();
        let mut acid = Mixture::new(298.0).unwrap();
        acid.add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        acid.add_substance(&registry, "destroy:strong_acid", 0.1)
            .unwrap();
        equilibrate_solution_equilibria(&registry, &mut acid).unwrap();
        assert!(acid.ph(&registry).unwrap().unwrap() < 2.0);

        let mut base = Mixture::new(298.0).unwrap();
        base.add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        base.add_substance(&registry, "destroy:hydroxide", 0.1)
            .unwrap();
        equilibrate_solution_equilibria(&registry, &mut base).unwrap();
        assert!(base.ph(&registry).unwrap().unwrap() > 12.0);
    }

    #[test]
    fn weak_acid_buffer_resists_added_proton_and_hydroxide() {
        let registry = acid_registry();
        let mut buffer = Mixture::new(298.0).unwrap();
        buffer
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        buffer
            .add_substance(&registry, "destroy:weak_acid", 0.1)
            .unwrap();
        buffer
            .add_substance(&registry, "destroy:weak_base", 0.1)
            .unwrap();
        equilibrate_solution_equilibria(&registry, &mut buffer).unwrap();
        let initial = buffer.ph(&registry).unwrap().unwrap();

        let mut acidified = buffer.clone();
        acidified
            .add_substance(&registry, "destroy:proton", 0.01)
            .unwrap();
        equilibrate_solution_equilibria(&registry, &mut acidified).unwrap();
        let acidified_ph = acidified.ph(&registry).unwrap().unwrap();

        let mut basified = buffer;
        basified
            .add_substance(&registry, "destroy:hydroxide", 0.01)
            .unwrap();
        equilibrate_solution_equilibria(&registry, &mut basified).unwrap();
        let basified_ph = basified.ph(&registry).unwrap().unwrap();

        assert!(acidified_ph < initial);
        assert!(basified_ph > initial);
        assert!(initial > 4.0 && initial < 6.0, "buffer pH was {initial}");
        assert!(initial - acidified_ph < 1.0);
        assert!(basified_ph - initial < 1.0);
    }

    #[test]
    fn ionic_strength_lowers_ion_activity() {
        let registry = acid_registry();
        let mut dilute = Mixture::new(298.0).unwrap();
        dilute
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        dilute
            .add_substance(&registry, "destroy:proton", 0.01)
            .unwrap();

        let mut salty = dilute.clone();
        salty
            .add_substance(&registry, "destroy:salt_cation", 0.5)
            .unwrap();
        salty
            .add_substance(&registry, "destroy:salt_anion", 0.5)
            .unwrap();

        assert!(
            salty.aqueous_ionic_strength(&registry).unwrap()
                > dilute.aqueous_ionic_strength(&registry).unwrap()
        );
        assert!(
            salty
                .activity_of(&registry, &"destroy:proton".into(), MixturePhase::Aqueous)
                .unwrap()
                < dilute
                    .activity_of(&registry, &"destroy:proton".into(), MixturePhase::Aqueous)
                    .unwrap()
        );
    }

    #[test]
    fn invalid_acid_base_pair_fails_registry_build() {
        let error = ChemistryRegistryBuilder::new()
            .substance(aqueous_substance("destroy:proton", 1, 1.0))
            .substance(aqueous_substance("destroy:bad_acid", 0, 10.0))
            .substance(aqueous_substance("destroy:bad_base", 0, 9.0))
            .acid_base_pair(AcidBaseSpec::new(
                "destroy:bad_pair",
                "destroy:bad_acid",
                "destroy:bad_base",
                4.0,
            ))
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::ChargeNotConserved { .. }));
    }

    #[test]
    fn invalid_equilibrium_fails_registry_build() {
        let error = ChemistryRegistryBuilder::new()
            .substance(aqueous_substance("destroy:a", 0, 10.0))
            .substance(aqueous_substance("destroy:b", 0, 11.0))
            .equilibrium(EquilibriumSpec::new(
                "destroy:bad_equilibrium",
                [(SubstanceId::from("destroy:a"), 1, MixturePhase::Aqueous)],
                [(SubstanceId::from("destroy:b"), 1, MixturePhase::Aqueous)],
                0.0,
            ))
            .build()
            .unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidReaction { .. }));
    }
}
