use super::error::{ChemistryError, ChemistryResult};
use super::mixture::MixturePhase;
use super::solution::EquilibriumSpec;
use super::substance::{
    LiquidPhasePreference, SolventRole, Substance, SubstanceId, SubstancePhaseProperties,
    SubstanceTagId,
};

#[derive(Debug, Clone, PartialEq)]
pub struct ComplexLigand {
    pub substance_id: SubstanceId,
    pub count: u32,
}

impl ComplexLigand {
    pub fn new(substance_id: impl Into<SubstanceId>, count: u32) -> Self {
        Self {
            substance_id: substance_id.into(),
            count,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ComplexSpec {
    pub id: SubstanceId,
    pub central_ion: SubstanceId,
    pub ligands: Vec<ComplexLigand>,
    pub charge: i32,
    pub formation_constant: f64,
    pub phase: MixturePhase,
    pub translation_key: Option<String>,
    pub color_argb: u32,
    pub tags: Vec<SubstanceTagId>,
}

impl ComplexSpec {
    pub fn new(
        id: impl Into<SubstanceId>,
        central_ion: impl Into<SubstanceId>,
        ligands: impl IntoIterator<Item = ComplexLigand>,
        charge: i32,
        formation_constant: f64,
    ) -> Self {
        Self {
            id: id.into(),
            central_ion: central_ion.into(),
            ligands: ligands.into_iter().collect(),
            charge,
            formation_constant,
            phase: MixturePhase::Aqueous,
            translation_key: None,
            color_argb: 0x20FF_FFFF,
            tags: Vec::new(),
        }
    }

    pub fn with_phase(mut self, phase: MixturePhase) -> Self {
        self.phase = phase;
        self
    }

    pub fn with_translation_key(mut self, translation_key: impl Into<String>) -> Self {
        self.translation_key = Some(translation_key.into());
        self
    }

    pub fn with_color_argb(mut self, color_argb: u32) -> Self {
        self.color_argb = color_argb;
        self
    }

    pub fn with_tags(mut self, tags: impl IntoIterator<Item = SubstanceTagId>) -> Self {
        self.tags = tags.into_iter().collect();
        self
    }

    pub fn validate_shape(&self) -> ChemistryResult<()> {
        if self.id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: "<complex>".to_string(),
                reason: "complex id must not be empty".to_string(),
            });
        }
        if self.central_ion.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: self.id.to_string(),
                reason: "complex central ion must not be empty".to_string(),
            });
        }
        if self.ligands.is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: self.id.to_string(),
                reason: "complex must contain at least one ligand".to_string(),
            });
        }
        if !self.formation_constant.is_finite() || self.formation_constant <= 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: format!("{}.formation", self.id),
                reason: "complex formation constant must be positive and finite".to_string(),
            });
        }
        if self.phase == MixturePhase::Gas || self.phase == MixturePhase::Solid {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: format!("{}.formation", self.id),
                reason: "complex formation equilibrium must use a liquid phase".to_string(),
            });
        }
        for ligand in &self.ligands {
            if ligand.count == 0 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: format!("{}.formation", self.id),
                    reason: "complex ligand count must be greater than zero".to_string(),
                });
            }
            if ligand.substance_id.as_str().trim().is_empty() {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: "complex ligand id must not be empty".to_string(),
                });
            }
        }
        Ok(())
    }

    pub(crate) fn to_equilibrium(&self) -> EquilibriumSpec {
        let reactants = std::iter::once((self.central_ion.clone(), 1, self.phase)).chain(
            self.ligands
                .iter()
                .map(|ligand| (ligand.substance_id.clone(), ligand.count, self.phase)),
        );
        EquilibriumSpec::new(
            format!("{}.formation", self.id),
            reactants,
            [(self.id.clone(), 1, self.phase)],
            self.formation_constant,
        )
    }

    pub(crate) fn to_substance(
        &self,
        molar_mass_grams: f64,
        charge: i32,
    ) -> ChemistryResult<Substance> {
        if charge != self.charge {
            return Err(ChemistryError::ChargeNotConserved {
                reaction_id: format!("{}.formation", self.id),
                reactants: charge,
                products: self.charge,
            });
        }
        Ok(Substance::new(
            self.id.clone(),
            self.charge,
            molar_mass_grams,
            1000.0,
            f64::MAX,
            100.0,
            20_000.0,
        )
        .with_catalog_metadata(
            None,
            self.translation_key.clone(),
            self.color_argb,
            self.tags.clone(),
        )
        .with_phase_properties(SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: Some(10.0),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: false,
            solvent_role: SolventRole::NotSolvent,
        }))
    }
}
