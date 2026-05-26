use std::fmt::{Display, Formatter};

use super::error::{ChemistryError, ChemistryResult};
use super::functional_group::{find_functional_groups, FunctionalGroup};
use super::molecule::MolecularStructure;
use super::redox::{
    assign_oxidation_states, explicit_oxidation_assignment, ExplicitOxidationState, RedoxRole,
};

const MOLECULAR_MASS_TOLERANCE_GRAMS_PER_MOL: f64 = 1.0e-6;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SubstanceId(String);

impl SubstanceId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: value,
                reason: "id must not be empty".to_string(),
            });
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for SubstanceId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for SubstanceId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SubstanceTagId(String);

impl SubstanceTagId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: "<tag>".to_string(),
                reason: "tag id must not be empty".to_string(),
            });
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for SubstanceTagId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for SubstanceTagId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

#[derive(Debug, Clone)]
pub struct Substance {
    pub id: SubstanceId,
    pub charge: i32,
    pub molar_mass_grams: f64,
    pub liquid_density_grams_per_bucket: f64,
    pub solid_density_grams_per_bucket: f64,
    pub melting_point_kelvin: f64,
    pub boiling_point_kelvin: f64,
    pub molar_heat_capacity_j_per_mol_kelvin: f64,
    pub fusion_heat_j_per_mol: f64,
    pub latent_heat_j_per_mol: f64,
    pub structure_code: Option<String>,
    pub molecular_structure: Option<MolecularStructure>,
    pub functional_groups: Vec<FunctionalGroup>,
    pub translation_key: Option<String>,
    pub color_argb: u32,
    pub tags: Vec<SubstanceTagId>,
    pub phase_properties: SubstancePhaseBehavior,
    pub redox_roles: Vec<RedoxRole>,
    pub explicit_oxidation_states: Vec<ExplicitOxidationState>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SubstanceAggregateState {
    Solid,
    Liquid,
    Gas,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum LiquidPhasePreference {
    Aqueous,
    Organic,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SolventRole {
    NotSolvent,
    KnownSolvent,
    ConservativePredictedSolvent,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SubstancePhaseBehavior {
    pub preferred_liquid_phase: LiquidPhasePreference,
    pub aqueous_solubility_mol_per_bucket: Option<f64>,
    pub organic_solubility_mol_per_bucket: Option<f64>,
    pub can_precipitate: bool,
    pub can_form_liquid_phase: bool,
    pub solvent_role: SolventRole,
}

pub type SubstancePhaseProperties = SubstancePhaseBehavior;

impl SubstancePhaseBehavior {
    pub fn aqueous_unlimited() -> Self {
        Self {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: None,
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: false,
            can_form_liquid_phase: true,
            solvent_role: SolventRole::NotSolvent,
        }
    }

    pub fn aqueous_solvent() -> Self {
        Self {
            solvent_role: SolventRole::KnownSolvent,
            ..Self::aqueous_unlimited()
        }
    }

    pub fn organic_unlimited(aqueous_solubility_mol_per_bucket: f64) -> Self {
        Self {
            preferred_liquid_phase: LiquidPhasePreference::Organic,
            aqueous_solubility_mol_per_bucket: Some(aqueous_solubility_mol_per_bucket),
            organic_solubility_mol_per_bucket: None,
            can_precipitate: false,
            can_form_liquid_phase: true,
            solvent_role: SolventRole::KnownSolvent,
        }
    }

    pub fn organic_solute(aqueous_solubility_mol_per_bucket: f64) -> Self {
        Self {
            preferred_liquid_phase: LiquidPhasePreference::Organic,
            aqueous_solubility_mol_per_bucket: Some(aqueous_solubility_mol_per_bucket),
            organic_solubility_mol_per_bucket: None,
            can_precipitate: false,
            can_form_liquid_phase: false,
            solvent_role: SolventRole::NotSolvent,
        }
    }

    pub fn validate(&self, substance_id: &SubstanceId, charge: i32) -> ChemistryResult<()> {
        validate_solubility_limit(
            substance_id,
            "aqueous solubility",
            self.aqueous_solubility_mol_per_bucket,
        )?;
        validate_solubility_limit(
            substance_id,
            "organic solubility",
            self.organic_solubility_mol_per_bucket,
        )?;
        if charge != 0 && self.preferred_liquid_phase != LiquidPhasePreference::Aqueous {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "charged substances must prefer the aqueous phase".to_string(),
            });
        }
        if charge != 0 && self.solvent_role != SolventRole::NotSolvent {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "charged substances must not be solvents".to_string(),
            });
        }
        if self.solvent_role != SolventRole::NotSolvent && !self.can_form_liquid_phase {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: "solvents must be able to form a liquid phase".to_string(),
            });
        }
        Ok(())
    }
}

impl Substance {
    pub fn new(
        id: impl Into<SubstanceId>,
        charge: i32,
        molar_mass_grams: f64,
        liquid_density_grams_per_bucket: f64,
        boiling_point_kelvin: f64,
        molar_heat_capacity_j_per_mol_kelvin: f64,
        latent_heat_j_per_mol: f64,
    ) -> Self {
        let id = id.into();
        let phase_properties = if id.as_str() == "destroy:water" || id.as_str() == "water" {
            SubstancePhaseProperties::aqueous_solvent()
        } else if charge == 0 {
            SubstancePhaseProperties::organic_unlimited(0.05)
        } else {
            SubstancePhaseProperties {
                preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                aqueous_solubility_mol_per_bucket: Some(10.0),
                organic_solubility_mol_per_bucket: Some(0.0),
                can_precipitate: true,
                can_form_liquid_phase: false,
                solvent_role: SolventRole::NotSolvent,
            }
        };
        Self {
            id,
            charge,
            molar_mass_grams,
            liquid_density_grams_per_bucket,
            solid_density_grams_per_bucket: liquid_density_grams_per_bucket,
            melting_point_kelvin: 0.0,
            boiling_point_kelvin,
            molar_heat_capacity_j_per_mol_kelvin,
            fusion_heat_j_per_mol: 0.0,
            latent_heat_j_per_mol,
            structure_code: None,
            molecular_structure: None,
            functional_groups: Vec::new(),
            translation_key: None,
            color_argb: 0x20FF_FFFF,
            tags: Vec::new(),
            phase_properties,
            redox_roles: Vec::new(),
            explicit_oxidation_states: Vec::new(),
        }
    }

    pub fn with_catalog_metadata(
        mut self,
        structure_code: Option<String>,
        translation_key: Option<String>,
        color_argb: u32,
        tags: Vec<SubstanceTagId>,
    ) -> Self {
        self.structure_code = structure_code;
        self.translation_key = translation_key;
        self.color_argb = color_argb;
        self.tags = tags;
        self
    }

    pub fn with_phase_properties(mut self, phase_properties: SubstancePhaseProperties) -> Self {
        self.phase_properties = phase_properties;
        self
    }

    pub fn with_solvent_role(mut self, solvent_role: SolventRole) -> Self {
        self.phase_properties.solvent_role = solvent_role;
        if solvent_role != SolventRole::NotSolvent {
            self.phase_properties.can_form_liquid_phase = true;
        }
        self
    }

    pub fn with_melting_point_kelvin(mut self, melting_point_kelvin: f64) -> Self {
        self.melting_point_kelvin = melting_point_kelvin;
        self
    }

    pub fn with_fusion_heat_j_per_mol(mut self, fusion_heat_j_per_mol: f64) -> Self {
        self.fusion_heat_j_per_mol = fusion_heat_j_per_mol;
        self
    }

    pub fn with_solid_density_grams_per_bucket(
        mut self,
        solid_density_grams_per_bucket: f64,
    ) -> Self {
        self.solid_density_grams_per_bucket = solid_density_grams_per_bucket;
        self
    }

    pub fn with_redox_roles(mut self, redox_roles: Vec<RedoxRole>) -> Self {
        self.redox_roles = redox_roles;
        self
    }

    pub fn with_explicit_oxidation_states(
        mut self,
        explicit_oxidation_states: Vec<ExplicitOxidationState>,
    ) -> Self {
        self.explicit_oxidation_states = explicit_oxidation_states;
        self
    }

    pub fn aggregate_state_at(
        &self,
        temperature_kelvin: f64,
    ) -> ChemistryResult<SubstanceAggregateState> {
        if !temperature_kelvin.is_finite() || temperature_kelvin < 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: self.id.to_string(),
                reason: "temperature must be non-negative and finite".to_string(),
            });
        }
        if temperature_kelvin > self.boiling_point_kelvin {
            Ok(SubstanceAggregateState::Gas)
        } else if temperature_kelvin < self.melting_point_kelvin {
            Ok(SubstanceAggregateState::Solid)
        } else {
            Ok(SubstanceAggregateState::Liquid)
        }
    }

    pub fn with_molecular_structure(mut self, molecular_structure: MolecularStructure) -> Self {
        self.functional_groups = find_functional_groups(&molecular_structure);
        self.molecular_structure = Some(molecular_structure);
        self
    }

    pub fn validate(&self) -> ChemistryResult<()> {
        let id = self.id.to_string();
        if !self.molar_mass_grams.is_finite() || self.molar_mass_grams <= 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "molar mass must be positive and finite".to_string(),
            });
        }
        if !self.liquid_density_grams_per_bucket.is_finite()
            || self.liquid_density_grams_per_bucket <= 0.0
        {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "liquid density must be positive and finite".to_string(),
            });
        }
        if !self.solid_density_grams_per_bucket.is_finite()
            || self.solid_density_grams_per_bucket <= 0.0
        {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "solid density must be positive and finite".to_string(),
            });
        }
        if !self.melting_point_kelvin.is_finite() || self.melting_point_kelvin < 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "melting point must be non-negative and finite".to_string(),
            });
        }
        if !self.boiling_point_kelvin.is_finite() || self.boiling_point_kelvin < 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "boiling point must be non-negative and finite".to_string(),
            });
        }
        if self.melting_point_kelvin > self.boiling_point_kelvin {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "melting point must not be above boiling point".to_string(),
            });
        }
        if !self.molar_heat_capacity_j_per_mol_kelvin.is_finite()
            || self.molar_heat_capacity_j_per_mol_kelvin < 0.0
        {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "molar heat capacity must be non-negative and finite".to_string(),
            });
        }
        if !self.latent_heat_j_per_mol.is_finite() || self.latent_heat_j_per_mol < 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "latent heat must be non-negative and finite".to_string(),
            });
        }
        if !self.fusion_heat_j_per_mol.is_finite() || self.fusion_heat_j_per_mol < 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "fusion heat must be non-negative and finite".to_string(),
            });
        }
        self.phase_properties.validate(&self.id, self.charge)?;
        if !self.explicit_oxidation_states.is_empty() {
            explicit_oxidation_assignment(&self.id, self.charge, &self.explicit_oxidation_states)?;
        }
        if let Some(structure) = &self.molecular_structure {
            let summary =
                structure
                    .summary()
                    .map_err(|error| ChemistryError::InvalidSubstance {
                        substance_id: self.id.to_string(),
                        reason: format!("invalid molecular structure: {error}"),
                    })?;
            if summary.charge != self.charge {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: format!(
                        "molecular structure charge {} does not match substance charge {}",
                        summary.charge, self.charge
                    ),
                });
            }
            if (summary.molar_mass_grams - self.molar_mass_grams).abs()
                > MOLECULAR_MASS_TOLERANCE_GRAMS_PER_MOL
            {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: format!(
                        "molecular structure mass {} does not match substance mass {}",
                        summary.molar_mass_grams, self.molar_mass_grams
                    ),
                });
            }
            let assignment = assign_oxidation_states(structure).map_err(|error| {
                ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: format!("invalid oxidation state assignment: {error}"),
                }
            })?;
            if assignment.total_charge != self.charge {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: format!(
                        "oxidation states sum to charge {}, expected {}",
                        assignment.total_charge, self.charge
                    ),
                });
            }
        }
        Ok(())
    }
}

fn validate_solubility_limit(
    substance_id: &SubstanceId,
    name: &str,
    value: Option<f64>,
) -> ChemistryResult<()> {
    if let Some(value) = value {
        if !value.is_finite() || value < 0.0 {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: substance_id.to_string(),
                reason: format!("{name} must be non-negative and finite"),
            });
        }
    }
    Ok(())
}
