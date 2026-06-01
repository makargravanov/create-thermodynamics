use std::fmt::{Display, Formatter};

use super::error::{ChemistryError, ChemistryResult};
use super::functional_group::{find_functional_groups, FunctionalGroup};
use super::molecule::MolecularStructure;
use super::redox::{
    assign_oxidation_states, explicit_oxidation_assignment, ExplicitOxidationState, RedoxRole,
};

const MOLECULAR_MASS_TOLERANCE_GRAMS_PER_MOL: f64 = 1.0e-6;
const GAS_CONSTANT_J_PER_MOL_KELVIN: f64 = 8.314_462_618_153_24;
const NORMAL_BOILING_PRESSURE_PASCAL: f64 = 101_325.0;

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
    pub representation: SubstanceRepresentation,
    pub charge: i32,
    pub molar_mass_grams: f64,
    pub liquid_density_grams_per_bucket: f64,
    pub solid_density_grams_per_bucket: f64,
    pub melting_point_kelvin: f64,
    pub boiling_point_kelvin: f64,
    pub molar_heat_capacity_j_per_mol_kelvin: f64,
    pub fusion_heat_j_per_mol: f64,
    pub latent_heat_j_per_mol: f64,
    pub critical_temperature_kelvin: Option<f64>,
    pub critical_pressure_pascal: Option<f64>,
    pub acentric_factor: Option<f64>,
    pub vapor_pressure_model: Option<VaporPressureModel>,
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

#[derive(Debug, Clone, PartialEq)]
pub enum SubstanceRepresentation {
    Molecular,
    Ion {
        parent_element: Option<String>,
    },
    IonicSolid {
        formula_units: Vec<MaterialFormulaUnit>,
    },
    Metal {
        element_symbol: String,
    },
    Oxide {
        formula_units: Vec<MaterialFormulaUnit>,
    },
    Hydrate {
        formula_units: Vec<MaterialFormulaUnit>,
        water_count: u32,
    },
    SurfaceMaterial {
        active_site: Option<String>,
    },
    UnspecifiedMaterial {
        reason: String,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaterialFormulaUnit {
    pub substance_id: SubstanceId,
    pub coefficient: u32,
}

impl MaterialFormulaUnit {
    pub fn new(substance_id: impl Into<SubstanceId>, coefficient: u32) -> Self {
        Self {
            substance_id: substance_id.into(),
            coefficient,
        }
    }
}

impl SubstanceRepresentation {
    pub fn validate(&self, substance: &Substance) -> ChemistryResult<()> {
        match self {
            SubstanceRepresentation::Molecular => Ok(()),
            SubstanceRepresentation::Ion { parent_element } => {
                if substance.charge == 0 {
                    return invalid_representation(
                        substance,
                        "ion representation requires non-zero charge",
                    );
                }
                validate_optional_non_empty_text(
                    substance,
                    parent_element.as_deref(),
                    "ion parent element",
                )
            }
            SubstanceRepresentation::IonicSolid { formula_units } => {
                if substance.charge != 0 {
                    return invalid_representation(substance, "ionic solids must be neutral");
                }
                if substance.molecular_structure.is_some() {
                    return invalid_representation(
                        substance,
                        "ionic solids must not pretend to be molecular graphs",
                    );
                }
                if !substance.phase_properties.can_precipitate {
                    return invalid_representation(
                        substance,
                        "ionic solids must be able to exist as a solid phase",
                    );
                }
                if substance.phase_properties.solvent_role != SolventRole::NotSolvent {
                    return invalid_representation(substance, "ionic solids must not be solvents");
                }
                validate_formula_units(substance, formula_units)
            }
            SubstanceRepresentation::Metal { element_symbol } => {
                if substance.charge != 0 {
                    return invalid_representation(substance, "metal materials must be neutral");
                }
                validate_non_empty_text(substance, element_symbol, "metal element symbol")?;
                if substance.phase_properties.solvent_role != SolventRole::NotSolvent {
                    return invalid_representation(
                        substance,
                        "metal materials must not be solvents",
                    );
                }
                Ok(())
            }
            SubstanceRepresentation::Oxide { formula_units } => {
                if substance.charge != 0 {
                    return invalid_representation(substance, "oxide materials must be neutral");
                }
                if substance.phase_properties.solvent_role != SolventRole::NotSolvent {
                    return invalid_representation(
                        substance,
                        "oxide materials must not be solvents",
                    );
                }
                validate_formula_units(substance, formula_units)
            }
            SubstanceRepresentation::Hydrate {
                formula_units,
                water_count,
            } => {
                if substance.charge != 0 {
                    return invalid_representation(substance, "hydrates must be neutral");
                }
                if *water_count == 0 {
                    return invalid_representation(substance, "hydrates must contain water");
                }
                validate_formula_units(substance, formula_units)
            }
            SubstanceRepresentation::SurfaceMaterial { active_site } => {
                if substance.charge != 0 {
                    return invalid_representation(substance, "surface materials must be neutral");
                }
                validate_optional_non_empty_text(
                    substance,
                    active_site.as_deref(),
                    "surface active site",
                )?;
                if substance.phase_properties.solvent_role != SolventRole::NotSolvent {
                    return invalid_representation(
                        substance,
                        "surface materials must not be solvents",
                    );
                }
                Ok(())
            }
            SubstanceRepresentation::UnspecifiedMaterial { reason } => {
                validate_non_empty_text(substance, reason, "unspecified material reason")
            }
        }?;
        if matches!(
            self,
            SubstanceRepresentation::IonicSolid { .. }
                | SubstanceRepresentation::Metal { .. }
                | SubstanceRepresentation::Oxide { .. }
                | SubstanceRepresentation::Hydrate { .. }
                | SubstanceRepresentation::SurfaceMaterial { .. }
        ) && substance.boiling_point_kelvin.is_finite()
            && substance.boiling_point_kelvin < substance.melting_point_kelvin
        {
            return invalid_representation(
                substance,
                "material boiling point must not be below its melting point",
            );
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum VaporPressureModel {
    ClausiusClapeyron {
        reference_temperature_kelvin: f64,
        reference_pressure_pascal: f64,
        enthalpy_j_per_mol: f64,
    },
    Log10PressurePascalAntoine {
        a: f64,
        b_kelvin: f64,
        c_kelvin: f64,
        min_temperature_kelvin: Option<f64>,
        max_temperature_kelvin: Option<f64>,
    },
}

impl VaporPressureModel {
    pub fn pressure_pascal(
        &self,
        substance_id: &SubstanceId,
        temperature_kelvin: f64,
    ) -> ChemistryResult<f64> {
        self.validate(substance_id)?;
        validate_positive_finite(substance_id, temperature_kelvin, "temperature")?;
        let pressure = match self {
            VaporPressureModel::ClausiusClapeyron {
                reference_temperature_kelvin,
                reference_pressure_pascal,
                enthalpy_j_per_mol,
            } => clausius_clapeyron_pressure_pascal(
                substance_id,
                *reference_temperature_kelvin,
                *reference_pressure_pascal,
                *enthalpy_j_per_mol,
                temperature_kelvin,
            )?,
            VaporPressureModel::Log10PressurePascalAntoine {
                a,
                b_kelvin,
                c_kelvin,
                min_temperature_kelvin,
                max_temperature_kelvin,
            } => {
                if min_temperature_kelvin.is_some_and(|min| temperature_kelvin < min)
                    || max_temperature_kelvin.is_some_and(|max| temperature_kelvin > max)
                {
                    return Err(ChemistryError::InvalidSubstance {
                        substance_id: substance_id.to_string(),
                        reason: format!(
                            "temperature {temperature_kelvin} K is outside Antoine model range"
                        ),
                    });
                }
                let denominator = temperature_kelvin + c_kelvin;
                if !denominator.is_finite() || denominator <= 0.0 {
                    return Err(ChemistryError::InvalidSubstance {
                        substance_id: substance_id.to_string(),
                        reason: "Antoine denominator must be positive at requested temperature"
                            .to_string(),
                    });
                }
                10.0_f64.powf(a - b_kelvin / denominator)
            }
        };
        validate_vapor_pressure_result(substance_id, pressure)?;
        Ok(pressure)
    }

    pub fn validate(&self, substance_id: &SubstanceId) -> ChemistryResult<()> {
        match self {
            VaporPressureModel::ClausiusClapeyron {
                reference_temperature_kelvin,
                reference_pressure_pascal,
                enthalpy_j_per_mol,
            } => {
                validate_positive_finite(
                    substance_id,
                    *reference_temperature_kelvin,
                    "vapor pressure reference temperature",
                )?;
                validate_positive_finite(
                    substance_id,
                    *reference_pressure_pascal,
                    "vapor pressure reference pressure",
                )?;
                validate_positive_finite(
                    substance_id,
                    *enthalpy_j_per_mol,
                    "vapor pressure enthalpy",
                )
            }
            VaporPressureModel::Log10PressurePascalAntoine {
                a,
                b_kelvin,
                c_kelvin,
                min_temperature_kelvin,
                max_temperature_kelvin,
            } => {
                if !a.is_finite() {
                    return Err(ChemistryError::InvalidSubstance {
                        substance_id: substance_id.to_string(),
                        reason: "Antoine coefficient a must be finite".to_string(),
                    });
                }
                validate_positive_finite(substance_id, *b_kelvin, "Antoine coefficient b")?;
                if !c_kelvin.is_finite() {
                    return Err(ChemistryError::InvalidSubstance {
                        substance_id: substance_id.to_string(),
                        reason: "Antoine coefficient c must be finite".to_string(),
                    });
                }
                validate_optional_positive_temperature(
                    substance_id,
                    *min_temperature_kelvin,
                    "minimum Antoine temperature",
                )?;
                validate_optional_positive_temperature(
                    substance_id,
                    *max_temperature_kelvin,
                    "maximum Antoine temperature",
                )?;
                if let (Some(min), Some(max)) = (min_temperature_kelvin, max_temperature_kelvin) {
                    if min > max {
                        return Err(ChemistryError::InvalidSubstance {
                            substance_id: substance_id.to_string(),
                            reason: "minimum Antoine temperature must not exceed maximum Antoine temperature"
                                .to_string(),
                        });
                    }
                    if min + c_kelvin <= 0.0 {
                        return Err(ChemistryError::InvalidSubstance {
                            substance_id: substance_id.to_string(),
                            reason: "Antoine denominator must remain positive in its valid range"
                                .to_string(),
                        });
                    }
                }
                Ok(())
            }
        }
    }
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
            representation: if charge == 0 {
                SubstanceRepresentation::Molecular
            } else {
                SubstanceRepresentation::Ion {
                    parent_element: None,
                }
            },
            charge,
            molar_mass_grams,
            liquid_density_grams_per_bucket,
            solid_density_grams_per_bucket: liquid_density_grams_per_bucket,
            melting_point_kelvin: 0.0,
            boiling_point_kelvin,
            molar_heat_capacity_j_per_mol_kelvin,
            fusion_heat_j_per_mol: 0.0,
            latent_heat_j_per_mol,
            critical_temperature_kelvin: None,
            critical_pressure_pascal: None,
            acentric_factor: None,
            vapor_pressure_model: None,
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

    pub fn with_representation(mut self, representation: SubstanceRepresentation) -> Self {
        self.representation = representation;
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

    pub fn with_critical_point(
        mut self,
        critical_temperature_kelvin: f64,
        critical_pressure_pascal: f64,
    ) -> Self {
        self.critical_temperature_kelvin = Some(critical_temperature_kelvin);
        self.critical_pressure_pascal = Some(critical_pressure_pascal);
        self
    }

    pub fn with_acentric_factor(mut self, acentric_factor: f64) -> Self {
        self.acentric_factor = Some(acentric_factor);
        self
    }

    pub fn with_vapor_pressure_model(mut self, vapor_pressure_model: VaporPressureModel) -> Self {
        self.vapor_pressure_model = Some(vapor_pressure_model);
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

    pub fn vapor_pressure_pascal(&self, temperature_kelvin: f64) -> ChemistryResult<Option<f64>> {
        validate_positive_finite(&self.id, temperature_kelvin, "temperature")?;
        if self.charge != 0 {
            return Ok(None);
        }
        if let Some(critical_temperature) = self.critical_temperature_kelvin {
            if temperature_kelvin >= critical_temperature {
                return Ok(None);
            }
        }
        let pressure = if let Some(model) = &self.vapor_pressure_model {
            model.pressure_pascal(&self.id, temperature_kelvin)?
        } else {
            if !self.boiling_point_kelvin.is_finite()
                || self.boiling_point_kelvin <= 0.0
                || !self.latent_heat_j_per_mol.is_finite()
                || self.latent_heat_j_per_mol <= 0.0
            {
                return Ok(None);
            }
            clausius_clapeyron_pressure_pascal(
                &self.id,
                self.boiling_point_kelvin,
                NORMAL_BOILING_PRESSURE_PASCAL,
                self.latent_heat_j_per_mol,
                temperature_kelvin,
            )?
        };
        validate_vapor_pressure_result(&self.id, pressure)?;
        Ok(Some(pressure))
    }

    pub fn with_molecular_structure(mut self, molecular_structure: MolecularStructure) -> Self {
        self.functional_groups = find_functional_groups(&molecular_structure);
        self.molecular_structure = Some(molecular_structure);
        self
    }

    pub fn validate(&self) -> ChemistryResult<()> {
        let id = self.id.to_string();
        if self.id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidSubstance {
                substance_id: id,
                reason: "id must not be empty".to_string(),
            });
        }
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
        match (
            self.critical_temperature_kelvin,
            self.critical_pressure_pascal,
        ) {
            (Some(temperature), Some(pressure)) => {
                validate_positive_finite(&self.id, temperature, "critical temperature")?;
                validate_positive_finite(&self.id, pressure, "critical pressure")?;
                if temperature <= self.boiling_point_kelvin && self.boiling_point_kelvin.is_finite()
                {
                    return Err(ChemistryError::InvalidSubstance {
                        substance_id: self.id.to_string(),
                        reason: "critical temperature must be above normal boiling point"
                            .to_string(),
                    });
                }
            }
            (None, None) => {}
            _ => {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: "critical temperature and pressure must be specified together"
                        .to_string(),
                });
            }
        }
        if let Some(acentric_factor) = self.acentric_factor {
            if !acentric_factor.is_finite() {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: "acentric factor must be finite".to_string(),
                });
            }
            if self.critical_temperature_kelvin.is_none() {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: self.id.to_string(),
                    reason: "acentric factor requires a critical point".to_string(),
                });
            }
        }
        if let Some(model) = &self.vapor_pressure_model {
            model.validate(&self.id)?;
        }
        self.phase_properties.validate(&self.id, self.charge)?;
        self.representation.validate(self)?;
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

fn validate_formula_units(
    substance: &Substance,
    formula_units: &[MaterialFormulaUnit],
) -> ChemistryResult<()> {
    if formula_units.is_empty() {
        return invalid_representation(substance, "material formula must not be empty");
    }
    for unit in formula_units {
        if unit.substance_id.as_str().trim().is_empty() {
            return invalid_representation(
                substance,
                "material formula component id must not be empty",
            );
        }
        if unit.coefficient == 0 {
            return invalid_representation(
                substance,
                "material formula coefficients must be greater than zero",
            );
        }
    }
    Ok(())
}

fn validate_non_empty_text(substance: &Substance, value: &str, name: &str) -> ChemistryResult<()> {
    if value.trim().is_empty() {
        return invalid_representation(substance, &format!("{name} must not be empty"));
    }
    Ok(())
}

fn validate_optional_non_empty_text(
    substance: &Substance,
    value: Option<&str>,
    name: &str,
) -> ChemistryResult<()> {
    if let Some(value) = value {
        validate_non_empty_text(substance, value, name)?;
    }
    Ok(())
}

fn invalid_representation<T>(substance: &Substance, reason: &str) -> ChemistryResult<T> {
    Err(ChemistryError::InvalidSubstance {
        substance_id: substance.id.to_string(),
        reason: reason.to_string(),
    })
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

fn validate_positive_finite(
    substance_id: &SubstanceId,
    value: f64,
    name: &str,
) -> ChemistryResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: substance_id.to_string(),
            reason: format!("{name} must be positive and finite"),
        });
    }
    Ok(())
}

fn validate_optional_positive_temperature(
    substance_id: &SubstanceId,
    value: Option<f64>,
    name: &str,
) -> ChemistryResult<()> {
    if let Some(value) = value {
        validate_positive_finite(substance_id, value, name)?;
    }
    Ok(())
}

fn clausius_clapeyron_pressure_pascal(
    substance_id: &SubstanceId,
    reference_temperature_kelvin: f64,
    reference_pressure_pascal: f64,
    enthalpy_j_per_mol: f64,
    temperature_kelvin: f64,
) -> ChemistryResult<f64> {
    let exponent = -enthalpy_j_per_mol / GAS_CONSTANT_J_PER_MOL_KELVIN
        * (1.0 / temperature_kelvin - 1.0 / reference_temperature_kelvin);
    let pressure = reference_pressure_pascal * exponent.exp();
    if !pressure.is_finite() || pressure < 0.0 {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: substance_id.to_string(),
            reason: "vapor pressure calculation produced an invalid value".to_string(),
        });
    }
    Ok(pressure)
}

fn validate_vapor_pressure_result(
    substance_id: &SubstanceId,
    pressure_pascal: f64,
) -> ChemistryResult<()> {
    if !pressure_pascal.is_finite() || pressure_pascal < 0.0 {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: substance_id.to_string(),
            reason: "vapor pressure must be non-negative and finite".to_string(),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn valid_substance() -> Substance {
        Substance::new("destroy:test", 0, 44.0, 44_000.0, 250.0, 60.0, 10_000.0)
    }

    #[test]
    fn gas_liquid_critical_properties_validate_as_a_pair() {
        let valid = valid_substance().with_critical_point(304.1, 7_377_000.0);
        assert!(valid.validate().is_ok());

        let mut incomplete = valid_substance();
        incomplete.critical_temperature_kelvin = Some(304.1);
        let error = incomplete.validate().unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }

    #[test]
    fn critical_temperature_must_be_above_normal_boiling_point() {
        let substance = valid_substance().with_critical_point(240.0, 7_377_000.0);
        let error = substance.validate().unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }

    #[test]
    fn acentric_factor_requires_critical_point() {
        let substance = valid_substance().with_acentric_factor(0.225);
        let error = substance.validate().unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }

    #[test]
    fn ionic_solid_representation_requires_neutral_precipitating_material() {
        let valid = valid_substance()
            .with_phase_properties(SubstancePhaseProperties {
                preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                aqueous_solubility_mol_per_bucket: Some(0.01),
                organic_solubility_mol_per_bucket: Some(0.0),
                can_precipitate: true,
                can_form_liquid_phase: false,
                solvent_role: SolventRole::NotSolvent,
            })
            .with_representation(SubstanceRepresentation::IonicSolid {
                formula_units: vec![
                    MaterialFormulaUnit::new("destroy:sodium_ion", 1),
                    MaterialFormulaUnit::new("destroy:chloride", 1),
                ],
            });
        assert!(valid.validate().is_ok());

        let charged = Substance::new(
            "destroy:bad_salt",
            1,
            58.5,
            2_160_000.0,
            f64::MAX,
            50.0,
            0.0,
        )
        .with_phase_properties(SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: Some(0.01),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: false,
            solvent_role: SolventRole::NotSolvent,
        })
        .with_representation(SubstanceRepresentation::IonicSolid {
            formula_units: vec![MaterialFormulaUnit::new("destroy:sodium_ion", 1)],
        });
        assert!(matches!(
            charged.validate().unwrap_err(),
            ChemistryError::InvalidSubstance { .. }
        ));
    }

    #[test]
    fn material_representation_rejects_solvent_metal() {
        let metal = valid_substance()
            .with_solvent_role(SolventRole::KnownSolvent)
            .with_representation(SubstanceRepresentation::Metal {
                element_symbol: "Fe".to_string(),
            });

        assert!(matches!(
            metal.validate().unwrap_err(),
            ChemistryError::InvalidSubstance { .. }
        ));
    }

    #[test]
    fn vapor_pressure_model_rejects_invalid_numbers() {
        let substance =
            valid_substance().with_vapor_pressure_model(VaporPressureModel::ClausiusClapeyron {
                reference_temperature_kelvin: 250.0,
                reference_pressure_pascal: 101_325.0,
                enthalpy_j_per_mol: f64::NAN,
            });
        let error = substance.validate().unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }

    #[test]
    fn antoine_range_must_be_ordered_and_non_singular() {
        let reversed = valid_substance().with_vapor_pressure_model(
            VaporPressureModel::Log10PressurePascalAntoine {
                a: 8.0,
                b_kelvin: 1_500.0,
                c_kelvin: -200.0,
                min_temperature_kelvin: Some(350.0),
                max_temperature_kelvin: Some(300.0),
            },
        );
        assert!(matches!(
            reversed.validate().unwrap_err(),
            ChemistryError::InvalidSubstance { .. }
        ));

        let singular = valid_substance().with_vapor_pressure_model(
            VaporPressureModel::Log10PressurePascalAntoine {
                a: 8.0,
                b_kelvin: 1_500.0,
                c_kelvin: -400.0,
                min_temperature_kelvin: Some(300.0),
                max_temperature_kelvin: Some(350.0),
            },
        );
        assert!(matches!(
            singular.validate().unwrap_err(),
            ChemistryError::InvalidSubstance { .. }
        ));
    }

    #[test]
    fn implicit_vapor_pressure_uses_normal_boiling_point() {
        let substance = valid_substance();

        let at_boiling = substance
            .vapor_pressure_pascal(substance.boiling_point_kelvin)
            .unwrap()
            .unwrap();
        let below_boiling = substance
            .vapor_pressure_pascal(substance.boiling_point_kelvin - 25.0)
            .unwrap()
            .unwrap();
        let above_boiling = substance
            .vapor_pressure_pascal(substance.boiling_point_kelvin + 25.0)
            .unwrap()
            .unwrap();

        assert!((at_boiling / NORMAL_BOILING_PRESSURE_PASCAL - 1.0).abs() < 1.0e-12);
        assert!(below_boiling < NORMAL_BOILING_PRESSURE_PASCAL);
        assert!(above_boiling > NORMAL_BOILING_PRESSURE_PASCAL);
    }

    #[test]
    fn explicit_clausius_vapor_pressure_uses_reference_point() {
        let substance =
            valid_substance().with_vapor_pressure_model(VaporPressureModel::ClausiusClapeyron {
                reference_temperature_kelvin: 300.0,
                reference_pressure_pascal: 50_000.0,
                enthalpy_j_per_mol: 20_000.0,
            });

        let reference = substance.vapor_pressure_pascal(300.0).unwrap().unwrap();
        let warmer = substance.vapor_pressure_pascal(320.0).unwrap().unwrap();

        assert!((reference - 50_000.0).abs() < 1.0e-9);
        assert!(warmer > reference);
    }

    #[test]
    fn antoine_vapor_pressure_respects_temperature_range() {
        let substance = valid_substance().with_vapor_pressure_model(
            VaporPressureModel::Log10PressurePascalAntoine {
                a: 5.0,
                b_kelvin: 100.0,
                c_kelvin: 0.0,
                min_temperature_kelvin: Some(250.0),
                max_temperature_kelvin: Some(350.0),
            },
        );

        let pressure = substance.vapor_pressure_pascal(300.0).unwrap().unwrap();
        let error = substance.vapor_pressure_pascal(200.0).unwrap_err();

        assert!((pressure - 10.0_f64.powf(5.0 - 100.0 / 300.0)).abs() < 1.0e-9);
        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }

    #[test]
    fn vapor_pressure_is_absent_for_ions_and_supercritical_temperature() {
        let ion = Substance::new("destroy:ion", 1, 20.0, 20_000.0, f64::MAX, 50.0, 20_000.0);
        let supercritical = valid_substance().with_critical_point(300.0, 7_000_000.0);

        assert_eq!(ion.vapor_pressure_pascal(298.0).unwrap(), None);
        assert_eq!(supercritical.vapor_pressure_pascal(300.0).unwrap(), None);
    }

    #[test]
    fn vapor_pressure_rejects_non_positive_temperature() {
        let substance = valid_substance();
        let error = substance.vapor_pressure_pascal(0.0).unwrap_err();

        assert!(matches!(error, ChemistryError::InvalidSubstance { .. }));
    }
}
