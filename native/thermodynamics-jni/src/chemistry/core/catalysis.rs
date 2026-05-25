use std::collections::BTreeMap;
use std::fmt::{Display, Formatter};

use super::error::{ChemistryError, ChemistryResult};
use super::mixture::MixturePhase;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CatalystSurfaceId(String);

impl CatalystSurfaceId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<surface>".to_string(),
                reason: "surface id must not be empty".to_string(),
            });
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for CatalystSurfaceId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for CatalystSurfaceId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

impl From<String> for CatalystSurfaceId {
    fn from(value: String) -> Self {
        Self(value)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SurfaceSiteId(String);

impl SurfaceSiteId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<surface-site>".to_string(),
                reason: "surface site id must not be empty".to_string(),
            });
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for SurfaceSiteId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for SurfaceSiteId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

impl From<String> for SurfaceSiteId {
    fn from(value: String) -> Self {
        Self(value)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CatalystSurfaceSpec {
    pub id: CatalystSurfaceId,
    pub site_id: SurfaceSiteId,
    pub accessible_phases: Vec<MixturePhase>,
    pub molar_mass_grams: Option<f64>,
    pub charge: Option<i32>,
    pub unchecked_mass_reason: Option<String>,
}

impl CatalystSurfaceSpec {
    pub fn unchecked(
        id: impl Into<CatalystSurfaceId>,
        unchecked_mass_reason: impl Into<String>,
    ) -> Self {
        Self {
            id: id.into(),
            site_id: "default".into(),
            accessible_phases: vec![
                MixturePhase::Aqueous,
                MixturePhase::Organic,
                MixturePhase::Gas,
            ],
            molar_mass_grams: None,
            charge: None,
            unchecked_mass_reason: Some(unchecked_mass_reason.into()),
        }
    }

    pub fn chemical(id: impl Into<CatalystSurfaceId>, molar_mass_grams: f64, charge: i32) -> Self {
        Self {
            id: id.into(),
            site_id: "default".into(),
            accessible_phases: vec![
                MixturePhase::Aqueous,
                MixturePhase::Organic,
                MixturePhase::Gas,
            ],
            molar_mass_grams: Some(molar_mass_grams),
            charge: Some(charge),
            unchecked_mass_reason: None,
        }
    }

    pub fn with_site(mut self, site_id: impl Into<SurfaceSiteId>) -> Self {
        self.site_id = site_id.into();
        self
    }

    pub fn with_accessible_phases(
        mut self,
        phases: impl IntoIterator<Item = MixturePhase>,
    ) -> Self {
        self.accessible_phases = phases.into_iter().collect();
        self
    }

    pub fn validate(&self) -> ChemistryResult<()> {
        if self.id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<surface>".to_string(),
                reason: "surface id must not be empty".to_string(),
            });
        }
        if self.site_id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: self.id.to_string(),
                reason: "surface site id must not be empty".to_string(),
            });
        }
        if self.accessible_phases.is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: self.id.to_string(),
                reason: "surface must expose at least one phase".to_string(),
            });
        }
        match (self.molar_mass_grams, self.charge) {
            (Some(mass), Some(_)) if mass.is_finite() && mass >= 0.0 => Ok(()),
            (None, None) if self.unchecked_mass_reason.is_some() => Ok(()),
            (Some(_), Some(_)) => Err(ChemistryError::InvalidReaction {
                reaction_id: self.id.to_string(),
                reason: "surface mass must be non-negative and finite".to_string(),
            }),
            _ => Err(ChemistryError::InvalidReaction {
                reaction_id: self.id.to_string(),
                reason:
                    "surface must either provide both mass and charge or an unchecked mass reason"
                        .to_string(),
            }),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceRequirement {
    pub surface_id: CatalystSurfaceId,
    pub site_id: SurfaceSiteId,
    pub sites_per_reaction: f64,
    pub phases: Vec<MixturePhase>,
}

impl SurfaceRequirement {
    pub fn new(surface_id: impl Into<CatalystSurfaceId>, sites_per_reaction: f64) -> Self {
        Self {
            surface_id: surface_id.into(),
            site_id: "default".into(),
            sites_per_reaction,
            phases: vec![
                MixturePhase::Aqueous,
                MixturePhase::Organic,
                MixturePhase::Gas,
            ],
        }
    }

    pub fn with_site(mut self, site_id: impl Into<SurfaceSiteId>) -> Self {
        self.site_id = site_id.into();
        self
    }

    pub fn with_phases(mut self, phases: impl IntoIterator<Item = MixturePhase>) -> Self {
        self.phases = phases.into_iter().collect();
        self
    }

    pub fn validate(&self, reaction_id: &str) -> ChemistryResult<()> {
        if self.surface_id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "surface requirement id must not be empty".to_string(),
            });
        }
        if self.site_id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "surface requirement site id must not be empty".to_string(),
            });
        }
        if !self.sites_per_reaction.is_finite() || self.sites_per_reaction <= 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "surface requirement sites must be positive and finite".to_string(),
            });
        }
        if self.phases.is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "surface requirement must expose at least one phase".to_string(),
            });
        }
        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceStep {
    Adsorb {
        surface_id: CatalystSurfaceId,
        site_id: SurfaceSiteId,
        sites_per_reaction: f64,
    },
    Desorb {
        surface_id: CatalystSurfaceId,
        site_id: SurfaceSiteId,
        sites_per_reaction: f64,
    },
    Poison {
        surface_id: CatalystSurfaceId,
        site_id: SurfaceSiteId,
        sites_per_reaction: f64,
    },
    Restore {
        surface_id: CatalystSurfaceId,
        sites_per_reaction: f64,
    },
}

impl SurfaceStep {
    pub fn validate(&self, reaction_id: &str) -> ChemistryResult<()> {
        let (surface_id, site_id, sites) = match self {
            SurfaceStep::Adsorb {
                surface_id,
                site_id,
                sites_per_reaction,
            }
            | SurfaceStep::Desorb {
                surface_id,
                site_id,
                sites_per_reaction,
            }
            | SurfaceStep::Poison {
                surface_id,
                site_id,
                sites_per_reaction,
            } => (surface_id, Some(site_id), *sites_per_reaction),
            SurfaceStep::Restore {
                surface_id,
                sites_per_reaction,
            } => (surface_id, None, *sites_per_reaction),
        };
        if surface_id.as_str().trim().is_empty() {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "surface step id must not be empty".to_string(),
            });
        }
        if site_id
            .map(|site| site.as_str().trim().is_empty())
            .unwrap_or(false)
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "surface step site id must not be empty".to_string(),
            });
        }
        if !sites.is_finite() || sites <= 0.0 {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: reaction_id.to_string(),
                reason: "surface step sites must be positive and finite".to_string(),
            });
        }
        Ok(())
    }

    pub fn surface_id(&self) -> &CatalystSurfaceId {
        match self {
            SurfaceStep::Adsorb { surface_id, .. }
            | SurfaceStep::Desorb { surface_id, .. }
            | SurfaceStep::Poison { surface_id, .. }
            | SurfaceStep::Restore { surface_id, .. } => surface_id,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CatalystSurfaceState {
    pub total_sites_mol_per_bucket: f64,
    pub occupied_sites: BTreeMap<SurfaceSiteId, f64>,
    pub poisoned_sites_mol_per_bucket: f64,
    pub temperature_kelvin: f64,
    pub accessible_phases: Vec<MixturePhase>,
}

impl CatalystSurfaceState {
    pub fn new(total_sites_mol_per_bucket: f64) -> ChemistryResult<Self> {
        if !total_sites_mol_per_bucket.is_finite() || total_sites_mol_per_bucket < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "surface site amount must be non-negative and finite".to_string(),
            ));
        }
        Ok(Self {
            total_sites_mol_per_bucket,
            occupied_sites: BTreeMap::new(),
            poisoned_sites_mol_per_bucket: 0.0,
            temperature_kelvin: 298.0,
            accessible_phases: vec![
                MixturePhase::Aqueous,
                MixturePhase::Organic,
                MixturePhase::Gas,
            ],
        })
    }

    pub fn with_accessible_phases(
        mut self,
        phases: impl IntoIterator<Item = MixturePhase>,
    ) -> ChemistryResult<Self> {
        let phases = phases.into_iter().collect::<Vec<_>>();
        if phases.is_empty() {
            return Err(ChemistryError::InvalidMixtureState(
                "surface must expose at least one phase".to_string(),
            ));
        }
        self.accessible_phases = phases;
        Ok(self)
    }

    pub fn free_sites(&self) -> f64 {
        let occupied = self.occupied_sites.values().sum::<f64>();
        (self.total_sites_mol_per_bucket - occupied - self.poisoned_sites_mol_per_bucket).max(0.0)
    }

    pub fn occupied_sites(&self) -> f64 {
        self.occupied_sites.values().sum()
    }

    pub fn occupied_site(&self, site_id: &SurfaceSiteId) -> f64 {
        self.occupied_sites.get(site_id).copied().unwrap_or(0.0)
    }

    pub fn poisoned_sites(&self) -> f64 {
        self.poisoned_sites_mol_per_bucket
    }

    pub fn validate(&self) -> ChemistryResult<()> {
        if !self.total_sites_mol_per_bucket.is_finite() || self.total_sites_mol_per_bucket < 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "surface site amount must be non-negative and finite".to_string(),
            ));
        }
        if !self.poisoned_sites_mol_per_bucket.is_finite()
            || self.poisoned_sites_mol_per_bucket < 0.0
        {
            return Err(ChemistryError::InvalidMixtureState(
                "poisoned surface site amount must be non-negative and finite".to_string(),
            ));
        }
        if !self.temperature_kelvin.is_finite() || self.temperature_kelvin <= 0.0 {
            return Err(ChemistryError::InvalidMixtureState(
                "surface temperature must be positive and finite".to_string(),
            ));
        }
        if self.accessible_phases.is_empty() {
            return Err(ChemistryError::InvalidMixtureState(
                "surface must expose at least one phase".to_string(),
            ));
        }
        for amount in self.occupied_sites.values() {
            if !amount.is_finite() || *amount < 0.0 {
                return Err(ChemistryError::InvalidMixtureState(
                    "occupied surface site amount must be non-negative and finite".to_string(),
                ));
            }
        }
        if self.occupied_sites() + self.poisoned_sites_mol_per_bucket
            > self.total_sites_mol_per_bucket + super::mixture::TRACE_CONCENTRATION_MOL_PER_BUCKET
        {
            return Err(ChemistryError::InvalidMixtureState(
                "surface cannot have more occupied and poisoned sites than total sites".to_string(),
            ));
        }
        Ok(())
    }
}
