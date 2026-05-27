//! Central selectivity engine for evaluating reaction favorability
//!
//! Coordinates all selectivity rules and provides unified API for generators.

use crate::chemistry::reactive_site::ReactiveSiteKind;
use super::types::*;
use super::nucleophilic_substitution::evaluate_sn2;
use super::elimination::{evaluate_e2, evaluate_e1};
use super::esterification::evaluate_fischer_esterification;
use super::carbonyl_addition::{evaluate_carbonyl_addition, NucleophileStrength};

/// Central engine for selectivity evaluation
pub struct SelectivityEngine;

impl SelectivityEngine {
    /// Evaluate SN2 substitution with E2 competition
    pub fn sn2_with_competition(
        site: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> SelectivityResult {
        evaluate_sn2(site, context)
    }
    
    /// Evaluate E2 elimination (direct)
    pub fn e2_elimination(
        site: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> ReactivityScore {
        evaluate_e2(site, context)
    }
    
    /// Evaluate E1 elimination (usually competing with SN1)
    pub fn e1_elimination(
        site: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> ReactivityScore {
        evaluate_e1(site, context)
    }
    
    /// Evaluate Fischer esterification with elimination competition
    pub fn fischer_esterification(
        acid: &SiteDescriptor,
        alcohol: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> SelectivityResult {
        evaluate_fischer_esterification(acid, alcohol, context)
    }
    
    /// Evaluate carbonyl addition (generic)
    pub fn carbonyl_addition(
        carbonyl: &SiteDescriptor,
        nucleophile_strength: NucleophileStrength,
        context: &SelectivityContext,
    ) -> ReactivityScore {
        evaluate_carbonyl_addition(carbonyl, nucleophile_strength, context)
    }
    
    /// Evaluate carbonyl addition with NaBH4 (moderate nucleophile)
    pub fn borohydride_reduction(
        carbonyl: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> ReactivityScore {
        evaluate_carbonyl_addition(carbonyl, NucleophileStrength::Strong, context)
    }
    
    /// Evaluate carbonyl addition with Grignard/organolithium (very strong nucleophile)
    pub fn organometallic_addition(
        carbonyl: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> ReactivityScore {
        evaluate_carbonyl_addition(carbonyl, NucleophileStrength::VeryStrong, context)
    }
    
    /// Quick check: will SN2 occur (vs elimination)?
    pub fn sn2_favored(site: &SiteDescriptor, context: &SelectivityContext) -> bool {
        let result = evaluate_sn2(site, context);
        matches!(
            result.recommendation,
            SelectivityRecommendation::Exclusive | SelectivityRecommendation::Preferred
        )
    }
    
    /// Quick check: will esterification occur (vs elimination)?
    pub fn esterification_favored(
        acid: &SiteDescriptor,
        alcohol: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> bool {
        let result = evaluate_fischer_esterification(acid, alcohol, context);
        !matches!(
            result.recommendation,
            SelectivityRecommendation::Suppressed | SelectivityRecommendation::None
        )
    }
    
    /// Quick check: will carbonyl addition occur faster to aldehyde than ketone?
    pub fn aldehyde_preferred_over_ketone(
        aldehyde: &SiteDescriptor,
        ketone: &SiteDescriptor,
        nucleophile: NucleophileStrength,
        context: &SelectivityContext,
    ) -> bool {
        let ald_score = evaluate_carbonyl_addition(aldehyde, nucleophile, context);
        let ket_score = evaluate_carbonyl_addition(ketone, nucleophile, context);
        ald_score.value > ket_score.value * 2.0  // At least 2:1 selectivity
    }
    
    /// Get activation energy modification from selectivity
    /// 
    /// Returns delta to add to base activation energy (negative = faster)
    pub fn activation_energy_delta(score: &ReactivityScore) -> f64 {
        score.activation_delta
    }
    
    /// Get pre-exponential factor modification from selectivity
    pub fn pre_exponential_multiplier(score: &ReactivityScore) -> f64 {
        score.pre_exp_multiplier
    }
}

/// Builder for creating site descriptors from molecular analysis
pub struct SiteDescriptorBuilder;

impl SiteDescriptorBuilder {
    /// Build descriptor from basic parameters
    /// 
    /// This is a simplified version that can be expanded later with
    /// full molecular structure analysis
    pub fn build(
        site_kind: ReactiveSiteKind,
        degree: SubstitutionDegree,
        edg_count: u32,
        ewg_count: u32,
        bulky_substituents: u32,
        has_beta_hydrogen: bool,
        resonance: bool,
        aromatic: bool,
    ) -> SiteDescriptor {
        let electronics = ElectronicEnvironment {
            electron_donating_groups: edg_count,
            electron_withdrawing_groups: ewg_count,
            resonance_stabilization: resonance,
            aromatic,
        };
        
        SiteDescriptor::new(
            site_kind,
            degree,
            electronics,
            bulky_substituents,
            has_beta_hydrogen,
        )
    }
    
    /// Create simple primary site (no substituents)
    pub fn primary_alcohol() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Primary,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create simple secondary site
    pub fn secondary_alcohol() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Secondary,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create simple tertiary site
    pub fn tertiary_alcohol() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Tertiary,
            0, 0, 0, false, false, false,
        )
    }
    
    /// Create benzylic site
    pub fn benzylic_alcohol() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Benzylic,
            0, 0, 0, true, true, false,
        )
    }
    
    /// Create simple aldehyde
    pub fn aldehyde() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Aldehyde,
            SubstitutionDegree::Primary,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create simple ketone
    pub fn ketone() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Ketone,
            SubstitutionDegree::Secondary,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create aromatic aldehyde (benzaldehyde-like)
    pub fn aromatic_aldehyde() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Aldehyde,
            SubstitutionDegree::Primary,
            0, 0, 0, false, true, true,
        )
    }
    
    /// Create primary halide
    pub fn primary_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Primary,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create secondary halide
    pub fn secondary_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Secondary,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create tertiary halide
    pub fn tertiary_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Tertiary,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create benzylic halide
    pub fn benzylic_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Benzylic,
            0, 0, 0, true, true, false,
        )
    }
    
    /// Create carboxylic acid
    pub fn carboxylic_acid() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::CarboxylicAcid,
            SubstitutionDegree::Primary,
            0, 0, 0, false, false, false,
        )
    }
    
    // Integration methods for creating descriptors from actual site data
    
    /// Create descriptor from alcohol site
    pub(crate) fn from_alcohol_site(site: &crate::chemistry::organic::centers::AlcoholSite) -> SiteDescriptor {
        let degree = if site.degree >= 3 {
            SubstitutionDegree::Tertiary
        } else if site.degree == 2 {
            SubstitutionDegree::Secondary
        } else {
            SubstitutionDegree::Primary
        };
        
        Self::build(
            ReactiveSiteKind::Alcohol,
            degree,
            0,  // EDG count - to be determined from structure analysis
            0,  // EWG count
            0,  // bulky substituents
            true,  // has beta hydrogen (simplified, should check structure)
            false,  // resonance
            false,  // aromatic
        )
    }
    
    /// Create descriptor from carbonyl site
    pub(crate) fn from_carbonyl_site(site: &crate::chemistry::organic::centers::CarbonylSite) -> SiteDescriptor {
        let degree = if site.is_ketone {
            SubstitutionDegree::Secondary
        } else {
            SubstitutionDegree::Primary
        };
        
        Self::build(
            if site.is_ketone { ReactiveSiteKind::Ketone } else { ReactiveSiteKind::Aldehyde },
            degree,
            0, 0, 0, true, false, false,
        )
    }
    
    /// Create descriptor from carboxylic acid site
    pub(crate) fn from_carboxylic_acid_site(_site: &crate::chemistry::organic::centers::CarboxylicAcidSite) -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::CarboxylicAcid,
            SubstitutionDegree::Primary,  // R-COOH, carbon is attached to R and =O
            0, 0, 0, false, false, false,
        )
    }
    
    /// Create descriptor from halide site
    pub(crate) fn from_halide_site(site: &crate::chemistry::organic::centers::HalideSite) -> SiteDescriptor {
        let degree = if site.degree >= 3 {
            SubstitutionDegree::Tertiary
        } else if site.degree == 2 {
            SubstitutionDegree::Secondary
        } else {
            SubstitutionDegree::Primary
        };
        
        Self::build(
            ReactiveSiteKind::Halide,
            degree,
            0, 0, 0, true, false, false,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn selectivity_engine_sn2_evaluation() {
        let site = SiteDescriptorBuilder::primary_halide();
        let context = SelectivityContext::default();
        let result = SelectivityEngine::sn2_with_competition(&site, &context);
        
        assert!(result.primary.value > 0.5);
        assert!(matches!(result.recommendation, 
            SelectivityRecommendation::Exclusive | SelectivityRecommendation::Preferred));
    }
    
    #[test]
    fn selectivity_engine_esterification_suppressed_for_tertiary() {
        let acid = SiteDescriptorBuilder::carboxylic_acid();
        let alcohol = SiteDescriptorBuilder::tertiary_alcohol();
        
        let context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 380.0,
            ph: Some(2.0),
        };
        
        let result = SelectivityEngine::fischer_esterification(&acid, &alcohol, &context);
        
        assert!(matches!(result.recommendation,
            SelectivityRecommendation::Suppressed | SelectivityRecommendation::None | 
            SelectivityRecommendation::Mixed));
    }
    
    #[test]
    fn aldehyde_vs_ketone_selectivity() {
        let aldehyde = SiteDescriptorBuilder::aldehyde();
        let ketone = SiteDescriptorBuilder::ketone();
        
        let context = SelectivityContext::default();
        
        assert!(SelectivityEngine::aldehyde_preferred_over_ketone(
            &aldehyde, &ketone, NucleophileStrength::Moderate, &context
        ));
    }
    
    #[test]
    fn site_descriptor_builder_convenience() {
        let primary = SiteDescriptorBuilder::primary_alcohol();
        assert!(matches!(primary.degree, SubstitutionDegree::Primary));
        assert!(primary.has_beta_hydrogen);
        
        let tertiary = SiteDescriptorBuilder::tertiary_alcohol();
        assert!(matches!(tertiary.degree, SubstitutionDegree::Tertiary));
        
        let benzaldehyde = SiteDescriptorBuilder::aromatic_aldehyde();
        assert!(benzaldehyde.electronics.aromatic);
    }
    
    #[test]
    fn builder_factory_methods() {
        // All factory methods should create valid descriptors
        let sites = vec![
            SiteDescriptorBuilder::primary_alcohol(),
            SiteDescriptorBuilder::secondary_alcohol(),
            SiteDescriptorBuilder::tertiary_alcohol(),
            SiteDescriptorBuilder::benzylic_alcohol(),
            SiteDescriptorBuilder::aldehyde(),
            SiteDescriptorBuilder::ketone(),
            SiteDescriptorBuilder::aromatic_aldehyde(),
            SiteDescriptorBuilder::primary_halide(),
            SiteDescriptorBuilder::secondary_halide(),
            SiteDescriptorBuilder::tertiary_halide(),
            SiteDescriptorBuilder::benzylic_halide(),
            SiteDescriptorBuilder::carboxylic_acid(),
        ];
        
        for site in &sites {
            assert!(site.steric_score >= 0.0 && site.steric_score <= 1.0);
            assert!(site.steric_accessibility() > 0.0);
        }
    }
}


