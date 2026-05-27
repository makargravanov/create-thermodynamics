//! Core types for kinetic selectivity system

use std::collections::HashMap;

/// Substitution degree at a reactive center
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubstitutionDegree {
    /// Primary: attached to ≤1 carbon (e.g., CH3OH, R-CH2-OH, R-CHO)
    Primary,
    /// Secondary: attached to 2 carbons (e.g., R2CH-OH, R2C=O ketone)
    Secondary,
    /// Tertiary: attached to 3 carbons (e.g., R3C-OH)
    Tertiary,
    /// Benzylic: attached to benzene ring (Ph-CH2-X)
    Benzylic,
    /// Allylic: attached to alkene carbon (C=C-CH2-X)
    Allylic,
}

impl SubstitutionDegree {
    /// Base steric score from degree alone (0.0 = unhindered, 1.0 = blocked)
    pub fn base_steric_score(&self) -> f64 {
        match self {
            SubstitutionDegree::Primary => 0.0,
            SubstitutionDegree::Secondary => 0.3,
            SubstitutionDegree::Tertiary => 0.6,
            SubstitutionDegree::Benzylic => 0.2,  // phenyl is flat, less hindrance
            SubstitutionDegree::Allylic => 0.25,
        }
    }
}

/// Electronic environment around a reactive site
#[derive(Debug, Clone, Default)]
pub struct ElectronicEnvironment {
    /// Count of electron-donating groups (+I, +M effects) in conjugation/position
    pub electron_donating_groups: u32,
    /// Count of electron-withdrawing groups (-I, -M effects)
    pub electron_withdrawing_groups: u32,
    /// Whether site is resonance-stabilized (conjugated)
    pub resonance_stabilization: bool,
    /// Whether site is part of aromatic system
    pub aromatic: bool,
}

impl ElectronicEnvironment {
    /// Net electronic effect: positive = activating, negative = deactivating
    pub fn net_effect(&self) -> f64 {
        let edg_factor = self.electron_donating_groups as f64 * 0.1;
        let ewg_factor = self.electron_withdrawing_groups as f64 * 0.1;
        let resonance_bonus = if self.resonance_stabilization { 0.15 } else { 0.0 };
        edg_factor - ewg_factor + resonance_bonus
    }
}

/// Complete descriptor for selectivity evaluation
#[derive(Debug, Clone)]
pub struct SiteDescriptor {
    /// Type of reactive site
    pub site_kind: crate::chemistry::reactive_site::ReactiveSiteKind,
    /// Substitution degree at the reactive center
    pub degree: SubstitutionDegree,
    /// Electronic environment
    pub electronics: ElectronicEnvironment,
    /// Steric hindrance score (0.0 = open, 1.0 = fully hindered)
    pub steric_score: f64,
    /// Whether β-hydrogen is available (for elimination reactions)
    pub has_beta_hydrogen: bool,
}

impl SiteDescriptor {
    /// Create descriptor with steric score calculated from degree
    pub fn new(
        site_kind: crate::chemistry::reactive_site::ReactiveSiteKind,
        degree: SubstitutionDegree,
        electronics: ElectronicEnvironment,
        bulky_substituents: u32,
        has_beta_hydrogen: bool,
    ) -> Self {
        let base_steric = degree.base_steric_score();
        let additional_steric = (bulky_substituents as f64 * 0.1).min(0.4);
        let steric_score = (base_steric + additional_steric).min(1.0);
        
        Self {
            site_kind,
            degree,
            electronics,
            steric_score,
            has_beta_hydrogen,
        }
    }
    
    /// Steric accessibility factor (1.0 = fully accessible, 0.0 = blocked)
    pub fn steric_accessibility(&self) -> f64 {
        (1.0 - self.steric_score).max(0.1)
    }
}

/// Selectivity evaluation conditions
/// 
/// This is separate from simulation::ReactionContext which tracks
/// actual reaction state (external reactants, catalysts, etc.)
#[derive(Debug, Clone)]
pub struct SelectivityContext {
    /// Temperature in Kelvin
    pub temperature: f64,
    /// pH if applicable (None for non-aqueous/neutral)
    pub ph: Option<f64>,
    /// Solvent type
    pub solvent_type: SolventType,
}

impl Default for SelectivityContext {
    fn default() -> Self {
        Self {
            temperature: 298.15,  // 25°C
            ph: None,
            solvent_type: SolventType::Neutral,
        }
    }
}

impl SelectivityContext {
    /// Create context for specific temperature
    pub fn at_temperature(kelvin: f64) -> Self {
        Self {
            temperature: kelvin,
            ..Default::default()
        }
    }
    
    /// Check if conditions are acidic (pH < 6)
    pub fn is_acidic(&self) -> bool {
        self.ph.map(|p| p < 6.0).unwrap_or(false)
    }
    
    /// Check if conditions are basic (pH > 8)
    pub fn is_basic(&self) -> bool {
        self.ph.map(|p| p > 8.0).unwrap_or(false)
    }
    
    /// Check if temperature is high (> 80°C)
    pub fn is_high_temperature(&self) -> bool {
        self.temperature > 353.15  // 80°C
    }
    
    /// Check if temperature is very high (> 150°C)
    pub fn is_very_high_temperature(&self) -> bool {
        self.temperature > 423.15  // 150°C
    }
}

/// Legacy alias for compatibility during transition
/// 
/// TODO: Remove this after all code is updated to use SelectivityContext
pub type ReactionContext = SelectivityContext;

/// Solvent classification for selectivity
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolventType {
    /// Protic polar (water, alcohols) - favors SN1, E1
    Protic,
    /// Aprotic polar (DMSO, DMF) - favors SN2
    AproticPolar,
    /// Non-polar (hexane, toluene)
    NonPolar,
    /// Basic conditions
    Basic,
    /// Acidic conditions
    Acidic,
    /// Neutral/default
    Neutral,
}

/// Type of organic reaction mechanism
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactionType {
    /// SN2 bimolecular substitution
    SN2,
    /// SN1 unimolecular substitution
    SN1,
    /// E2 bimolecular elimination
    E2,
    /// E1 unimolecular elimination
    E1,
    /// Acid-catalyzed esterification (Fischer)
    FischerEsterification,
    /// Nucleophilic addition to carbonyl
    CarbonylAddition,
    /// Electrophilic addition to alkene/alkyne
    ElectrophilicAddition,
}

/// Reactivity score with metadata
#[derive(Debug, Clone)]
pub struct ReactivityScore {
    /// Primary score (0.0 to 1.0+, >1.0 means faster than reference)
    pub value: f64,
    /// Activation energy modification in kJ/mol (negative = faster)
    pub activation_delta: f64,
    /// Pre-exponential factor multiplier
    pub pre_exp_multiplier: f64,
    /// Competing mechanisms with their scores
    pub competing: Vec<(CompetingMechanism, f64)>,
    /// Reasoning for this score
    pub reason: String,
}

impl ReactivityScore {
    /// Create new score with default competing mechanisms
    pub fn new(value: f64, reason: impl Into<String>) -> Self {
        Self {
            value: value.max(0.0),
            activation_delta: Self::value_to_activation_delta(value),
            pre_exp_multiplier: 1.0,
            competing: Vec::new(),
            reason: reason.into(),
        }
    }
    
    /// Create score from explicit activation energy delta
    pub fn with_activation_delta(delta: f64, reason: impl Into<String>) -> Self {
        Self {
            value: Self::activation_delta_to_value(delta),
            activation_delta: delta,
            pre_exp_multiplier: 1.0,
            competing: Vec::new(),
            reason: reason.into(),
        }
    }
    
    /// Add competing mechanism
    pub fn with_competing(mut self, mechanism: CompetingMechanism, score: f64) -> Self {
        self.competing.push((mechanism, score));
        self
    }
    
    /// Set pre-exponential multiplier
    pub fn with_pre_exp_multiplier(mut self, multiplier: f64) -> Self {
        self.pre_exp_multiplier = multiplier;
        self
    }
    
    /// Convert relative value to activation energy delta
    /// Uses approximation: ΔEa ≈ -RT ln(k/k0) at 298K
    fn value_to_activation_delta(value: f64) -> f64 {
        const R: f64 = 8.314;  // J/(mol·K)
        const T: f64 = 298.15; // K
        if value <= 0.0 {
            return 50.0; // Very slow
        }
        -R * T * value.ln() / 1000.0  // Convert to kJ/mol
    }
    
    /// Convert activation energy delta to relative value
    fn activation_delta_to_value(delta: f64) -> f64 {
        const R: f64 = 8.314;
        const T: f64 = 298.15;
        (-delta * 1000.0 / (R * T)).exp()
    }
}

/// Competing reaction mechanisms
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompetingMechanism {
    None,
    SN1,
    SN2,
    E1,
    E2,
    /// General elimination (E1/E2 unspecified)
    Elimination,
    /// General substitution
    Substitution,
    Rearrangement,
}

/// Selectivity evaluation recommendation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectivityRecommendation {
    /// Reaction proceeds exclusively (ratio > 10:1)
    Exclusive,
    /// Reaction is preferred but mixture expected (ratio 3:1 to 10:1)
    Preferred,
    /// Significant mixture of products (ratio 1:3 to 3:1)
    Mixed,
    /// Reaction is suppressed (ratio < 1:3)
    Suppressed,
    /// Reaction does not occur (score effectively zero)
    None,
}

/// Complete result of selectivity evaluation
#[derive(Debug, Clone)]
pub struct SelectivityResult {
    /// Primary mechanism score
    pub primary: ReactivityScore,
    /// All competing mechanisms with scores
    pub all_scores: HashMap<ReactionType, ReactivityScore>,
    /// Final recommendation
    pub recommendation: SelectivityRecommendation,
    /// Dominant competing mechanism (if any)
    pub dominant_competitor: Option<CompetingMechanism>,
}

impl SelectivityResult {
    /// Create result from primary score alone (no competition)
    pub fn exclusive(score: ReactivityScore) -> Self {
        let mut all_scores = HashMap::new();
        all_scores.insert(ReactionType::SN2, score.clone());
        
        Self {
            primary: score,
            all_scores,
            recommendation: SelectivityRecommendation::Exclusive,
            dominant_competitor: None,
        }
    }
    
    /// Determine recommendation from competing scores
    pub fn from_scores(
        primary_type: ReactionType,
        primary: ReactivityScore,
        competitors: Vec<(ReactionType, ReactivityScore)>,
    ) -> Self {
        let mut all_scores = HashMap::new();
        all_scores.insert(primary_type, primary.clone());
        
        let mut max_competitor_score = 0.0;
        let mut dominant_competitor = None;
        
        for (mech, score) in &competitors {
            all_scores.insert(*mech, score.clone());
            if score.value > max_competitor_score {
                max_competitor_score = score.value;
                dominant_competitor = CompetingMechanism::from_reaction_type(*mech);
            }
        }
        
        let ratio = if max_competitor_score > 0.0 {
            primary.value / max_competitor_score
        } else {
            f64::INFINITY
        };
        
        let recommendation = if primary.value < 0.01 {
            SelectivityRecommendation::None
        } else if ratio >= 10.0 {
            SelectivityRecommendation::Exclusive
        } else if ratio >= 3.0 {
            SelectivityRecommendation::Preferred
        } else if ratio >= 0.5 {
            SelectivityRecommendation::Mixed
        } else if ratio >= 0.1 {
            SelectivityRecommendation::Suppressed
        } else {
            SelectivityRecommendation::None
        };
        
        Self {
            primary,
            all_scores,
            recommendation,
            dominant_competitor,
        }
    }
}

impl CompetingMechanism {
    fn from_reaction_type(rt: ReactionType) -> Option<Self> {
        match rt {
            ReactionType::SN1 => Some(Self::SN1),
            ReactionType::SN2 => Some(Self::SN2),
            ReactionType::E1 => Some(Self::E1),
            ReactionType::E2 => Some(Self::E2),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn substitution_degree_base_scores() {
        assert!(SubstitutionDegree::Primary.base_steric_score() < SubstitutionDegree::Secondary.base_steric_score());
        assert!(SubstitutionDegree::Secondary.base_steric_score() < SubstitutionDegree::Tertiary.base_steric_score());
        assert!(SubstitutionDegree::Benzylic.base_steric_score() < SubstitutionDegree::Allylic.base_steric_score());
    }
    
    #[test]
    fn site_descriptor_steric_calculation() {
        let desc = SiteDescriptor::new(
            crate::chemistry::reactive_site::ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Secondary,
            ElectronicEnvironment::default(),
            2,  // two bulky groups
            true,
        );
        // Base 0.3 + 0.2 = 0.5
        assert!((desc.steric_score - 0.5).abs() < 0.01);
        assert!((desc.steric_accessibility() - 0.5).abs() < 0.01);
    }
    
    #[test]
    fn reactivity_score_activation_conversion() {
        let score = ReactivityScore::new(1.0, "reference");
        assert!(score.activation_delta.abs() < 0.1); // ~0 for k=k0
        
        let fast = ReactivityScore::new(2.0, "2x faster");
        assert!(fast.activation_delta < 0.0); // negative = faster
        
        let slow = ReactivityScore::new(0.5, "half speed");
        assert!(slow.activation_delta > 0.0); // positive = slower
    }
    
    #[test]
    fn selectivity_result_ratios() {
        let primary = ReactivityScore::new(1.0, "primary");
        let competitor = ReactivityScore::new(0.1, "competitor");
        
        let result = SelectivityResult::from_scores(
            ReactionType::SN2,
            primary,
            vec![(ReactionType::E2, competitor)],
        );
        
        assert_eq!(result.recommendation, SelectivityRecommendation::Exclusive);
    }
}
