//! SN2/SN1 selectivity rules with E2 competition

use crate::chemistry::selectivity::types::*;

/// Evaluate SN2 reactivity with competing E2 mechanism
/// 
/// Key principles:
/// - SN2: 1° > benzylic/allylic > 2° > 3° (3° essentially zero)
/// - E2 competes at high temperature and for 2°/3° substrates
/// - Basic conditions favor E2
/// - Steric hindrance dramatically reduces SN2
pub fn evaluate_sn2(
    site: &SiteDescriptor,
    context: &SelectivityContext,
) -> SelectivityResult {
    let base_score = match site.degree {
        SubstitutionDegree::Primary => 1.0,
        SubstitutionDegree::Benzylic => 0.8,  // resonance stabilization helps SN2
        SubstitutionDegree::Allylic => 0.7,
        SubstitutionDegree::Secondary => 0.25,
        SubstitutionDegree::Tertiary => 0.005,  // effectively SN1 only
    };
    
    // Steric factor (critical for SN2)
    let steric_factor = site.steric_accessibility().powi(3);  // cubic penalty
    
    // Electronic effects (minor for SN2)
    let electronic_factor = 1.0 + site.electronics.net_effect() * 0.2;
    
    // Temperature: SN2 decreases at high temp (E2 competition)
    let temp_factor = if context.is_very_high_temperature() {
        0.7
    } else if context.is_high_temperature() {
        0.85
    } else {
        1.0
    };
    
    // Solvent effects
    let solvent_factor = match context.solvent_type {
        SolventType::AproticPolar => 1.5,  // SN2 loves polar aprotic
        SolventType::Protic => 0.6,        // protic slows SN2
        SolventType::Basic => 0.8,         // basic may start E2
        _ => 1.0,
    };
    
    let sn2_score = base_score * steric_factor * electronic_factor * temp_factor * solvent_factor;
    
    // Build primary score
    let mut sn2_reactivity = ReactivityScore::new(
        sn2_score,
        format!(
            "SN2: {:?}, steric={:.2}, electronic={:.2}",
            site.degree, site.steric_score, electronic_factor
        ),
    );
    
    // Adjust activation energy based on sterics
    let steric_ea_delta = site.steric_score * 20.0;  // up to +20 kJ/mol for hindered sites
    sn2_reactivity.activation_delta += steric_ea_delta;
    
    // Pre-exponential penalty for sterics
    sn2_reactivity.pre_exp_multiplier = (1.0 - site.steric_score * 0.5).max(0.3);
    
    // Evaluate E2 competition (only if β-hydrogen available)
    let competitors = if site.has_beta_hydrogen {
        let e2_score = evaluate_e2_competition(site, context);
        vec![(ReactionType::E2, e2_score)]
    } else {
        vec![]
    };
    
    SelectivityResult::from_scores(ReactionType::SN2, sn2_reactivity, competitors)
}

/// Quick E2 evaluation for competition purposes
fn evaluate_e2_competition(site: &SiteDescriptor, context: &SelectivityContext) -> ReactivityScore {
    let base_score = match site.degree {
        SubstitutionDegree::Tertiary => 1.0,
        SubstitutionDegree::Secondary => 0.5,
        SubstitutionDegree::Allylic => 0.4,
        SubstitutionDegree::Benzylic => 0.2,
        SubstitutionDegree::Primary => 0.05,  // E2 rare for 1°
    };
    
    // Temperature strongly favors E2
    let temp_bonus = if context.is_very_high_temperature() {
        2.0
    } else if context.is_high_temperature() {
        1.5
    } else {
        1.0
    };
    
    // Basic conditions strongly favor E2
    let base_bonus = if context.is_basic() {
        2.5
    } else if context.solvent_type == SolventType::Basic {
        2.0
    } else {
        1.0
    };
    
    let score = base_score * temp_bonus * base_bonus;
    
    let mut reactivity = ReactivityScore::new(
        score,
        format!("E2 competition: {:?}, temp_bonus={:.1}, base_bonus={:.1}",
            site.degree, temp_bonus, base_bonus),
    );
    
    // E2 has lower activation energy than SN2 for 3°
    if matches!(site.degree, SubstitutionDegree::Tertiary | SubstitutionDegree::Secondary) {
        reactivity.activation_delta -= 8.0;  // -8 kJ/mol advantage
    }
    
    reactivity
}

/// Evaluate SN1 (only relevant for 2° and 3°, especially benzylic/allylic)
pub fn evaluate_sn1(site: &SiteDescriptor, context: &SelectivityContext) -> ReactivityScore {
    let base_score = match site.degree {
        SubstitutionDegree::Tertiary => 1.0,
        SubstitutionDegree::Secondary => 0.2,
        SubstitutionDegree::Benzylic => 0.9,   // resonance stabilization!
        SubstitutionDegree::Allylic => 0.8,     // allylic cation stabilized
        SubstitutionDegree::Primary => 0.001, // essentially never
    };
    
    // Protic solvents favor SN1
    let solvent_factor = match context.solvent_type {
        SolventType::Protic => 1.5,
        SolventType::Acidic => 1.3,
        SolventType::AproticPolar => 0.5,
        SolventType::Basic => 0.1,  // base kills carbocations
        _ => 1.0,
    };
    
    // Electronic effects: EDG stabilize carbocation
    let edg_bonus = site.electronics.electron_donating_groups as f64 * 0.15;
    let ewg_penalty = site.electronics.electron_withdrawing_groups as f64 * 0.2;
    let electronic_factor = 1.0 + edg_bonus - ewg_penalty;
    
    // Resonance stabilization (already in degree for benzylic/allylic, but check environment)
    let resonance_bonus = if site.electronics.resonance_stabilization {
        0.3
    } else {
        0.0
    };
    
    let score = base_score * solvent_factor * (electronic_factor + resonance_bonus);
    
    let mut reactivity = ReactivityScore::new(
        score,
        format!(
            "SN1: {:?}, protic_factor={:.1}, EDG={}",
            site.degree, solvent_factor, 
            site.electronics.electron_donating_groups
        ),
    );
    
    // SN1 has very temperature-dependent activation (entropy driven)
    if context.is_high_temperature() {
        reactivity.activation_delta -= 5.0;
    }
    
    reactivity
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::reactive_site::ReactiveSiteKind;
    
    fn test_site(degree: SubstitutionDegree, steric: f64) -> SiteDescriptor {
        SiteDescriptor {
            site_kind: ReactiveSiteKind::Halide,
            degree,
            electronics: ElectronicEnvironment::default(),
            steric_score: steric,
            has_beta_hydrogen: true,
        }
    }
    
    #[test]
    fn sn2_primary_vs_tertiary() {
        let context = SelectivityContext::default();
        let primary = test_site(SubstitutionDegree::Primary, 0.0);
        let tertiary = test_site(SubstitutionDegree::Tertiary, 0.6);
        
        let r1 = evaluate_sn2(&primary, &context);
        let r2 = evaluate_sn2(&tertiary, &context);
        
        // Primary halides react well with SN2
        assert!(r1.primary.value > 0.5);
        // Tertiary halides have negligible SN2 rate (base rate 0.005)
        assert!(r2.primary.value < 0.1);
        // For tertiary with competing E2, SN2 is essentially None (not just suppressed)
        assert!(matches!(r2.recommendation, 
            SelectivityRecommendation::Suppressed | SelectivityRecommendation::None));
    }
    
    #[test]
    fn sn2_steric_hindrance_effect() {
        let context = SelectivityContext::default();
        let unhindered = test_site(SubstitutionDegree::Secondary, 0.1);
        let hindered = test_site(SubstitutionDegree::Secondary, 0.7);
        
        let r1 = evaluate_sn2(&unhindered, &context);
        let r2 = evaluate_sn2(&hindered, &context);
        
        assert!(r1.primary.value > r2.primary.value);
        // Hindered should have higher activation energy
        assert!(r1.primary.activation_delta < r2.primary.activation_delta);
    }
    
    #[test]
    fn temperature_favors_e2_competition() {
        let context_cold = SelectivityContext::at_temperature(298.0);
        let context_hot = SelectivityContext::at_temperature(400.0);
        let secondary = test_site(SubstitutionDegree::Secondary, 0.3);
        
        let _r_cold = evaluate_sn2(&secondary, &context_cold);
        let r_hot = evaluate_sn2(&secondary, &context_hot);
        
        // Hot conditions should have stronger E2 competition
        assert!(r_hot.dominant_competitor.is_some());
    }
    
    #[test]
    fn benzylic_sn1_favored() {
        let context = SelectivityContext {
            solvent_type: SolventType::Protic,
            ..Default::default()
        };
        let benzylic = test_site(SubstitutionDegree::Benzylic, 0.2);
        
        let sn1 = evaluate_sn1(&benzylic, &context);
        assert!(sn1.value > 0.5);  // Benzylic SN1 is viable
    }
}


