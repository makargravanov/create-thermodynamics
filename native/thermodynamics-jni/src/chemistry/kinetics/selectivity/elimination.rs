//! E1/E2 elimination selectivity with competition handling

use crate::chemistry::selectivity::types::*;

/// Evaluate E2 elimination reactivity
/// 
/// Key principles:
/// - E2: 3° > 2° > 1° (strong base required for 1°)
/// - Requires β-hydrogen
/// - Basic conditions strongly favor E2
/// - High temperature favors E2
/// - Steric bulk at β-carbon affects Zaitsev/Hofmann selectivity
pub fn evaluate_e2(
    site: &SiteDescriptor,
    context: &SelectivityContext,
) -> ReactivityScore {
    // Check prerequisite: must have beta-hydrogen
    if !site.has_beta_hydrogen {
        return ReactivityScore::new(0.0, "No β-hydrogen available for elimination");
    }
    
    let base_score = match site.degree {
        SubstitutionDegree::Tertiary => 1.0,
        SubstitutionDegree::Secondary => 0.6,
        SubstitutionDegree::Allylic => 0.5,
        SubstitutionDegree::Benzylic => 0.4,
        SubstitutionDegree::Primary => 0.15,  // needs strong base
    };
    
    // Temperature strongly favors elimination (higher entropy)
    let temp_factor = if context.is_very_high_temperature() {
        2.5
    } else if context.is_high_temperature() {
        1.8
    } else {
        1.0
    };
    
    // Basic conditions strongly favor E2
    let base_factor = if context.is_basic() {
        3.0
    } else if context.solvent_type == SolventType::Basic {
        2.5
    } else if context.ph.map(|p| p > 10.0).unwrap_or(false) {
        2.0
    } else {
        0.3  // E2 needs base
    };
    
    // Solvent effects
    let solvent_factor = match context.solvent_type {
        SolventType::AproticPolar => 1.4,  // polar aprotic favors elimination
        SolventType::Protic => 0.8,         // protic slows E2
        _ => 1.0,
    };
    
    // Steric effects: bulky substrate favors elimination over substitution
    // (SN2 is more sterically demanding)
    let steric_bonus = site.steric_score * 0.5;
    
    let score = base_score * temp_factor * base_factor * solvent_factor * (1.0 + steric_bonus);
    
    let mut reactivity = ReactivityScore::new(
        score,
        format!(
            "E2: {:?}, temp={:.1}, base={:.1}, has_β-H={}",
            site.degree, temp_factor, base_factor, site.has_beta_hydrogen
        ),
    );
    
    // E2 activation energy decreases with temperature (entropy driven)
    let temp_ea_bonus = if context.is_very_high_temperature() {
        -12.0  // -12 kJ/mol at high T
    } else if context.is_high_temperature() {
        -6.0
    } else {
        0.0
    };
    reactivity.activation_delta += temp_ea_bonus;
    
    reactivity
}

/// Evaluate E1 elimination (competes with SN1)
/// 
/// E1 is generally preferred when:
/// - Protic solvents
/// - Tertiary substrates with stable carbocations
/// - High temperature favors E1 over SN1
pub fn evaluate_e1(site: &SiteDescriptor, context: &SelectivityContext) -> ReactivityScore {
    let base_score = match site.degree {
        SubstitutionDegree::Tertiary => 0.7,
        SubstitutionDegree::Secondary => 0.15,
        SubstitutionDegree::Benzylic => 0.6,   // benzylic cation stable
        SubstitutionDegree::Allylic => 0.5,  // allylic cation stable
        SubstitutionDegree::Primary => 0.001, // never E1 for 1°
    };
    
    // Must have beta-hydrogen
    if !site.has_beta_hydrogen {
        return ReactivityScore::new(0.0, "E1 requires β-hydrogen");
    }
    
    // Temperature: E1 more entropy-driven than SN1
    let temp_factor = if context.is_very_high_temperature() {
        1.6
    } else if context.is_high_temperature() {
        1.3
    } else {
        1.0
    };
    
    // Solvent: E1 loves protic like SN1
    let solvent_factor = match context.solvent_type {
        SolventType::Protic => 1.4,
        SolventType::Acidic => 1.2,
        SolventType::Basic => 0.05,  // base destroys E1
        _ => 1.0,
    };
    
    // Electronic: carbocation stability
    let edg_stabilization = site.electronics.electron_donating_groups as f64 * 0.2;
    let resonance_bonus = if site.electronics.resonance_stabilization {
        0.4
    } else {
        0.0
    };
    let electronic_factor = 1.0 + edg_stabilization + resonance_bonus;
    
    let score = base_score * temp_factor * solvent_factor * electronic_factor;
    
    let mut reactivity = ReactivityScore::new(
        score,
        format!(
            "E1: {:?}, protic={:.1}, EDG={}, resonance={}",
            site.degree, solvent_factor, 
            site.electronics.electron_donating_groups,
            site.electronics.resonance_stabilization
        ),
    );
    
    // High temperature favors E1 over SN1 by entropy
    if context.is_high_temperature() {
        reactivity.activation_delta -= 4.0;
    }
    
    reactivity
}

/// Determine Zaitsev vs Hofmann preference for E2
/// 
/// Returns ratio of Zaitsev (more substituted alkene) to Hofmann (less substituted)
pub fn zaitsev_hofmann_ratio(
    site: &SiteDescriptor,
    context: &SelectivityContext,
) -> f64 {
    // Base factor: bulky bases favor Hofmann
    let base_factor = match context.solvent_type {
        SolventType::Basic => 0.5,  // strong bulky base -> Hofmann
        _ => 1.0,
    };
    
    // Substrate factor: steric hindrance at β-carbon favors Hofmann
    let steric_factor = 1.0 - site.steric_score * 0.5;
    
    // Temperature: high T slightly favors more stable (Zaitsev)
    let temp_factor = if context.is_high_temperature() { 1.2 } else { 1.0 };
    
    // Base ratio (Zaitsev preferred under normal conditions)
    let base_ratio = 4.0;
    
    base_ratio * base_factor * steric_factor * temp_factor
}

/// Determine E2 vs E1 preference for elimination
/// 
/// E2 is generally preferred when:
/// - Strong base present
/// - Good leaving group
/// - Primary/secondary substrates
/// E1 is preferred when:
/// - Protic solvents
/// - Tertiary substrates with stable carbocations
pub fn elimination_mechanism_preference(
    site: &SiteDescriptor,
    context: &SelectivityContext,
) -> (f64, f64) {  // (E2_score, E1_score)
    let e2 = evaluate_e2(site, context);
    let e1 = evaluate_e1(site, context);
    (e2.value, e1.value)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::reactive_site::ReactiveSiteKind;
    
    fn test_site(degree: SubstitutionDegree, has_beta: bool) -> SiteDescriptor {
        SiteDescriptor {
            site_kind: ReactiveSiteKind::Alcohol,
            degree,
            electronics: ElectronicEnvironment::default(),
            steric_score: 0.3,
            has_beta_hydrogen: has_beta,
        }
    }
    
    #[test]
    fn e2_requires_beta_hydrogen() {
        let context = SelectivityContext {
            solvent_type: SolventType::Basic,
            temperature: 350.0,
            ph: Some(12.0),
        };
        
        let with_beta = test_site(SubstitutionDegree::Secondary, true);
        let no_beta = test_site(SubstitutionDegree::Secondary, false);
        
        let r1 = evaluate_e2(&with_beta, &context);
        let r2 = evaluate_e2(&no_beta, &context);
        
        assert!(r1.value > 0.5);
        assert!(r2.value < 0.01);
    }
    
    #[test]
    fn e2_tertiary_vs_primary() {
        let context = SelectivityContext {
            solvent_type: SolventType::Basic,
            temperature: 350.0,
            ph: Some(12.0),
        };
        
        let tertiary = test_site(SubstitutionDegree::Tertiary, true);
        let primary = test_site(SubstitutionDegree::Primary, true);
        
        let r3 = evaluate_e2(&tertiary, &context);
        let r1 = evaluate_e2(&primary, &context);
        
        assert!(r3.value > r1.value);
        assert!(r3.value > 2.0);  // Strong base + heat = good E2
    }
    
    #[test]
    fn basic_conditions_required_for_e2() {
        let context_neutral = SelectivityContext::default();
        let context_basic = SelectivityContext {
            solvent_type: SolventType::Basic,
            ph: Some(12.0),
            ..Default::default()
        };
        
        let secondary = test_site(SubstitutionDegree::Secondary, true);
        
        let r_neutral = evaluate_e2(&secondary, &context_neutral);
        let r_basic = evaluate_e2(&secondary, &context_basic);
        
        assert!(r_basic.value > r_neutral.value * 5.0);
    }
    
    #[test]
    fn e1_tertiary_in_protic() {
        let context = SelectivityContext {
            solvent_type: SolventType::Protic,
            temperature: 350.0,
            ..Default::default()
        };
        
        let tertiary = test_site(SubstitutionDegree::Tertiary, true);
        
        let r = evaluate_e1(&tertiary, &context);
        assert!(r.value > 0.5);
    }
    
    #[test]
    fn temperature_favors_elimination() {
        let context_cold = SelectivityContext::at_temperature(298.0);
        let context_hot = SelectivityContext::at_temperature(400.0);
        
        let secondary = test_site(SubstitutionDegree::Secondary, true);
        
        let r_cold = evaluate_e2(&secondary, &context_cold);
        let r_hot = evaluate_e2(&secondary, &context_hot);
        
        assert!(r_hot.value > r_cold.value);
    }
    
    #[test]
    fn zaitsev_favored_over_hofmann() {
        let context = SelectivityContext {
            solvent_type: SolventType::Basic,
            temperature: 350.0,
            ph: Some(12.0),
        };
        
        let secondary = test_site(SubstitutionDegree::Secondary, true);
        let ratio = zaitsev_hofmann_ratio(&secondary, &context);
        
        // Zaitsev should be preferred (ratio > 1)
        assert!(ratio > 1.0);
        // But not exclusively (real world ~4:1 to 10:1)
        assert!(ratio < 20.0);
    }
}


