//! Carbonyl addition selectivity (aldehydes vs ketones)

use crate::chemistry::selectivity::types::*;

/// Evaluate nucleophilic addition to carbonyl selectivity
///
/// Key principles:
/// - Aldehydes more reactive than ketones (steric + electronic)
/// - Electronic: EDG on carbonyl carbon decrease reactivity (stabilize C=O)
/// - Steric hindrance crucial
/// - Aromatic aldehydes less reactive (resonance deactivation)
/// - Cyclic ketones more reactive than acyclic (angle strain)
pub fn evaluate_carbonyl_addition(
    carbonyl: &SiteDescriptor,
    nucleophile_strength: NucleophileStrength,
    context: &SelectivityContext,
) -> ReactivityScore {
    // Base reactivity: aldehyde vs ketone
    let carbonyl_factor = if carbonyl.degree == SubstitutionDegree::Primary {
        1.0 // aldehyde R-CHO
    } else {
        0.3 // ketone R2C=O (sterica + electronic effect)
    };

    // Steric factor: major contributor for ketones
    // Aldehydes: H is small; Ketones: two R groups
    let steric_factor = if carbonyl.degree == SubstitutionDegree::Primary {
        0.95 // aldehydes slightly hindered
    } else {
        (0.4 - carbonyl.steric_score * 0.2).max(0.15)
    };

    // Electronic effects: EDG destabilize transition state
    // EWG activate carbonyl carbon
    let electronic_factor = if carbonyl.electronics.electron_donating_groups > 0 {
        1.0 - carbonyl.electronics.electron_donating_groups as f64 * 0.15
    } else if carbonyl.electronics.electron_withdrawing_groups > 0 {
        1.0 + carbonyl.electronics.electron_withdrawing_groups as f64 * 0.2
    } else {
        1.0
    };

    // Aromatic deactivation (resonance with carbonyl)
    let aromatic_penalty = if carbonyl.electronics.aromatic {
        0.6 // aromatic aldehydes (benzaldehyde) less reactive
    } else {
        1.0
    };

    // Nucleophile strength factor
    let nucleophile_factor = match nucleophile_strength {
        NucleophileStrength::VeryStrong => 2.0, // RMgX, RLi
        NucleophileStrength::Strong => 1.5,     // NaBH4
        NucleophileStrength::Moderate => 1.0,   // NaOH/H2O addition
        NucleophileStrength::Weak => 0.5,       // water, alcohols
    };

    // Solvent effects
    let solvent_factor = match context.solvent_type {
        SolventType::AproticPolar => 1.3, // polar aprotic accelerates
        SolventType::Protic => 0.8,       // protic slows (H-bonding to nucleophile)
        _ => 1.0,
    };

    let score = carbonyl_factor
        * steric_factor
        * electronic_factor
        * aromatic_penalty
        * nucleophile_factor
        * solvent_factor;

    let mut reactivity = ReactivityScore::new(
        score,
        format!(
            "Carbonyl add: {:?}, steric={:.2}, EDG={}, aromatic={}",
            carbonyl.degree,
            steric_factor,
            carbonyl.electronics.electron_donating_groups,
            carbonyl.electronics.aromatic
        ),
    );

    // Steric activation energy penalty (major for ketones)
    let steric_ea = match carbonyl.degree {
        SubstitutionDegree::Primary => carbonyl.steric_score * 5.0, // 0-5 kJ/mol
        SubstitutionDegree::Secondary => carbonyl.steric_score * 15.0, // 0-15 kJ/mol
        _ => carbonyl.steric_score * 10.0,
    };
    reactivity.activation_delta += steric_ea;

    // Strong nucleophiles reduce steric sensitivity
    let nucleophile_steric_relief = match nucleophile_strength {
        NucleophileStrength::VeryStrong => 0.5,
        NucleophileStrength::Strong => 0.7,
        _ => 1.0,
    };
    reactivity.activation_delta *= nucleophile_steric_relief;

    reactivity
}

/// Compare reactivity of aldehyde vs ketone
///
/// Returns relative rates: (aldehyde_rate, ketone_rate)
pub fn aldehyde_vs_ketone(
    aldehyde: &SiteDescriptor,
    ketone: &SiteDescriptor,
    nucleophile: NucleophileStrength,
    context: &SelectivityContext,
) -> (f64, f64) {
    let ald_score = evaluate_carbonyl_addition(aldehyde, nucleophile, context);
    let ket_score = evaluate_carbonyl_addition(ketone, nucleophile, context);
    (ald_score.value, ket_score.value)
}

/// Evaluate selectivity when both aldehyde and ketone present
///
/// Returns aldehyde selectivity factor (higher = more selective for aldehyde)
pub fn chemoselectivity_aldehyde_over_ketone(
    aldehyde: &SiteDescriptor,
    ketone: &SiteDescriptor,
    nucleophile: NucleophileStrength,
    context: &SelectivityContext,
) -> f64 {
    let (ald_rate, ket_rate) = aldehyde_vs_ketone(aldehyde, ketone, nucleophile, context);
    if ket_rate == 0.0 {
        return f64::INFINITY;
    }
    ald_rate / ket_rate
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::reactive_site::ReactiveSiteKind;

    fn aldehyde_site(aromatic: bool, steric: f64) -> SiteDescriptor {
        let mut electronics = ElectronicEnvironment::default();
        electronics.aromatic = aromatic;

        SiteDescriptor {
            site_kind: ReactiveSiteKind::Aldehyde,
            degree: SubstitutionDegree::Primary,
            electronics,
            steric_score: steric,
            has_beta_hydrogen: true,
        }
    }

    fn ketone_site(edg: u32, ewg: u32, steric: f64) -> SiteDescriptor {
        let electronics = ElectronicEnvironment {
            electron_donating_groups: edg,
            electron_withdrawing_groups: ewg,
            resonance_stabilization: false,
            aromatic: false,
        };

        SiteDescriptor {
            site_kind: ReactiveSiteKind::Ketone,
            degree: SubstitutionDegree::Secondary,
            electronics,
            steric_score: steric,
            has_beta_hydrogen: true,
        }
    }

    #[test]
    fn aldehyde_more_reactive_than_ketone() {
        let context = SelectivityContext::default();
        let acetaldehyde = aldehyde_site(false, 0.1);
        let acetone = ketone_site(0, 0, 0.1);

        let ald =
            evaluate_carbonyl_addition(&acetaldehyde, NucleophileStrength::Moderate, &context);
        let ket = evaluate_carbonyl_addition(&acetone, NucleophileStrength::Moderate, &context);

        assert!(ald.value > ket.value);
        assert!(ald.value > 0.5);
        assert!(ket.value < 0.4);
    }

    #[test]
    fn aromatic_aldehyde_deactivated() {
        let context = SelectivityContext::default();
        let benzaldehyde = aldehyde_site(true, 0.1);
        let acetaldehyde = aldehyde_site(false, 0.1);

        let arom =
            evaluate_carbonyl_addition(&benzaldehyde, NucleophileStrength::Moderate, &context);
        let aliph =
            evaluate_carbonyl_addition(&acetaldehyde, NucleophileStrength::Moderate, &context);

        assert!(aliph.value > arom.value);
        // Aromatic penalty is substantial
        assert!(arom.value < aliph.value * 0.8);
    }

    #[test]
    fn edg_deactivate_ketone() {
        let context = SelectivityContext::default();
        let plain = ketone_site(0, 0, 0.1);
        let activated = ketone_site(2, 0, 0.1); // Two EDG

        let r1 = evaluate_carbonyl_addition(&plain, NucleophileStrength::Moderate, &context);
        let r2 = evaluate_carbonyl_addition(&activated, NucleophileStrength::Moderate, &context);

        assert!(r1.value > r2.value);
    }

    #[test]
    fn ewg_activate_ketone() {
        let context = SelectivityContext::default();
        let plain = ketone_site(0, 0, 0.1);
        let activated = ketone_site(0, 2, 0.1); // Two EWG

        let r1 = evaluate_carbonyl_addition(&plain, NucleophileStrength::Moderate, &context);
        let r2 = evaluate_carbonyl_addition(&activated, NucleophileStrength::Moderate, &context);

        assert!(r2.value > r1.value);
    }

    #[test]
    fn steric_hindrance_ketone() {
        let context = SelectivityContext::default();
        let unhindered = ketone_site(0, 0, 0.1);
        let hindered = ketone_site(0, 0, 0.6); // Sterically hindered ketone

        let r1 = evaluate_carbonyl_addition(&unhindered, NucleophileStrength::Moderate, &context);
        let r2 = evaluate_carbonyl_addition(&hindered, NucleophileStrength::Moderate, &context);

        // Unhindered should be faster than hindered
        assert!(r1.value > r2.value);
        // The difference includes both steric factor and Ea penalty
        // steric_score=0.1: factor 0.38, Ea penalty 1.5
        // steric_score=0.6: factor 0.28, Ea penalty 9.0
        assert!(r1.value > r2.value * 1.2);
    }

    #[test]
    fn nucleophile_strength_effect() {
        let context = SelectivityContext::default();
        let acetone = ketone_site(0, 0, 0.3);

        let weak = evaluate_carbonyl_addition(&acetone, NucleophileStrength::Weak, &context);
        let strong =
            evaluate_carbonyl_addition(&acetone, NucleophileStrength::VeryStrong, &context);

        assert!(strong.value > weak.value);
    }

    #[test]
    fn chemoselectivity_strong() {
        let context = SelectivityContext::default();
        let acetaldehyde = aldehyde_site(false, 0.1);
        let acetone = ketone_site(0, 0, 0.1);

        let selectivity = chemoselectivity_aldehyde_over_ketone(
            &acetaldehyde,
            &acetone,
            NucleophileStrength::Moderate,
            &context,
        );

        // Should be quite selective for aldehyde
        assert!(selectivity > 2.0);
    }
}
