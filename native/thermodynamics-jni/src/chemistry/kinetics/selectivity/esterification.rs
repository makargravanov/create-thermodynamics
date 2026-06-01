//! Fischer esterification selectivity with competition from elimination

use crate::chemistry::selectivity::elimination::evaluate_e2;
use crate::chemistry::selectivity::types::*;

/// Evaluate Fischer esterification (acid-catalyzed) selectivity
///
/// Key principles:
/// - 1° alcohols esterify well
/// - 2° alcohols esterify slower, compete with elimination at high T
/// - 3° alcohols DO NOT esterify (elimination dominates)
/// - Steric hindrance in acid also matters
/// - High temperature and strong acid favor elimination over esterification for 2°/3°
pub fn evaluate_fischer_esterification(
    acid: &SiteDescriptor,    // carboxylic acid
    alcohol: &SiteDescriptor, // alcohol
    context: &SelectivityContext,
) -> SelectivityResult {
    // Check: must be carboxylic acid (though site descriptor abstracts this)
    // and must have acidic conditions
    let acid_factor = match acid.degree {
        SubstitutionDegree::Primary => 1.0,   // acetic, propionic
        SubstitutionDegree::Secondary => 0.8, // isobutyric
        SubstitutionDegree::Tertiary => 0.4,  // pivalic - very hindered
        _ => 0.9,
    };

    // Main factor: alcohol degree
    let alcohol_factor = match alcohol.degree {
        SubstitutionDegree::Primary => 1.0,    // excellent
        SubstitutionDegree::Benzylic => 0.95,  // benzylic alcohols work well
        SubstitutionDegree::Allylic => 0.9,    // allylic also OK
        SubstitutionDegree::Secondary => 0.35, // slow, competitive
        SubstitutionDegree::Tertiary => 0.02,  // essentially no esterification
    };

    // Steric factor: multiply acid and alcohol hindrance
    let steric_factor = (1.0 - acid.steric_score * 0.5) * (1.0 - alcohol.steric_score * 0.7);

    // Electronic effects (minor for esterification)
    let electronic_factor = 1.0 + 0.1 * alcohol.electronics.net_effect();

    // Temperature effects: high T bad for 2°/3° (elimination competes)
    let temp_factor = match alcohol.degree {
        SubstitutionDegree::Primary
        | SubstitutionDegree::Benzylic
        | SubstitutionDegree::Allylic => {
            if context.is_very_high_temperature() {
                0.7
            } else {
                1.0
            }
        }
        SubstitutionDegree::Secondary => {
            if context.is_high_temperature() {
                0.5
            } else if context.is_very_high_temperature() {
                0.2
            } else {
                0.9
            }
        }
        SubstitutionDegree::Tertiary => 0.1, // always suppressed
    };

    let score = acid_factor * alcohol_factor * steric_factor * electronic_factor * temp_factor;

    let mut esterification = ReactivityScore::new(
        score,
        format!(
            "Fischer: acid={:?} ({:.2}), alcohol={:?} ({:.2}), T_factor={:.2}",
            acid.degree, acid_factor, alcohol.degree, alcohol_factor, temp_factor
        ),
    );

    // Activation energy adjustments
    let steric_ea = (acid.steric_score + alcohol.steric_score) * 10.0;
    esterification.activation_delta += steric_ea;

    // Pre-exponential factor: steric hindrance reduces collision efficiency
    esterification.pre_exp_multiplier =
        (1.0 - (acid.steric_score + alcohol.steric_score) * 0.3).max(0.5);

    // Evaluate competition from elimination (only for 2° and 3° alcohols)
    let competitors = if matches!(
        alcohol.degree,
        SubstitutionDegree::Secondary | SubstitutionDegree::Tertiary
    ) {
        // For esterification, we need acidic conditions
        // But alcohol can undergo E2 with acid as proton source
        let mut elimination_context = context.clone();
        elimination_context.solvent_type = SolventType::Acidic;

        let e2_score = evaluate_e2(alcohol, &elimination_context);
        vec![(ReactionType::E2, e2_score)]
    } else {
        vec![]
    };

    SelectivityResult::from_scores(
        ReactionType::FischerEsterification,
        esterification,
        competitors,
    )
}

/// Evaluate esterification rate for mixed alcohol systems
///
/// Returns relative rates for competing alcohols
pub fn esterification_competition(
    acid: &SiteDescriptor,
    alcohols: &[&SiteDescriptor],
    context: &SelectivityContext,
) -> Vec<(usize, f64)> {
    // (index, relative_rate)
    let scores: Vec<f64> = alcohols
        .iter()
        .map(|alc| {
            let result = evaluate_fischer_esterification(acid, alc, context);
            result.primary.value
        })
        .collect();

    let max_score = scores.iter().cloned().fold(0.0, f64::max);
    if max_score == 0.0 {
        return vec![];
    }

    scores
        .into_iter()
        .enumerate()
        .map(|(i, s)| (i, s / max_score))
        .collect()
}

/// Quick check: will this alcohol esterify under these conditions?
pub fn will_esterify(alcohol: &SiteDescriptor, context: &SelectivityContext) -> bool {
    // Primary always yes (with acid)
    if matches!(
        alcohol.degree,
        SubstitutionDegree::Primary | SubstitutionDegree::Benzylic
    ) {
        return true;
    }

    // Secondary: yes at moderate T, check competition
    if matches!(alcohol.degree, SubstitutionDegree::Secondary) {
        return !context.is_very_high_temperature();
    }

    // Tertiary: essentially no
    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::reactive_site::ReactiveSiteKind;

    fn alcohol_site(degree: SubstitutionDegree, steric: f64) -> SiteDescriptor {
        SiteDescriptor {
            site_kind: ReactiveSiteKind::Alcohol,
            degree,
            electronics: ElectronicEnvironment::default(),
            steric_score: steric,
            has_beta_hydrogen: true,
        }
    }

    fn acid_site(degree: SubstitutionDegree, steric: f64) -> SiteDescriptor {
        SiteDescriptor {
            site_kind: ReactiveSiteKind::CarboxylicAcid,
            degree,
            electronics: ElectronicEnvironment::default(),
            steric_score: steric,
            has_beta_hydrogen: false, // acids don't eliminate
        }
    }

    #[test]
    fn primary_alcohol_esterifies_well() {
        let context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 373.0, // 100°C
            ph: Some(2.0),
            ..Default::default()
        };

        let acetic = acid_site(SubstitutionDegree::Primary, 0.0);
        let ethanol = alcohol_site(SubstitutionDegree::Primary, 0.1);

        let result = evaluate_fischer_esterification(&acetic, &ethanol, &context);

        assert!(result.primary.value > 0.5);
        assert_eq!(result.recommendation, SelectivityRecommendation::Exclusive);
    }

    #[test]
    fn tertiary_alcohol_suppressed() {
        let context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 350.0,
            ph: Some(2.0),
            ..Default::default()
        };

        let acetic = acid_site(SubstitutionDegree::Primary, 0.0);
        let t_butanol = alcohol_site(SubstitutionDegree::Tertiary, 0.6);

        let result = evaluate_fischer_esterification(&acetic, &t_butanol, &context);

        assert!(result.primary.value < 0.1);
        assert!(result.dominant_competitor.is_some()); // elimination wins
    }

    #[test]
    fn secondary_temperature_dependence() {
        let cold_context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 330.0, // ~60°C
            ph: Some(2.0),
            ..Default::default()
        };

        let hot_context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 420.0, // ~150°C
            ph: Some(2.0),
            ..Default::default()
        };

        let acetic = acid_site(SubstitutionDegree::Primary, 0.0);
        let isopropanol = alcohol_site(SubstitutionDegree::Secondary, 0.3);

        let cold = evaluate_fischer_esterification(&acetic, &isopropanol, &cold_context);
        let hot = evaluate_fischer_esterification(&acetic, &isopropanol, &hot_context);

        // Cold: esterification preferred
        assert!(cold.primary.value > 0.2);
        // Hot: elimination competes strongly
        assert!(hot.primary.value < cold.primary.value);
        assert!(
            hot.dominant_competitor.is_some()
                || hot.recommendation == SelectivityRecommendation::Mixed
        );
    }

    #[test]
    fn steric_acid_effect() {
        let context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 373.0,
            ph: Some(2.0),
            ..Default::default()
        };

        let pivalic = acid_site(SubstitutionDegree::Tertiary, 0.6); // (CH3)3CCOOH
        let acetic = acid_site(SubstitutionDegree::Primary, 0.0);
        let ethanol = alcohol_site(SubstitutionDegree::Primary, 0.1);

        let with_pivalic = evaluate_fischer_esterification(&pivalic, &ethanol, &context);
        let with_acetic = evaluate_fischer_esterification(&acetic, &ethanol, &context);

        assert!(with_acetic.primary.value > with_pivalic.primary.value);
    }

    #[test]
    fn competition_multiple_alcohols() {
        let context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 350.0,
            ph: Some(2.0),
            ..Default::default()
        };

        let acetic = acid_site(SubstitutionDegree::Primary, 0.0);
        let ethanol = alcohol_site(SubstitutionDegree::Primary, 0.1);
        let isopropanol = alcohol_site(SubstitutionDegree::Secondary, 0.3);
        let t_butanol = alcohol_site(SubstitutionDegree::Tertiary, 0.6);

        let alcohols = vec![&ethanol, &isopropanol, &t_butanol];
        let rates = esterification_competition(&acetic, &alcohols, &context);

        // Primary should be fastest
        assert!(rates[0].1 > rates[1].1); // ethanol > isopropanol
        assert!(rates[1].1 > rates[2].1); // isopropanol > t-butanol (which is ~0)
    }

    #[test]
    fn benzylic_esterifies_well() {
        let context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 350.0,
            ph: Some(2.0),
            ..Default::default()
        };

        let acetic = acid_site(SubstitutionDegree::Primary, 0.0);
        let benzyl = alcohol_site(SubstitutionDegree::Benzylic, 0.2);

        let result = evaluate_fischer_esterification(&acetic, &benzyl, &context);

        // Should esterify nearly as well as primary
        assert!(result.primary.value > 0.6);
    }
}
