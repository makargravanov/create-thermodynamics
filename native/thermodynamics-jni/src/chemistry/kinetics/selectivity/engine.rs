//! Central selectivity engine for evaluating reaction favorability
//!
//! Coordinates all selectivity rules and provides unified API for generators.

use super::carbonyl_addition::evaluate_carbonyl_addition;
use super::elimination::{evaluate_e1, evaluate_e2};
use super::esterification::evaluate_fischer_esterification;
use super::nucleophilic_substitution::evaluate_sn2;
use super::types::*;
use crate::chemistry::mixture::TRACE_CONCENTRATION_MOL_PER_BUCKET;
use crate::chemistry::molecule::{bond_order_matches, MolecularStructure};
use crate::chemistry::reactive_site::ReactiveSiteKind;

/// Central engine for selectivity evaluation
pub struct SelectivityEngine;

impl SelectivityEngine {
    pub fn evaluate_profile(
        profile: &SelectivityProfile,
        context: &SelectivityContext,
    ) -> SelectivityRuntimeEffect {
        let score = match profile.mechanism {
            ReactionType::SN2 => evaluate_sn2(&profile.primary_site, context).primary,
            ReactionType::SN1 => {
                super::nucleophilic_substitution::evaluate_sn1(&profile.primary_site, context)
            }
            ReactionType::E2 => evaluate_e2(&profile.primary_site, context),
            ReactionType::E1 => evaluate_e1(&profile.primary_site, context),
            ReactionType::FischerEsterification => {
                let Some(alcohol) = profile.secondary_site.as_ref() else {
                    return missing_secondary_site(profile);
                };
                evaluate_fischer_esterification(&profile.primary_site, alcohol, context).primary
            }
            ReactionType::CarbonylAddition => {
                let nucleophile_strength = profile
                    .nucleophile_strength
                    .unwrap_or(NucleophileStrength::Moderate);
                let mut score = evaluate_carbonyl_addition(
                    &profile.primary_site,
                    nucleophile_strength,
                    context,
                );
                apply_inorganic_environment_to_carbonyl_score(
                    &mut score,
                    nucleophile_strength,
                    context,
                );
                score
            }
            ReactionType::CarbonylReduction => {
                let mut score = evaluate_carbonyl_addition(
                    &profile.primary_site,
                    profile
                        .nucleophile_strength
                        .unwrap_or(NucleophileStrength::Strong),
                    context,
                );
                if context.is_acidic() {
                    score.value *= 0.2;
                    score.activation_delta += 8.0;
                    score.reason = format!(
                        "hydride is quenched in strongly acidic medium; {}",
                        score.reason
                    );
                } else if context.is_oxidizing() {
                    let redox_penalty = redox_competition_penalty(context);
                    score.value *= redox_penalty;
                    score.activation_delta += 6.0 + (1.0 - redox_penalty) * 8.0;
                    score.reason =
                        format!("oxidizing medium competes with reduction; {}", score.reason);
                } else {
                    score.reason = format!("hydride carbonyl reduction: {}", score.reason);
                }
                score
            }
            ReactionType::OrganicOxidation => {
                let mut score = ReactivityScore::new(1.0, "organic oxidation");
                if context.is_oxidizing() {
                    score.value *= 1.5 + context.oxidizing_strength.min(4.0) * 0.5;
                    score.activation_delta -= (context.oxidizing_strength.min(4.0) * 1.5).max(1.0);
                    score.reason = "oxidizing medium favors organic oxidation".to_string();
                }
                if context.is_reducing() {
                    let penalty = redox_competition_penalty(context);
                    score.value *= penalty;
                    score.activation_delta += 5.0 + (1.0 - penalty) * 8.0;
                    score.reason = "reducing medium competes with organic oxidation".to_string();
                }
                if context.is_water_rich()
                    && matches!(
                        profile.primary_site.site_kind,
                        ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Alcohol
                    )
                {
                    score.value *= 1.15;
                    score.activation_delta -= 0.8;
                }
                if context.medium.is_supercritical()
                    && matches!(context.solvent_type, SolventType::NonPolar)
                {
                    score.value *= 0.85;
                    score.activation_delta += 1.0;
                }
                score
            }
            ReactionType::WittigOlefination
            | ReactionType::HornerWadsworthEmmonsOlefination
            | ReactionType::JuliaOlefination => {
                let Some(olefination_reagent) = profile.secondary_site.as_ref() else {
                    return missing_secondary_site(profile);
                };
                let mut score = evaluate_carbonyl_addition(
                    &profile.primary_site,
                    profile
                        .nucleophile_strength
                        .unwrap_or(NucleophileStrength::VeryStrong),
                    context,
                );
                score.value *= olefination_reagent.steric_accessibility().max(0.2);
                if olefination_reagent.electronics.resonance_stabilization {
                    score.value *= 0.85;
                    score.activation_delta += 1.5;
                }
                score.reason = format!("carbonyl olefination: {}", score.reason);
                score
            }
            ReactionType::ElectrophilicAddition => ReactivityScore::new(
                1.0,
                "electrophilic addition has no specialized selectivity profile yet",
            ),
            ReactionType::SilylEtherFormation
            | ReactionType::SilylEtherCleavage
            | ReactionType::AcetalFormation
            | ReactionType::AcetalHydrolysis
            | ReactionType::CarbamateFormation
            | ReactionType::CarbamateCleavage
            | ReactionType::EsterProtection
            | ReactionType::EsterHydrolysis => evaluate_protecting_group_profile(profile, context),
            ReactionType::AlphaHalogenation
            | ReactionType::AldolAddition
            | ReactionType::AldolDehydration
            | ReactionType::EnamineFormation
            | ReactionType::EnolateAlkylation
            | ReactionType::MichaelAddition
            | ReactionType::ClaisenCondensation
            | ReactionType::PhosphoniumSaltFormation
            | ReactionType::PhosphoniumYlideFormation => {
                evaluate_alpha_carbon_profile(profile, context)
            }
            ReactionType::Lactonization
            | ReactionType::Lactamization
            | ReactionType::HeterocycleCondensation => {
                evaluate_cyclization_profile(profile, context)
            }
            ReactionType::DielsAlder => {
                let mut score = ReactivityScore::new(1.0, "Diels-Alder [4+2] cycloaddition");
                score.value *= profile.primary_site.steric_accessibility().max(0.2);
                if let Some(dienophile) = profile.secondary_site.as_ref() {
                    // Electron-poor dienophiles accelerate the normal-demand cycloaddition.
                    if dienophile.electronics.electron_withdrawing_groups >= 1 {
                        score.value *= 1.0
                            + 0.4 * dienophile.electronics.electron_withdrawing_groups as f64;
                        score.activation_delta -= 3.0
                            * dienophile.electronics.electron_withdrawing_groups.min(3) as f64;
                        score.reason =
                            "electron-poor dienophile accelerates Diels-Alder".to_string();
                    }
                    score.value *= dienophile.steric_accessibility().max(0.2);
                }
                if context.is_high_temperature() {
                    score.value *= 1.2;
                    score.activation_delta -= 2.0;
                }
                score
            }
            ReactionType::NAlkylation => evaluate_sn2(&profile.primary_site, context).primary,
        };
        let recommendation = if matches!(profile.mechanism, ReactionType::SN2) {
            evaluate_sn2(&profile.primary_site, context).recommendation
        } else if matches!(profile.mechanism, ReactionType::FischerEsterification) {
            match profile.secondary_site.as_ref() {
                Some(alcohol) => {
                    evaluate_fischer_esterification(&profile.primary_site, alcohol, context)
                        .recommendation
                }
                None => SelectivityRecommendation::None,
            }
        } else if score.value < 0.01 {
            SelectivityRecommendation::None
        } else if score.value < 0.1 {
            SelectivityRecommendation::Suppressed
        } else {
            SelectivityRecommendation::Preferred
        };
        let suppressed = matches!(
            profile.suppression_policy,
            SelectivitySuppressionPolicy::SuppressWhenDisfavored
        ) && matches!(
            recommendation,
            SelectivityRecommendation::Suppressed | SelectivityRecommendation::None
        );
        SelectivityRuntimeEffect {
            rate_multiplier: if suppressed {
                0.0
            } else {
                score.value.max(0.0)
            },
            activation_delta_kj_per_mol: score.activation_delta,
            pre_exp_multiplier: score.pre_exp_multiplier,
            suppressed,
            reason: score.reason,
        }
    }

    /// Evaluate SN2 substitution with E2 competition
    pub fn sn2_with_competition(
        site: &SiteDescriptor,
        context: &SelectivityContext,
    ) -> SelectivityResult {
        evaluate_sn2(site, context)
    }

    /// Evaluate E2 elimination (direct)
    pub fn e2_elimination(site: &SiteDescriptor, context: &SelectivityContext) -> ReactivityScore {
        evaluate_e2(site, context)
    }

    /// Evaluate E1 elimination (usually competing with SN1)
    pub fn e1_elimination(site: &SiteDescriptor, context: &SelectivityContext) -> ReactivityScore {
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
        ald_score.value > ket_score.value * 2.0 // At least 2:1 selectivity
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

fn apply_inorganic_environment_to_carbonyl_score(
    score: &mut ReactivityScore,
    nucleophile_strength: NucleophileStrength,
    context: &SelectivityContext,
) {
    if nucleophile_strength != NucleophileStrength::VeryStrong {
        return;
    }
    let mut reasons = Vec::new();
    if context.is_water_rich() {
        score.value *= 0.02;
        score.activation_delta += 24.0;
        reasons.push("water quenches strongly basic organometallic nucleophiles");
    } else if context.water_activity > 0.02 {
        score.value *= 0.2;
        score.activation_delta += 10.0;
        reasons.push("trace water competes with organometallic addition");
    }
    if context.is_oxygen_rich() || context.is_oxidizing() {
        score.value *= 0.25;
        score.activation_delta += 8.0;
        reasons.push("oxidizing medium consumes strongly reducing organometallic reagent");
    }
    if context.has_free_complexable_metal() {
        let metal_penalty =
            (1.0 / (1.0 + context.free_complexable_metal_activity * 20.0)).clamp(0.05, 1.0);
        score.value *= metal_penalty;
        score.activation_delta += (1.0 - metal_penalty) * 12.0;
        reasons.push("free metal ions coordinate or quench the organometallic reagent");
    }
    if !reasons.is_empty() {
        score.reason = format!("{}; {}", score.reason, reasons.join("; "));
    }
}

fn redox_competition_penalty(context: &SelectivityContext) -> f64 {
    let total = context.oxidizing_strength + context.reducing_strength;
    if total <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
        return 1.0;
    }
    let oxidizing_fraction = context.oxidizing_strength / total;
    (1.0 - oxidizing_fraction * 0.9).clamp(0.05, 1.0)
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
            0,
            0,
            0,
            true,
            false,
            false,
        )
    }

    /// Create simple secondary site
    pub fn secondary_alcohol() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Secondary,
            0,
            0,
            0,
            true,
            false,
            false,
        )
    }

    /// Create simple tertiary site
    pub fn tertiary_alcohol() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Tertiary,
            0,
            0,
            0,
            false,
            false,
            false,
        )
    }

    /// Create benzylic site
    pub fn benzylic_alcohol() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Benzylic,
            0,
            0,
            0,
            true,
            true,
            false,
        )
    }

    /// Create simple aldehyde
    pub fn aldehyde() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Aldehyde,
            SubstitutionDegree::Primary,
            0,
            0,
            0,
            true,
            false,
            false,
        )
    }

    /// Create simple ketone
    pub fn ketone() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Ketone,
            SubstitutionDegree::Secondary,
            0,
            0,
            0,
            true,
            false,
            false,
        )
    }

    /// Create aromatic aldehyde (benzaldehyde-like)
    pub fn aromatic_aldehyde() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Aldehyde,
            SubstitutionDegree::Primary,
            0,
            0,
            0,
            false,
            true,
            true,
        )
    }

    /// Create primary halide
    pub fn primary_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Primary,
            0,
            0,
            0,
            true,
            false,
            false,
        )
    }

    /// Create secondary halide
    pub fn secondary_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Secondary,
            0,
            0,
            0,
            true,
            false,
            false,
        )
    }

    /// Create tertiary halide
    pub fn tertiary_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Tertiary,
            0,
            0,
            0,
            true,
            false,
            false,
        )
    }

    /// Create benzylic halide
    pub fn benzylic_halide() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Halide,
            SubstitutionDegree::Benzylic,
            0,
            0,
            0,
            true,
            true,
            false,
        )
    }

    /// Create carboxylic acid
    pub fn carboxylic_acid() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::CarboxylicAcid,
            SubstitutionDegree::Primary,
            0,
            0,
            0,
            false,
            false,
            false,
        )
    }

    pub fn silyl_ether() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::SilylEther,
            SubstitutionDegree::Primary,
            0,
            0,
            1,
            false,
            false,
            false,
        )
    }

    pub fn acetal() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::Acetal,
            SubstitutionDegree::Secondary,
            0,
            1,
            1,
            false,
            false,
            false,
        )
    }

    pub fn boc_carbamate() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::BocCarbamate,
            SubstitutionDegree::Secondary,
            0,
            1,
            2,
            false,
            true,
            false,
        )
    }

    pub fn cbz_carbamate() -> SiteDescriptor {
        Self::build(
            ReactiveSiteKind::CbzCarbamate,
            SubstitutionDegree::Secondary,
            0,
            1,
            2,
            false,
            true,
            true,
        )
    }

    // Integration methods for creating descriptors from actual site data

    /// Create descriptor from alcohol site
    pub(crate) fn from_alcohol_site(
        site: &crate::chemistry::organic::centers::AlcoholSite,
    ) -> SiteDescriptor {
        descriptor_from_carbon(
            site.participant.structure,
            ReactiveSiteKind::Alcohol,
            site.carbon,
            site.degree,
        )
    }

    /// Create descriptor from carbonyl site
    pub(crate) fn from_carbonyl_site(
        site: &crate::chemistry::organic::centers::CarbonylSite,
    ) -> SiteDescriptor {
        let mut descriptor = descriptor_from_carbon(
            site.participant.structure,
            if site.is_ketone {
                ReactiveSiteKind::Ketone
            } else {
                ReactiveSiteKind::Aldehyde
            },
            site.carbon,
            if site.is_ketone { 2 } else { 1 },
        );
        if carbonyl_is_conjugated_to_aromatic(site.participant.structure, site.carbon, site.oxygen)
        {
            descriptor.electronics.aromatic = true;
            descriptor.electronics.resonance_stabilization = true;
        }
        descriptor
    }

    /// Create descriptor from carboxylic acid site
    pub(crate) fn from_carboxylic_acid_site(
        site: &crate::chemistry::organic::centers::CarboxylicAcidSite,
    ) -> SiteDescriptor {
        descriptor_from_carbon(
            site.participant.structure,
            ReactiveSiteKind::CarboxylicAcid,
            site.carbon,
            1,
        )
    }

    pub(crate) fn from_unsaturated_bond_site(
        site: &crate::chemistry::organic::centers::UnsaturatedBondSite,
    ) -> SiteDescriptor {
        descriptor_from_carbon(
            site.participant.structure,
            site.participant.site.kind.clone(),
            site.high_degree_carbon,
            site.participant
                .structure
                .carbon_degree(site.high_degree_carbon),
        )
    }

    /// Create descriptor from halide site
    pub(crate) fn from_halide_site(
        site: &crate::chemistry::organic::centers::HalideSite,
    ) -> SiteDescriptor {
        descriptor_from_carbon(
            site.participant.structure,
            ReactiveSiteKind::Halide,
            site.carbon,
            site.degree,
        )
    }

    pub(crate) fn from_amine_site(
        site: &crate::chemistry::organic::centers::AmineSite,
    ) -> SiteDescriptor {
        let degree = match site.participant.structure.carbon_degree(site.nitrogen) {
            0 | 1 => SubstitutionDegree::Primary,
            2 => SubstitutionDegree::Secondary,
            _ => SubstitutionDegree::Tertiary,
        };
        SiteDescriptor::new(
            site.participant.site.kind.clone(),
            degree,
            electronic_environment(site.participant.structure, site.nitrogen),
            bulky_substituent_count(site.participant.structure, site.nitrogen),
            has_beta_hydrogen(site.participant.structure, site.nitrogen),
        )
    }

    pub(crate) fn from_phosphonium_salt_site(
        site: &crate::chemistry::organic::centers::PhosphoniumSaltSite,
    ) -> SiteDescriptor {
        descriptor_from_carbon(
            site.participant.structure,
            ReactiveSiteKind::PhosphoniumSalt,
            site.alpha_carbon,
            site.participant.structure.carbon_degree(site.alpha_carbon),
        )
    }

    pub(crate) fn from_thiol_site(
        site: &crate::chemistry::organic::centers::ThiolSite,
    ) -> SiteDescriptor {
        let degree = match site.participant.structure.carbon_degree(site.sulfur) {
            0 | 1 => SubstitutionDegree::Primary,
            2 => SubstitutionDegree::Secondary,
            _ => SubstitutionDegree::Tertiary,
        };
        SiteDescriptor::new(
            ReactiveSiteKind::Thiol,
            degree,
            electronic_environment(site.participant.structure, site.sulfur),
            bulky_substituent_count(site.participant.structure, site.sulfur),
            has_beta_hydrogen(site.participant.structure, site.sulfur),
        )
    }

    pub(crate) fn from_phosphorus_ylide_site(
        site: &crate::chemistry::organic::centers::PhosphorusYlideSite,
    ) -> SiteDescriptor {
        descriptor_from_carbon(
            site.participant.structure,
            ReactiveSiteKind::PhosphorusYlide,
            site.alpha_carbon,
            site.participant.structure.carbon_degree(site.alpha_carbon),
        )
    }

    pub(crate) fn from_phosphonate_carbanion_site(
        site: &crate::chemistry::organic::centers::PhosphonateCarbanionSite,
    ) -> SiteDescriptor {
        let mut descriptor = descriptor_from_carbon(
            site.participant.structure,
            ReactiveSiteKind::PhosphonateCarbanion,
            site.alpha_carbon,
            site.participant.structure.carbon_degree(site.alpha_carbon),
        );
        descriptor.electronics.electron_withdrawing_groups += 1;
        descriptor.electronics.resonance_stabilization = true;
        descriptor
    }

    pub(crate) fn from_sulfone_carbanion_site(
        site: &crate::chemistry::organic::centers::SulfoneCarbanionSite,
    ) -> SiteDescriptor {
        let mut descriptor = descriptor_from_carbon(
            site.participant.structure,
            ReactiveSiteKind::SulfoneCarbanion,
            site.alpha_carbon,
            site.participant.structure.carbon_degree(site.alpha_carbon),
        );
        descriptor.electronics.electron_withdrawing_groups += 2;
        descriptor.electronics.resonance_stabilization = true;
        descriptor
    }
}

fn missing_secondary_site(profile: &SelectivityProfile) -> SelectivityRuntimeEffect {
    SelectivityRuntimeEffect {
        rate_multiplier: 0.0,
        activation_delta_kj_per_mol: 50.0,
        pre_exp_multiplier: 1.0,
        suppressed: true,
        reason: format!(
            "{:?} selectivity profile is missing a secondary site",
            profile.mechanism
        ),
    }
}

fn evaluate_alpha_carbon_profile(
    profile: &SelectivityProfile,
    context: &SelectivityContext,
) -> ReactivityScore {
    let mut value = profile.primary_site.steric_accessibility();
    let mut activation_delta = 0.0;
    let mut reason = format!("{:?} alpha-carbon selectivity", profile.mechanism);

    if profile.primary_site.electronics.electron_withdrawing_groups >= 2 {
        value *= 2.0;
        activation_delta -= 4.0;
        reason.push_str("; activated alpha carbon");
    } else if profile.primary_site.electronics.electron_withdrawing_groups == 0 {
        value *= 0.6;
        activation_delta += 3.0;
        reason.push_str("; weakly acidic alpha carbon");
    }

    match profile.mechanism {
        ReactionType::AlphaHalogenation => {
            if context.is_acidic() || context.is_basic() {
                value *= 1.4;
                activation_delta -= 2.0;
            } else {
                value *= 0.4;
                activation_delta += 8.0;
            }
        }
        ReactionType::AldolAddition
        | ReactionType::EnolateAlkylation
        | ReactionType::MichaelAddition
        | ReactionType::ClaisenCondensation => {
            if context.is_basic() {
                value *= 1.8;
                activation_delta -= 4.0;
            } else {
                value *= 0.35;
                activation_delta += 10.0;
            }
            if context.is_water_rich() {
                value *= 0.65;
                activation_delta += 2.0;
                reason.push_str("; water competes with enolate chemistry");
            }
        }
        ReactionType::AldolDehydration => {
            if context.is_acidic() || context.is_high_temperature() {
                value *= 1.5;
                activation_delta -= 3.0;
            } else {
                value *= 0.3;
                activation_delta += 12.0;
            }
        }
        ReactionType::EnamineFormation => {
            if context.is_acidic() {
                value *= 1.2;
                activation_delta -= 2.0;
            }
            if context.is_basic() {
                value *= 0.4;
                activation_delta += 8.0;
            }
        }
        _ => {}
    }

    ReactivityScore::with_activation_delta(activation_delta, reason)
        .with_pre_exp_multiplier(value.max(0.05))
}

fn evaluate_protecting_group_profile(
    profile: &SelectivityProfile,
    context: &SelectivityContext,
) -> ReactivityScore {
    let mut value: f64 = 1.0;
    let mut activation_delta = 0.0;
    let mut reasons = Vec::new();

    match profile.mechanism {
        ReactionType::SilylEtherFormation => {
            if context.is_water_poor() {
                value *= 2.0;
                activation_delta -= 4.0;
                reasons.push("dry medium favors silyl ether formation");
            } else {
                value *= 0.15;
                activation_delta += 12.0;
                reasons.push("water suppresses silyl ether formation");
            }
            if context.is_basic() {
                value *= 1.4;
                activation_delta -= 2.0;
                reasons.push("basic medium scavenges acid during silylation");
            }
            if context.is_acidic() {
                value *= 0.25;
                activation_delta += 8.0;
                reasons.push("acidic medium destabilizes silyl ether formation");
            }
        }
        ReactionType::SilylEtherCleavage => {
            if context.has_fluoride() {
                value *= 4.0;
                activation_delta -= 10.0;
                reasons.push("fluoride accelerates silyl ether cleavage");
            } else {
                value *= 0.02;
                activation_delta += 22.0;
                reasons.push("silyl ether cleavage lacks fluoride");
            }
            if context.is_water_rich() || context.is_acidic() {
                value *= 1.25;
                activation_delta -= 1.5;
                reasons.push("protic conditions help proton transfer after cleavage");
            }
        }
        ReactionType::AcetalFormation => {
            if context.is_acidic() {
                value *= 2.2;
                activation_delta -= 5.0;
                reasons.push("acid catalyzes acetal formation");
            } else {
                value *= 0.08;
                activation_delta += 16.0;
                reasons.push("acetal formation lacks acid catalysis");
            }
            if context.is_water_poor() {
                value *= 2.0;
                activation_delta -= 4.0;
                reasons.push("dry medium shifts acetal equilibrium toward product");
            } else {
                value *= 0.12;
                activation_delta += 14.0;
                reasons.push("water suppresses acetal formation");
            }
        }
        ReactionType::AcetalHydrolysis => {
            if context.is_acidic() {
                value *= 2.0;
                activation_delta -= 4.0;
                reasons.push("acid catalyzes acetal hydrolysis");
            } else {
                value *= 0.05;
                activation_delta += 18.0;
                reasons.push("acetal hydrolysis lacks acid catalysis");
            }
            if context.is_water_rich() {
                value *= 2.5;
                activation_delta -= 5.0;
                reasons.push("water-rich medium drives acetal hydrolysis");
            } else {
                value *= 0.15;
                activation_delta += 12.0;
                reasons.push("dry medium suppresses acetal hydrolysis");
            }
        }
        ReactionType::CarbamateFormation => {
            if context.is_basic() {
                value *= 1.8;
                activation_delta -= 4.0;
                reasons.push("basic medium favors carbamate protection");
            } else {
                value *= 0.4;
                activation_delta += 6.0;
                reasons.push("carbamate protection lacks basic acid scavenging");
            }
            if context.is_water_poor() {
                value *= 1.5;
                activation_delta -= 2.5;
                reasons.push("water-poor medium favors carbamate formation");
            } else {
                value *= 0.35;
                activation_delta += 7.0;
                reasons.push("water competes with carbamate formation");
            }
        }
        ReactionType::CarbamateCleavage => {
            let acid_path = context.is_acidic() && context.is_water_rich();
            let hydrogenolysis_path = context.has_hydrogen()
                && context.palladium_available
                && context.has_available_surface();
            if acid_path {
                value *= 2.5;
                activation_delta -= 7.0;
                reasons.push("acidic water-rich medium cleaves acid-labile carbamates");
            }
            if hydrogenolysis_path {
                value *= 3.0;
                activation_delta -= 9.0;
                reasons.push("hydrogen and palladium enable carbamate hydrogenolysis");
            }
            if !acid_path && !hydrogenolysis_path {
                value *= 0.03;
                activation_delta += 24.0;
                reasons
                    .push("carbamate cleavage lacks acid hydrolysis or hydrogenolysis conditions");
            }
        }
        ReactionType::EsterProtection => {
            if context.is_acidic() {
                value *= 1.5;
                activation_delta -= 3.0;
                reasons.push("acid catalyzes ester protection");
            }
            if context.is_water_poor() {
                value *= 1.8;
                activation_delta -= 3.5;
                reasons.push("dry medium shifts esterification toward ester");
            } else {
                value *= 0.25;
                activation_delta += 9.0;
                reasons.push("water suppresses ester protection");
            }
        }
        ReactionType::EsterHydrolysis => {
            if context.is_acidic() || context.is_basic() {
                value *= 1.8;
                activation_delta -= 4.0;
                reasons.push("acid or base catalyzes ester hydrolysis");
            } else {
                value *= 0.2;
                activation_delta += 10.0;
                reasons.push("ester hydrolysis lacks acid or base catalysis");
            }
            if context.is_water_rich() {
                value *= 1.8;
                activation_delta -= 3.5;
                reasons.push("water-rich medium favors ester hydrolysis");
            } else {
                value *= 0.25;
                activation_delta += 8.0;
                reasons.push("dry medium suppresses ester hydrolysis");
            }
        }
        _ => {}
    }

    ReactivityScore::with_activation_delta(activation_delta, reasons.join("; "))
        .with_pre_exp_multiplier(value.max(0.01))
}

/// Medium/condition selectivity for ring-closing condensations. The ring-size
/// strain term (Baldwin's rules) is applied separately by the generator via
/// `ring_closure_activation_penalty_kj_per_mol`, because it depends on the
/// concrete atoms being bonded rather than on the reaction medium.
fn evaluate_cyclization_profile(
    profile: &SelectivityProfile,
    context: &SelectivityContext,
) -> ReactivityScore {
    let mut value: f64 = 1.0;
    let mut activation_delta = 0.0;
    let mut reasons = Vec::new();

    match profile.mechanism {
        ReactionType::Lactonization => {
            if context.is_acidic() {
                value *= 1.6;
                activation_delta -= 3.0;
                reasons.push("acid catalyzes lactonization");
            }
            if context.is_water_poor() {
                value *= 1.8;
                activation_delta -= 3.5;
                reasons.push("dry medium drives the dehydrative ring closure");
            } else {
                value *= 0.3;
                activation_delta += 9.0;
                reasons.push("water reverses lactonization");
            }
        }
        ReactionType::Lactamization => {
            if context.is_high_temperature() {
                value *= 1.5;
                activation_delta -= 3.0;
                reasons.push("heat drives amide ring closure");
            }
            if context.is_water_poor() {
                value *= 1.6;
                activation_delta -= 3.0;
                reasons.push("dry medium favors lactam formation");
            } else {
                value *= 0.4;
                activation_delta += 7.0;
                reasons.push("water competes with lactam formation");
            }
        }
        ReactionType::HeterocycleCondensation => {
            if context.is_acidic() {
                value *= 2.0;
                activation_delta -= 4.0;
                reasons.push("acid catalyzes the heterocyclic condensation");
            } else {
                value *= 0.2;
                activation_delta += 12.0;
                reasons.push("heterocyclic condensation lacks acid catalysis");
            }
            if context.is_water_poor() {
                value *= 1.8;
                activation_delta -= 3.5;
                reasons.push("dry medium drives the dehydrative aromatization");
            } else {
                value *= 0.4;
                activation_delta += 6.0;
                reasons.push("water suppresses the dehydrative aromatization");
            }
        }
        _ => {}
    }

    ReactivityScore::with_activation_delta(activation_delta, reasons.join("; "))
        .with_pre_exp_multiplier(value.max(0.01))
}

fn descriptor_from_carbon(
    structure: &MolecularStructure,
    kind: ReactiveSiteKind,
    carbon: usize,
    carbon_degree: usize,
) -> SiteDescriptor {
    let degree = if is_benzylic(structure, carbon) {
        SubstitutionDegree::Benzylic
    } else if is_allylic(structure, carbon) {
        SubstitutionDegree::Allylic
    } else if carbon_degree >= 3 {
        SubstitutionDegree::Tertiary
    } else if carbon_degree == 2 {
        SubstitutionDegree::Secondary
    } else {
        SubstitutionDegree::Primary
    };
    let electronics = electronic_environment(structure, carbon);
    let bulky = bulky_substituent_count(structure, carbon);
    SiteDescriptor::new(
        kind,
        degree,
        electronics,
        bulky,
        has_beta_hydrogen(structure, carbon),
    )
}

fn electronic_environment(structure: &MolecularStructure, center: usize) -> ElectronicEnvironment {
    let mut env = ElectronicEnvironment::default();
    env.aromatic = is_aromatic_carbon(structure, center);
    env.resonance_stabilization =
        env.aromatic || is_benzylic(structure, center) || is_allylic(structure, center);
    let mut visited = std::collections::BTreeSet::new();
    let mut queue = std::collections::VecDeque::new();
    visited.insert(center);
    queue.push_back((center, 0usize));
    while let Some((atom, distance)) = queue.pop_front() {
        if distance >= 3 {
            continue;
        }
        for (neighbor, order) in structure.neighbors(atom) {
            if !visited.insert(neighbor) {
                continue;
            }
            if neighbor != center {
                match electronic_atom_effect(structure, neighbor, order) {
                    ElectronicAtomEffect::Donating => env.electron_donating_groups += 1,
                    ElectronicAtomEffect::Withdrawing => env.electron_withdrawing_groups += 1,
                    ElectronicAtomEffect::Neutral => {}
                }
            }
            queue.push_back((neighbor, distance + 1));
        }
    }
    env
}

enum ElectronicAtomEffect {
    Donating,
    Withdrawing,
    Neutral,
}

fn electronic_atom_effect(
    structure: &MolecularStructure,
    atom: usize,
    bond_order: f64,
) -> ElectronicAtomEffect {
    let element = structure.atoms[atom].element.as_str();
    if structure.atoms[atom].charge > 0.1 {
        return ElectronicAtomEffect::Withdrawing;
    }
    match element {
        "O" | "N" | "S" if bond_order_matches(bond_order, 1.0) => ElectronicAtomEffect::Donating,
        "O" | "N" | "F" | "Cl" | "Br" | "I" => ElectronicAtomEffect::Withdrawing,
        "C" if has_double_bonded_heteroatom(structure, atom)
            || has_triple_bonded_nitrogen(structure, atom) =>
        {
            ElectronicAtomEffect::Withdrawing
        }
        "C" => ElectronicAtomEffect::Donating,
        _ => ElectronicAtomEffect::Neutral,
    }
}

fn bulky_substituent_count(structure: &MolecularStructure, center: usize) -> u32 {
    structure
        .neighbors(center)
        .into_iter()
        .filter(|(neighbor, _)| {
            let atom = &structure.atoms[*neighbor];
            matches!(atom.element.as_str(), "Br" | "I")
                || (atom.element == "C"
                    && structure
                        .neighbors(*neighbor)
                        .iter()
                        .filter(|(n, _)| *n != center && structure.atoms[*n].element != "H")
                        .count()
                        >= 2)
        })
        .count() as u32
}

fn has_beta_hydrogen(structure: &MolecularStructure, center: usize) -> bool {
    structure
        .neighbors(center)
        .into_iter()
        .any(|(neighbor, order)| {
            neighbor != center
                && bond_order_matches(order, 1.0)
                && structure.atoms[neighbor].element == "C"
                && structure.neighbors(neighbor).iter().any(|(n, h_order)| {
                    structure.atoms[*n].element == "H" && bond_order_matches(*h_order, 1.0)
                })
        })
}

fn is_benzylic(structure: &MolecularStructure, atom: usize) -> bool {
    structure
        .neighbors(atom)
        .into_iter()
        .any(|(neighbor, order)| {
            bond_order_matches(order, 1.0) && is_aromatic_carbon(structure, neighbor)
        })
}

fn is_allylic(structure: &MolecularStructure, atom: usize) -> bool {
    structure
        .neighbors(atom)
        .into_iter()
        .any(|(neighbor, order)| {
            bond_order_matches(order, 1.0)
                && structure
                    .neighbors(neighbor)
                    .into_iter()
                    .any(|(other, other_order)| {
                        other != atom
                            && structure.atoms[other].element == "C"
                            && bond_order_matches(other_order, 2.0)
                    })
        })
}

fn is_aromatic_carbon(structure: &MolecularStructure, atom: usize) -> bool {
    structure.atoms.get(atom).is_some_and(|a| a.element == "C")
        && structure
            .neighbors(atom)
            .iter()
            .filter(|(_, order)| bond_order_matches(*order, 1.5))
            .count()
            >= 2
}

fn carbonyl_is_conjugated_to_aromatic(
    structure: &MolecularStructure,
    carbon: usize,
    oxygen: usize,
) -> bool {
    structure
        .neighbors(carbon)
        .into_iter()
        .any(|(neighbor, order)| {
            neighbor != oxygen
                && bond_order_matches(order, 1.0)
                && is_aromatic_carbon(structure, neighbor)
        })
}

fn has_double_bonded_heteroatom(structure: &MolecularStructure, atom: usize) -> bool {
    structure.neighbors(atom).iter().any(|(neighbor, order)| {
        matches!(structure.atoms[*neighbor].element.as_str(), "O" | "N" | "S")
            && bond_order_matches(*order, 2.0)
    })
}

fn has_triple_bonded_nitrogen(structure: &MolecularStructure, atom: usize) -> bool {
    structure.neighbors(atom).iter().any(|(neighbor, order)| {
        structure.atoms[*neighbor].element == "N" && bond_order_matches(*order, 3.0)
    })
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
        assert!(matches!(
            result.recommendation,
            SelectivityRecommendation::Exclusive | SelectivityRecommendation::Preferred
        ));
    }

    #[test]
    fn selectivity_engine_esterification_suppressed_for_tertiary() {
        let acid = SiteDescriptorBuilder::carboxylic_acid();
        let alcohol = SiteDescriptorBuilder::tertiary_alcohol();

        let context = SelectivityContext {
            solvent_type: SolventType::Acidic,
            temperature: 380.0,
            ph: Some(2.0),
            ..Default::default()
        };

        let result = SelectivityEngine::fischer_esterification(&acid, &alcohol, &context);

        assert!(matches!(
            result.recommendation,
            SelectivityRecommendation::Suppressed
                | SelectivityRecommendation::None
                | SelectivityRecommendation::Mixed
        ));
    }

    #[test]
    fn aldehyde_vs_ketone_selectivity() {
        let aldehyde = SiteDescriptorBuilder::aldehyde();
        let ketone = SiteDescriptorBuilder::ketone();

        let context = SelectivityContext::default();

        assert!(SelectivityEngine::aldehyde_preferred_over_ketone(
            &aldehyde,
            &ketone,
            NucleophileStrength::Moderate,
            &context
        ));
    }

    #[test]
    fn strongly_basic_carbonyl_addition_uses_inorganic_mixture_state() {
        let profile = SelectivityProfile::new(
            ReactionType::CarbonylAddition,
            SiteDescriptorBuilder::ketone(),
        )
        .with_nucleophile_strength(NucleophileStrength::VeryStrong)
        .never_suppress();
        let dry = SelectivityContext {
            water_activity: 0.0,
            oxygen_activity: 0.0,
            oxidizing_strength: 0.0,
            free_complexable_metal_activity: 0.0,
            ..Default::default()
        };
        let wet_metal_oxidizing = SelectivityContext {
            water_activity: 0.8,
            oxygen_activity: 0.2,
            oxidizing_strength: 0.3,
            free_complexable_metal_activity: 0.05,
            ..Default::default()
        };

        let dry_effect = SelectivityEngine::evaluate_profile(&profile, &dry);
        let wet_effect = SelectivityEngine::evaluate_profile(&profile, &wet_metal_oxidizing);

        assert!(dry_effect.rate_multiplier > wet_effect.rate_multiplier);
        assert!(dry_effect.activation_delta_kj_per_mol < wet_effect.activation_delta_kj_per_mol);
    }

    #[test]
    fn protecting_group_profiles_respond_to_water_acid_and_fluoride() {
        let silylation = SelectivityProfile::new(
            ReactionType::SilylEtherFormation,
            SiteDescriptorBuilder::primary_alcohol(),
        )
        .never_suppress();
        let dry_basic = SelectivityContext {
            solvent_type: SolventType::Basic,
            ph: Some(10.0),
            water_activity: 0.02,
            ..Default::default()
        };
        let wet_acidic = SelectivityContext {
            solvent_type: SolventType::Acidic,
            ph: Some(2.0),
            water_activity: 0.8,
            ..Default::default()
        };
        let dry_effect = SelectivityEngine::evaluate_profile(&silylation, &dry_basic);
        let wet_effect = SelectivityEngine::evaluate_profile(&silylation, &wet_acidic);
        assert!(dry_effect.rate_multiplier > wet_effect.rate_multiplier);
        assert!(dry_effect.activation_delta_kj_per_mol < wet_effect.activation_delta_kj_per_mol);

        let cleavage = SelectivityProfile::new(
            ReactionType::SilylEtherCleavage,
            SiteDescriptorBuilder::silyl_ether(),
        )
        .never_suppress();
        let no_fluoride =
            SelectivityEngine::evaluate_profile(&cleavage, &SelectivityContext::default());
        let with_fluoride = SelectivityEngine::evaluate_profile(
            &cleavage,
            &SelectivityContext {
                fluoride_mol_per_bucket: 0.1,
                water_activity: 0.5,
                ..Default::default()
            },
        );
        assert!(with_fluoride.rate_multiplier > no_fluoride.rate_multiplier);
        assert!(
            with_fluoride.activation_delta_kj_per_mol < no_fluoride.activation_delta_kj_per_mol
        );
    }

    #[test]
    fn carbamate_cleavage_profile_sees_acid_hydrolysis_and_hydrogenolysis() {
        let profile = SelectivityProfile::new(
            ReactionType::CarbamateCleavage,
            SiteDescriptorBuilder::boc_carbamate(),
        )
        .never_suppress();
        let neutral = SelectivityEngine::evaluate_profile(&profile, &SelectivityContext::default());
        let acid_water = SelectivityEngine::evaluate_profile(
            &profile,
            &SelectivityContext {
                solvent_type: SolventType::Acidic,
                ph: Some(2.0),
                water_activity: 0.8,
                ..Default::default()
            },
        );
        let hydrogenolysis = SelectivityEngine::evaluate_profile(
            &profile,
            &SelectivityContext {
                hydrogen_mol_per_bucket: 0.5,
                palladium_available: true,
                available_surface_sites_mol_per_bucket: 0.1,
                ..Default::default()
            },
        );
        assert!(acid_water.rate_multiplier > neutral.rate_multiplier);
        assert!(hydrogenolysis.rate_multiplier > neutral.rate_multiplier);
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

    #[test]
    fn graph_descriptor_detects_benzylic_and_beta_hydrogen() {
        let structure = crate::chemistry::frowns::parse_frowns(
            "destroy:graph:atoms=C.C.C.C.C.C.C.Cl.H.H.H.H.H.H.H;\
             bonds=0-a-1,1-a-2,2-a-3,3-a-4,4-a-5,5-a-0,\
                   0-s-6,6-s-7,6-s-8,6-s-9,1-s-10,2-s-11,3-s-12,4-s-13,5-s-14",
        )
        .unwrap();
        let descriptor = super::descriptor_from_carbon(
            &structure,
            ReactiveSiteKind::Halide,
            6,
            structure.carbon_degree(6),
        );

        assert_eq!(descriptor.degree, SubstitutionDegree::Benzylic);
        assert!(descriptor.electronics.resonance_stabilization);

        let alkyl_chloride = crate::chemistry::frowns::parse_frowns("CCCCl").unwrap();
        let alkyl_descriptor = super::descriptor_from_carbon(
            &alkyl_chloride,
            ReactiveSiteKind::Halide,
            2,
            alkyl_chloride.carbon_degree(2),
        );
        assert!(alkyl_descriptor.has_beta_hydrogen);
    }
}
