//! Central selectivity engine for evaluating reaction favorability
//!
//! Coordinates all selectivity rules and provides unified API for generators.

use super::carbonyl_addition::evaluate_carbonyl_addition;
use super::elimination::{evaluate_e1, evaluate_e2};
use super::esterification::evaluate_fischer_esterification;
use super::nucleophilic_substitution::evaluate_sn2;
use super::types::*;
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
            ReactionType::CarbonylAddition => evaluate_carbonyl_addition(
                &profile.primary_site,
                profile
                    .nucleophile_strength
                    .unwrap_or(NucleophileStrength::Moderate),
                context,
            ),
            ReactionType::ElectrophilicAddition => ReactivityScore::new(
                1.0,
                "electrophilic addition has no specialized selectivity profile yet",
            ),
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
