use super::super::centers::{
    HalideSite, NucleophilicPhosphorusSite, PhosphineSite, PhosphonateCarbanionSite,
    PhosphoniumSaltSite, PhosphorusYlideSite, SulfoneCarbanionSite, YlideStability,
};
use super::super::resolver::DerivedSubstanceResolver;
use super::addition::{expand_stereo_product_distribution_with_parameters, StereoProductVariant};
use super::common::*;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::kinetics::ReactionChannel;
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor};
use crate::chemistry::reaction::{Reaction, StoichiometricTerm};
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityContext, SelectivityProfile},
    NucleophileStrength,
};

pub(crate) fn generate_phosphonium_salt_formation(
    halide_site: &super::super::centers::HalideSite<'_>,
    phosphine_site: &PhosphineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    if halide_site.degree >= 3 {
        return Ok(None);
    }
    let halide_desc = SiteDescriptorBuilder::from_halide_site(halide_site);
    let substance = halide_site.participant.substance;
    let halide_structure = halide_site.participant.structure;
    let phosphine = phosphine_site.participant.substance;
    let phosphine_structure = phosphine_site.participant.structure;
    let carbon = halide_site.carbon;
    let halogen = halide_site.halogen;
    let halide_ion = halide_ion(
        halide_structure,
        halogen,
        "phosphonium_salt_formation",
        &halide_site.participant,
    )?;

    let mut halide_editor = MolecularEditor::new(halide_structure);
    let mapping = halide_editor.remove_atoms(&[halogen])?;
    let halide_carbon = mapped_atom(&mapping, carbon, "halide carbon")?;
    let halide_fragment = halide_editor.finish()?;

    let mut phosphine_editor = MolecularEditor::new(phosphine_structure);
    let phosphorus = phosphine_site.phosphorus;
    phosphine_editor.add_group(phosphorus, &halide_fragment, halide_carbon, 1.0)?;
    phosphine_editor.replace_atom(phosphorus, "P", 1.0)?;
    let product = resolver.resolve(phosphine_editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "phosphonium_salt_formation",
            &halide_site.participant,
            &phosphine_site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant(phosphine.id.clone(), 1, 1)
        .product(product, 1)
        .product(halide_ion, 1)
        .activation_energy_kj_per_mol(22.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::PhosphoniumSaltFormation, halide_desc)
                .with_nucleophile_strength(NucleophileStrength::VeryStrong),
        )
        .build(),
    ))
}

pub(crate) fn generate_phosphonium_ylide_formation(
    site: &PhosphoniumSaltSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let alpha_hydrogen = match site.alpha_hydrogens.first().copied() {
        Some(index) => index,
        None => return Ok(None),
    };
    let phosphorus_desc = SiteDescriptorBuilder::from_phosphonium_salt_site(site);
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[alpha_hydrogen])?;
    let phosphorus = mapped_atom(&mapping, site.phosphorus, "phosphonium phosphorus")?;
    let alpha_carbon = mapped_atom(&mapping, site.alpha_carbon, "phosphonium alpha carbon")?;
    editor.replace_atom(phosphorus, "P", 1.0)?;
    editor.replace_atom(alpha_carbon, "C", -1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "phosphonium_ylide_formation",
            &site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant("destroy:ethoxide", 1, 1)
        .product(product, 1)
        .product("destroy:ethanol", 1)
        .activation_energy_kj_per_mol(28.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::PhosphoniumYlideFormation, phosphorus_desc)
                .with_nucleophile_strength(NucleophileStrength::Strong),
        )
        .build(),
    ))
}

pub(crate) fn generate_wittig_olefination(
    ylide_site: &PhosphorusYlideSite<'_>,
    carbonyl_site: &super::super::centers::CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(carbonyl_site);
    let ylide_desc = SiteDescriptorBuilder::from_phosphorus_ylide_site(ylide_site);
    let carbonyl = carbonyl_site.participant.substance;
    let carbonyl_structure = carbonyl_site.participant.structure;
    let ylide = ylide_site.participant.substance;
    let ylide_structure = ylide_site.participant.structure;
    let carbonyl_carbon = carbonyl_site.carbon;
    let carbonyl_oxygen = carbonyl_site.oxygen;
    let phosphorus = ylide_site.phosphorus;
    let alpha_carbon = ylide_site.alpha_carbon;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_oxygen])?;
    let carbonyl_carbon = mapped_atom(&carbonyl_mapping, carbonyl_carbon, "carbonyl carbon")?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let (first_fragment, first_map, second_fragment, second_map) =
        MolecularEditor::split_at_bond(ylide_structure, phosphorus, alpha_carbon)?;
    let (p_fragment, p_map, alpha_fragment, alpha_map) = if first_map[phosphorus].is_some() {
        (first_fragment, first_map, second_fragment, second_map)
    } else {
        (second_fragment, second_map, first_fragment, first_map)
    };
    let p_index = mapped_atom(&p_map, phosphorus, "ylide phosphorus")?;
    let alpha_index = mapped_atom(&alpha_map, alpha_carbon, "ylide alpha carbon")?;

    let mut p_editor = MolecularEditor::new(&p_fragment);
    p_editor.add_atom(p_index, "O", 0.0, 2.0)?;
    p_editor.replace_atom(p_index, "P", 0.0)?;
    let phosphine_oxide = resolver.resolve(p_editor.finish()?)?;

    let mut alpha_editor = MolecularEditor::new(&alpha_fragment);
    alpha_editor.replace_atom(alpha_index, "C", 0.0)?;
    let alpha_fragment = alpha_editor.finish()?;

    let alkene_structure = MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &alpha_fragment,
        alpha_index,
        2.0,
    )?;
    let joined_alpha_carbon = carbonyl_fragment.atom_count() + alpha_index;
    let mut alkene_editor = MolecularEditor::new(&alkene_structure);
    let has_stereo = alkene_editor
        .mark_double_bond_stereo_mixture_if_valid(carbonyl_carbon, joined_alpha_carbon)?;
    let mut variants = if has_stereo {
        expand_stereo_product_distribution_with_parameters(
            alkene_editor.finish()?,
            "wittig".to_string(),
            0.0,
            1.0,
        )?
    } else {
        vec![StereoProductVariant {
            structure: alkene_editor.finish()?,
            channel_suffix: "single".to_string(),
            activation_delta_kj_per_mol: 0.0,
            pre_exponential_factor_multiplier: 1.0,
        }]
    };

    if variants.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_pair_site_reaction_id(
                "wittig_olefination",
                &ylide_site.participant,
                &carbonyl_site.participant,
            ),
            reason: "wittig olefination produced no product variants".to_string(),
        });
    }

    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "wittig_olefination",
        &ylide_site.participant,
        &carbonyl_site.participant,
    ))
    .reactant(ylide.id.clone(), 1, 1)
    .reactant(carbonyl.id.clone(), 1, 1)
    .activation_energy_kj_per_mol(24.0);

    if variants.len() == 1 {
        builder = builder
            .product(resolver.resolve(variants.remove(0).structure)?, 1)
            .product(phosphine_oxide, 1)
            .selectivity_profile(
                SelectivityProfile::new(ReactionType::WittigOlefination, carbonyl_desc)
                    .with_secondary_site(ylide_desc)
                    .with_nucleophile_strength(match ylide_site.stability {
                        YlideStability::Unstabilized => NucleophileStrength::VeryStrong,
                        YlideStability::SemiStabilized => NucleophileStrength::Strong,
                        YlideStability::Stabilized => NucleophileStrength::Strong,
                    }),
            );
        return Ok(Some(builder.build()));
    }

    for variant in variants {
        let mut activation_delta = variant.activation_delta_kj_per_mol;
        let mut pre_exp = 10_000.0 * variant.pre_exponential_factor_multiplier;
        let suffix = variant.channel_suffix.as_str();
        match (ylide_site.stability, suffix) {
            (YlideStability::Unstabilized, suffix) if suffix.ends_with("db_z") => {
                activation_delta -= 2.0;
                pre_exp *= 1.25;
            }
            (YlideStability::Unstabilized, _) => {
                activation_delta += 1.0;
                pre_exp *= 0.85;
            }
            (YlideStability::SemiStabilized, suffix) if suffix.ends_with("db_z") => {
                activation_delta -= 0.75;
                pre_exp *= 1.1;
            }
            (YlideStability::SemiStabilized, _) => {
                activation_delta += 0.5;
                pre_exp *= 0.95;
            }
            (YlideStability::Stabilized, suffix) if suffix.ends_with("db_e") => {
                activation_delta -= 1.5;
                pre_exp *= 1.15;
            }
            (YlideStability::Stabilized, _) => {
                activation_delta += 1.0;
                pre_exp *= 0.9;
            }
        }
        builder = builder.channel(
            ReactionChannel::new(
                format!("wittig_olefination:{}", variant.channel_suffix),
                [
                    StoichiometricTerm::new(resolver.resolve(variant.structure)?, 1),
                    StoichiometricTerm::new(phosphine_oxide.clone(), 1),
                ],
                24.0 + activation_delta,
            )
            .with_pre_exponential_factor(pre_exp)
            .with_selectivity_profile(
                SelectivityProfile::new(ReactionType::WittigOlefination, carbonyl_desc.clone())
                    .with_secondary_site(ylide_desc.clone())
                    .with_nucleophile_strength(match ylide_site.stability {
                        YlideStability::Unstabilized => NucleophileStrength::VeryStrong,
                        YlideStability::SemiStabilized => NucleophileStrength::Strong,
                        YlideStability::Stabilized => NucleophileStrength::Strong,
                    }),
            ),
        );
    }

    Ok(Some(builder.build()))
}

pub(crate) fn generate_horner_wadsworth_emmons_olefination(
    phosphonate_site: &PhosphonateCarbanionSite<'_>,
    carbonyl_site: &super::super::centers::CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(carbonyl_site);
    let phosphonate_desc = SiteDescriptorBuilder::from_phosphonate_carbanion_site(phosphonate_site);
    let carbonyl = carbonyl_site.participant.substance;
    let phosphonate = phosphonate_site.participant.substance;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_site.participant.structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_site.oxygen])?;
    let carbonyl_carbon = mapped_atom(&carbonyl_mapping, carbonyl_site.carbon, "carbonyl carbon")?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let (first_fragment, first_map, second_fragment, second_map) = MolecularEditor::split_at_bond(
        phosphonate_site.participant.structure,
        phosphonate_site.phosphorus,
        phosphonate_site.alpha_carbon,
    )?;
    let (p_fragment, p_map, alpha_fragment, alpha_map) =
        if first_map[phosphonate_site.phosphorus].is_some() {
            (first_fragment, first_map, second_fragment, second_map)
        } else {
            (second_fragment, second_map, first_fragment, first_map)
        };
    let p_index = mapped_atom(
        &p_map,
        phosphonate_site.phosphorus,
        "phosphonate phosphorus",
    )?;
    let alpha_index = mapped_atom(
        &alpha_map,
        phosphonate_site.alpha_carbon,
        "phosphonate alpha carbon",
    )?;

    let mut phosphate_editor = MolecularEditor::new(&p_fragment);
    phosphate_editor.add_atom(p_index, "O", -1.0, 1.0)?;
    let phosphate = resolver.resolve(phosphate_editor.finish()?)?;

    let mut alpha_editor = MolecularEditor::new(&alpha_fragment);
    alpha_editor.replace_atom(alpha_index, "C", 0.0)?;
    let alpha_fragment = alpha_editor.finish()?;

    let alkene_structure = MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &alpha_fragment,
        alpha_index,
        2.0,
    )?;
    let joined_alpha_carbon = carbonyl_fragment.atom_count() + alpha_index;
    build_olefination_reaction(
        "horner_wadsworth_emmons_olefination",
        ReactionType::HornerWadsworthEmmonsOlefination,
        phosphonate.id.clone(),
        carbonyl.id.clone(),
        &phosphonate_site.participant,
        &carbonyl_site.participant,
        alkene_structure,
        carbonyl_carbon,
        joined_alpha_carbon,
        phosphate,
        carbonyl_desc,
        phosphonate_desc,
        EPreference::StrongE,
        resolver,
    )
}

pub(crate) fn generate_julia_olefination(
    sulfone_site: &SulfoneCarbanionSite<'_>,
    carbonyl_site: &super::super::centers::CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(carbonyl_site);
    let sulfone_desc = SiteDescriptorBuilder::from_sulfone_carbanion_site(sulfone_site);
    let carbonyl = carbonyl_site.participant.substance;
    let sulfone = sulfone_site.participant.substance;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_site.participant.structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_site.oxygen])?;
    let carbonyl_carbon = mapped_atom(&carbonyl_mapping, carbonyl_site.carbon, "carbonyl carbon")?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let (first_fragment, first_map, second_fragment, second_map) = MolecularEditor::split_at_bond(
        sulfone_site.participant.structure,
        sulfone_site.sulfur,
        sulfone_site.alpha_carbon,
    )?;
    let (s_fragment, s_map, alpha_fragment, alpha_map) = if first_map[sulfone_site.sulfur].is_some()
    {
        (first_fragment, first_map, second_fragment, second_map)
    } else {
        (second_fragment, second_map, first_fragment, first_map)
    };
    let s_index = mapped_atom(&s_map, sulfone_site.sulfur, "sulfone sulfur")?;
    let alpha_index = mapped_atom(
        &alpha_map,
        sulfone_site.alpha_carbon,
        "sulfone alpha carbon",
    )?;

    let mut sulfonate_editor = MolecularEditor::new(&s_fragment);
    sulfonate_editor.add_atom(s_index, "O", -1.0, 1.0)?;
    let sulfonate = resolver.resolve(sulfonate_editor.finish()?)?;

    let mut alpha_editor = MolecularEditor::new(&alpha_fragment);
    alpha_editor.replace_atom(alpha_index, "C", 0.0)?;
    let alpha_fragment = alpha_editor.finish()?;

    let alkene_structure = MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &alpha_fragment,
        alpha_index,
        2.0,
    )?;
    let joined_alpha_carbon = carbonyl_fragment.atom_count() + alpha_index;
    build_olefination_reaction(
        "julia_olefination",
        ReactionType::JuliaOlefination,
        sulfone.id.clone(),
        carbonyl.id.clone(),
        &sulfone_site.participant,
        &carbonyl_site.participant,
        alkene_structure,
        carbonyl_carbon,
        joined_alpha_carbon,
        sulfonate,
        carbonyl_desc,
        sulfone_desc,
        EPreference::ModerateE,
        resolver,
    )
}

#[derive(Clone, Copy)]
enum EPreference {
    ModerateE,
    StrongE,
}

#[allow(clippy::too_many_arguments)]
fn build_olefination_reaction(
    prefix: &str,
    reaction_type: ReactionType,
    reagent_id: crate::chemistry::substance::SubstanceId,
    carbonyl_id: crate::chemistry::substance::SubstanceId,
    reagent_participant: &super::super::space::SiteParticipant<'_>,
    carbonyl_participant: &super::super::space::SiteParticipant<'_>,
    alkene_structure: crate::chemistry::molecule::MolecularStructure,
    first_alkene_carbon: usize,
    second_alkene_carbon: usize,
    inorganic_product: crate::chemistry::substance::SubstanceId,
    carbonyl_desc: crate::chemistry::selectivity::types::SiteDescriptor,
    reagent_desc: crate::chemistry::selectivity::types::SiteDescriptor,
    preference: EPreference,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let mut alkene_editor = MolecularEditor::new(&alkene_structure);
    let has_stereo = alkene_editor
        .mark_double_bond_stereo_mixture_if_valid(first_alkene_carbon, second_alkene_carbon)?;
    let mut variants = if has_stereo {
        expand_stereo_product_distribution_with_parameters(
            alkene_editor.finish()?,
            prefix.to_string(),
            0.0,
            1.0,
        )?
    } else {
        vec![StereoProductVariant {
            structure: alkene_editor.finish()?,
            channel_suffix: "single".to_string(),
            activation_delta_kj_per_mol: 0.0,
            pre_exponential_factor_multiplier: 1.0,
        }]
    };

    if variants.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_pair_site_reaction_id(
                prefix,
                reagent_participant,
                carbonyl_participant,
            ),
            reason: "olefination produced no product variants".to_string(),
        });
    }

    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        prefix,
        reagent_participant,
        carbonyl_participant,
    ))
    .reactant(reagent_id, 1, 1)
    .reactant(carbonyl_id, 1, 1)
    .activation_energy_kj_per_mol(22.0);

    let profile = || {
        SelectivityProfile::new(reaction_type, carbonyl_desc.clone())
            .with_secondary_site(reagent_desc.clone())
            .with_nucleophile_strength(NucleophileStrength::Strong)
    };

    if variants.len() == 1 {
        builder = builder
            .product(resolver.resolve(variants.remove(0).structure)?, 1)
            .product(inorganic_product, 1)
            .selectivity_profile(profile());
        return Ok(Some(builder.build()));
    }

    for variant in variants {
        let suffix = variant.channel_suffix.as_str();
        let mut activation_delta = variant.activation_delta_kj_per_mol;
        let mut pre_exp = 10_000.0 * variant.pre_exponential_factor_multiplier;
        match (preference, suffix) {
            (EPreference::StrongE, suffix) if suffix.ends_with("db_e") => {
                activation_delta -= 2.0;
                pre_exp *= 1.25;
            }
            (EPreference::StrongE, _) => {
                activation_delta += 1.5;
                pre_exp *= 0.85;
            }
            (EPreference::ModerateE, suffix) if suffix.ends_with("db_e") => {
                activation_delta -= 1.0;
                pre_exp *= 1.1;
            }
            (EPreference::ModerateE, _) => {
                activation_delta += 0.75;
                pre_exp *= 0.9;
            }
        }
        builder = builder.channel(
            ReactionChannel::new(
                format!("{prefix}:{}", variant.channel_suffix),
                [
                    StoichiometricTerm::new(resolver.resolve(variant.structure)?, 1),
                    StoichiometricTerm::new(inorganic_product.clone(), 1),
                ],
                22.0 + activation_delta,
            )
            .with_pre_exponential_factor(pre_exp)
            .with_selectivity_profile(profile()),
        );
    }
    Ok(Some(builder.build()))
}

/// Nucleophilic substitution at phosphorus: a phosphine with a P-H bond
/// (PH3, R-PH2, or R2PH) attacks an alkyl halide, replacing one P-H with a
/// P-C bond and releasing HX. Mechanistically identical to amine alkylation
/// (SN2 at the halide carbon by the phosphorus lone pair), but the nucleophile
/// is phosphorus rather than nitrogen.
pub(crate) fn generate_nucleophilic_phosphorus_alkylation(
    halide_site: &HalideSite<'_>,
    phosphorus_site: &NucleophilicPhosphorusSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if halide_site.degree >= 3 {
        return Ok(None);
    }
    let halide_structure = halide_site.participant.structure;
    if halide_structure
        .neighbors(halide_site.carbon)
        .into_iter()
        .any(|(n, order)| {
            bond_order_matches(order, 2.0) && halide_structure.atoms[n].element == "C"
        })
    {
        return Ok(None);
    }

    let halide = halide_site.participant.substance;
    let phosphine = phosphorus_site.participant.substance;
    let phosphine_structure = phosphorus_site.participant.structure;
    let halide_carbon = halide_site.carbon;
    let halogen = halide_site.halogen;
    let phosphorus = phosphorus_site.phosphorus;
    let phosphorus_hydrogen = *phosphorus_site.hydrogens.first().ok_or_else(|| {
        ChemistryError::InvalidReaction {
            reaction_id: generated_pair_site_reaction_id(
                "nucleophilic_phosphorus_alkylation",
                &halide_site.participant,
                &phosphorus_site.participant,
            ),
            reason: "nucleophilic phosphorus has no P-H bond".to_string(),
        }
    })?;

    let mut halide_editor = MolecularEditor::new(halide_structure);
    let halide_mapping = halide_editor.remove_atoms(&[halogen])?;
    let halide_carbon = mapped_atom(&halide_mapping, halide_carbon, "halide carbon")?;
    let halide_fragment = halide_editor.finish()?;

    let mut phosphine_editor = MolecularEditor::new(phosphine_structure);
    let phosphine_mapping = phosphine_editor.remove_atoms(&[phosphorus_hydrogen])?;
    let phosphorus_atom =
        mapped_atom(&phosphine_mapping, phosphorus, "nucleophilic phosphorus")?;
    phosphine_editor.add_group(phosphorus_atom, &halide_fragment, halide_carbon, 1.0)?;
    let product = resolver.resolve(phosphine_editor.finish()?)?;

    let halide_desc = SiteDescriptorBuilder::from_halide_site(halide_site);
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "nucleophilic_phosphorus_alkylation",
            &halide_site.participant,
            &phosphorus_site.participant,
        ))
        .reactant(halide.id.clone(), 1, 1)
        .reactant(phosphine.id.clone(), 1, 1)
        .product(product, 1)
        .product(
            halide_ion(
                halide_structure,
                halogen,
                "nucleophilic_phosphorus_alkylation",
                &halide_site.participant,
            )?,
            1,
        )
        .product("destroy:proton", 1)
        .activation_energy_kj_per_mol(28.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::SN2, halide_desc)
                .with_nucleophile_strength(NucleophileStrength::VeryStrong),
        )
        .build(),
    ))
}

/// Hydrolysis of white phosphorus (P4) under basic conditions produces phosphine
/// (PH3). This is the industrial route and the missing entry point for phosphorus
/// into the organic reaction network.
///
/// Disproportionation: P4 + 6 H2O → PH3 + 3 H3PO2 (hypophosphorous acid)
/// Both products are charge-neutral, so charge is conserved.
pub(crate) fn generate_p4_hydrolysis(
    substance: &crate::chemistry::substance::Substance,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let structure = match &substance.molecular_structure {
        Some(s) => s,
        None => return Ok(None),
    };
    let phosphorus_atoms: Vec<usize> = (0..structure.atoms.len())
        .filter(|&i| structure.atoms[i].element == "P")
        .collect();
    if phosphorus_atoms.len() != 4 {
        return Ok(None);
    }
    for &p in &phosphorus_atoms {
        let p_neighbors: usize = structure
            .neighbors(p)
            .into_iter()
            .filter(|(n, order)| {
                structure.atoms[*n].element == "P" && bond_order_matches(*order, 1.0)
            })
            .count();
        if p_neighbors != 3 {
            return Ok(None);
        }
    }

    let phosphine = resolver.resolve(crate::chemistry::frowns::parse_frowns("P")?)?;
    let hypophosphorous_acid =
        resolver.resolve(crate::chemistry::frowns::parse_frowns("OP(O)")?)?;

    Ok(Some(
        Reaction::builder(format!("p4_hydrolysis_{}", substance.id))
            .reactant(substance.id.clone(), 1, 1)
            .reactant("destroy:water", 6, 1)
            .product(phosphine, 1)
            .product(hypophosphorous_acid, 3)
            .condition(
                crate::chemistry::condition::ReactionCondition::new(
                    "P4 hydrolysis requires moderate temperature",
                )
                .min_temperature_kelvin(300.0),
            )
            .activation_energy_kj_per_mol(60.0)
            .build(),
    ))
}
