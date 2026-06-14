use std::collections::{BTreeMap, BTreeSet};

use super::catalog::GeneratedOrganicCatalog;
use super::generators::*;
use super::resolver::DerivedSubstanceResolver;
use super::space::{GenerationScope, OrganicGenerationSpace, SiteParticipant};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::reactive_site::ReactiveSiteKind;
use crate::chemistry::registry::{ChemistryRegistry, ChemistryRegistryBuilder};
use crate::chemistry::selectivity::types::SelectivityContext;
use crate::chemistry::substance::{Substance, SubstanceId};

pub fn destroy_registry_with_generated_reactions_builder(
) -> ChemistryResult<ChemistryRegistryBuilder> {
    let base_registry = crate::chemistry::destroy_registry_builder()?.build()?;
    let generated = generate_organic_reactions(&base_registry)?;
    let mut builder = ChemistryRegistryBuilder::from_registry(&base_registry);
    for substance in generated.substances {
        builder = builder.substance(substance);
    }
    for reaction in generated.reactions {
        builder = builder.reaction(reaction);
    }
    Ok(builder)
}

pub(crate) fn generate_organic_reactions(
    registry: &ChemistryRegistry,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let scope = GenerationScope::all(registry);
    let space = OrganicGenerationSpace::new(registry.substances(), &scope)?;
    generate_organic_reactions_with_space(&space, None, &SelectivityContext::default())
}

#[cfg(test)]
pub(crate) fn generate_organic_reactions_for_substances(
    substances: &[&Substance],
    seeds: &BTreeSet<SubstanceId>,
    scope: &BTreeSet<SubstanceId>,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let scope = GenerationScope::from_substances(scope);
    let space = OrganicGenerationSpace::new(substances.iter().copied(), &scope)?;
    generate_organic_reactions_with_space(&space, Some(seeds), &SelectivityContext::default())
}

fn generate_organic_reactions_with_space(
    space: &OrganicGenerationSpace<'_>,
    seeds: Option<&BTreeSet<SubstanceId>>,
    context: &SelectivityContext,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let canonical_to_id = canonical_to_id_from_substances(space.all_substances.iter().copied())?;
    let seed_ids = seeds.cloned().unwrap_or_else(|| {
        space
            .all_substances
            .iter()
            .map(|substance| substance.id.clone())
            .collect()
    });
    generate_organic_reactions_for_seed_substances(space, &seed_ids, canonical_to_id, context)
}

pub(crate) fn generate_organic_reactions_for_seed_participants<'a>(
    space: &OrganicGenerationSpace<'a>,
    seed_participants: impl IntoIterator<Item = SiteParticipant<'a>>,
    canonical_to_id: BTreeMap<String, SubstanceId>,
    context: &SelectivityContext,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let mut resolver = DerivedSubstanceResolver::new_from_canonical_to_id(canonical_to_id);
    let mut reactions = Vec::new();
    let mut reaction_ids = BTreeSet::new();

    for participant in seed_participants {
        match participant.site.kind {
            ReactiveSiteKind::Halide => {
                let site = participant.clone().halide_site()?;
                if let Some(reaction) =
                    generate_halide_hydroxide_substitution(&site, &mut resolver, context)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) =
                    generate_halide_ammonia_substitution(&site, &mut resolver, context)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) =
                    generate_halide_cyanide_substitution(&site, &mut resolver, context)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_halide_dehydrohalogenation(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                for metal in [
                    OrganometallicFormationMetal::Magnesium,
                    OrganometallicFormationMetal::Lithium,
                ] {
                    if let Some(reaction) =
                        generate_organometallic_formation(&site, metal, &mut resolver)?
                    {
                        push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::Alcohol => {
                let site = participant.clone().alcohol_site()?;
                for reaction in generate_alcohol_oxidations(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_alcohol_dehydration(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_keto_enol_tautomerization(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                let reaction = generate_thionyl_chloride_substitution(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_alcohol_silyl_protection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Alkoxide => {
                let site = participant.clone().alkoxide_site()?;
                let reaction = generate_alkoxide_protonation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Nitrile => {
                let site = participant.clone().nitrile_site()?;
                let reaction = generate_nitrile_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_nitrile_hydrogenation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Nitro => {
                let site = participant.clone().nitro_site()?;
                let reaction = generate_nitro_hydrogenation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Oxime => {
                let site = participant.clone().oxime_site()?;
                for reaction in generate_beckmann_rearrangements(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::AcylChloride => {
                let site = participant.clone().acyl_chloride_site()?;
                let reaction = generate_acyl_chloride_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::AcidAnhydride => {
                let site = participant.clone().acid_anhydride_site()?;
                for reaction in generate_anhydride_hydrolysis(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::CarboxylicAcid => {
                let site = participant.clone().carboxylic_acid_site()?;
                let reaction = generate_acyl_chloride_formation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => {
                let site = participant.clone().carbonyl_site()?;
                for reaction in generate_baeyer_villiger_rearrangements(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                for reaction in generate_aldehyde_oxidations(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                let reaction =
                    generate_cyanide_nucleophilic_addition(&site, &mut resolver, context)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction =
                    generate_borohydride_carbonyl_reduction(&site, &mut resolver, context)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_wolff_kishner_reduction(&site, &mut resolver, context)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Amide => {
                let site = participant.clone().amide_site()?;
                if let Some(reaction) = generate_amide_hydrolysis(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Ester => {
                let site = participant.clone().ester_site()?;
                let reaction = generate_ester_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_lah_ester_reduction(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::PrimaryAmine => {
                let site = participant.clone().amine_site()?;
                let reaction = generate_amine_phosgenation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_amine_boc_protection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_amine_cbz_protection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::NonTertiaryAmine => {
                let site = participant.clone().amine_site()?;
                let reaction = generate_cyanamide_addition(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_amine_boc_protection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_amine_cbz_protection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::PhosphoniumSalt => {
                let site = participant.clone().phosphonium_salt_site()?;
                if let Some(reaction) = generate_phosphonium_ylide_formation(&site, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::PhosphorusYlide => {
                let ylide_site = participant.clone().phosphorus_ylide_site()?;
                for carbonyl_kind in carbonyl_site_kinds() {
                    for carbonyl in space.sites_of(&carbonyl_kind) {
                        let carbonyl_site = carbonyl.carbonyl_site()?;
                        if let Some(reaction) =
                            generate_wittig_olefination(&ylide_site, &carbonyl_site, &mut resolver)?
                        {
                            push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                        }
                    }
                }
            }
            ReactiveSiteKind::PhosphonateCarbanion => {
                let phosphonate_site = participant.clone().phosphonate_carbanion_site()?;
                for carbonyl_kind in carbonyl_site_kinds() {
                    for carbonyl in space.sites_of(&carbonyl_kind) {
                        let carbonyl_site = carbonyl.carbonyl_site()?;
                        if let Some(reaction) = generate_horner_wadsworth_emmons_olefination(
                            &phosphonate_site,
                            &carbonyl_site,
                            &mut resolver,
                        )? {
                            push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                        }
                    }
                }
            }
            ReactiveSiteKind::SulfoneCarbanion => {
                let sulfone_site = participant.clone().sulfone_carbanion_site()?;
                for carbonyl_kind in carbonyl_site_kinds() {
                    for carbonyl in space.sites_of(&carbonyl_kind) {
                        let carbonyl_site = carbonyl.carbonyl_site()?;
                        if let Some(reaction) = generate_julia_olefination(
                            &sulfone_site,
                            &carbonyl_site,
                            &mut resolver,
                        )? {
                            push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                        }
                    }
                }
            }
            ReactiveSiteKind::Isocyanate => {
                let site = participant.clone().isocyanate_site()?;
                let reaction = generate_isocyanate_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_isocyanate_ammonolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Borane => {
                let site = participant.clone().borane_site()?;
                let reaction = generate_borane_oxidation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::BorateEster => {
                let site = participant.clone().borate_ester_site()?;
                let reaction = generate_borate_ester_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Alkene => {
                let site = participant.clone().unsaturated_bond_site()?;
                if let Some(reaction) = generate_alkene_photoisomerization(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_chain_growth_polymerization(&site, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                for reaction in generate_retro_diels_alder(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_alkene_epoxidation(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                for spec in electrophilic_addition_specs(false) {
                    let reaction = match generate_electrophilic_addition(&site, spec, &mut resolver)
                    {
                        Ok(reaction) => reaction,
                        Err(error) if is_unknown_stereo_distribution(&error) => continue,
                        Err(error) => return Err(error),
                    };
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                for enol in space.sites_of(&ReactiveSiteKind::Enol) {
                    let center = enol.alpha_carbon_center()?;
                    if let Some(reaction) =
                        generate_michael_addition(&center, &site, &mut resolver)?
                    {
                        push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                    }
                }
                // Diels–Alder: this alkene is the diene seed; pair it with every
                // other alkene (on a different molecule) as the dienophile.
                for other in space.sites_of(&ReactiveSiteKind::Alkene) {
                    let dienophile = other.unsaturated_bond_site()?;
                    if let Some(reaction) = generate_diels_alder(&site, &dienophile, &mut resolver)?
                    {
                        push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::Alkyne => {
                let site = participant.clone().unsaturated_bond_site()?;
                for spec in electrophilic_addition_specs(true) {
                    let reaction = match generate_electrophilic_addition(&site, spec, &mut resolver)
                    {
                        Ok(reaction) => reaction,
                        Err(error) if is_unknown_stereo_distribution(&error) => continue,
                        Err(error) => return Err(error),
                    };
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            // Protecting group sites - generate deprotection reactions
            ReactiveSiteKind::SilylEther => {
                let site = participant.clone().silyl_ether_center()?;
                let reaction = generate_silyl_ether_deprotection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Acetal | ReactiveSiteKind::Ketal => {
                let site = participant.clone().acetal_center()?;
                let reaction = generate_acetal_deprotection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::BocCarbamate => {
                let site = participant.clone().boc_carbamate_center()?;
                let reaction = generate_boc_deprotection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::CbzCarbamate => {
                let site = participant.clone().cbz_carbamate_center()?;
                let reaction = generate_cbz_deprotection(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Sulfide => {
                let site = participant.clone().sulfide_site()?;
                if let Some(reaction) =
                    generate_sulfide_oxidation_to_sulfoxide(&site, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Sulfoxide => {
                let site = participant.clone().sulfoxide_site()?;
                if let Some(reaction) =
                    generate_sulfoxide_oxidation_to_sulfone(&site, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            _ => {}
        }

        generate_pair_reactions_for_seed(
            participant,
            space,
            &mut resolver,
            &mut reactions,
            &mut reaction_ids,
            context,
        )?;
    }

    Ok(GeneratedOrganicCatalog {
        substances: resolver.substances,
        reactions,
    })
}

pub(crate) fn generate_organic_reactions_for_seed_substances<'a>(
    space: &OrganicGenerationSpace<'a>,
    seeds: &BTreeSet<SubstanceId>,
    canonical_to_id: BTreeMap<String, SubstanceId>,
    context: &SelectivityContext,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let seed_participants = space
        .site_participants()
        .filter(|participant| participant.is_seed(Some(seeds)))
        .filter(|participant| is_generator_seed_site(&participant.site.kind))
        .collect::<Vec<_>>();
    let mut generated = generate_organic_reactions_for_seed_participants(
        space,
        seed_participants,
        canonical_to_id,
        context,
    )?;
    let canonical_to_id = canonical_to_id_from_generated(space, &generated)?;
    let mut resolver = DerivedSubstanceResolver::new_from_canonical_to_id(canonical_to_id);
    let mut reaction_ids = generated
        .reactions
        .iter()
        .map(|reaction| reaction.id.to_string())
        .collect::<BTreeSet<_>>();
    for substance in &space.all_substances {
        if seeds.contains(&substance.id) {
            if let Some(reaction) = generate_complete_combustion(substance)? {
                push_unique_reaction(&mut generated.reactions, &mut reaction_ids, reaction)?;
            }
            let mut radical_halogens = Vec::new();
            if space.contains_substance("destroy:chlorine") {
                radical_halogens.push(RadicalHalogen::Chlorine);
            }
            if space.contains_substance("destroy:bromine") {
                radical_halogens.push(RadicalHalogen::Bromine);
            }
            for reaction in
                generate_radical_halogenations(substance, &mut resolver, &radical_halogens)?
            {
                push_unique_reaction(&mut generated.reactions, &mut reaction_ids, reaction)?;
            }
            for reaction in generate_cracking(substance, &mut resolver)? {
                push_unique_reaction(&mut generated.reactions, &mut reaction_ids, reaction)?;
            }
            for reaction in generate_pyrolysis(substance, &mut resolver)? {
                push_unique_reaction(&mut generated.reactions, &mut reaction_ids, reaction)?;
            }
            // Bimolecular dehydrogenative coupling joins this light hydrocarbon to
            // every other (including itself: 2 CH4 -> C2H6 + H2). It is a GROWTH
            // reaction — its alkane product is itself a valid coupling substrate — but
            // generation is always bounded by a step limit (max_iterations), so the
            // growth cascade terminates at the requested depth. There is no unbounded
            // enumeration mode for it to diverge in.
            for partner in &space.all_substances {
                for reaction in
                    generate_dehydrogenative_coupling(substance, partner, &mut resolver)?
                {
                    push_unique_reaction(&mut generated.reactions, &mut reaction_ids, reaction)?;
                }
            }
            // Step-growth polycondensation pairs this diacid with every other
            // difunctional comonomer (diol/diamine). The generator self-gates to
            // clean difunctional pairs, so most pairings return None cheaply.
            for comonomer in &space.all_substances {
                if let Some(reaction) =
                    generate_polycondensation(substance, comonomer, &mut resolver)?
                {
                    push_unique_reaction(&mut generated.reactions, &mut reaction_ids, reaction)?;
                }
            }
            // White phosphorus hydrolysis: P4 + OH⁻ + H2O → PH3 + hypophosphite
            if substance.id.as_str() == "destroy:white_phosphorus" {
                if let Some(reaction) = generate_p4_hydrolysis(substance, &mut resolver)? {
                    push_unique_reaction(&mut generated.reactions, &mut reaction_ids, reaction)?;
                }
            }
        }
    }
    generate_site_reactions_for_seed_participants(
        space,
        space
            .site_participants()
            .filter(|participant| participant.is_seed(Some(seeds)))
            .collect::<Vec<_>>(),
        &mut resolver,
        &mut generated.reactions,
        &mut reaction_ids,
        context,
    )?;
    generated.substances.extend(resolver.substances);
    Ok(generated)
}

fn is_generator_seed_site(kind: &ReactiveSiteKind) -> bool {
    matches!(
        kind,
        ReactiveSiteKind::Halide
            | ReactiveSiteKind::Alcohol
            | ReactiveSiteKind::Alkoxide
            | ReactiveSiteKind::Nitrile
            | ReactiveSiteKind::Nitro
            | ReactiveSiteKind::Oxime
            | ReactiveSiteKind::AcylChloride
            | ReactiveSiteKind::AcidAnhydride
            | ReactiveSiteKind::CarboxylicAcid
            | ReactiveSiteKind::Aldehyde
            | ReactiveSiteKind::Ketone
            | ReactiveSiteKind::Carbonyl
            | ReactiveSiteKind::Ester
            | ReactiveSiteKind::Amide
            | ReactiveSiteKind::PrimaryAmine
            | ReactiveSiteKind::NonTertiaryAmine
            | ReactiveSiteKind::NucleophilicPhosphorus
            | ReactiveSiteKind::Phosphine
            | ReactiveSiteKind::PhosphoniumSalt
            | ReactiveSiteKind::PhosphorusYlide
            | ReactiveSiteKind::SilylEther
            | ReactiveSiteKind::Acetal
            | ReactiveSiteKind::Ketal
            | ReactiveSiteKind::BocCarbamate
            | ReactiveSiteKind::CbzCarbamate
            | ReactiveSiteKind::PhosphonateCarbanion
            | ReactiveSiteKind::SulfoneCarbanion
            | ReactiveSiteKind::Sulfide
            | ReactiveSiteKind::Sulfoxide
            | ReactiveSiteKind::Isocyanate
            | ReactiveSiteKind::Borane
            | ReactiveSiteKind::BorateEster
            | ReactiveSiteKind::Alkene
            | ReactiveSiteKind::Alkyne
            | ReactiveSiteKind::ArylHalide
            | ReactiveSiteKind::Enol
            | ReactiveSiteKind::Enolate
    )
}

fn canonical_to_id_from_generated(
    space: &OrganicGenerationSpace<'_>,
    generated: &GeneratedOrganicCatalog,
) -> ChemistryResult<BTreeMap<String, SubstanceId>> {
    let mut canonical_to_id =
        canonical_to_id_from_substances(space.all_substances.iter().copied())?;
    for substance in &generated.substances {
        if let Some(structure) = &substance.molecular_structure {
            canonical_to_id
                .entry(crate::chemistry::frowns::write_frowns(structure)?)
                .or_insert_with(|| substance.id.clone());
        }
    }
    Ok(canonical_to_id)
}

fn is_unknown_stereo_distribution(error: &ChemistryError) -> bool {
    matches!(
        error,
        ChemistryError::InvalidReaction { reason, .. }
            if reason.contains("stereo distribution")
                || reason.contains("stereo mixture has no quantitative distribution")
    )
}

fn canonical_to_id_from_substances<'a>(
    substances: impl IntoIterator<Item = &'a Substance>,
) -> ChemistryResult<BTreeMap<String, SubstanceId>> {
    let mut canonical_to_id = BTreeMap::new();
    for substance in substances {
        if let Some(structure) = &substance.molecular_structure {
            canonical_to_id
                .entry(crate::chemistry::frowns::write_frowns(structure)?)
                .or_insert_with(|| substance.id.clone());
        }
    }
    Ok(canonical_to_id)
}

fn generate_pair_reactions_for_seed(
    seed: SiteParticipant<'_>,
    space: &OrganicGenerationSpace<'_>,
    resolver: &mut DerivedSubstanceResolver,
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
    context: &SelectivityContext,
) -> ChemistryResult<()> {
    match seed.site.kind {
        ReactiveSiteKind::CarboxylicAcid => {
            let acid_site = seed.clone().carboxylic_acid_site()?;
            for alcohol in space.sites_of(&ReactiveSiteKind::Alcohol) {
                let alcohol_site = alcohol.alcohol_site()?;
                if let Some(reaction) = generate_carboxylic_acid_esterification(
                    &acid_site,
                    &alcohol_site,
                    resolver,
                    context,
                )? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                // Intramolecular closure to a lactone when the alcohol is on the
                // same molecule (self-gated by substance id and ring size).
                if let Some(reaction) = generate_lactonization(&acid_site, &alcohol_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            // Condensation with a second carboxylic acid to an anhydride (incl.
            // self-condensation to the symmetric anhydride). The generator folds
            // the acid pair into a canonical order, so seeding from either acid
            // yields one reaction that push_unique_reaction then collapses.
            for other_acid in space.sites_of(&ReactiveSiteKind::CarboxylicAcid) {
                let other_acid_site = other_acid.carboxylic_acid_site()?;
                if let Some(reaction) =
                    generate_acid_anhydride_formation(&acid_site, &other_acid_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            // Intramolecular closure to a lactam when an amine is on the same molecule.
            for amine_kind in [
                ReactiveSiteKind::PrimaryAmine,
                ReactiveSiteKind::NonTertiaryAmine,
            ] {
                for amine in space.sites_of(&amine_kind) {
                    let amine_site = amine.amine_site()?;
                    if let Some(reaction) =
                        generate_lactamization(&acid_site, &amine_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                    // Intermolecular condensation to an open-chain amide when the
                    // amine is on a different molecule (self-gated by substance id).
                    if let Some(reaction) = generate_amidation(&acid_site, &amine_site, resolver)? {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
        }
        ReactiveSiteKind::Alcohol => {
            let alcohol_site = seed.clone().alcohol_site()?;
            for acid in space.sites_of(&ReactiveSiteKind::CarboxylicAcid) {
                let acid_site = acid.carboxylic_acid_site()?;
                if let Some(reaction) = generate_carboxylic_acid_esterification(
                    &acid_site,
                    &alcohol_site,
                    resolver,
                    context,
                )? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for acyl_chloride in space.sites_of(&ReactiveSiteKind::AcylChloride) {
                let acyl_chloride_site = acyl_chloride.acyl_chloride_site()?;
                let reaction = generate_acyl_chloride_esterification(
                    &acyl_chloride_site,
                    &alcohol_site,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for anhydride in space.sites_of(&ReactiveSiteKind::AcidAnhydride) {
                let anhydride_site = anhydride.acid_anhydride_site()?;
                for reaction in
                    generate_anhydride_alcohol_acylation(&anhydride_site, &alcohol_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for carbonyl_kind in carbonyl_site_kinds() {
                for carbonyl in space.sites_of(&carbonyl_kind) {
                    let carbonyl_site = carbonyl.carbonyl_site()?;
                    let reaction = generate_acetal_formation(
                        &carbonyl_site,
                        &alcohol_site,
                        resolver,
                        context,
                    )?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::AcylChloride => {
            let acyl_chloride_site = seed.clone().acyl_chloride_site()?;
            for alcohol in space.sites_of(&ReactiveSiteKind::Alcohol) {
                let alcohol_site = alcohol.alcohol_site()?;
                let reaction = generate_acyl_chloride_esterification(
                    &acyl_chloride_site,
                    &alcohol_site,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for amine_kind in [
                ReactiveSiteKind::PrimaryAmine,
                ReactiveSiteKind::NonTertiaryAmine,
            ] {
                for amine in space.sites_of(&amine_kind) {
                    let amine_site = amine.amine_site()?;
                    let reaction = generate_acyl_chloride_amidation(
                        &acyl_chloride_site,
                        &amine_site,
                        resolver,
                    )?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for thiol in space.sites_of(&ReactiveSiteKind::Thiol) {
                let thiol_site = thiol.thiol_site()?;
                let reaction = generate_acyl_chloride_thioesterification(
                    &acyl_chloride_site,
                    &thiol_site,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for aromatic in space.sites_of(&ReactiveSiteKind::AromaticRing) {
                if let Some(reaction) = generate_fc_acylation(aromatic, seed.clone(), resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::AcidAnhydride => {
            let anhydride_site = seed.clone().acid_anhydride_site()?;
            for alcohol in space.sites_of(&ReactiveSiteKind::Alcohol) {
                let alcohol_site = alcohol.alcohol_site()?;
                for reaction in
                    generate_anhydride_alcohol_acylation(&anhydride_site, &alcohol_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for amine_kind in [
                ReactiveSiteKind::PrimaryAmine,
                ReactiveSiteKind::NonTertiaryAmine,
            ] {
                for amine in space.sites_of(&amine_kind) {
                    let amine_site = amine.amine_site()?;
                    for reaction in
                        generate_anhydride_amine_acylation(&anhydride_site, &amine_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            for thiol in space.sites_of(&ReactiveSiteKind::Thiol) {
                let thiol_site = thiol.thiol_site()?;
                for reaction in
                    generate_anhydride_thiol_acylation(&anhydride_site, &thiol_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::Halide => {
            let halide_site = seed.clone().halide_site()?;
            for metal in [
                OrganometallicFormationMetal::Magnesium,
                OrganometallicFormationMetal::Lithium,
            ] {
                if let Some(reaction) =
                    generate_organometallic_formation(&halide_site, metal, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for enol in space.sites_of(&ReactiveSiteKind::Enol) {
                let center = enol.alpha_carbon_center()?;
                if let Some(reaction) =
                    generate_enolate_alkylation(&center, &halide_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for amine in space.sites_of(&ReactiveSiteKind::NonTertiaryAmine) {
                let amine_site = amine.amine_site()?;
                if let Some(reaction) = generate_halide_amine_substitution(
                    &halide_site,
                    &amine_site,
                    resolver,
                    context,
                )? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            // Amide/imide/lactam N-H are alkylated too (weak nucleophile, needs
            // base). Carried as a dedicated AmideNitrogen kind so only this path
            // sees it. This is the step that methylates xanthine to caffeine.
            for amide_nitrogen in space.sites_of(&ReactiveSiteKind::AmideNitrogen) {
                let amide_nitrogen_site = amide_nitrogen.amine_site()?;
                if let Some(reaction) =
                    generate_amide_n_alkylation(&halide_site, &amide_nitrogen_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            // Intramolecular N-alkylation: an amine on the same molecule as this
            // halide closes a saturated N-heterocycle (self-gated by substance id
            // and ring size).
            for amine_kind in [
                ReactiveSiteKind::PrimaryAmine,
                ReactiveSiteKind::NonTertiaryAmine,
            ] {
                for amine in space.sites_of(&amine_kind) {
                    let amine_site = amine.amine_site()?;
                    if let Some(reaction) =
                        generate_intramolecular_n_alkylation(&amine_site, &halide_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            for aromatic in space.sites_of(&ReactiveSiteKind::AromaticRing) {
                if let Some(reaction) = generate_fc_alkylation(aromatic, seed.clone(), resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::Phosphine => {
            let phosphine_site = seed.clone().phosphine_site()?;
            for halide in space.sites_of(&ReactiveSiteKind::Halide) {
                let halide_site = halide.halide_site()?;
                if let Some(reaction) = generate_phosphonium_salt_formation(
                    &halide_site,
                    &phosphine_site,
                    resolver,
                    context,
                )? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::NucleophilicPhosphorus => {
            let phosphorus_site = seed.clone().nucleophilic_phosphorus_site()?;
            for halide in space.sites_of(&ReactiveSiteKind::Halide) {
                let halide_site = halide.halide_site()?;
                if let Some(reaction) = generate_nucleophilic_phosphorus_alkylation(
                    &halide_site,
                    &phosphorus_site,
                    resolver,
                )? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::NonTertiaryAmine => {
            let amine_site = seed.clone().amine_site()?;
            for halide in space.sites_of(&ReactiveSiteKind::Halide) {
                let halide_site = halide.halide_site()?;
                if let Some(reaction) = generate_halide_amine_substitution(
                    &halide_site,
                    &amine_site,
                    resolver,
                    context,
                )? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for carbonyl_kind in carbonyl_site_kinds() {
                for carbonyl in space.sites_of(&carbonyl_kind) {
                    let carbonyl_site = carbonyl.clone().carbonyl_site()?;
                    for alpha in space
                        .sites_of(&ReactiveSiteKind::Enol)
                        .filter(|site| site.substance.id == carbonyl_site.participant.substance.id)
                    {
                        let alpha_center = alpha.alpha_carbon_center()?;
                        if let Some(reaction) = generate_enamine_formation(
                            &carbonyl_site,
                            &amine_site,
                            &alpha_center,
                            resolver,
                        )? {
                            push_unique_reaction(reactions, reaction_ids, reaction)?;
                        }
                    }
                }
            }
            for acyl_chloride in space.sites_of(&ReactiveSiteKind::AcylChloride) {
                let acyl_chloride_site = acyl_chloride.acyl_chloride_site()?;
                let reaction =
                    generate_acyl_chloride_amidation(&acyl_chloride_site, &amine_site, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for anhydride in space.sites_of(&ReactiveSiteKind::AcidAnhydride) {
                let anhydride_site = anhydride.acid_anhydride_site()?;
                for reaction in
                    generate_anhydride_amine_acylation(&anhydride_site, &amine_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for isocyanate in space.sites_of(&ReactiveSiteKind::Isocyanate) {
                let isocyanate_site = isocyanate.isocyanate_site()?;
                if let Some(reaction) =
                    generate_isocyanate_amine_addition(&isocyanate_site, &amine_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for donor in space.sites_of(&ReactiveSiteKind::FormylationDonor) {
                let donor_site = donor.formylation_donor_center()?;
                if let Some(reaction) =
                    generate_amine_formylation(&amine_site, &donor_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::Isocyanate => {
            let isocyanate_site = seed.clone().isocyanate_site()?;
            for amine in space.sites_of(&ReactiveSiteKind::NonTertiaryAmine) {
                let amine_site = amine.amine_site()?;
                if let Some(reaction) =
                    generate_isocyanate_amine_addition(&isocyanate_site, &amine_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::FormylationDonor => {
            let donor_site = seed.clone().formylation_donor_center()?;
            for amine in space.sites_of(&ReactiveSiteKind::NonTertiaryAmine) {
                let amine_site = amine.amine_site()?;
                if let Some(reaction) =
                    generate_amine_formylation(&amine_site, &donor_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::Ester => {
            let ester_site = seed.clone().ester_site()?;
            for enol in space.sites_of(&ReactiveSiteKind::Enol) {
                let center = enol.alpha_carbon_center()?;
                if let Some(reaction) =
                    generate_claisen_condensation(&center, &ester_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => {
            let carbonyl_site = seed.clone().carbonyl_site()?;
            for carbonyl_kind in carbonyl_site_kinds() {
                for other in space
                    .sites_of(&carbonyl_kind)
                    .filter(|site| site.substance.id == carbonyl_site.participant.substance.id)
                {
                    let other_site = other.carbonyl_site()?;
                    if let Some(reaction) =
                        generate_paal_knorr_furan(&carbonyl_site, &other_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                    for amine in space.sites_of(&ReactiveSiteKind::PrimaryAmine) {
                        let amine_site = amine.amine_site()?;
                        if let Some(reaction) = generate_paal_knorr_pyrrole(
                            &carbonyl_site,
                            &other_site,
                            &amine_site,
                            resolver,
                        )? {
                            push_unique_reaction(reactions, reaction_ids, reaction)?;
                        }
                    }
                    for thiol in space.sites_of(&ReactiveSiteKind::Thiol) {
                        let thiol_site = thiol.thiol_site()?;
                        if let Some(reaction) = generate_paal_knorr_thiophene(
                            &carbonyl_site,
                            &other_site,
                            &thiol_site,
                            resolver,
                        )? {
                            push_unique_reaction(reactions, reaction_ids, reaction)?;
                        }
                    }
                }
            }
            for alcohol in space.sites_of(&ReactiveSiteKind::Alcohol) {
                let alcohol_site = alcohol.alcohol_site()?;
                let reaction =
                    generate_acetal_formation(&carbonyl_site, &alcohol_site, resolver, context)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for amine in space.sites_of(&ReactiveSiteKind::PrimaryAmine) {
                let amine_site = amine.amine_site()?;
                let reaction =
                    generate_imine_formation(&carbonyl_site, &amine_site, resolver, context)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for bis_nucleophile in space.sites_of(&ReactiveSiteKind::BisNucleophile) {
                let nucleophile_site = bis_nucleophile.bis_nucleophile_center()?;
                if let Some(reaction) = generate_hydrazone_formation(
                    &carbonyl_site,
                    &nucleophile_site,
                    resolver,
                    context,
                )? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for alpha in space
                .sites_of(&ReactiveSiteKind::Enol)
                .filter(|site| site.substance.id == carbonyl_site.participant.substance.id)
            {
                let alpha_center = alpha.alpha_carbon_center()?;
                for amine in space.sites_of(&ReactiveSiteKind::NonTertiaryAmine) {
                    let amine_site = amine.amine_site()?;
                    if let Some(reaction) = generate_enamine_formation(
                        &carbonyl_site,
                        &amine_site,
                        &alpha_center,
                        resolver,
                    )? {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
        }
        ReactiveSiteKind::BisNucleophile => {
            let nucleophile_site = seed.clone().bis_nucleophile_center()?;
            for carbonyl_kind in carbonyl_site_kinds() {
                for carbonyl in space.sites_of(&carbonyl_kind) {
                    let carbonyl_site = carbonyl.carbonyl_site()?;
                    if let Some(reaction) = generate_hydrazone_formation(
                        &carbonyl_site,
                        &nucleophile_site,
                        resolver,
                        context,
                    )? {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
        }
        ReactiveSiteKind::PrimaryAmine => {
            let amine_site = seed.clone().amine_site()?;
            for carbonyl_kind in carbonyl_site_kinds() {
                for carbonyl in space.sites_of(&carbonyl_kind) {
                    let carbonyl_site = carbonyl.carbonyl_site()?;
                    let reaction =
                        generate_imine_formation(&carbonyl_site, &amine_site, resolver, context)?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            // Intramolecular dehydrative closure to a cyclic amidine (imidazole)
            // when an amide is on the same molecule — self-gated by substance id
            // and ring size. Fuses onto an existing ring (purine) for free.
            for amide in space.sites_of(&ReactiveSiteKind::Amide) {
                let amide_site = amide.amide_site()?;
                if let Some(reaction) =
                    generate_amidine_cyclization(&amine_site, &amide_site, resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::AromaticRing => {
            for halide in space.sites_of(&ReactiveSiteKind::Halide) {
                if let Some(reaction) = generate_fc_alkylation(seed.clone(), halide, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            for acyl in space.sites_of(&ReactiveSiteKind::AcylChloride) {
                if let Some(reaction) = generate_fc_acylation(seed.clone(), acyl, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        _ => {}
    }
    Ok(())
}

fn carbonyl_site_kinds() -> [ReactiveSiteKind; 3] {
    [
        ReactiveSiteKind::Aldehyde,
        ReactiveSiteKind::Ketone,
        ReactiveSiteKind::Carbonyl,
    ]
}

fn generate_site_reactions_for_seed_participants<'a>(
    space: &OrganicGenerationSpace<'a>,
    seed_sites: impl IntoIterator<Item = SiteParticipant<'a>>,
    resolver: &mut DerivedSubstanceResolver,
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
    context: &SelectivityContext,
) -> ChemistryResult<()> {
    for seed in seed_sites {
        match seed.site.kind {
            ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => {
                let carbonyl_site = seed.clone().carbonyl_site()?;
                for reaction in generate_baeyer_villiger_rearrangements(&carbonyl_site, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                for organometallic_kind in [
                    ReactiveSiteKind::Organomagnesium,
                    ReactiveSiteKind::Organolithium,
                    ReactiveSiteKind::Organocopper,
                ] {
                    for organometallic in space.sites_of(&organometallic_kind) {
                        let reaction = generate_organometallic_carbonyl_addition(
                            seed.clone(),
                            organometallic,
                            resolver,
                            context,
                        )?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
                for enol in space.sites_of(&ReactiveSiteKind::Enol) {
                    let reaction = generate_aldol_addition(enol, seed.clone(), resolver)?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Oxime => {
                let oxime_site = seed.clone().oxime_site()?;
                for reaction in generate_beckmann_rearrangements(&oxime_site, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Organomagnesium
            | ReactiveSiteKind::Organolithium
            | ReactiveSiteKind::Organocopper => {
                for carbonyl_kind in [
                    ReactiveSiteKind::Aldehyde,
                    ReactiveSiteKind::Ketone,
                    ReactiveSiteKind::Carbonyl,
                ] {
                    for carbonyl in space.sites_of(&carbonyl_kind) {
                        let reaction = generate_organometallic_carbonyl_addition(
                            carbonyl,
                            seed.clone(),
                            resolver,
                            context,
                        )?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
                for nitrile in space.sites_of(&ReactiveSiteKind::Nitrile) {
                    let nitrile_site = nitrile.nitrile_site()?;
                    let reaction = generate_organometallic_nitrile_addition(
                        &nitrile_site,
                        seed.clone(),
                        resolver,
                    )?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                for epoxide in space.sites_of(&ReactiveSiteKind::Epoxide) {
                    let reaction =
                        generate_organometallic_epoxide_opening(epoxide, seed.clone(), resolver)?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Enol => {
                let center = seed.clone().alpha_carbon_center()?;
                for reaction in generate_alpha_halogenation(&center, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_aldol_dehydration(&center, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                for carbonyl_kind in [
                    ReactiveSiteKind::Aldehyde,
                    ReactiveSiteKind::Ketone,
                    ReactiveSiteKind::Carbonyl,
                ] {
                    for carbonyl in space.sites_of(&carbonyl_kind) {
                        let reaction = generate_aldol_addition(seed.clone(), carbonyl, resolver)?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
                for halide in space.sites_of(&ReactiveSiteKind::Halide) {
                    let halide_site = halide.halide_site()?;
                    if let Some(reaction) =
                        generate_enolate_alkylation(&center, &halide_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
                for alkene in space.sites_of(&ReactiveSiteKind::Alkene) {
                    let alkene_site = alkene.unsaturated_bond_site()?;
                    if let Some(reaction) =
                        generate_michael_addition(&center, &alkene_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
                for ester in space.sites_of(&ReactiveSiteKind::Ester) {
                    let ester_site = ester.ester_site()?;
                    if let Some(reaction) =
                        generate_claisen_condensation(&center, &ester_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::Ester => {
                let ester_site = seed.clone().ester_site()?;
                for enol in space.sites_of(&ReactiveSiteKind::Enol) {
                    let center = enol.alpha_carbon_center()?;
                    if let Some(reaction) =
                        generate_claisen_condensation(&center, &ester_site, resolver)?
                    {
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::AromaticRing => {
                if let Some(reaction) = generate_aromatic_nitration(seed.clone(), resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_aromatic_chlorination(seed.clone(), resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_aromatic_bromination(seed.clone(), resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_aromatic_sulfonation(seed, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Epoxide => {
                let reaction = generate_epoxide_hydrolysis(seed.clone(), resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
                for organometallic_kind in [
                    ReactiveSiteKind::Organomagnesium,
                    ReactiveSiteKind::Organolithium,
                    ReactiveSiteKind::Organocopper,
                ] {
                    for organometallic in space.sites_of(&organometallic_kind) {
                        let reaction = generate_organometallic_epoxide_opening(
                            seed.clone(),
                            organometallic,
                            resolver,
                        )?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::Nitrile => {
                let nitrile_site = seed.clone().nitrile_site()?;
                for organometallic_kind in [
                    ReactiveSiteKind::Organomagnesium,
                    ReactiveSiteKind::Organolithium,
                    ReactiveSiteKind::Organocopper,
                ] {
                    for organometallic in space.sites_of(&organometallic_kind) {
                        let reaction = generate_organometallic_nitrile_addition(
                            &nitrile_site,
                            organometallic,
                            resolver,
                        )?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::ArylHalide => {
                if let Some(reaction) =
                    generate_aryl_halide_hydroxide_substitution(seed.clone(), resolver)?
                {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_aryl_halide_ammonia_substitution(seed, resolver)? {
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            _ => {}
        }
    }
    Ok(())
}

fn push_unique_reaction(
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
    reaction: Reaction,
) -> ChemistryResult<()> {
    let id = reaction.id.to_string();
    if reaction_ids.insert(id.clone()) {
        reactions.push(reaction);
    }
    Ok(())
}
