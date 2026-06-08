use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{MolecularEditor, MolecularStructure};
use crate::chemistry::organic::space::SiteParticipant;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{
        ReactionType, SelectivityContext, SelectivityProfile, SiteDescriptor, SubstitutionDegree,
    },
};
use crate::chemistry::substance::SubstanceId;
use std::collections::BTreeSet;

#[derive(Clone, Copy)]
enum AcylDonorKind {
    AcylChloride,
    AcidAnhydride,
}

struct AcylDonor<'a> {
    participant: &'a SiteParticipant<'a>,
    structure: &'a MolecularStructure,
    acyl_carbon: usize,
    leaving_atoms: Vec<usize>,
    leaving_attachment: Option<usize>,
    leaving_product_static: Option<&'static str>,
    leaving_product_prefix: &'static str,
    kind: AcylDonorKind,
}

#[derive(Clone, Copy)]
enum AcylNucleophileKind {
    Alcohol,
    Amine,
    Thiol,
}

struct AcylNucleophile<'a> {
    participant: &'a SiteParticipant<'a>,
    structure: &'a MolecularStructure,
    attack_atom: usize,
    proton: usize,
    kind: AcylNucleophileKind,
}

pub(crate) fn generate_carboxylic_acid_esterification(
    acid_site: &CarboxylicAcidSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    let acid_desc = SiteDescriptorBuilder::from_carboxylic_acid_site(acid_site);
    let alcohol_desc = SiteDescriptorBuilder::from_alcohol_site(alcohol_site);

    let base_ea = 25.0;

    let acid = acid_site.participant.substance;
    let acid_structure = acid_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let alcohol_structure = alcohol_site.participant.structure;
    let acid_carbon = acid_site.carbon;
    let acid_hydroxyl_oxygen = acid_site.hydroxyl_oxygen;
    let acid_proton = acid_site.hydroxyl_hydrogen;
    let alcohol_oxygen = alcohol_site.oxygen;
    let alcohol_proton = alcohol_site.hydrogen;

    let mut acid_editor = MolecularEditor::new(acid_structure);
    let acid_mapping = acid_editor.remove_atoms(&[acid_proton, acid_hydroxyl_oxygen])?;
    let acid_carbon = mapped_atom(&acid_mapping, acid_carbon, "acid carbon")?;
    let acid_fragment = acid_editor.finish()?;

    let mut alcohol_editor = MolecularEditor::new(alcohol_structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol_proton])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol_oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let product_structure = MolecularEditor::join_structures(
        &acid_fragment,
        acid_carbon,
        &alcohol_fragment,
        alcohol_oxygen,
        1.0,
    )?;
    let product = resolver.resolve(product_structure)?;
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "carboxylic_acid_esterification",
            &acid_site.participant,
            &alcohol_site.participant,
        ))
        .reactant(acid.id.clone(), 1, 1)
        .reactant(alcohol.id.clone(), 1, 0)
        .catalyst_order("destroy:sulfuric_acid", 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .activation_energy_kj_per_mol(base_ea)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::EsterProtection, acid_desc)
                .with_secondary_site(alcohol_desc)
                .never_suppress(),
        )
        .build(),
    ))
}

/// Intermolecular amidation (Fischer-type): a carboxylic acid and an amine on
/// SEPARATE molecules condense to an amide, expelling water. Mirrors
/// `generate_carboxylic_acid_esterification` but the nucleophile is the amine
/// nitrogen rather than an alcohol oxygen, so the new bond is acyl-C–N. The
/// intramolecular (same-molecule) case is handled separately by
/// `generate_lactamization`, which closes a ring instead.
pub(crate) fn generate_amidation(
    acid_site: &CarboxylicAcidSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Same-molecule pairs are a lactam closure, not an intermolecular condensation.
    if acid_site.participant.substance.id == amine_site.participant.substance.id {
        return Ok(None);
    }
    let acid_desc = SiteDescriptorBuilder::from_carboxylic_acid_site(acid_site);
    let amine_desc = SiteDescriptorBuilder::from_amine_site(amine_site);

    let acid = acid_site.participant.substance;
    let amine = amine_site.participant.substance;
    let acid_structure = acid_site.participant.structure;
    let amine_structure = amine_site.participant.structure;

    // Acid sheds its -OH (oxygen + proton); the amine sheds one N-H. Together they
    // leave as water, and the acyl carbon bonds to the now-freed amine nitrogen.
    let mut acid_editor = MolecularEditor::new(acid_structure);
    let acid_mapping =
        acid_editor.remove_atoms(&[acid_site.hydroxyl_hydrogen, acid_site.hydroxyl_oxygen])?;
    let acid_carbon = mapped_atom(&acid_mapping, acid_site.carbon, "acid carbon")?;
    let acid_fragment = acid_editor.finish()?;

    let mut amine_editor = MolecularEditor::new(amine_structure);
    let amine_mapping = amine_editor.remove_atoms(&[amine_site.hydrogens[0]])?;
    let amine_nitrogen = mapped_atom(&amine_mapping, amine_site.nitrogen, "amine nitrogen")?;
    let amine_fragment = amine_editor.finish()?;

    let product_structure = MolecularEditor::join_structures(
        &acid_fragment,
        acid_carbon,
        &amine_fragment,
        amine_nitrogen,
        1.0,
    )?;
    let product = resolver.resolve(product_structure)?;
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "amidation",
            &acid_site.participant,
            &amine_site.participant,
        ))
        .reactant(acid.id.clone(), 1, 1)
        .reactant(amine.id.clone(), 1, 0)
        .product(product, 1)
        .product("destroy:water", 1)
        // Amide condensation needs heat and a dry medium (water reverses it); the
        // Lactamization selectivity arm models exactly that pull, so reuse it
        // rather than adding a near-duplicate ReactionType.
        .condition(
            ReactionCondition::new("amidation is driven by heat and a water-poor medium")
                .max_water_activity(0.5),
        )
        .activation_energy_kj_per_mol(45.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::Lactamization, acid_desc)
                .with_secondary_site(amine_desc),
        )
        .build(),
    ))
}

pub(crate) fn generate_acyl_chloride_formation(
    site: &CarboxylicAcidSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let hydroxyl_oxygen = site.hydroxyl_oxygen;
    let proton = site.hydroxyl_hydrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydroxyl_oxygen, proton])?;
    let carbon = mapped_atom(&mapping, carbon, "carboxylic acid carbon")?;
    editor.add_atom(carbon, "Cl", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "acyl_chloride_formation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:phosgene", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .product("destroy:carbon_dioxide", 1)
    .build())
}

pub(crate) fn generate_acyl_chloride_hydrolysis(
    site: &AcylChlorideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    generate_acyl_donor_hydrolysis(&AcylDonor::from_acyl_chloride(site), resolver)
}

pub(crate) fn generate_acyl_chloride_esterification(
    acyl_chloride_site: &AcylChlorideSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    generate_acyl_transfer(
        &AcylDonor::from_acyl_chloride(acyl_chloride_site),
        &AcylNucleophile::from_alcohol(alcohol_site),
        resolver,
    )
}

pub(crate) fn generate_acyl_chloride_amidation(
    acyl_chloride_site: &AcylChlorideSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    generate_acyl_transfer(
        &AcylDonor::from_acyl_chloride(acyl_chloride_site),
        &AcylNucleophile::from_amine(amine_site),
        resolver,
    )
}

pub(crate) fn generate_acyl_chloride_thioesterification(
    acyl_chloride_site: &AcylChlorideSite<'_>,
    thiol_site: &ThiolSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    generate_acyl_transfer(
        &AcylDonor::from_acyl_chloride(acyl_chloride_site),
        &AcylNucleophile::from_thiol(thiol_site),
        resolver,
    )
}

/// Intermolecular acid anhydride formation: two carboxylic acids condense,
/// expelling water. One acid keeps its carbonyl and its hydroxyl oxygen, which
/// becomes the bridge; the other sheds its entire -OH (oxygen + proton) and
/// donates its acyl carbon to that bridge oxygen. This is the dehydrative mirror
/// of esterification, with a second acid playing the alcohol role. The reaction
/// is symmetric in its two acid partners, so the id folds the substance ids in a
/// canonical order and the donor/bridge roles are assigned by that same order —
/// seeding from either partner yields one identical reaction, which
/// `push_unique_reaction` then collapses. Self-condensation (one acid with
/// itself, 2 RCOOH -> (RCO)2O) is allowed and produces the symmetric anhydride.
pub(crate) fn generate_acid_anhydride_formation(
    first_site: &CarboxylicAcidSite<'_>,
    second_site: &CarboxylicAcidSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Canonical ordering: the bridge-oxygen donor is the lexicographically
    // smaller substance id (ties broken by acyl-carbon index for two sites on
    // one molecule), so the symmetric pair maps to a single reaction.
    let first_key = (first_site.participant.substance.id.as_str(), first_site.carbon);
    let second_key = (
        second_site.participant.substance.id.as_str(),
        second_site.carbon,
    );
    if first_key == second_key {
        return Ok(None);
    }
    let (bridge_site, acyl_site) = if first_key <= second_key {
        (first_site, second_site)
    } else {
        (second_site, first_site)
    };

    // Bridge acid: drop only the hydroxyl proton, keeping its oxygen as the
    // bridge. Acyl acid: drop its whole -OH so its acyl carbon bonds the bridge.
    let mut bridge_editor = MolecularEditor::new(bridge_site.participant.structure);
    let bridge_mapping = bridge_editor.remove_atoms(&[bridge_site.hydroxyl_hydrogen])?;
    let bridge_oxygen = mapped_atom(
        &bridge_mapping,
        bridge_site.hydroxyl_oxygen,
        "anhydride bridge oxygen",
    )?;
    let bridge_fragment = bridge_editor.finish()?;

    let mut acyl_editor = MolecularEditor::new(acyl_site.participant.structure);
    let acyl_mapping =
        acyl_editor.remove_atoms(&[acyl_site.hydroxyl_oxygen, acyl_site.hydroxyl_hydrogen])?;
    let acyl_carbon = mapped_atom(&acyl_mapping, acyl_site.carbon, "anhydride acyl carbon")?;
    let acyl_fragment = acyl_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &bridge_fragment,
        bridge_oxygen,
        &acyl_fragment,
        acyl_carbon,
        1.0,
    )?)?;

    let self_condensation =
        bridge_site.participant.substance.id == acyl_site.participant.substance.id;
    let acyl_coefficient = if self_condensation { 2 } else { 1 };
    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "acid_anhydride_formation",
        &bridge_site.participant,
        &acyl_site.participant,
    ))
    .product(product, 1)
    .product("destroy:water", 1)
    .catalyst_order("destroy:sulfuric_acid", 1)
    .condition(
        ReactionCondition::new(
            "anhydride formation is an equilibrium dehydration favored only under \
             acid catalysis and water-poor conditions",
        )
        .acidity(AcidityCondition::Acidic)
        .max_water_activity(0.5),
    )
    .activation_energy_kj_per_mol(35.0)
    .selectivity_profile(SelectivityProfile::new(
        ReactionType::AcylSubstitution,
        SiteDescriptorBuilder::from_carboxylic_acid_site(bridge_site),
    ));
    builder = if self_condensation {
        builder.reactant(bridge_site.participant.substance.id.clone(), acyl_coefficient, 1)
    } else {
        builder
            .reactant(bridge_site.participant.substance.id.clone(), 1, 1)
            .reactant(acyl_site.participant.substance.id.clone(), 1, 1)
    };
    Ok(Some(builder.build()))
}

pub(crate) fn generate_anhydride_hydrolysis(
    site: &AcidAnhydrideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    let donors = AcylDonor::from_anhydride(site);
    let mut reactions = Vec::new();
    let mut ids = BTreeSet::new();
    for donor in donors {
        let reaction = generate_acyl_donor_hydrolysis(&donor, resolver)?;
        if ids.insert(reaction.id.to_string()) {
            reactions.push(reaction);
        }
    }
    Ok(reactions)
}

pub(crate) fn generate_anhydride_alcohol_acylation(
    site: &AcidAnhydrideSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    generate_anhydride_acyl_transfer(site, &AcylNucleophile::from_alcohol(alcohol_site), resolver)
}

pub(crate) fn generate_anhydride_amine_acylation(
    site: &AcidAnhydrideSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    generate_anhydride_acyl_transfer(site, &AcylNucleophile::from_amine(amine_site), resolver)
}

pub(crate) fn generate_anhydride_thiol_acylation(
    site: &AcidAnhydrideSite<'_>,
    thiol_site: &ThiolSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    generate_anhydride_acyl_transfer(site, &AcylNucleophile::from_thiol(thiol_site), resolver)
}

pub(crate) fn generate_ester_hydrolysis(
    site: &EsterSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let alkoxy_oxygen = site.alkoxy_oxygen;

    let alkoxy_branch = ester_alkoxy_branch(structure, alkoxy_oxygen, carbon)?;

    let mut acid_editor = MolecularEditor::new(structure);
    let acid_mapping = acid_editor.remove_atoms(&alkoxy_branch)?;
    let acid_carbon = mapped_atom(&acid_mapping, carbon, "ester carbonyl carbon")?;
    add_hydroxyl(&mut acid_editor, acid_carbon)?;
    let acid = resolver.resolve(acid_editor.finish()?)?;

    let mut alcohol_editor = MolecularEditor::new(structure);
    let keep = alkoxy_branch
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let remove = (0..structure.atoms.len())
        .filter(|atom| !keep.contains(atom))
        .collect::<Vec<_>>();
    let alcohol_mapping = alcohol_editor.remove_atoms(&remove)?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alkoxy_oxygen, "ester alkoxy oxygen")?;
    alcohol_editor.add_atom(alcohol_oxygen, "H", 0.0, 1.0)?;
    let alcohol = resolver.resolve(alcohol_editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "ester_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(acid, 1)
    .product(alcohol, 1)
    .condition(
        ReactionCondition::new("ester hydrolysis requires acidic, water-rich conditions")
            .acidity(AcidityCondition::Acidic)
            .min_water_activity(0.35),
    )
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::EsterHydrolysis,
            SiteDescriptorBuilder::carboxylic_acid(),
        )
        .never_suppress(),
    )
    .activation_energy_kj_per_mol(42.0)
    .build())
}

pub(crate) fn generate_lah_ester_reduction(
    site: &EsterSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let carbonyl_oxygen = site.carbonyl_oxygen;
    let alkoxy_oxygen = site.alkoxy_oxygen;
    let alkoxy_branch = ester_alkoxy_branch(structure, alkoxy_oxygen, carbon)?;

    let mut acyl_editor = MolecularEditor::new(structure);
    let acyl_mapping = acyl_editor.remove_atoms(&alkoxy_branch)?;
    let carbon = mapped_atom(&acyl_mapping, carbon, "ester carbonyl carbon")?;
    let carbonyl_oxygen = mapped_atom(&acyl_mapping, carbonyl_oxygen, "ester carbonyl oxygen")?;
    acyl_editor.set_bond_order(carbon, carbonyl_oxygen, 1.0)?;
    acyl_editor.add_atom(carbonyl_oxygen, "H", 0.0, 1.0)?;
    acyl_editor.add_atom(carbon, "H", 0.0, 1.0)?;
    acyl_editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let acyl_alcohol = resolver.resolve(acyl_editor.finish()?)?;

    let mut alkoxy_editor = MolecularEditor::new(structure);
    let keep = alkoxy_branch
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let remove = (0..structure.atoms.len())
        .filter(|atom| !keep.contains(atom))
        .collect::<Vec<_>>();
    let alkoxy_mapping = alkoxy_editor.remove_atoms(&remove)?;
    let alcohol_oxygen = mapped_atom(&alkoxy_mapping, alkoxy_oxygen, "ester alkoxy oxygen")?;
    alkoxy_editor.add_atom(alcohol_oxygen, "H", 0.0, 1.0)?;
    let alkoxy_alcohol = resolver.resolve(alkoxy_editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "lah_ester_reduction",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .chemical_external_reactant(
        "lithium aluminium hydride hydride/proton equivalents",
        1.0,
        4.04,
        0,
    )
    .product(acyl_alcohol, 1)
    .product(alkoxy_alcohol, 1)
    .condition(
        ReactionCondition::new("LAH ester reduction requires dry aprotic conditions")
            .max_water_activity(0.02),
    )
    .activation_energy_kj_per_mol(18.0)
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::CarbonylReduction,
            SiteDescriptorBuilder::carboxylic_acid(),
        )
        .with_nucleophile_strength(crate::chemistry::selectivity::NucleophileStrength::VeryStrong)
        .never_suppress(),
    )
    .build())
}

pub(crate) fn generate_amide_hydrolysis(
    site: &AmideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if amide_leaving_branch(site.participant.structure, site.nitrogen, site.carbon)?.is_none() {
        return Ok(None);
    }
    generate_acyl_leaving_group_hydrolysis(
        site.participant.substance.id.as_str(),
        &site.participant,
        site.carbon,
        site.nitrogen,
        "amide_hydrolysis",
        "amide nitrogen",
        resolver,
    )
    .map(Some)
}

fn generate_anhydride_acyl_transfer(
    site: &AcidAnhydrideSite<'_>,
    nucleophile: &AcylNucleophile<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    let donors = AcylDonor::from_anhydride(site);
    let mut reactions = Vec::new();
    let mut ids = BTreeSet::new();
    for donor in donors {
        let reaction = generate_acyl_transfer(&donor, nucleophile, resolver)?;
        if ids.insert(reaction.id.to_string()) {
            reactions.push(reaction);
        }
    }
    Ok(reactions)
}

fn generate_acyl_transfer(
    donor: &AcylDonor<'_>,
    nucleophile: &AcylNucleophile<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let mut acyl_editor = MolecularEditor::new(donor.structure);
    let acyl_mapping = acyl_editor.remove_atoms(&donor.leaving_atoms)?;
    let acyl_carbon = mapped_atom(&acyl_mapping, donor.acyl_carbon, "acyl carbon")?;
    let acyl_fragment = acyl_editor.finish()?;

    let mut nucleophile_editor = MolecularEditor::new(nucleophile.structure);
    let nucleophile_mapping = nucleophile_editor.remove_atoms(&[nucleophile.proton])?;
    let attack_atom = mapped_atom(
        &nucleophile_mapping,
        nucleophile.attack_atom,
        "nucleophile atom",
    )?;
    let nucleophile_fragment = nucleophile_editor.finish()?;

    let acylated = resolver.resolve(MolecularEditor::join_structures(
        &acyl_fragment,
        acyl_carbon,
        &nucleophile_fragment,
        attack_atom,
        1.0,
    )?)?;

    let byproduct = donor.leaving_product(resolver)?;
    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        acyl_transfer_id_prefix(donor, nucleophile),
        donor.participant,
        nucleophile.participant,
    ))
    .reactant(donor.participant.substance.id.clone(), 1, 1)
    .reactant(nucleophile.participant.substance.id.clone(), 1, 1)
    .product(acylated, 1)
    .product(byproduct, 1)
    .activation_energy_kj_per_mol(acyl_transfer_activation_energy(donor, nucleophile))
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::AcylSubstitution,
            SiteDescriptorBuilder::carboxylic_acid(),
        )
        .with_secondary_site(nucleophile_descriptor(nucleophile))
        .never_suppress(),
    );
    if matches!(donor.kind, AcylDonorKind::AcidAnhydride) {
        builder = builder.condition(
            ReactionCondition::new("anhydride acyl transfer is suppressed by strongly wet media")
                .max_water_activity(0.8),
        );
    }
    Ok(builder.build())
}

fn generate_acyl_donor_hydrolysis(
    donor: &AcylDonor<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let mut editor = MolecularEditor::new(donor.structure);
    let mapping = editor.remove_atoms(&donor.leaving_atoms)?;
    let carbon = mapped_atom(&mapping, donor.acyl_carbon, "acyl donor carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let acid = resolver.resolve(editor.finish()?)?;
    let byproduct = donor.leaving_product(resolver)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        donor_hydrolysis_id_prefix(donor),
        donor.participant,
    ))
    .reactant(donor.participant.substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product(acid, 1)
    .product(byproduct, 1)
    .condition(
        ReactionCondition::new("acyl derivative hydrolysis requires available water")
            .min_water_activity(0.1),
    )
    .activation_energy_kj_per_mol(match donor.kind {
        AcylDonorKind::AcylChloride => 12.0,
        AcylDonorKind::AcidAnhydride => 18.0,
    })
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::AcylSubstitution,
            SiteDescriptorBuilder::carboxylic_acid(),
        )
        .never_suppress(),
    )
    .build())
}

fn generate_acyl_leaving_group_hydrolysis(
    substance_id: &str,
    participant: &SiteParticipant<'_>,
    acyl_carbon: usize,
    leaving_atom: usize,
    id_prefix: &'static str,
    leaving_label: &'static str,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let structure = participant.structure;
    let leaving_branch = leaving_branch(structure, leaving_atom, acyl_carbon, id_prefix)?
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id(id_prefix, participant),
            reason: format!("{leaving_label} is part of a ring; use ring-opening chemistry"),
        })?;

    let mut acid_editor = MolecularEditor::new(structure);
    let acid_mapping = acid_editor.remove_atoms(&leaving_branch)?;
    let acid_carbon = mapped_atom(&acid_mapping, acyl_carbon, "acyl carbon")?;
    add_hydroxyl(&mut acid_editor, acid_carbon)?;
    let acid = resolver.resolve(acid_editor.finish()?)?;

    let mut leaving_editor = MolecularEditor::new(structure);
    let keep = leaving_branch.iter().copied().collect::<BTreeSet<_>>();
    let remove = (0..structure.atoms.len())
        .filter(|atom| !keep.contains(atom))
        .collect::<Vec<_>>();
    let leaving_mapping = leaving_editor.remove_atoms(&remove)?;
    let leaving_atom = mapped_atom(&leaving_mapping, leaving_atom, leaving_label)?;
    leaving_editor.add_atom(leaving_atom, "H", 0.0, 1.0)?;
    let leaving_product = resolver.resolve(leaving_editor.finish()?)?;

    Ok(
        Reaction::builder(generated_site_reaction_id(id_prefix, participant))
            .reactant(substance_id, 1, 1)
            .reactant("destroy:water", 1, 1)
            .catalyst_order("destroy:proton", 1)
            .product(acid, 1)
            .product(leaving_product, 1)
            .condition(
                ReactionCondition::new(
                    "acyl leaving-group hydrolysis requires acidic, water-rich conditions",
                )
                .acidity(AcidityCondition::Acidic)
                .min_water_activity(0.35),
            )
            .selectivity_profile(
                SelectivityProfile::new(
                    ReactionType::AcylSubstitution,
                    SiteDescriptorBuilder::carboxylic_acid(),
                )
                .never_suppress(),
            )
            .activation_energy_kj_per_mol(42.0)
            .build(),
    )
}

impl<'a> AcylDonor<'a> {
    fn from_acyl_chloride(site: &'a AcylChlorideSite<'a>) -> Self {
        Self {
            participant: &site.participant,
            structure: site.participant.structure,
            acyl_carbon: site.carbon,
            leaving_atoms: vec![site.chlorine],
            leaving_attachment: Some(site.chlorine),
            leaving_product_static: Some("destroy:hydrochloric_acid"),
            leaving_product_prefix: "acyl_chloride_leaving_group",
            kind: AcylDonorKind::AcylChloride,
        }
    }

    fn from_anhydride(site: &'a AcidAnhydrideSite<'a>) -> [Self; 2] {
        [
            Self {
                participant: &site.participant,
                structure: site.participant.structure,
                acyl_carbon: site.carbon_a,
                leaving_atoms: anhydride_leaving_atoms(site, site.carbon_b, site.oxygen_b),
                leaving_attachment: Some(site.bridge_oxygen),
                leaving_product_static: None,
                leaving_product_prefix: "anhydride_leaving_acid",
                kind: AcylDonorKind::AcidAnhydride,
            },
            Self {
                participant: &site.participant,
                structure: site.participant.structure,
                acyl_carbon: site.carbon_b,
                leaving_atoms: anhydride_leaving_atoms(site, site.carbon_a, site.oxygen_a),
                leaving_attachment: Some(site.bridge_oxygen),
                leaving_product_static: None,
                leaving_product_prefix: "anhydride_leaving_acid",
                kind: AcylDonorKind::AcidAnhydride,
            },
        ]
    }

    fn leaving_product(
        &self,
        resolver: &mut DerivedSubstanceResolver,
    ) -> ChemistryResult<SubstanceId> {
        if let Some(static_id) = self.leaving_product_static {
            return Ok(SubstanceId::from(static_id));
        }
        let leaving_attachment =
            self.leaving_attachment
                .ok_or_else(|| ChemistryError::InvalidReaction {
                    reaction_id: generated_site_reaction_id(
                        self.leaving_product_prefix,
                        self.participant,
                    ),
                    reason: "dynamic acyl leaving group has no attachment atom".to_string(),
                })?;
        let mut editor = MolecularEditor::new(self.structure);
        let keep = self.leaving_atoms.iter().copied().collect::<BTreeSet<_>>();
        let remove = (0..self.structure.atoms.len())
            .filter(|atom| !keep.contains(atom))
            .collect::<Vec<_>>();
        let mapping = editor.remove_atoms(&remove)?;
        let attachment = mapped_atom(&mapping, leaving_attachment, "acyl leaving group atom")?;
        editor.add_atom(attachment, "H", 0.0, 1.0)?;
        resolver.resolve(editor.finish()?)
    }
}

impl<'a> AcylNucleophile<'a> {
    fn from_alcohol(site: &'a AlcoholSite<'a>) -> Self {
        Self {
            participant: &site.participant,
            structure: site.participant.structure,
            attack_atom: site.oxygen,
            proton: site.hydrogen,
            kind: AcylNucleophileKind::Alcohol,
        }
    }

    fn from_amine(site: &'a AmineSite<'a>) -> Self {
        Self {
            participant: &site.participant,
            structure: site.participant.structure,
            attack_atom: site.nitrogen,
            proton: site.hydrogens[0],
            kind: AcylNucleophileKind::Amine,
        }
    }

    fn from_thiol(site: &'a ThiolSite<'a>) -> Self {
        Self {
            participant: &site.participant,
            structure: site.participant.structure,
            attack_atom: site.sulfur,
            proton: site.hydrogens[0],
            kind: AcylNucleophileKind::Thiol,
        }
    }
}

fn anhydride_leaving_atoms(
    site: &AcidAnhydrideSite<'_>,
    leaving_carbon: usize,
    leaving_oxygen: usize,
) -> Vec<usize> {
    let mut atoms =
        acyl_substituent_branch(site.participant.structure, leaving_carbon, leaving_oxygen);
    atoms.push(leaving_carbon);
    atoms.push(leaving_oxygen);
    atoms.push(site.bridge_oxygen);
    atoms.sort_unstable();
    atoms.dedup();
    atoms
}

fn acyl_substituent_branch(
    structure: &MolecularStructure,
    acyl_carbon: usize,
    carbonyl_oxygen: usize,
) -> Vec<usize> {
    let mut branch = Vec::new();
    for (neighbor, order) in structure.neighbors(acyl_carbon) {
        if neighbor == carbonyl_oxygen
            || structure.atoms[neighbor].element == "O"
            || !crate::chemistry::molecule::bond_order_matches(order, 1.0)
        {
            continue;
        }
        let mut stack = vec![neighbor];
        let mut visited = vec![false; structure.atoms.len()];
        visited[acyl_carbon] = true;
        while let Some(atom) = stack.pop() {
            if visited[atom] {
                continue;
            }
            visited[atom] = true;
            branch.push(atom);
            for (next, _) in structure.neighbors(atom) {
                if !visited[next] {
                    stack.push(next);
                }
            }
        }
    }
    branch
}

fn acyl_transfer_id_prefix(
    donor: &AcylDonor<'_>,
    nucleophile: &AcylNucleophile<'_>,
) -> &'static str {
    match (donor.kind, nucleophile.kind) {
        (AcylDonorKind::AcylChloride, AcylNucleophileKind::Alcohol) => {
            "acyl_chloride_esterification"
        }
        (AcylDonorKind::AcylChloride, AcylNucleophileKind::Amine) => "acyl_chloride_amidation",
        (AcylDonorKind::AcylChloride, AcylNucleophileKind::Thiol) => {
            "acyl_chloride_thioesterification"
        }
        (AcylDonorKind::AcidAnhydride, AcylNucleophileKind::Alcohol) => {
            "anhydride_alcohol_acylation"
        }
        (AcylDonorKind::AcidAnhydride, AcylNucleophileKind::Amine) => "anhydride_amine_acylation",
        (AcylDonorKind::AcidAnhydride, AcylNucleophileKind::Thiol) => "anhydride_thiol_acylation",
    }
}

fn donor_hydrolysis_id_prefix(donor: &AcylDonor<'_>) -> &'static str {
    match donor.kind {
        AcylDonorKind::AcylChloride => "acyl_chloride_hydrolysis",
        AcylDonorKind::AcidAnhydride => "anhydride_hydrolysis",
    }
}

fn acyl_transfer_activation_energy(
    donor: &AcylDonor<'_>,
    nucleophile: &AcylNucleophile<'_>,
) -> f64 {
    let donor_term = match donor.kind {
        AcylDonorKind::AcylChloride => 10.0,
        AcylDonorKind::AcidAnhydride => 18.0,
    };
    let nucleophile_term = match nucleophile.kind {
        AcylNucleophileKind::Amine => -2.0,
        AcylNucleophileKind::Alcohol => 4.0,
        AcylNucleophileKind::Thiol => 2.0,
    };
    donor_term + nucleophile_term
}

fn nucleophile_descriptor(nucleophile: &AcylNucleophile<'_>) -> SiteDescriptor {
    match nucleophile.kind {
        AcylNucleophileKind::Alcohol => SiteDescriptorBuilder::build(
            crate::chemistry::reactive_site::ReactiveSiteKind::Alcohol,
            degree_from_usize(
                nucleophile
                    .participant
                    .site
                    .substitution_degree
                    .unwrap_or(1),
            ),
            0,
            0,
            0,
            true,
            false,
            false,
        ),
        AcylNucleophileKind::Amine => SiteDescriptorBuilder::build(
            crate::chemistry::reactive_site::ReactiveSiteKind::NonTertiaryAmine,
            SubstitutionDegree::Primary,
            1,
            0,
            0,
            false,
            false,
            false,
        ),
        AcylNucleophileKind::Thiol => SiteDescriptorBuilder::build(
            crate::chemistry::reactive_site::ReactiveSiteKind::Thiol,
            SubstitutionDegree::Primary,
            0,
            0,
            0,
            false,
            false,
            false,
        ),
    }
}

fn degree_from_usize(degree: usize) -> SubstitutionDegree {
    match degree {
        0 | 1 => SubstitutionDegree::Primary,
        2 => SubstitutionDegree::Secondary,
        _ => SubstitutionDegree::Tertiary,
    }
}

fn ester_alkoxy_branch(
    structure: &MolecularStructure,
    alkoxy_oxygen: usize,
    carbonyl_carbon: usize,
) -> ChemistryResult<Vec<usize>> {
    leaving_branch(
        structure,
        alkoxy_oxygen,
        carbonyl_carbon,
        "ester_hydrolysis",
    )?
    .ok_or_else(|| ChemistryError::InvalidReaction {
        reaction_id: "ester_hydrolysis".to_string(),
        reason: "ester alkoxy oxygen is part of a ring; use lactone ring-opening chemistry"
            .to_string(),
    })
}

fn amide_leaving_branch(
    structure: &MolecularStructure,
    nitrogen: usize,
    carbonyl_carbon: usize,
) -> ChemistryResult<Option<Vec<usize>>> {
    leaving_branch(structure, nitrogen, carbonyl_carbon, "amide_hydrolysis")
}

fn leaving_branch(
    structure: &MolecularStructure,
    leaving_atom: usize,
    carbonyl_carbon: usize,
    reaction_id: &str,
) -> ChemistryResult<Option<Vec<usize>>> {
    if leaving_atom >= structure.atoms.len() || carbonyl_carbon >= structure.atoms.len() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: "acyl leaving group references an atom outside the structure".to_string(),
        });
    }
    let mut stack = vec![leaving_atom];
    let mut visited = vec![false; structure.atoms.len()];
    visited[carbonyl_carbon] = true;
    while let Some(atom) = stack.pop() {
        if visited[atom] {
            continue;
        }
        visited[atom] = true;
        for (neighbor, _) in structure.neighbors(atom) {
            if !visited[neighbor] {
                stack.push(neighbor);
            }
        }
    }
    let branch = visited
        .into_iter()
        .enumerate()
        .filter_map(|(atom, seen)| (seen && atom != carbonyl_carbon).then_some(atom))
        .collect::<Vec<_>>();
    if !branch.contains(&leaving_atom) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: reaction_id.to_string(),
            reason: "acyl leaving branch does not contain the leaving atom".to_string(),
        });
    }
    if structure
        .neighbors(carbonyl_carbon)
        .into_iter()
        .any(|(neighbor, _)| neighbor != leaving_atom && branch.contains(&neighbor))
    {
        return Ok(None);
    }
    Ok(Some(branch))
}
