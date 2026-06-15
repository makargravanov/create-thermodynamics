//! Generators for protecting group reactions
//!
//! This module implements chemical reactions for protection and deprotection
//! of functional groups used in organic synthesis.
//!
//! NOTE: This module currently covers TMS alcohol protection, acetals, Boc, and
//! Cbz amine protection.

use super::super::centers::{
    AcetalCenter, AlcoholSite, AmineSite, BocCarbamateCenter, CbzCarbamateCenter,
    ChloroformateSite, SilylEtherCenter,
};
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};
use crate::chemistry::substance::SubstanceId;

/// Add TMS (trimethylsilyl) protecting group to an alcohol
///
/// Reaction: R-OH + Me3SiCl -> R-OSiMe3 + HCl
/// Note: In real synthesis, a base is used to scavenge HCl, but for simplicity
/// we model this as the direct reaction.
pub(crate) fn generate_alcohol_silyl_protection(
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = alcohol_site.participant.substance;
    let structure = alcohol_site.participant.structure;
    let oxygen = alcohol_site.oxygen;
    let hydrogen = alcohol_site.hydrogen;

    let mut editor = MolecularEditor::new(structure);
    // Remove the hydroxyl hydrogen
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let oxygen = mapped_atom(&mapping, oxygen, "alcohol oxygen")?;

    // Add silicon atom attached to oxygen
    let silicon = editor.add_atom(oxygen, "Si", 0.0, 1.0)?;

    // Add three methyl groups to silicon
    for _ in 0..3 {
        let methyl_carbon = editor.add_atom(silicon, "C", 0.0, 1.0)?;
        editor.add_atom(methyl_carbon, "H", 0.0, 1.0)?;
        editor.add_atom(methyl_carbon, "H", 0.0, 1.0)?;
        editor.add_atom(methyl_carbon, "H", 0.0, 1.0)?;
    }

    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "alcohol_silyl_protection",
        &alcohol_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:trimethylsilyl_chloride", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .condition(
        ReactionCondition::new("silyl protection requires dry conditions").max_water_activity(0.1),
    )
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::SilylEtherFormation,
            SiteDescriptorBuilder::from_alcohol_site(alcohol_site),
        )
        .never_suppress(),
    )
    .activation_energy_kj_per_mol(15.0)
    .build())
}

pub(crate) fn generate_alcohol_chloroformate_formation(
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let alcohol = alcohol_site.participant.substance;
    let phosgene_id = SubstanceId::from("destroy:phosgene");
    let phosgene =
        resolver
            .known_structure(&phosgene_id)
            .ok_or_else(|| ChemistryError::InvalidReaction {
                reaction_id: generated_site_reaction_id(
                    "alcohol_chloroformate_formation",
                    &alcohol_site.participant,
                ),
                reason: "required reagent 'destroy:phosgene' has no known molecular structure"
                    .to_string(),
            })?;
    let (phosgene_carbon, _phosgene_oxygen, leaving_chlorine) =
        phosgene_chloroformate_atoms(phosgene)?;

    let mut alcohol_editor = MolecularEditor::new(alcohol_site.participant.structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol_site.hydrogen])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol_site.oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let mut phosgene_editor = MolecularEditor::new(phosgene);
    let phosgene_mapping = phosgene_editor.remove_atoms(&[leaving_chlorine])?;
    let phosgene_carbon = mapped_atom(
        &phosgene_mapping,
        phosgene_carbon,
        "phosgene carbonyl carbon",
    )?;
    let phosgene_fragment = phosgene_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &alcohol_fragment,
        alcohol_oxygen,
        &phosgene_fragment,
        phosgene_carbon,
        1.0,
    )?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "alcohol_chloroformate_formation",
        &alcohol_site.participant,
    ))
    .reactant(alcohol.id.clone(), 1, 1)
    .reactant(phosgene_id, 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .condition(
        ReactionCondition::new("chloroformate formation requires dry alcoholysis of phosgene")
            .max_water_activity(0.05),
    )
    .activation_energy_kj_per_mol(18.0)
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::AcylSubstitution,
            SiteDescriptorBuilder::from_alcohol_site(alcohol_site),
        )
        .never_suppress(),
    )
    .build())
}

fn phosgene_chloroformate_atoms(
    structure: &MolecularStructure,
) -> ChemistryResult<(usize, usize, usize)> {
    for carbon in 0..structure.atoms.len() {
        if structure.atoms[carbon].element != "C" {
            continue;
        }
        let mut oxygen = None;
        let mut chlorines = Vec::new();
        for (neighbor, order) in structure.neighbors(carbon) {
            if structure.atoms[neighbor].element == "O" && bond_order_matches(order, 2.0) {
                oxygen = Some(neighbor);
            }
            if structure.atoms[neighbor].element == "Cl" && bond_order_matches(order, 1.0) {
                chlorines.push(neighbor);
            }
        }
        if let (Some(oxygen), Some(chlorine)) = (oxygen, chlorines.first().copied()) {
            if chlorines.len() == 2 {
                return Ok((carbon, oxygen, chlorine));
            }
        }
    }
    Err(ChemistryError::InvalidReaction {
        reaction_id: "alcohol_chloroformate_formation".to_string(),
        reason: "phosgene structure must contain C(=O)Cl2".to_string(),
    })
}

pub(crate) fn generate_chloroformate_alcohol_carbonate_formation(
    chloroformate: &ChloroformateSite<'_>,
    alcohol: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    generate_chloroformate_transfer(
        chloroformate,
        alcohol.participant.structure,
        alcohol.oxygen,
        alcohol.hydrogen,
        &alcohol.participant,
        "chloroformate_alcohol_carbonate_formation",
        ReactionType::AcylSubstitution,
        SiteDescriptorBuilder::from_alcohol_site(alcohol),
        resolver,
    )
}

pub(crate) fn generate_chloroformate_amine_carbamate_formation(
    chloroformate: &ChloroformateSite<'_>,
    amine: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let Some(&hydrogen) = amine.hydrogens.first() else {
        return Ok(None);
    };
    Ok(Some(generate_chloroformate_transfer(
        chloroformate,
        amine.participant.structure,
        amine.nitrogen,
        hydrogen,
        &amine.participant,
        "chloroformate_amine_carbamate_formation",
        ReactionType::AcylSubstitution,
        SiteDescriptorBuilder::from_amine_site(amine),
        resolver,
    )?))
}

fn generate_chloroformate_transfer(
    chloroformate: &ChloroformateSite<'_>,
    nucleophile_structure: &MolecularStructure,
    nucleophile_atom: usize,
    nucleophile_hydrogen: usize,
    nucleophile_participant: &crate::chemistry::organic::space::SiteParticipant<'_>,
    prefix: &'static str,
    reaction_type: ReactionType,
    nucleophile_descriptor: crate::chemistry::selectivity::types::SiteDescriptor,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let mut donor_editor = MolecularEditor::new(chloroformate.participant.structure);
    let donor_mapping = donor_editor.remove_atoms(&[chloroformate.chlorine])?;
    let carbon = mapped_atom(
        &donor_mapping,
        chloroformate.carbon,
        "chloroformate carbonyl carbon",
    )?;
    let donor_fragment = donor_editor.finish()?;

    let mut nucleophile_editor = MolecularEditor::new(nucleophile_structure);
    let nucleophile_mapping = nucleophile_editor.remove_atoms(&[nucleophile_hydrogen])?;
    let nucleophile_atom = mapped_atom(
        &nucleophile_mapping,
        nucleophile_atom,
        "chloroformate nucleophile atom",
    )?;
    let nucleophile_fragment = nucleophile_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &donor_fragment,
        carbon,
        &nucleophile_fragment,
        nucleophile_atom,
        1.0,
    )?)?;

    Ok(Reaction::builder(generated_pair_site_reaction_id(
        prefix,
        &chloroformate.participant,
        nucleophile_participant,
    ))
    .reactant(chloroformate.participant.substance.id.clone(), 1, 1)
    .reactant(nucleophile_participant.substance.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .condition(
        ReactionCondition::new("chloroformate transfer is suppressed by wet media")
            .max_water_activity(0.2),
    )
    .activation_energy_kj_per_mol(16.0)
    .selectivity_profile(
        SelectivityProfile::new(
            reaction_type,
            SiteDescriptorBuilder::build(
                crate::chemistry::reactive_site::ReactiveSiteKind::Chloroformate,
                crate::chemistry::selectivity::types::SubstitutionDegree::Primary,
                0,
                0,
                0,
                false,
                false,
                false,
            ),
        )
        .with_secondary_site(nucleophile_descriptor)
        .never_suppress(),
    )
    .build())
}

/// Remove TMS protecting group from a silyl ether
///
/// Reaction: R-OSiMe3 + F- + H+ -> R-OH + Me3SiF
/// Fluoride cleavage is the standard method for removing TMS groups.
pub(crate) fn generate_silyl_ether_deprotection(
    silyl_site: &SilylEtherCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let substance = silyl_site.participant.substance;
    let structure = silyl_site.participant.structure;
    let oxygen = silyl_site.oxygen;
    let silicon = silyl_site.silicon;

    let Some(mut atoms_to_remove) = trimethylsilyl_fragment_atoms(structure, oxygen, silicon)
    else {
        return Ok(None);
    };
    atoms_to_remove.push(silicon);
    atoms_to_remove.sort_unstable();
    atoms_to_remove.dedup();

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&atoms_to_remove)?;
    let oxygen = mapped_atom(&mapping, oxygen, "silyl ether oxygen")?;

    editor.add_atom(oxygen, "H", 0.0, 1.0)?;

    let product = resolver.resolve(editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "silyl_ether_deprotection",
            &silyl_site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant("destroy:fluoride", 1, 1) // F- from TBAF or similar
        .reactant("destroy:proton", 1, 1)
        .product(product, 1)
        .product("destroy:trimethylsilyl_fluoride", 1)
        .condition(ReactionCondition::new(
            "fluoride deprotection requires fluoride source",
        ))
        .selectivity_profile(
            SelectivityProfile::new(
                ReactionType::SilylEtherCleavage,
                SiteDescriptorBuilder::silyl_ether(),
            )
            .never_suppress(),
        )
        .activation_energy_kj_per_mol(20.0)
        .build(),
    ))
}

/// Hydrolyze an acetal or ketal back to the carbonyl compound and concrete alcohols.
///
/// Reaction: R2C(OR')2 + H2O -> R2C=O + 2 R'OH
pub(crate) fn generate_acetal_deprotection(
    acetal_site: &AcetalCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = acetal_site.participant.substance;
    let structure = acetal_site.participant.structure;
    let acetal_carbon = acetal_site.acetal_carbon;
    let oxygen_a = acetal_site.oxygen_a;
    let oxygen_b = acetal_site.oxygen_b;

    let branch_a = branch_atoms(structure, oxygen_a, acetal_carbon)?;
    let branch_b = branch_atoms(structure, oxygen_b, acetal_carbon)?;
    if branch_a.iter().any(|atom| branch_b.contains(atom)) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("acetal_deprotection", &acetal_site.participant),
            reason: "acetal oxygens belong to the same branch; cyclic acetals need a dedicated diol product path".to_string(),
        });
    }

    let alcohol_a = alcohol_fragment_product(structure, &branch_a, oxygen_a, resolver)?;
    let alcohol_b = alcohol_fragment_product(structure, &branch_b, oxygen_b, resolver)?;

    let mut atoms_to_remove = branch_a;
    atoms_to_remove.extend(branch_b);
    atoms_to_remove.sort_unstable();
    atoms_to_remove.dedup();

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&atoms_to_remove)?;
    let carbonyl_carbon = mapped_atom(&mapping, acetal_carbon, "acetal carbon")?;
    editor.add_atom(carbonyl_carbon, "O", 0.0, 2.0)?;
    let carbonyl = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "acetal_deprotection",
        &acetal_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(carbonyl, 1)
    .product(alcohol_a, 1)
    .product(alcohol_b, 1)
    .condition(
        ReactionCondition::new("acetal hydrolysis requires acidic, water-rich conditions")
            .acidity(AcidityCondition::Acidic)
            .min_water_activity(0.35),
    )
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::AcetalHydrolysis,
            SiteDescriptorBuilder::acetal(),
        )
        .never_suppress(),
    )
    .activation_energy_kj_per_mol(25.0)
    .build())
}

/// Protect an amine as a Boc carbamate.
///
/// Reaction: R2NH + Boc2O -> R2NCO2tBu + tBuOH + CO2
pub(crate) fn generate_amine_boc_protection(
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = amine_site.participant.substance;
    let structure = amine_site.participant.structure;
    let nitrogen = amine_site.nitrogen;
    let hydrogen =
        *amine_site
            .hydrogens
            .first()
            .ok_or_else(|| ChemistryError::InvalidReaction {
                reaction_id: generated_site_reaction_id(
                    "amine_boc_protection",
                    &amine_site.participant,
                ),
                reason: "Boc protection requires an explicit N-H bond".to_string(),
            })?;

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "amine nitrogen")?;
    add_boc_group(&mut editor, nitrogen)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "amine_boc_protection",
        &amine_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:di_tert_butyl_dicarbonate", 1, 1)
    .product(product, 1)
    .product("destroy:tert_butanol", 1)
    .product("destroy:carbon_dioxide", 1)
    .condition(
        ReactionCondition::new("Boc protection prefers basic, water-poor conditions")
            .acidity(AcidityCondition::Basic)
            .max_water_activity(0.2),
    )
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::CarbamateFormation,
            SiteDescriptorBuilder::from_amine_site(amine_site),
        )
        .never_suppress(),
    )
    .activation_energy_kj_per_mol(20.0)
    .build())
}

/// Protect an amine as a Cbz carbamate.
///
/// Reaction: R2NH + benzyl chloroformate -> R2NCO2CH2Ph + HCl
pub(crate) fn generate_amine_cbz_protection(
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = amine_site.participant.substance;
    let structure = amine_site.participant.structure;
    let nitrogen = amine_site.nitrogen;
    let hydrogen =
        *amine_site
            .hydrogens
            .first()
            .ok_or_else(|| ChemistryError::InvalidReaction {
                reaction_id: generated_site_reaction_id(
                    "amine_cbz_protection",
                    &amine_site.participant,
                ),
                reason: "Cbz protection requires an explicit N-H bond".to_string(),
            })?;

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "amine nitrogen")?;
    add_cbz_group(&mut editor, nitrogen)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "amine_cbz_protection",
        &amine_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:benzyl_chloroformate", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .condition(
        ReactionCondition::new("Cbz protection prefers basic, water-poor conditions")
            .acidity(AcidityCondition::Basic)
            .max_water_activity(0.2),
    )
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::CarbamateFormation,
            SiteDescriptorBuilder::from_amine_site(amine_site),
        )
        .never_suppress(),
    )
    .activation_energy_kj_per_mol(20.0)
    .build())
}

/// Hydrolyze an acid-labile Boc carbamate back to the amine.
///
/// Reaction: R2NCO2tBu + H2O -> R2NH + tBuOH + CO2
pub(crate) fn generate_boc_deprotection(
    boc_site: &BocCarbamateCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = boc_site.participant.substance;
    let structure = boc_site.participant.structure;
    let mut atoms_to_remove = vec![
        boc_site.carbonyl_carbon,
        boc_site.carbonyl_oxygen,
        boc_site.alkoxy_oxygen,
        boc_site.tert_butyl_carbon,
    ];
    for (methyl, order) in structure.neighbors(boc_site.tert_butyl_carbon) {
        if structure.atoms[methyl].element == "C" && bond_order_matches(order, 1.0) {
            atoms_to_remove.push(methyl);
            for (hydrogen, h_order) in structure.neighbors(methyl) {
                if structure.atoms[hydrogen].element == "H" && bond_order_matches(h_order, 1.0) {
                    atoms_to_remove.push(hydrogen);
                }
            }
        }
    }
    atoms_to_remove.sort_unstable();
    atoms_to_remove.dedup();

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&atoms_to_remove)?;
    let nitrogen = mapped_atom(&mapping, boc_site.nitrogen, "Boc carbamate nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let amine = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "boc_deprotection",
        &boc_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(amine, 1)
    .product("destroy:tert_butanol", 1)
    .product("destroy:carbon_dioxide", 1)
    .condition(
        ReactionCondition::new("Boc deprotection requires acidic, water-rich conditions")
            .acidity(AcidityCondition::Acidic)
            .min_water_activity(0.2),
    )
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::CarbamateCleavage,
            SiteDescriptorBuilder::boc_carbamate(),
        )
        .never_suppress(),
    )
    .activation_energy_kj_per_mol(18.0)
    .build())
}

/// Hydrogenolyze a Cbz carbamate back to the amine.
///
/// Reaction: R2NCO2CH2Ph + H2 -> R2NH + toluene + CO2
pub(crate) fn generate_cbz_deprotection(
    cbz_site: &CbzCarbamateCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = cbz_site.participant.substance;
    let structure = cbz_site.participant.structure;
    let mut atoms_to_remove =
        branch_atoms(structure, cbz_site.alkoxy_oxygen, cbz_site.carbonyl_carbon)?;
    atoms_to_remove.push(cbz_site.carbonyl_carbon);
    atoms_to_remove.push(cbz_site.carbonyl_oxygen);
    atoms_to_remove.sort_unstable();
    atoms_to_remove.dedup();

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&atoms_to_remove)?;
    let nitrogen = mapped_atom(&mapping, cbz_site.nitrogen, "Cbz carbamate nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let amine = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "cbz_deprotection",
        &cbz_site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen", 1, 1)
    .product(amine, 1)
    .product("destroy:toluene", 1)
    .product("destroy:carbon_dioxide", 1)
    .external_catalyst("forge:dusts/palladium", 1.0)
    .condition(ReactionCondition::new(
        "Cbz deprotection requires hydrogen and palladium catalyst",
    ))
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::CarbamateCleavage,
            SiteDescriptorBuilder::cbz_carbamate(),
        )
        .never_suppress(),
    )
    .activation_energy_kj_per_mol(12.0)
    .build())
}

fn add_boc_group(editor: &mut MolecularEditor, nitrogen: usize) -> ChemistryResult<()> {
    let carbonyl_carbon = editor.add_atom(nitrogen, "C", 0.0, 1.0)?;
    editor.add_atom(carbonyl_carbon, "O", 0.0, 2.0)?;
    let alkoxy_oxygen = editor.add_atom(carbonyl_carbon, "O", 0.0, 1.0)?;
    let tert_butyl_carbon = editor.add_atom(alkoxy_oxygen, "C", 0.0, 1.0)?;
    for _ in 0..3 {
        let methyl = editor.add_atom(tert_butyl_carbon, "C", 0.0, 1.0)?;
        for _ in 0..3 {
            editor.add_atom(methyl, "H", 0.0, 1.0)?;
        }
    }
    Ok(())
}

fn add_cbz_group(editor: &mut MolecularEditor, nitrogen: usize) -> ChemistryResult<()> {
    let carbonyl_carbon = editor.add_atom(nitrogen, "C", 0.0, 1.0)?;
    editor.add_atom(carbonyl_carbon, "O", 0.0, 2.0)?;
    let alkoxy_oxygen = editor.add_atom(carbonyl_carbon, "O", 0.0, 1.0)?;
    let benzyl_carbon = editor.add_atom(alkoxy_oxygen, "C", 0.0, 1.0)?;
    editor.add_atom(benzyl_carbon, "H", 0.0, 1.0)?;
    editor.add_atom(benzyl_carbon, "H", 0.0, 1.0)?;
    let mut ring = Vec::with_capacity(6);
    ring.push(editor.add_atom(benzyl_carbon, "C", 0.0, 1.0)?);
    for index in 1..6 {
        let carbon = editor.add_atom(ring[index - 1], "C", 0.0, 1.5)?;
        ring.push(carbon);
    }
    editor.add_bond(ring[5], ring[0], 1.5)?;
    for carbon in ring.iter().skip(1) {
        editor.add_atom(*carbon, "H", 0.0, 1.0)?;
    }
    Ok(())
}

fn trimethylsilyl_fragment_atoms(
    structure: &MolecularStructure,
    protected_oxygen: usize,
    silicon: usize,
) -> Option<Vec<usize>> {
    let substituents = structure
        .neighbors(silicon)
        .into_iter()
        .filter(|(neighbor, order)| {
            *neighbor != protected_oxygen && bond_order_matches(*order, 1.0)
        })
        .collect::<Vec<_>>();
    if substituents.len() != 3 {
        return None;
    }

    let mut atoms = Vec::new();
    for (carbon, _) in substituents {
        if structure.atoms[carbon].element != "C" {
            return None;
        }
        let mut hydrogens = Vec::new();
        for (neighbor, order) in structure.neighbors(carbon) {
            if neighbor == silicon {
                continue;
            }
            if structure.atoms[neighbor].element == "H" && bond_order_matches(order, 1.0) {
                hydrogens.push(neighbor);
            } else {
                return None;
            }
        }
        if hydrogens.len() != 3 {
            return None;
        }
        atoms.push(carbon);
        atoms.extend(hydrogens);
    }
    Some(atoms)
}

fn branch_atoms(
    structure: &MolecularStructure,
    start: usize,
    blocked: usize,
) -> ChemistryResult<Vec<usize>> {
    if start >= structure.atoms.len() || blocked >= structure.atoms.len() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "acetal_deprotection".to_string(),
            reason: "acetal branch references an atom outside the structure".to_string(),
        });
    }
    let mut stack = vec![start];
    let mut visited = vec![false; structure.atoms.len()];
    visited[blocked] = true;
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
    Ok(visited
        .into_iter()
        .enumerate()
        .filter_map(|(atom, seen)| (seen && atom != blocked).then_some(atom))
        .collect())
}

fn alcohol_fragment_product(
    structure: &MolecularStructure,
    atoms: &[usize],
    oxygen: usize,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<crate::chemistry::substance::SubstanceId> {
    let mut editor = MolecularEditor::new(structure);
    let keep = atoms
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let remove = (0..structure.atoms.len())
        .filter(|atom| !keep.contains(atom))
        .collect::<Vec<_>>();
    let mapping = editor.remove_atoms(&remove)?;
    let oxygen = mapped_atom(&mapping, oxygen, "acetal alcohol oxygen")?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    resolver.resolve(editor.finish()?)
}
