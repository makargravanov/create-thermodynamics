use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityContext, SelectivityProfile},
};

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
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let chlorine = site.chlorine;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[chlorine])?;
    let carbon = mapped_atom(&mapping, carbon, "acyl chloride carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "acyl_chloride_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
}

pub(crate) fn generate_acyl_chloride_esterification(
    acyl_chloride_site: &AcylChlorideSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let acyl_chloride = acyl_chloride_site.participant.substance;
    let acyl_chloride_structure = acyl_chloride_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let alcohol_structure = alcohol_site.participant.structure;
    let acyl_carbon = acyl_chloride_site.carbon;
    let chlorine = acyl_chloride_site.chlorine;
    let alcohol_oxygen = alcohol_site.oxygen;
    let alcohol_proton = alcohol_site.hydrogen;
    let mut acyl_editor = MolecularEditor::new(acyl_chloride_structure);
    let acyl_mapping = acyl_editor.remove_atoms(&[chlorine])?;
    let acyl_carbon = mapped_atom(&acyl_mapping, acyl_carbon, "acyl chloride carbon")?;
    let acyl_fragment = acyl_editor.finish()?;

    let mut alcohol_editor = MolecularEditor::new(alcohol_structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol_proton])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol_oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &acyl_fragment,
        acyl_carbon,
        &alcohol_fragment,
        alcohol_oxygen,
        1.0,
    )?)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "acyl_chloride_esterification",
        &acyl_chloride_site.participant,
        &alcohol_site.participant,
    ))
    .reactant(acyl_chloride.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
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
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let nitrogen = site.nitrogen;
    let hydrogens = &site.nitrogen_hydrogens;
    if hydrogens.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("amide_hydrolysis", &site.participant),
            reason: "unsubstituted amide must have exactly two explicit nitrogen hydrogens"
                .to_string(),
        });
    }
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen, hydrogens[0], hydrogens[1]])?;
    let carbon = mapped_atom(&mapping, carbon, "amide carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "amide_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(product, 1)
    .product("destroy:ammonia", 1)
    .build())
}

fn ester_alkoxy_branch(
    structure: &MolecularStructure,
    alkoxy_oxygen: usize,
    carbonyl_carbon: usize,
) -> ChemistryResult<Vec<usize>> {
    if alkoxy_oxygen >= structure.atoms.len() || carbonyl_carbon >= structure.atoms.len() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "ester_hydrolysis".to_string(),
            reason: "ester site references an atom outside the structure".to_string(),
        });
    }
    let mut stack = vec![alkoxy_oxygen];
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
    if !branch.contains(&alkoxy_oxygen) {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "ester_hydrolysis".to_string(),
            reason: "ester alkoxy branch does not contain the alkoxy oxygen".to_string(),
        });
    }
    Ok(branch)
}
