use super::super::centers::{CarbonylSite, OximeSite};
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};

pub(crate) fn generate_baeyer_villiger_rearrangements(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    if !site.is_ketone {
        return Ok(Vec::new());
    }
    let mut reactions = Vec::new();
    for migratory_carbon in migratory_carbon_neighbors(site.participant.structure, site.carbon) {
        let mut editor = MolecularEditor::new(site.participant.structure);
        editor.insert_bridging_atom(site.carbon, migratory_carbon, "O", 0.0)?;
        let product = resolver.resolve(editor.finish()?)?;
        reactions.push(
            Reaction::builder(format!(
                "{}/migrates_{}",
                generated_site_reaction_id("baeyer_villiger_rearrangement", &site.participant),
                migratory_carbon
            ))
            .reactant(site.participant.substance.id.clone(), 1, 1)
            .reactant("destroy:hydrogen_peroxide", 1, 1)
            .product(product, 1)
            .product("destroy:water", 1)
            .condition(
                ReactionCondition::new("Baeyer-Villiger oxygen insertion requires acidic oxidizing medium")
                    .acidity(AcidityCondition::Acidic),
            )
            .reactant_phase_access(site.participant.substance.id.clone(), [MixturePhase::Organic])
            .activation_energy_kj_per_mol(31.0 + migration_penalty(site.participant.structure, migratory_carbon))
            .selectivity_profile(SelectivityProfile::new(
                ReactionType::SkeletalRearrangement,
                SiteDescriptorBuilder::from_carbonyl_site(site),
            ))
            .build(),
        );
    }
    Ok(reactions)
}

pub(crate) fn generate_beckmann_rearrangements(
    site: &OximeSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    let mut reactions = Vec::new();
    for migratory_carbon in migratory_carbon_neighbors(site.participant.structure, site.carbon) {
        let mut editor = MolecularEditor::new(site.participant.structure);
        let mapping = editor.remove_atoms(&[site.oxygen, site.hydrogen])?;
        let carbon = mapped_atom(&mapping, site.carbon, "oxime carbon")?;
        let nitrogen = mapped_atom(&mapping, site.nitrogen, "oxime nitrogen")?;
        let migratory_carbon = mapped_atom(&mapping, migratory_carbon, "migrating carbon")?;
        editor.set_bond_order(carbon, nitrogen, 1.0)?;
        editor.move_bond_attachment(carbon, migratory_carbon, nitrogen, 1.0)?;
        editor.add_atom(carbon, "O", 0.0, 2.0)?;
        editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
        let product = resolver.resolve(editor.finish()?)?;
        reactions.push(
            Reaction::builder(format!(
                "{}/migrates_{}",
                generated_site_reaction_id("beckmann_rearrangement", &site.participant),
                migratory_carbon
            ))
            .reactant(site.participant.substance.id.clone(), 1, 1)
            .product(product, 1)
            .condition(
                ReactionCondition::new("Beckmann rearrangement requires acid activation of the oxime")
                    .acidity(AcidityCondition::Acidic),
            )
            .reactant_phase_access(site.participant.substance.id.clone(), [MixturePhase::Organic])
            .activation_energy_kj_per_mol(34.0 + migration_penalty(site.participant.structure, migratory_carbon))
            .selectivity_profile(SelectivityProfile::new(
                ReactionType::SkeletalRearrangement,
                SiteDescriptorBuilder::build(
                    site.participant.site.kind.clone(),
                    crate::chemistry::selectivity::types::SubstitutionDegree::Secondary,
                    0,
                    1,
                    site.participant.structure.carbon_degree(site.carbon) as u32,
                    false,
                    false,
                    false,
                ),
            ))
            .build(),
        );
    }
    Ok(reactions)
}

fn migratory_carbon_neighbors(
    structure: &crate::chemistry::molecule::MolecularStructure,
    center: usize,
) -> Vec<usize> {
    structure
        .neighbors(center)
        .into_iter()
        .filter_map(|(neighbor, order)| {
            (structure.atoms[neighbor].element == "C" && bond_order_matches(order, 1.0))
                .then_some(neighbor)
        })
        .collect()
}

fn migration_penalty(
    structure: &crate::chemistry::molecule::MolecularStructure,
    migratory_atom: usize,
) -> f64 {
    let carbon_degree = structure.carbon_degree(migratory_atom);
    if carbon_degree >= 3 {
        -2.0
    } else if carbon_degree == 2 {
        0.0
    } else {
        2.0
    }
}
