use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};

pub(crate) fn generate_hydrazone_aryl_annulation(
    center: &ArylHydrazoneCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    let mut reactions = Vec::new();
    for annulation_site in center.annulation_sites() {
        let center = &annulation_site.aryl_hydrazone;
        for sidechain_carbon in annulation_sidechain_carbons(center.participant.structure, center) {
            let original_sidechain_carbon = sidechain_carbon;
            let sidechain_hydrogens =
                center.participant.structure.bonded_atoms_by_element(sidechain_carbon, "H");
            if sidechain_hydrogens.len() < 2 {
                continue;
            }

            let mut editor = MolecularEditor::new(center.participant.structure);
            let mapping = editor.remove_atoms(&[
                center.imine_nitrogen,
                annulation_site.ortho_hydrogen,
                sidechain_hydrogens[0],
                sidechain_hydrogens[1],
            ])?;
            let terminal_nitrogen = mapped_atom(
                &mapping,
                center.terminal_nitrogen,
                "annulation terminal nitrogen",
            )?;
            let hydrazone_carbon =
                mapped_atom(&mapping, center.carbon, "annulation hydrazone carbon")?;
            let sidechain_carbon =
                mapped_atom(&mapping, sidechain_carbon, "annulation sidechain carbon")?;
            let ortho_atom = mapped_atom(&mapping, annulation_site.ortho_atom, "annulation ortho atom")?;
            editor.add_bond(terminal_nitrogen, hydrazone_carbon, 1.0)?;
            editor.set_bond_order(hydrazone_carbon, sidechain_carbon, 2.0)?;
            editor.add_bond(sidechain_carbon, ortho_atom, 1.0)?;
            let product = resolver.resolve(editor.finish()?)?;

            reactions.push(
                Reaction::builder(format!(
                    "{}/ortho_{}_side_{}",
                    generated_site_reaction_id(
                        "hydrazone_aryl_annulation",
                        &center.participant,
                    ),
                    annulation_site.ortho_atom,
                    original_sidechain_carbon
                ))
                .reactant(center.participant.substance.id.clone(), 1, 1)
                .product(product, 1)
                .product("destroy:ammonia", 1)
                .catalyst_order("destroy:proton", 1)
                .condition(
                    ReactionCondition::new(
                        "aryl hydrazone annulation requires acidic dehydrating conditions",
                    )
                    .acidity(AcidityCondition::Acidic)
                    .max_water_activity(0.25),
                )
                .activation_energy_kj_per_mol(58.0)
                .selectivity_profile(SelectivityProfile::new(
                    ReactionType::SkeletalRearrangement,
                    SiteDescriptorBuilder::build(
                        center.participant.site.kind.clone(),
                        crate::chemistry::selectivity::types::SubstitutionDegree::Secondary,
                        0,
                        1,
                        center
                            .participant
                            .structure
                            .carbon_degree(center.carbon) as u32,
                        true,
                        true,
                        true,
                    ),
                ))
                .build(),
            );
        }
    }
    Ok(reactions)
}

fn annulation_sidechain_carbons(
    structure: &MolecularStructure,
    center: &ArylHydrazoneCenter<'_>,
) -> Vec<usize> {
    structure
        .neighbors(center.carbon)
        .into_iter()
        .filter_map(|(neighbor, order)| {
            (neighbor != center.imine_nitrogen
                && structure.atoms[neighbor].element == "C"
                && bond_order_matches(order, 1.0))
            .then_some(neighbor)
        })
        .collect()
}
