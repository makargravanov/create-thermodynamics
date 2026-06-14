use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::substance::SubstanceId;

pub(crate) fn generate_bis_nucleophile_dicarbonyl_condensation(
    nucleophile: &BisNucleophileCenter<'_>,
    dicarbonyl: &DicarbonylElectrophileCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    if matches!(
        nucleophile.class,
        BisNucleophileClass::HydrazineLike | BisNucleophileClass::DiamineLike
    ) {
        return Ok(Vec::new());
    }
    let Some(activated) = dicarbonyl.activated_methylene_center() else {
        return Ok(Vec::new());
    };
    if activated.hydrogens.len() < 2 {
        return Ok(Vec::new());
    }
    let Some(first_n_hydrogen) = first_hydrogen_on(
        nucleophile.participant.structure,
        nucleophile.first_nucleophile,
    ) else {
        return Ok(Vec::new());
    };
    let Some(second_n_hydrogen) = first_hydrogen_on(
        nucleophile.participant.structure,
        nucleophile.second_nucleophile,
    ) else {
        return Ok(Vec::new());
    };

    let mut reactions = Vec::new();
    let mut product_ids = Vec::<SubstanceId>::new();
    for (orientation, first_nitrogen, second_nitrogen) in [
        (
            "forward",
            nucleophile.first_nucleophile,
            nucleophile.second_nucleophile,
        ),
        (
            "reverse",
            nucleophile.second_nucleophile,
            nucleophile.first_nucleophile,
        ),
    ] {
        let mut nucleophile_editor = MolecularEditor::new(nucleophile.participant.structure);
        let nucleophile_mapping =
            nucleophile_editor.remove_atoms(&[first_n_hydrogen, second_n_hydrogen])?;
        let first_nitrogen = mapped_atom(
            &nucleophile_mapping,
            first_nitrogen,
            "heterocycle first nucleophile nitrogen",
        )?;
        let second_nitrogen = mapped_atom(
            &nucleophile_mapping,
            second_nitrogen,
            "heterocycle second nucleophile nitrogen",
        )?;
        let nucleophile_fragment = nucleophile_editor.finish()?;

        let mut dicarbonyl_editor = MolecularEditor::new(dicarbonyl.participant.structure);
        let dicarbonyl_mapping = dicarbonyl_editor.remove_atoms(&[
            dicarbonyl.first_carbonyl_oxygen,
            dicarbonyl.second_carbonyl_oxygen,
            activated.hydrogens[0],
            activated.hydrogens[1],
        ])?;
        let first_carbonyl = mapped_atom(
            &dicarbonyl_mapping,
            dicarbonyl.first_carbonyl_carbon,
            "heterocycle first carbonyl carbon",
        )?;
        let second_carbonyl = mapped_atom(
            &dicarbonyl_mapping,
            dicarbonyl.second_carbonyl_carbon,
            "heterocycle second carbonyl carbon",
        )?;
        let bridge = mapped_atom(
            &dicarbonyl_mapping,
            activated.carbon,
            "heterocycle activated bridge carbon",
        )?;
        dicarbonyl_editor.set_bond_order(first_carbonyl, bridge, 2.0)?;
        dicarbonyl_editor.set_bond_order(second_carbonyl, bridge, 2.0)?;
        let dicarbonyl_fragment = dicarbonyl_editor.finish()?;

        let mut product_editor = MolecularEditor::new(&dicarbonyl_fragment);
        let nucleophile_mapping = product_editor.add_group_with_mapping(
            first_carbonyl,
            &nucleophile_fragment,
            first_nitrogen,
            1.0,
        )?;
        let second_nitrogen = nucleophile_mapping[second_nitrogen];
        product_editor.add_bond(second_carbonyl, second_nitrogen, 1.0)?;
        let product = resolver.resolve(product_editor.finish()?)?;
        if product_ids.contains(&product) {
            continue;
        }
        product_ids.push(product.clone());

        reactions.push(
            Reaction::builder(format!(
                "{}/{}",
                generated_pair_site_reaction_id(
                    "bis_nucleophile_dicarbonyl_condensation",
                    &nucleophile.participant,
                    &dicarbonyl.participant,
                ),
                orientation,
            ))
            .reactant(nucleophile.participant.substance.id.clone(), 1, 1)
            .reactant(dicarbonyl.participant.substance.id.clone(), 1, 1)
            .product(product, 1)
            .product("destroy:water", 2)
            .catalyst_order("destroy:proton", 1)
            .condition(
                ReactionCondition::new(
                    "bis-nucleophile dicarbonyl condensation requires acidic dehydration",
                )
                .acidity(AcidityCondition::Acidic)
                .max_water_activity(0.35),
            )
            .activation_energy_kj_per_mol(45.0)
            .build(),
        );
    }
    Ok(reactions)
}

fn first_hydrogen_on(
    structure: &crate::chemistry::molecule::MolecularStructure,
    atom: usize,
) -> Option<usize> {
    structure
        .bonded_atoms_by_element(atom, "H")
        .into_iter()
        .next()
}
