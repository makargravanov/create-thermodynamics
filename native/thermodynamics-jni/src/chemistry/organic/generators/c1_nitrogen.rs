use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::ReactionCondition;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;

pub(crate) fn generate_isocyanate_ammonolysis(
    site: &IsocyanateSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let mut editor = MolecularEditor::new(site.participant.structure);
    editor.set_bond_order(site.nitrogen, site.functional_carbon, 1.0)?;
    editor.add_atom(site.nitrogen, "H", 0.0, 1.0)?;
    let ammonia_nitrogen = editor.add_atom(site.functional_carbon, "N", 0.0, 1.0)?;
    editor.add_atom(ammonia_nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(ammonia_nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "isocyanate_ammonolysis",
        &site.participant,
    ))
    .reactant(site.participant.substance.id.clone(), 1, 1)
    .reactant("destroy:ammonia", 1, 1)
    .product(product, 1)
    .condition(
        ReactionCondition::new("isocyanate ammonolysis requires a dry medium")
            .max_water_activity(0.1),
    )
    .activation_energy_kj_per_mol(20.0)
    .build())
}
