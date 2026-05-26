use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;
pub(crate) fn generate_aldol_addition(
    enol: SiteParticipant<'_>,
    acceptor: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let (_, alpha_carbon) = enol_atoms(enol.structure, &enol.site)?;
    let alpha_hydrogen = first_bonded_hydrogen(enol.structure, alpha_carbon).ok_or_else(|| {
        ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("aldol_addition", &enol),
            reason: "aldol donor alpha carbon has no explicit hydrogen".to_string(),
        }
    })?;
    let (acceptor_carbon, acceptor_oxygen) =
        carbonyl_atoms_from_site(acceptor.structure, &acceptor.site, "aldol addition")?;

    let mut donor_editor = MolecularEditor::new(enol.structure);
    let donor_mapping = donor_editor.remove_atoms(&[alpha_hydrogen])?;
    let alpha_carbon = mapped_atom(&donor_mapping, alpha_carbon, "aldol alpha carbon")?;
    let donor_fragment = donor_editor.finish()?;

    let mut acceptor_editor = MolecularEditor::new(acceptor.structure);
    acceptor_editor.set_bond_order(acceptor_carbon, acceptor_oxygen, 1.0)?;
    acceptor_editor.add_atom(acceptor_oxygen, "H", 0.0, 1.0)?;
    let acceptor_fragment = acceptor_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &donor_fragment,
        alpha_carbon,
        &acceptor_fragment,
        acceptor_carbon,
        1.0,
    )?)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "aldol_addition",
        &enol,
        &acceptor,
    ))
    .reactant(enol.substance.id.clone(), 1, 1)
    .reactant(acceptor.substance.id.clone(), 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .condition(
        ReactionCondition::new("aldol addition requires basic carbonyl enolization")
            .acidity(AcidityCondition::Basic)
            .max_temperature_kelvin(323.15),
    )
    .build())
}
