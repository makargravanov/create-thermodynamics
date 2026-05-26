use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::kinetics::ReactionChannel;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::{Reaction, StoichiometricTerm};
pub(crate) fn generate_aromatic_nitration(
    aromatic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let mut variants = Vec::new();
    for carbon in aromatic_substitution_carbons(aromatic.structure, &aromatic.site) {
        let Some(hydrogen) = first_bonded_hydrogen(aromatic.structure, carbon) else {
            continue;
        };
        let mut editor = MolecularEditor::new(aromatic.structure);
        let mapping = editor.remove_atoms(&[hydrogen])?;
        let carbon = mapped_atom(&mapping, carbon, "aromatic nitration carbon")?;
        let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
        editor.add_atom(nitrogen, "O", 0.0, 1.5)?;
        editor.add_atom(nitrogen, "O", 0.0, 1.5)?;
        let product = resolver.resolve(editor.finish()?)?;
        variants.push((
            product,
            aromatic_activation_delta(aromatic.structure, carbon),
        ));
    }
    if variants.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("aromatic_nitration", &aromatic),
            reason: "aromatic nitration found no aromatic carbon with explicit hydrogen"
                .to_string(),
        });
    }
    let mut builder =
        Reaction::builder(generated_site_reaction_id("aromatic_nitration", &aromatic))
            .reactant(aromatic.substance.id.clone(), 1, 1)
            .reactant("destroy:nitric_acid", 1, 1)
            .catalyst_order("destroy:sulfuric_acid", 1)
            .condition(
                ReactionCondition::new("aromatic nitration requires strongly acidic conditions")
                    .acidity(AcidityCondition::Acidic)
                    .max_water_activity(0.65),
            );
    if variants.len() == 1 {
        builder = builder
            .product(variants[0].0.clone(), 1)
            .product("destroy:water", 1);
    } else {
        for (index, (product, activation_delta)) in variants.into_iter().enumerate() {
            builder = builder.channel(ReactionChannel::new(
                format!("aromatic_nitration:position_{index}"),
                [
                    StoichiometricTerm::new(product, 1),
                    StoichiometricTerm::new("destroy:water", 1),
                ],
                30.0 + activation_delta,
            ));
        }
    }
    Ok(builder.build())
}
