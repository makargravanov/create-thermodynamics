use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;

pub(crate) fn generate_alcohol_dehydration(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let alcohol_carbon = site.carbon;
    let oxygen = site.oxygen;
    let proton = site.hydrogen;
    let mut products = Vec::new();
    for (neighbor, order) in structure.neighbors(alcohol_carbon) {
        if structure.atoms[neighbor].element != "C"
            || !crate::chemistry::molecule::bond_order_matches(order, 1.0)
        {
            continue;
        }
        let Some(beta_hydrogen) = first_bonded_hydrogen(structure, neighbor) else {
            continue;
        };
        let mut editor = MolecularEditor::new(structure);
        let mapping = editor.remove_atoms(&[beta_hydrogen, oxygen, proton])?;
        let carbon = mapped_atom(&mapping, alcohol_carbon, "alcohol carbon")?;
        let neighbor = mapped_atom(&mapping, neighbor, "dehydration neighbor carbon")?;
        editor.set_bond_order(carbon, neighbor, 2.0)?;
        products.push(resolver.resolve(editor.finish()?)?);
    }
    if products.is_empty() {
        return Ok(None);
    }
    let mut builder = Reaction::builder(generated_site_reaction_id(
        "alcohol_dehydration",
        &site.participant,
    ))
    .reactant(substance.id.clone(), products.len() as u32, 1)
    .reactant("destroy:oleum", products.len() as u32, 1)
    .product("destroy:sulfuric_acid", (products.len() * 2) as u32);
    for product in products {
        builder = builder.product(product, 1);
    }
    Ok(Some(builder.build()))
}
