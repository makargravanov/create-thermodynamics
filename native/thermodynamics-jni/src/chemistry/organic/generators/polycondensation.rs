use super::super::resolver::DerivedSubstanceResolver;
use super::common::generated_pair_reaction_id;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::polymer::{step_growth_polymer_substance, StepGrowthLinkage};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::substance::{Substance, SubstanceTagId};

const WATER_ID: &str = "destroy:water";
/// Polycondensation is an equilibrium driven by removing water. Water activity is
/// handled by the general medium/equilibrium layer rather than by a binary gate
/// here.
const POLYCONDENSATION_ACTIVATION_ENERGY_KJ_PER_MOL: f64 = 55.0;

/// Generates an AA+BB step-growth polycondensation edge: a difunctional diacid
/// plus a difunctional diol (polyester) or diamine (polyamide) condense into one
/// repeat unit, expelling two waters. The edge `1 diacid + 1 comonomer →
/// 1 repeat-unit + 2 water` is mass-exact (see
/// [`step_growth_repeat_unit_structure`]). Returns `None` unless the pair is a
/// clean difunctional match.
pub(crate) fn generate_polycondensation(
    diacid: &Substance,
    comonomer: &Substance,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // A monomer condensing with itself is the AB / ring-closure case, handled
    // elsewhere; here we only pair two distinct difunctional molecules.
    if diacid.id == comonomer.id {
        return Ok(None);
    }
    if !is_concrete_neutral(diacid) || !is_concrete_neutral(comonomer) {
        return Ok(None);
    }
    let (Some(acid_structure), Some(comonomer_structure)) = (
        diacid.molecular_structure.as_ref(),
        comonomer.molecular_structure.as_ref(),
    ) else {
        return Ok(None);
    };
    let Some((polymer, linkage, repeat_count)) = step_growth_polymer_substance(
        &diacid.id,
        acid_structure,
        &comonomer.id,
        comonomer_structure,
    )?
    else {
        return Ok(None);
    };
    let product = resolver.resolve_substance(polymer)?;
    let prefix = match linkage {
        StepGrowthLinkage::Polyester => "polyesterification",
        StepGrowthLinkage::Polyamide => "polyamidation",
    };
    Ok(Some(
        Reaction::builder(generated_pair_reaction_id(prefix, diacid, comonomer))
            .reactant(diacid.id.clone(), repeat_count, 1)
            .reactant(comonomer.id.clone(), repeat_count, 0)
            .product(product, 1)
            .product(WATER_ID, repeat_count * 2)
            .activation_energy_kj_per_mol(POLYCONDENSATION_ACTIVATION_ENERGY_KJ_PER_MOL)
            .build(),
    ))
}

fn is_concrete_neutral(substance: &Substance) -> bool {
    substance.charge == 0
        && !substance
            .tags
            .iter()
            .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
}
