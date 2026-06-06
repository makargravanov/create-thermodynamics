use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::{would_form_ring_of_size, MolecularEditor};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};

/// Smallest ring (in atoms) we will close intramolecularly. 3-membered
/// lactones/lactams (α-lactones, aziridinones) are too strained to form under
/// ordinary conditions; closures below this are rejected outright.
const MIN_CLOSABLE_RING: usize = 4;
/// Largest ring we will bother enumerating. Beyond medium rings the entropic
/// penalty makes closure negligible versus intermolecular pathways, and
/// enumerating them would only spam improbable products.
const MAX_CLOSABLE_RING: usize = 16;

fn ring_size_is_closable(ring_size: usize) -> bool {
    (MIN_CLOSABLE_RING..=MAX_CLOSABLE_RING).contains(&ring_size)
}

/// Intramolecular esterification: a carboxylic acid and an alcohol on the SAME
/// molecule condense, expelling water and closing a lactone ring. The acid and
/// alcohol must be the same substance, and the ring that would form must be of
/// a closable size (Baldwin's-rules gate via `ring_closure_activation_penalty`).
pub(crate) fn generate_lactonization(
    acid_site: &CarboxylicAcidSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Both functional groups must live on one molecule for an intramolecular closure.
    if acid_site.participant.substance.id != alcohol_site.participant.substance.id {
        return Ok(None);
    }
    // The acid carbon and the alcohol oxygen must be distinct atoms; a hydroxy
    // acid where they coincide is not a ring precursor.
    if acid_site.carbon == alcohol_site.oxygen {
        return Ok(None);
    }
    let structure = acid_site.participant.structure;
    // Size of the ring the new C(acyl)-O(alkyl) bond would close.
    let Some(ring_size) = would_form_ring_of_size(structure, acid_site.carbon, alcohol_site.oxygen)
    else {
        return Ok(None);
    };
    if !ring_size_is_closable(ring_size) {
        return Ok(None);
    }

    let acid_desc = SiteDescriptorBuilder::from_carboxylic_acid_site(acid_site);
    let alcohol_desc = SiteDescriptorBuilder::from_alcohol_site(alcohol_site);
    let substance = acid_site.participant.substance;
    let ring_penalty = ring_closure_activation_penalty_kj_per_mol(ring_size);

    // Remove the acid -OH (oxygen + its proton) and the alcohol proton, then
    // bond the acyl carbon to the (now ester) alkyl oxygen, all in one editor so
    // the closure stays intramolecular.
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[
        acid_site.hydroxyl_oxygen,
        acid_site.hydroxyl_hydrogen,
        alcohol_site.hydrogen,
    ])?;
    let acyl_carbon = mapped_atom(&mapping, acid_site.carbon, "lactone acyl carbon")?;
    let alkyl_oxygen = mapped_atom(&mapping, alcohol_site.oxygen, "lactone alkyl oxygen")?;
    editor.add_bond(acyl_carbon, alkyl_oxygen, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_intramolecular_pair_site_reaction_id(
            &format!("lactonization_{ring_size}"),
            &acid_site.participant,
            &alcohol_site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .catalyst_order("destroy:sulfuric_acid", 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .activation_energy_kj_per_mol(25.0 + ring_penalty)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::Lactonization, acid_desc)
                .with_secondary_site(alcohol_desc),
        )
        .build(),
    ))
}

/// Intramolecular amidation: a carboxylic acid and an amine on the SAME
/// molecule condense, expelling water and closing a lactam ring. Mirrors
/// `generate_lactonization` but bonds the acyl carbon to the amine nitrogen.
pub(crate) fn generate_lactamization(
    acid_site: &CarboxylicAcidSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if acid_site.participant.substance.id != amine_site.participant.substance.id {
        return Ok(None);
    }
    if acid_site.carbon == amine_site.nitrogen {
        return Ok(None);
    }
    let Some(nitrogen_hydrogen) = amine_site.hydrogens.first().copied() else {
        return Ok(None);
    };
    let structure = acid_site.participant.structure;
    let Some(ring_size) =
        would_form_ring_of_size(structure, acid_site.carbon, amine_site.nitrogen)
    else {
        return Ok(None);
    };
    if !ring_size_is_closable(ring_size) {
        return Ok(None);
    }

    let acid_desc = SiteDescriptorBuilder::from_carboxylic_acid_site(acid_site);
    let amine_desc = SiteDescriptorBuilder::from_amine_site(amine_site);
    let substance = acid_site.participant.substance;
    let ring_penalty = ring_closure_activation_penalty_kj_per_mol(ring_size);

    // Remove the acid -OH and one N-H, then bond the acyl carbon to the nitrogen.
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[
        acid_site.hydroxyl_oxygen,
        acid_site.hydroxyl_hydrogen,
        nitrogen_hydrogen,
    ])?;
    let acyl_carbon = mapped_atom(&mapping, acid_site.carbon, "lactam acyl carbon")?;
    let amide_nitrogen = mapped_atom(&mapping, amine_site.nitrogen, "lactam nitrogen")?;
    editor.add_bond(acyl_carbon, amide_nitrogen, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_intramolecular_pair_site_reaction_id(
            &format!("lactamization_{ring_size}"),
            &acid_site.participant,
            &amine_site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .activation_energy_kj_per_mol(40.0 + ring_penalty)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::Lactamization, acid_desc)
                .with_secondary_site(amine_desc),
        )
        .build(),
    ))
}
