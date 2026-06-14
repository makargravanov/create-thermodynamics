use super::super::resolver::DerivedSubstanceResolver;
use super::common::{mapped_atom, ring_closure_activation_penalty_kj_per_mol, sanitize_id};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::{would_form_ring_of_size, MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::types::SelectivityProfile;
use crate::chemistry::substance::SubstanceId;

/// Smallest ring (in atoms) we will close intramolecularly. 3-membered
/// lactones/lactams (alpha-lactones, aziridinones) are too strained to form under
/// ordinary conditions; closures below this are rejected outright.
const MIN_CLOSABLE_RING: usize = 4;
/// Largest ring we will enumerate. Beyond medium rings the entropic penalty makes
/// closure negligible versus intermolecular pathways, and enumerating them would
/// mostly create improbable products.
const MAX_CLOSABLE_RING: usize = 16;

fn ring_size_is_closable(ring_size: usize) -> bool {
    (MIN_CLOSABLE_RING..=MAX_CLOSABLE_RING).contains(&ring_size)
}

/// Specification for a generic intramolecular ring closure: a nucleophilic atom
/// and an electrophilic atom on one molecule bond together, expelling a set of
/// leaving atoms and emitting byproducts. Lactonization, lactamization,
/// intramolecular N-alkylation and amidine cyclization are all this one graph
/// operation with different atoms, leaving groups and energetics.
pub(crate) struct IntramolecularClosure<'a> {
    /// Reaction-id prefix; the closed ring size is appended as `{prefix}_{size}`.
    pub(crate) prefix: &'static str,
    pub(crate) structure: &'a MolecularStructure,
    pub(crate) substance_id: SubstanceId,
    /// Atom that keeps its identity and gains the new bond.
    pub(crate) nucleophile: usize,
    /// Atom the nucleophile bonds to.
    pub(crate) electrophile: usize,
    pub(crate) nucleophile_label: &'static str,
    pub(crate) electrophile_label: &'static str,
    /// Bond order of the new ring-closing bond.
    pub(crate) closure_bond_order: f64,
    /// Atoms removed before the new bond forms.
    pub(crate) leaving_atoms: Vec<usize>,
    /// Small-molecule byproducts.
    pub(crate) byproducts: Vec<&'static str>,
    pub(crate) catalysts: Vec<&'static str>,
    pub(crate) base_activation_energy_kj_per_mol: f64,
    pub(crate) selectivity_profile: SelectivityProfile,
}

/// Performs the topology-and-bookkeeping half of every intramolecular closure:
/// validates the closure is a real ring of closable size, edits the molecule,
/// resolves the product and assembles the reaction. Mechanism-specific screening
/// belongs in thin wrappers that build [`IntramolecularClosure`].
pub(crate) fn close_intramolecular_ring(
    spec: IntramolecularClosure<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let structure = spec.structure;
    if spec.nucleophile == spec.electrophile
        || structure.are_bonded(spec.nucleophile, spec.electrophile)
    {
        return Ok(None);
    }
    let Some(ring_size) = would_form_ring_of_size(structure, spec.nucleophile, spec.electrophile)
    else {
        return Ok(None);
    };
    if !ring_size_is_closable(ring_size) {
        return Ok(None);
    }
    let ring_penalty = ring_closure_activation_penalty_kj_per_mol(ring_size);

    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&spec.leaving_atoms)?;
    let ring_nucleophile = mapped_atom(&mapping, spec.nucleophile, spec.nucleophile_label)?;
    let ring_electrophile = mapped_atom(&mapping, spec.electrophile, spec.electrophile_label)?;
    editor.add_bond(ring_nucleophile, ring_electrophile, spec.closure_bond_order)?;
    let product = resolver.resolve(editor.finish()?)?;

    let mut builder = Reaction::builder(format!(
        "{}_{}/{}/{}_{}",
        spec.prefix,
        ring_size,
        sanitize_id(spec.substance_id.as_str()),
        spec.nucleophile,
        spec.electrophile
    ))
    .reactant(spec.substance_id.clone(), 1, 1)
    .product(product, 1)
    .activation_energy_kj_per_mol(spec.base_activation_energy_kj_per_mol + ring_penalty)
    .selectivity_profile(spec.selectivity_profile);
    for byproduct in spec.byproducts {
        builder = builder.product(byproduct, 1);
    }
    for catalyst in spec.catalysts {
        builder = builder.catalyst_order(catalyst, 1);
    }
    Ok(Some(builder.build()))
}
