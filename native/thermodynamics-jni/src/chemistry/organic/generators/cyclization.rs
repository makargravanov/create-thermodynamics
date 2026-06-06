use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::{
    bond_order_matches, would_form_ring_of_size, MolecularEditor, MolecularStructure,
};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};
use crate::chemistry::substance::SubstanceId;

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

/// Specification for a generic intramolecular ring closure: a nucleophilic atom
/// and an electrophilic atom on ONE molecule bond together, expelling a set of
/// leaving atoms and emitting byproducts. Lactonization, lactamization and
/// intramolecular N-alkylation are all this one operation with different atoms,
/// leaving groups and energetics. The topology engine (`would_form_ring_of_size`)
/// is blind to whether the closing atoms already sit in another ring, so closures
/// that FUSE a new ring onto an existing one fall out of this core for free.
struct IntramolecularClosure<'a> {
    /// Reaction-id prefix; the closed ring size is appended as `{prefix}_{size}`.
    prefix: &'static str,
    structure: &'a MolecularStructure,
    substance_id: SubstanceId,
    /// Atom that keeps its identity and gains the new bond (e.g. alcohol O, amine N).
    nucleophile: usize,
    /// Atom the nucleophile bonds to (e.g. acyl C, alkyl-halide C).
    electrophile: usize,
    nucleophile_label: &'static str,
    electrophile_label: &'static str,
    /// Atoms removed before the new bond forms (leaving group + shed hydrogens).
    leaving_atoms: Vec<usize>,
    /// Small-molecule byproducts (water, halide ion, proton, …).
    byproducts: Vec<&'static str>,
    catalysts: Vec<&'static str>,
    base_activation_energy_kj_per_mol: f64,
    selectivity_profile: SelectivityProfile,
}

/// Performs the topology-and-bookkeeping half of every intramolecular closure:
/// validates the closure is a real ring of closable size, edits the molecule
/// (drop leaving atoms, bond nucleophile→electrophile, re-perceive aromaticity),
/// and assembles the reaction. Mechanism-specific screening (which atoms, which
/// leaving group, sp2 rejection, …) lives in the thin wrappers that build the
/// `IntramolecularClosure`. Returns `None` whenever the closure is not viable.
///
/// The reaction id is keyed on the canonical `(nucleophile, electrophile)` atom
/// pair, NOT on site atom-tokens: a single primary amine surfaces as three amine
/// sites, which would otherwise emit three distinct ids for ONE ring closure.
/// Pair-keying lets `push_unique_reaction` collapse the duplicates.
fn close_intramolecular_ring(
    spec: IntramolecularClosure<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let structure = spec.structure;
    if spec.nucleophile == spec.electrophile
        || structure.are_bonded(spec.nucleophile, spec.electrophile)
    {
        return Ok(None);
    }
    let Some(ring_size) =
        would_form_ring_of_size(structure, spec.nucleophile, spec.electrophile)
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
    editor.add_bond(ring_nucleophile, ring_electrophile, 1.0)?;
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


pub(crate) fn generate_lactonization(
    acid_site: &CarboxylicAcidSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Both functional groups must live on one molecule for an intramolecular closure.
    if acid_site.participant.substance.id != alcohol_site.participant.substance.id {
        return Ok(None);
    }
    let acid_desc = SiteDescriptorBuilder::from_carboxylic_acid_site(acid_site);
    let alcohol_desc = SiteDescriptorBuilder::from_alcohol_site(alcohol_site);
    // Nucleophile = alkyl oxygen, electrophile = acyl carbon. The acid loses its
    // -OH (oxygen + proton) and the alcohol loses its proton; water is expelled,
    // acid catalysed. The new acyl-C–O bond closes the lactone ring.
    close_intramolecular_ring(
        IntramolecularClosure {
            prefix: "lactonization",
            structure: acid_site.participant.structure,
            substance_id: acid_site.participant.substance.id.clone(),
            nucleophile: alcohol_site.oxygen,
            electrophile: acid_site.carbon,
            nucleophile_label: "lactone alkyl oxygen",
            electrophile_label: "lactone acyl carbon",
            leaving_atoms: vec![
                acid_site.hydroxyl_oxygen,
                acid_site.hydroxyl_hydrogen,
                alcohol_site.hydrogen,
            ],
            byproducts: vec!["destroy:water"],
            catalysts: vec!["destroy:sulfuric_acid"],
            base_activation_energy_kj_per_mol: 25.0,
            selectivity_profile: SelectivityProfile::new(
                ReactionType::Lactonization,
                acid_desc,
            )
            .with_secondary_site(alcohol_desc),
        },
        resolver,
    )
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
    // A primary amine surfaces as several amine sites; pick any one N-H to shed.
    let Some(nitrogen_hydrogen) = amine_site.hydrogens.first().copied() else {
        return Ok(None);
    };
    let acid_desc = SiteDescriptorBuilder::from_carboxylic_acid_site(acid_site);
    let amine_desc = SiteDescriptorBuilder::from_amine_site(amine_site);
    // Nucleophile = amine nitrogen, electrophile = acyl carbon. The acid loses its
    // -OH and the amine loses one N-H; water is expelled and the acyl-C–N bond
    // closes the lactam ring.
    close_intramolecular_ring(
        IntramolecularClosure {
            prefix: "lactamization",
            structure: acid_site.participant.structure,
            substance_id: acid_site.participant.substance.id.clone(),
            nucleophile: amine_site.nitrogen,
            electrophile: acid_site.carbon,
            nucleophile_label: "lactam nitrogen",
            electrophile_label: "lactam acyl carbon",
            leaving_atoms: vec![
                acid_site.hydroxyl_oxygen,
                acid_site.hydroxyl_hydrogen,
                nitrogen_hydrogen,
            ],
            byproducts: vec!["destroy:water"],
            catalysts: vec![],
            base_activation_energy_kj_per_mol: 40.0,
            selectivity_profile: SelectivityProfile::new(
                ReactionType::Lactamization,
                acid_desc,
            )
            .with_secondary_site(amine_desc),
        },
        resolver,
    )
}

/// Intramolecular N-alkylation: an amine nitrogen and an alkyl halide carbon on
/// the SAME molecule close to a saturated N-heterocycle (pyrrolidine, piperidine,
/// …), expelling HX. The nitrogen displaces the halide in an internal SN2, so the
/// new N–C bond forms while the C–X bond breaks. Mirrors the lactam closure but
/// the leaving group is a halide ion plus a proton rather than water.
pub(crate) fn generate_intramolecular_n_alkylation(
    amine_site: &AmineSite<'_>,
    halide_site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if amine_site.participant.substance.id != halide_site.participant.substance.id {
        return Ok(None);
    }
    // The nitrogen must not coincide with the halogen atom (a degenerate N-haloamine
    // topology); the generic core handles nitrogen == carbon and already-bonded.
    if amine_site.nitrogen == halide_site.halogen {
        return Ok(None);
    }
    let structure = amine_site.participant.structure;
    // Reject sp2 (vinyl) halide carbons: an internal SN2 needs back-side attack on
    // a tetrahedral sp3 carbon, geometrically impossible at a planar alkene carbon.
    // The functional-group detector tags vinyl halides as plain Halide sites, so
    // this generator must screen them out itself. (Mechanism-specific: an acyl
    // electrophile in the lactone/lactam closures legitimately IS sp2, so this
    // guard belongs here in the SN2 wrapper, not in the generic core.)
    if structure
        .neighbors(halide_site.carbon)
        .into_iter()
        .any(|(n, order)| bond_order_matches(order, 2.0) && structure.atoms[n].element == "C")
    {
        return Ok(None);
    }
    let Some(nitrogen_hydrogen) = amine_site.hydrogens.first().copied() else {
        return Ok(None);
    };
    let leaving_halide = halide_ion(
        structure,
        halide_site.halogen,
        "intramolecular_n_alkylation",
        &halide_site.participant,
    )?;
    let amine_desc = SiteDescriptorBuilder::from_amine_site(amine_site);
    let halide_desc = SiteDescriptorBuilder::from_halide_site(halide_site);
    // Nucleophile = amine nitrogen, electrophile = alkyl-halide carbon. The halide
    // ion and a proton leave; the new N–C bond closes the azacycle. Intramolecular
    // SN2 forming a favorable 5/6-ring beats the intermolecular halide+amine path
    // (base EA 25) thanks to effective molarity, so base EA is 18.
    close_intramolecular_ring(
        IntramolecularClosure {
            prefix: "intramolecular_n_alkylation",
            structure,
            substance_id: amine_site.participant.substance.id.clone(),
            nucleophile: amine_site.nitrogen,
            electrophile: halide_site.carbon,
            nucleophile_label: "azacycle nitrogen",
            electrophile_label: "azacycle carbon",
            leaving_atoms: vec![halide_site.halogen, nitrogen_hydrogen],
            byproducts: vec![leaving_halide, "destroy:proton"],
            catalysts: vec![],
            base_activation_energy_kj_per_mol: 18.0,
            selectivity_profile: SelectivityProfile::new(
                ReactionType::NAlkylation,
                amine_desc,
            )
            .with_secondary_site(halide_desc),
        },
        resolver,
    )
}

/// For two carbonyl carbons `first` and `second` on one molecule, find the two
/// bridging carbons of a 1,4-dicarbonyl (`first`–Cβ–Cγ–`second`). Returns
/// `(c_beta, c_gamma)` with each bridge carbon carrying at least one hydrogen to
/// shed on aromatisation, or `None` when the carbons are not 1,4-related.
fn one_four_dicarbonyl_bridge(
    structure: &MolecularStructure,
    first: usize,
    second: usize,
) -> Option<(usize, usize)> {
    for (c_beta, beta_order) in structure.neighbors(first) {
        if !bond_order_matches(beta_order, 1.0)
            || structure.atoms[c_beta].element != "C"
            || c_beta == second
            || bonded_hydrogens(structure, c_beta).is_empty()
        {
            continue;
        }
        for (c_gamma, gamma_order) in structure.neighbors(second) {
            if !bond_order_matches(gamma_order, 1.0)
                || structure.atoms[c_gamma].element != "C"
                || c_gamma == first
                || c_gamma == c_beta
                || bonded_hydrogens(structure, c_gamma).is_empty()
            {
                continue;
            }
            if structure
                .neighbors(c_beta)
                .iter()
                .any(|(neighbor, order)| *neighbor == c_gamma && bond_order_matches(*order, 1.0))
            {
                return Some((c_beta, c_gamma));
            }
        }
    }
    None
}

/// Paal–Knorr furan synthesis: a 1,4-dicarbonyl on one molecule closes to an
/// (aromatic) furan with loss of water. One carbonyl oxygen becomes the ring
/// oxygen bridging the two former carbonyl carbons; the other leaves as water,
/// taking one hydrogen from each bridging carbon. The Kekulé ring laid down here
/// is aromatised by `editor.finish()`.
pub(crate) fn generate_paal_knorr_furan(
    first: &CarbonylSite<'_>,
    second: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if first.participant.substance.id != second.participant.substance.id {
        return Ok(None);
    }
    // Deterministic ordering so each unordered carbonyl pair is enumerated once.
    if first.carbon >= second.carbon {
        return Ok(None);
    }
    let structure = first.participant.structure;
    let Some((c_beta, c_gamma)) =
        one_four_dicarbonyl_bridge(structure, first.carbon, second.carbon)
    else {
        return Ok(None);
    };
    // The ring the new O–C bond closes (ring O + four carbons) must be furan-sized.
    if would_form_ring_of_size(structure, first.oxygen, second.carbon) != Some(5) {
        return Ok(None);
    }
    let Some(beta_hydrogen) = first_bonded_hydrogen(structure, c_beta) else {
        return Ok(None);
    };
    let Some(gamma_hydrogen) = first_bonded_hydrogen(structure, c_gamma) else {
        return Ok(None);
    };

    let first_desc = SiteDescriptorBuilder::from_carbonyl_site(first);
    let second_desc = SiteDescriptorBuilder::from_carbonyl_site(second);
    let substance = first.participant.substance;

    // Remove the departing carbonyl oxygen and one H from each bridge carbon,
    // then lay down a Kekulé furan: O(ring)-Ca=Cb-Cc=Cd-O(ring).
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[second.oxygen, beta_hydrogen, gamma_hydrogen])?;
    let ring_oxygen = mapped_atom(&mapping, first.oxygen, "furan ring oxygen")?;
    let carbon_a = mapped_atom(&mapping, first.carbon, "furan carbon a")?;
    let carbon_b = mapped_atom(&mapping, c_beta, "furan carbon b")?;
    let carbon_c = mapped_atom(&mapping, c_gamma, "furan carbon c")?;
    let carbon_d = mapped_atom(&mapping, second.carbon, "furan carbon d")?;
    editor.set_bond_order(carbon_a, ring_oxygen, 1.0)?;
    editor.set_bond_order(carbon_a, carbon_b, 2.0)?;
    editor.set_bond_order(carbon_c, carbon_d, 2.0)?;
    editor.add_bond(ring_oxygen, carbon_d, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_intramolecular_pair_site_reaction_id(
            "paal_knorr_furan",
            &first.participant,
            &second.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .catalyst_order("destroy:proton", 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .condition(
            ReactionCondition::new("Paal–Knorr furan closure needs acidic, water-poor conditions")
                .acidity(AcidityCondition::Acidic)
                .max_water_activity(0.35),
        )
        .activation_energy_kj_per_mol(45.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::HeterocycleCondensation, first_desc)
                .with_secondary_site(second_desc),
        )
        .build(),
    ))
}

/// Paal–Knorr pyrrole synthesis: a 1,4-dicarbonyl condenses with a primary amine
/// to an (aromatic) pyrrole, expelling two waters. The amine nitrogen bridges the
/// two former carbonyl carbons (taking the place of furan's ring oxygen); both
/// carbonyl oxygens leave as water, each pulling a hydrogen — two from the amine
/// N-H pair and one from each bridging carbon. Unlike furan closure this is
/// intermolecular, so the carbonyl and amine may be different substances.
pub(crate) fn generate_paal_knorr_pyrrole(
    first: &CarbonylSite<'_>,
    second: &CarbonylSite<'_>,
    amine: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // The two carbonyls must share one molecule; the amine is a separate donor.
    if first.participant.substance.id != second.participant.substance.id {
        return Ok(None);
    }
    if first.carbon >= second.carbon {
        return Ok(None);
    }
    if amine.hydrogens.len() < 2 {
        return Ok(None);
    }
    // The amine must be a separate molecule: this generator models the
    // intermolecular condensation and registers two reactant terms. An
    // amino-diketone closing on its own nitrogen would be a different
    // (intramolecular) mechanism and would also corrupt the stoichiometry by
    // listing one substance twice.
    if amine.participant.substance.id == first.participant.substance.id {
        return Ok(None);
    }
    let structure = first.participant.structure;
    let Some((c_beta, c_gamma)) =
        one_four_dicarbonyl_bridge(structure, first.carbon, second.carbon)
    else {
        return Ok(None);
    };
    let Some(beta_hydrogen) = first_bonded_hydrogen(structure, c_beta) else {
        return Ok(None);
    };
    let Some(gamma_hydrogen) = first_bonded_hydrogen(structure, c_gamma) else {
        return Ok(None);
    };

    let first_desc = SiteDescriptorBuilder::from_carbonyl_site(first);
    let amine_desc = SiteDescriptorBuilder::from_amine_site(amine);
    let carbonyl_substance = first.participant.substance;
    let amine_substance = amine.participant.substance;

    // Build the bridging-nitrogen fragment: strip two N-H so the nitrogen has the
    // open valences to bond to both ring carbons. Left unvalidated (under-valenced)
    // on purpose — the combined `finish()` validates the assembled pyrrole.
    let mut amine_editor = MolecularEditor::new(amine.participant.structure);
    let amine_mapping = amine_editor.remove_atoms(&[amine.hydrogens[0], amine.hydrogens[1]])?;
    let amine_nitrogen = mapped_atom(&amine_mapping, amine.nitrogen, "pyrrole nitrogen")?;
    let amine_fragment = amine_editor.structure();

    // Strip both carbonyl oxygens and one H from each bridge carbon, splice in the
    // nitrogen, then lay down the Kekulé pyrrole N-Ca=Cb-Cc=Cd-N.
    let mut editor = MolecularEditor::new(structure);
    let mapping =
        editor.remove_atoms(&[first.oxygen, second.oxygen, beta_hydrogen, gamma_hydrogen])?;
    let carbon_a = mapped_atom(&mapping, first.carbon, "pyrrole carbon a")?;
    let carbon_b = mapped_atom(&mapping, c_beta, "pyrrole carbon b")?;
    let carbon_c = mapped_atom(&mapping, c_gamma, "pyrrole carbon c")?;
    let carbon_d = mapped_atom(&mapping, second.carbon, "pyrrole carbon d")?;
    let ring_nitrogen = editor.add_group(carbon_a, &amine_fragment, amine_nitrogen, 1.0)?;
    editor.add_bond(ring_nitrogen, carbon_d, 1.0)?;
    editor.set_bond_order(carbon_a, carbon_b, 2.0)?;
    editor.set_bond_order(carbon_c, carbon_d, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_triple_site_reaction_id(
            "paal_knorr_pyrrole",
            &first.participant,
            &second.participant,
            &amine.participant,
        ))
        .reactant(carbonyl_substance.id.clone(), 1, 1)
        .reactant(amine_substance.id.clone(), 1, 1)
        .catalyst_order("destroy:proton", 1)
        .product(product, 1)
        .product("destroy:water", 2)
        .condition(
            ReactionCondition::new("Paal–Knorr pyrrole closure needs acidic, water-poor conditions")
                .acidity(AcidityCondition::Acidic)
                .max_water_activity(0.35),
        )
        .activation_energy_kj_per_mol(50.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::HeterocycleCondensation, first_desc)
                .with_secondary_site(amine_desc),
        )
        .build(),
    ))
}

/// Paal–Knorr thiophene synthesis: a 1,4-dicarbonyl condenses with a sulfur
/// donor (hydrogen sulfide or a thiol bearing two S–H) to an (aromatic)
/// thiophene, expelling two waters. Mirrors `generate_paal_knorr_pyrrole` but the
/// bridging heteroatom is sulfur. Both carbonyl oxygens leave as water; the donor
/// sheds its two S–H protons and each bridge carbon loses one hydrogen.
pub(crate) fn generate_paal_knorr_thiophene(
    first: &CarbonylSite<'_>,
    second: &CarbonylSite<'_>,
    thiol: &ThiolSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if first.participant.substance.id != second.participant.substance.id {
        return Ok(None);
    }
    if first.carbon >= second.carbon {
        return Ok(None);
    }
    // The donor must give up two S–H to bridge both ring carbons (H2S, or a
    // dithiol-like center with two sulfur hydrogens).
    if thiol.hydrogens.len() < 2 {
        return Ok(None);
    }
    // The sulfur donor must be a separate molecule (see the pyrrole note): this is
    // the intermolecular condensation and lists two reactant terms.
    if thiol.participant.substance.id == first.participant.substance.id {
        return Ok(None);
    }
    let structure = first.participant.structure;
    let Some((c_beta, c_gamma)) =
        one_four_dicarbonyl_bridge(structure, first.carbon, second.carbon)
    else {
        return Ok(None);
    };
    let Some(beta_hydrogen) = first_bonded_hydrogen(structure, c_beta) else {
        return Ok(None);
    };
    let Some(gamma_hydrogen) = first_bonded_hydrogen(structure, c_gamma) else {
        return Ok(None);
    };

    let first_desc = SiteDescriptorBuilder::from_carbonyl_site(first);
    let thiol_desc = SiteDescriptorBuilder::from_thiol_site(thiol);
    let carbonyl_substance = first.participant.substance;
    let thiol_substance = thiol.participant.substance;

    // Under-valenced sulfur fragment: strip two S–H so sulfur can bond to both
    // ring carbons. The combined `finish()` validates the assembled thiophene.
    let mut thiol_editor = MolecularEditor::new(thiol.participant.structure);
    let thiol_mapping = thiol_editor.remove_atoms(&[thiol.hydrogens[0], thiol.hydrogens[1]])?;
    let donor_sulfur = mapped_atom(&thiol_mapping, thiol.sulfur, "thiophene sulfur")?;
    let sulfur_fragment = thiol_editor.structure();

    let mut editor = MolecularEditor::new(structure);
    let mapping =
        editor.remove_atoms(&[first.oxygen, second.oxygen, beta_hydrogen, gamma_hydrogen])?;
    let carbon_a = mapped_atom(&mapping, first.carbon, "thiophene carbon a")?;
    let carbon_b = mapped_atom(&mapping, c_beta, "thiophene carbon b")?;
    let carbon_c = mapped_atom(&mapping, c_gamma, "thiophene carbon c")?;
    let carbon_d = mapped_atom(&mapping, second.carbon, "thiophene carbon d")?;
    let ring_sulfur = editor.add_group(carbon_a, &sulfur_fragment, donor_sulfur, 1.0)?;
    editor.add_bond(ring_sulfur, carbon_d, 1.0)?;
    editor.set_bond_order(carbon_a, carbon_b, 2.0)?;
    editor.set_bond_order(carbon_c, carbon_d, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_triple_site_reaction_id(
            "paal_knorr_thiophene",
            &first.participant,
            &second.participant,
            &thiol.participant,
        ))
        .reactant(carbonyl_substance.id.clone(), 1, 1)
        .reactant(thiol_substance.id.clone(), 1, 1)
        .catalyst_order("destroy:proton", 1)
        .product(product, 1)
        .product("destroy:water", 2)
        .condition(
            ReactionCondition::new(
                "Paal–Knorr thiophene closure needs acidic, water-poor conditions",
            )
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.35),
        )
        .activation_energy_kj_per_mol(55.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::HeterocycleCondensation, first_desc)
                .with_secondary_site(thiol_desc),
        )
        .build(),
    ))
}

/// The four carbons of a conjugated diene C1=C2–C3=C4, seeded from one of its two
/// double bonds. `seed` provides the C1=C2 double bond; this finds the C2–C3
/// single bond and the conjugated C3=C4 double bond. Returns `(c1, c2, c3, c4)`
/// where c1 and c4 are the terminal carbons that form the new σ bonds in a
/// Diels–Alder cycloaddition. Returns `None` if the seed is not part of a
/// conjugated diene, or if it is the higher-indexed of the two double bonds (so a
/// given diene is enumerated once, from its lower-indexed double bond only).
fn conjugated_diene_carbons(
    structure: &MolecularStructure,
    seed: &UnsaturatedBondSite<'_>,
) -> Option<(usize, usize, usize, usize)> {
    if seed.is_alkyne {
        return None;
    }
    let p = seed.high_degree_carbon;
    let q = seed.low_degree_carbon;
    // Try each seed carbon as the inner carbon (c2) bonded to the second alkene.
    for (c1, c2) in [(p, q), (q, p)] {
        for (c3, c2_c3_order) in structure.neighbors(c2) {
            if c3 == c1
                || structure.atoms[c3].element != "C"
                || !bond_order_matches(c2_c3_order, 1.0)
            {
                continue;
            }
            // c3 must carry a second C=C double bond to a carbon c4 (not back to c2).
            let Some(c4) = structure.neighbors(c3).into_iter().find_map(|(n, order)| {
                (n != c2 && structure.atoms[n].element == "C" && bond_order_matches(order, 2.0))
                    .then_some(n)
            }) else {
                continue;
            };
            // Enumerate each diene once: seed must be the lower-indexed double
            // bond. Skip only this (c3,c4) pair — a branched c2 may have another
            // neighbor that yields a valid, non-duplicate diene. (Role-swapped
            // diene/dienophile pairs are intentionally NOT suppressed here; they
            // give distinct regiochemical products with distinct reaction ids.)
            if c1.min(c2) > c3.min(c4) {
                continue;
            }
            return Some((c1, c2, c3, c4));
        }
    }
    None
}

/// Diels–Alder [4+2] cycloaddition: a conjugated diene (C1=C2–C3=C4 on the
/// `diene` molecule) and a dienophile alkene (Ca=Cb on the `dienophile` molecule)
/// combine into a cyclohexene. Two new σ bonds form (C1–Ca and C4–Cb), the
/// diene's two double bonds collapse to a single C2=C3 double bond, and the
/// dienophile's C=C becomes a single bond. No atoms are lost — the six ring
/// carbons are C1,C2,C3,C4 (from the diene) and Ca,Cb (from the dienophile).
pub(crate) fn generate_diels_alder(
    diene: &UnsaturatedBondSite<'_>,
    dienophile: &UnsaturatedBondSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Same-species pairs ARE allowed: this models intermolecular homodimerization
    // (e.g. two butadiene molecules → 4-vinylcyclohexene). add_group always splices
    // a fresh copy of the dienophile structure, so the product is always a two-
    // molecule adduct — this generator is structurally incapable of a true
    // intramolecular Diels–Alder, which would bond two parts of one molecule
    // without copying and is a distinct mechanism left unmodelled.
    if dienophile.is_alkyne {
        return Ok(None);
    }
    let Some((c1, _c2, _c3, c4)) =
        conjugated_diene_carbons(diene.participant.structure, diene)
    else {
        return Ok(None);
    };

    let diene_desc = SiteDescriptorBuilder::from_unsaturated_bond_site(diene);
    let dienophile_desc = SiteDescriptorBuilder::from_unsaturated_bond_site(dienophile);
    let diene_substance = diene.participant.substance;
    let dienophile_substance = dienophile.participant.substance;
    let ca = dienophile.high_degree_carbon;
    let cb = dienophile.low_degree_carbon;

    // Splice the dienophile onto the diene: the C1–Ca σ bond is created by
    // add_group, which appends the dienophile's atoms (offset) and returns Ca's
    // new index. Atom indices on the diene side are unchanged by add_group.
    let mut editor = MolecularEditor::new(diene.participant.structure);
    let ca_joined = editor.add_group(c1, dienophile.participant.structure, ca, 1.0)?;
    // add_group appends the dienophile's atoms contiguously at a fixed offset and
    // returns Ca's new index, so every original dienophile index i maps to
    // i + offset. Cb is therefore at cb + offset.
    let offset = ca_joined - ca;
    let cb_joined = cb + offset;
    // Close the second σ bond C4–Cb.
    editor.add_bond(c4, cb_joined, 1.0)?;
    // Collapse the diene: C1=C2 and C3=C4 become single, the central C2–C3
    // becomes the lone ring double bond. set_bond_order rewrites existing bonds.
    editor.set_bond_order(c1, _c2, 1.0)?;
    editor.set_bond_order(_c2, _c3, 2.0)?;
    editor.set_bond_order(_c3, c4, 1.0)?;
    // The dienophile's C=C is now a saturated ring bond.
    editor.set_bond_order(ca_joined, cb_joined, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    let same_species = diene_substance.id == dienophile_substance.id;
    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "diels_alder",
        &diene.participant,
        &dienophile.participant,
    ));
    builder = if same_species {
        // Homodimerization (e.g. 2 butadiene → 4-vinylcyclohexene): one species
        // consumed twice, so a single coefficient-2, second-order term rather than
        // two duplicate terms (which would corrupt the stoichiometry).
        builder.reactant(diene_substance.id.clone(), 2, 2)
    } else {
        builder
            .reactant(diene_substance.id.clone(), 1, 1)
            .reactant(dienophile_substance.id.clone(), 1, 1)
    };

    Ok(Some(
        builder
        .product(product, 1)
        // No hard temperature gate: a Diels–Alder is thermally promoted, but there
        // is no distinct physical cutoff below which it stops — low temperature just
        // makes it slow. The Arrhenius dependence on this barrier already suppresses
        // the rate in the cold, so a binary min-temperature ban would be wrong.
        // Flat, substituent-independent barrier sized for an unactivated pair
        // (butadiene + ethylene computes ~100 kJ/mol); electron-poor dienophiles
        // react faster but that demand-lowering is not modelled per-substituent.
        .activation_energy_kj_per_mol(95.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::DielsAlder, diene_desc)
                .with_secondary_site(dienophile_desc),
        )
        .build(),
    ))
}
