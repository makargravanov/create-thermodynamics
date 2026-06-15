use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use super::ring_closure::{close_intramolecular_ring, IntramolecularClosure};
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::kinetics::{
    ChannelConditionEffect, LightBand, ReactionChannel, ReactionChannelMode,
};
use crate::chemistry::molecule::{
    bond_order_matches, would_form_ring_of_size, MolecularBond, MolecularEditor,
    MolecularStructure, StereoDescriptor,
};
use crate::chemistry::reaction::{Reaction, StoichiometricTerm};
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};

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
            closure_bond_order: 1.0,
            leaving_atoms: vec![
                acid_site.hydroxyl_oxygen,
                acid_site.hydroxyl_hydrogen,
                alcohol_site.hydrogen,
            ],
            byproducts: vec!["destroy:water"],
            catalysts: vec!["destroy:sulfuric_acid"],
            base_activation_energy_kj_per_mol: 25.0,
            selectivity_profile: SelectivityProfile::new(ReactionType::Lactonization, acid_desc)
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
            closure_bond_order: 1.0,
            leaving_atoms: vec![
                acid_site.hydroxyl_oxygen,
                acid_site.hydroxyl_hydrogen,
                nitrogen_hydrogen,
            ],
            byproducts: vec!["destroy:water"],
            catalysts: vec![],
            base_activation_energy_kj_per_mol: 40.0,
            selectivity_profile: SelectivityProfile::new(ReactionType::Lactamization, acid_desc)
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
            closure_bond_order: 1.0,
            leaving_atoms: vec![halide_site.halogen, nitrogen_hydrogen],
            byproducts: vec![leaving_halide, "destroy:proton"],
            catalysts: vec![],
            base_activation_energy_kj_per_mol: 18.0,
            selectivity_profile: SelectivityProfile::new(ReactionType::NAlkylation, amine_desc)
                .with_secondary_site(halide_desc),
        },
        resolver,
    )
}

/// Intramolecular dehydrative amidine cyclization: a primary amine and an amide
/// on the SAME molecule close to a cyclic amidine, expelling water. The amine
/// nitrogen attacks the amide carbon, forming a C=N double bond while the amide
/// carbonyl oxygen leaves (with both amine hydrogens) as water. This is the
/// generic Traube/imidazole-forming step: applied to a 4,5-diaminopyrimidinedione
/// bearing a 5-formamido group it closes the imidazole of a purine — and because
/// `would_form_ring_of_size` is blind to existing rings, the fused bicycle (and
/// its aromatisation in `finish()`) falls out of the shared core for free, with
/// no purine/xanthine-specific code. Equally closes a free-standing imidazoline
/// or a benzimidazole from an ortho-formamido aniline.
pub(crate) fn generate_amidine_cyclization(
    amine_site: &AmineSite<'_>,
    amide_site: &AmideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Both groups must live on one molecule for an intramolecular closure.
    if amine_site.participant.substance.id != amide_site.participant.substance.id {
        return Ok(None);
    }
    // The amine must shed two N-H to become the imine-type (=N-) ring nitrogen.
    if amine_site.hydrogens.len() < 2 {
        return Ok(None);
    }
    let amine_desc = SiteDescriptorBuilder::from_amine_site(amine_site);
    let amide_desc = SiteDescriptorBuilder::from_amide_site(amide_site);
    // Nucleophile = amine nitrogen, electrophile = amide carbon. The amide oxygen
    // and both amine hydrogens leave as water; the new N=C double bond closes the
    // (aromatic, once Hückel is satisfied) amidine ring. Acid-catalysed, water-poor.
    close_intramolecular_ring(
        IntramolecularClosure {
            prefix: "amidine_cyclization",
            structure: amine_site.participant.structure,
            substance_id: amine_site.participant.substance.id.clone(),
            nucleophile: amine_site.nitrogen,
            electrophile: amide_site.carbon,
            nucleophile_label: "amidine imine nitrogen",
            electrophile_label: "amidine carbon",
            closure_bond_order: 2.0,
            leaving_atoms: vec![
                amide_site.carbonyl_oxygen,
                amine_site.hydrogens[0],
                amine_site.hydrogens[1],
            ],
            byproducts: vec!["destroy:water"],
            catalysts: vec!["destroy:proton"],
            base_activation_energy_kj_per_mol: 50.0,
            selectivity_profile: SelectivityProfile::new(
                ReactionType::HeterocycleCondensation,
                amine_desc,
            )
            .with_secondary_site(amide_desc),
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
fn generate_paal_knorr_external_heteroatom(
    first: &CarbonylSite<'_>,
    second: &CarbonylSite<'_>,
    donor_participant: &SiteParticipant<'_>,
    donor_atom: usize,
    donor_hydrogens: &[usize],
    donor_element: &'static str,
    donor_label: &'static str,
    reaction_prefix: &'static str,
    activation_energy_kj_per_mol: f64,
    donor_desc: crate::chemistry::selectivity::types::SiteDescriptor,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if first.participant.substance.id != second.participant.substance.id {
        return Ok(None);
    }
    if first.carbon >= second.carbon {
        return Ok(None);
    }
    if donor_hydrogens.len() < 2 {
        return Ok(None);
    }
    if donor_participant.substance.id == first.participant.substance.id {
        return Ok(None);
    }
    if donor_participant.structure.atoms[donor_atom].element != donor_element {
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
    let carbonyl_substance = first.participant.substance;
    let donor_substance = donor_participant.substance;

    let mut donor_editor = MolecularEditor::new(donor_participant.structure);
    let donor_mapping = donor_editor.remove_atoms(&[donor_hydrogens[0], donor_hydrogens[1]])?;
    let donor_atom = mapped_atom(&donor_mapping, donor_atom, donor_label)?;
    let donor_fragment = donor_editor.structure();

    let mut editor = MolecularEditor::new(structure);
    let mapping =
        editor.remove_atoms(&[first.oxygen, second.oxygen, beta_hydrogen, gamma_hydrogen])?;
    let carbon_a = mapped_atom(&mapping, first.carbon, "Paal-Knorr carbon a")?;
    let carbon_b = mapped_atom(&mapping, c_beta, "Paal-Knorr carbon b")?;
    let carbon_c = mapped_atom(&mapping, c_gamma, "Paal-Knorr carbon c")?;
    let carbon_d = mapped_atom(&mapping, second.carbon, "Paal-Knorr carbon d")?;
    let ring_heteroatom = editor.add_group(carbon_a, &donor_fragment, donor_atom, 1.0)?;
    editor.add_bond(ring_heteroatom, carbon_d, 1.0)?;
    editor.set_bond_order(carbon_a, carbon_b, 2.0)?;
    editor.set_bond_order(carbon_c, carbon_d, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_triple_site_reaction_id(
            reaction_prefix,
            &first.participant,
            &second.participant,
            donor_participant,
        ))
        .reactant(carbonyl_substance.id.clone(), 1, 1)
        .reactant(donor_substance.id.clone(), 1, 1)
        .catalyst_order("destroy:proton", 1)
        .product(product, 1)
        .product("destroy:water", 2)
        .condition(
            ReactionCondition::new(
                "Paal-Knorr heterocycle closure needs acidic, water-poor conditions",
            )
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.35),
        )
        .activation_energy_kj_per_mol(activation_energy_kj_per_mol)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::HeterocycleCondensation, first_desc)
                .with_secondary_site(donor_desc),
        )
        .build(),
    ))
}

pub(crate) fn generate_paal_knorr_pyrrole(
    first: &CarbonylSite<'_>,
    second: &CarbonylSite<'_>,
    amine: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    generate_paal_knorr_external_heteroatom(
        first,
        second,
        &amine.participant,
        amine.nitrogen,
        &amine.hydrogens,
        "N",
        "Paal-Knorr donor nitrogen",
        "paal_knorr_pyrrole",
        50.0,
        SiteDescriptorBuilder::from_amine_site(amine),
        resolver,
    )
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
    generate_paal_knorr_external_heteroatom(
        first,
        second,
        &thiol.participant,
        thiol.sulfur,
        &thiol.hydrogens,
        "S",
        "Paal-Knorr donor sulfur",
        "paal_knorr_thiophene",
        55.0,
        SiteDescriptorBuilder::from_thiol_site(thiol),
        resolver,
    )
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
    let Some((c1, _c2, _c3, c4)) = conjugated_diene_carbons(diene.participant.structure, diene)
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
        builder.reactant(diene_substance.id.clone(), 1, 1).reactant(
            dienophile_substance.id.clone(),
            1,
            1,
        )
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

pub(crate) fn generate_alkene_photocycloaddition(
    first_alkene: &UnsaturatedBondSite<'_>,
    second_alkene: &UnsaturatedBondSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if first_alkene.is_alkyne || second_alkene.is_alkyne {
        return Ok(None);
    }
    if !photocycloaddition_pair_is_canonical(first_alkene, second_alkene) {
        return Ok(None);
    }

    let first_substance = first_alkene.participant.substance;
    let second_substance = second_alkene.participant.substance;
    let a1 = first_alkene.high_degree_carbon;
    let a2 = first_alkene.low_degree_carbon;
    let b1 = second_alkene.high_degree_carbon;
    let b2 = second_alkene.low_degree_carbon;

    let mut editor = MolecularEditor::new(first_alkene.participant.structure);
    let b1_joined = editor.add_group(a1, second_alkene.participant.structure, b1, 1.0)?;
    let offset = b1_joined - b1;
    let b2_joined = b2 + offset;
    editor.add_bond(a2, b2_joined, 1.0)?;
    editor.set_bond_order(a1, a2, 1.0)?;
    editor.set_bond_order(b1_joined, b2_joined, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    let same_species = first_substance.id == second_substance.id;
    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "alkene_photocycloaddition",
        &first_alkene.participant,
        &second_alkene.participant,
    ));
    builder = if same_species {
        builder.reactant(first_substance.id.clone(), 2, 2)
    } else {
        builder.reactant(first_substance.id.clone(), 1, 1).reactant(
            second_substance.id.clone(),
            1,
            1,
        )
    };

    let first_desc = SiteDescriptorBuilder::from_unsaturated_bond_site(first_alkene);
    let second_desc = SiteDescriptorBuilder::from_unsaturated_bond_site(second_alkene);
    Ok(Some(
        builder
            .product(product, 1)
            .activation_energy_kj_per_mol(105.0)
            .requires_uv()
            .selectivity_profile(
                SelectivityProfile::new(ReactionType::PhotochemicalCycloaddition, first_desc)
                    .with_secondary_site(second_desc)
                    .never_suppress(),
            )
            .build(),
    ))
}

fn photocycloaddition_pair_is_canonical(
    first: &UnsaturatedBondSite<'_>,
    second: &UnsaturatedBondSite<'_>,
) -> bool {
    let first_key = alkene_pair_key(first);
    let second_key = alkene_pair_key(second);
    first_key <= second_key
}

fn alkene_pair_key(site: &UnsaturatedBondSite<'_>) -> (String, usize, usize) {
    (
        site.participant.substance.id.to_string(),
        site.high_degree_carbon.min(site.low_degree_carbon),
        site.high_degree_carbon.max(site.low_degree_carbon),
    )
}

pub(crate) fn generate_retro_diels_alder(
    alkene: &UnsaturatedBondSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    if alkene.is_alkyne {
        return Ok(Vec::new());
    }
    let structure = alkene.participant.structure;
    let c2 = alkene.high_degree_carbon;
    let c3 = alkene.low_degree_carbon;
    let mut reactions = Vec::new();
    for (c1, c2_c1_order) in structure.neighbors(c2) {
        if c1 == c3 || structure.atoms[c1].element != "C" || !bond_order_matches(c2_c1_order, 1.0) {
            continue;
        }
        for (c4, c3_c4_order) in structure.neighbors(c3) {
            if c4 == c2
                || c4 == c1
                || structure.atoms[c4].element != "C"
                || !bond_order_matches(c3_c4_order, 1.0)
            {
                continue;
            }
            for (ca, c1_ca_order) in structure.neighbors(c1) {
                if [c1, c2, c3, c4].contains(&ca)
                    || structure.atoms[ca].element != "C"
                    || !bond_order_matches(c1_ca_order, 1.0)
                {
                    continue;
                }
                for (cb, ca_cb_order) in structure.neighbors(ca) {
                    if [c1, c2, c3, c4, ca].contains(&cb)
                        || structure.atoms[cb].element != "C"
                        || !bond_order_matches(ca_cb_order, 1.0)
                    {
                        continue;
                    }
                    if !structure
                        .neighbors(cb)
                        .iter()
                        .any(|(n, order)| *n == c4 && bond_order_matches(*order, 1.0))
                    {
                        continue;
                    }
                    let Some((diene, dienophile)) =
                        retro_diels_alder_fragments(structure, [c1, c2, c3, c4], [ca, cb])
                    else {
                        continue;
                    };
                    let diene_id = resolver.resolve(diene)?;
                    let dienophile_id = resolver.resolve(dienophile)?;
                    let alkene_desc = SiteDescriptorBuilder::from_unsaturated_bond_site(alkene);
                    reactions.push(
                        Reaction::builder(format!(
                            "retro_diels_alder/{}/{c1}_{c2}_{c3}_{c4}_{ca}_{cb}",
                            sanitize_id(alkene.participant.substance.id.as_str())
                        ))
                        .reactant(alkene.participant.substance.id.clone(), 1, 1)
                        .product(diene_id, 1)
                        .product(dienophile_id, 1)
                        .activation_energy_kj_per_mol(125.0)
                        .selectivity_profile(
                            SelectivityProfile::new(ReactionType::RetroDielsAlder, alkene_desc)
                                .never_suppress(),
                        )
                        .build(),
                    );
                }
            }
        }
    }
    Ok(reactions)
}

pub(crate) fn generate_alkene_photoisomerization(
    alkene: &UnsaturatedBondSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if alkene.is_alkyne {
        return Ok(None);
    }
    let structure = alkene.participant.structure;
    let first = alkene.high_degree_carbon;
    let second = alkene.low_degree_carbon;
    let Some(first_substituent) = stereo_substituent(structure, first, second) else {
        return Ok(None);
    };
    let Some(second_substituent) = stereo_substituent(structure, second, first) else {
        return Ok(None);
    };
    if stereo_side_is_symmetric(structure, first, second)
        || stereo_side_is_symmetric(structure, second, first)
    {
        return Ok(None);
    }

    let mut variants = Vec::new();
    for (descriptor, suffix, activation_delta) in [
        (StereoDescriptor::E, "e", 0.0),
        (StereoDescriptor::Z, "z", 4.0),
    ] {
        let mut editor = MolecularEditor::new(structure);
        editor.set_double_bond_stereo(
            first,
            second,
            first_substituent,
            second_substituent,
            descriptor,
        )?;
        let product = resolver.resolve(editor.finish()?)?;
        variants.push((suffix, product, activation_delta));
    }
    if variants[0].1 == variants[1].1 {
        return Ok(None);
    }

    let alkene_desc = SiteDescriptorBuilder::from_unsaturated_bond_site(alkene);
    let mut builder = Reaction::builder(generated_site_reaction_id(
        "alkene_photoisomerization",
        &alkene.participant,
    ))
    .reactant(alkene.participant.substance.id.clone(), 1, 1)
    .activation_energy_kj_per_mol(65.0)
    .selectivity_profile(
        SelectivityProfile::new(
            ReactionType::PhotochemicalIsomerization,
            alkene_desc.clone(),
        )
        .never_suppress(),
    );
    for (suffix, product, activation_delta) in variants {
        builder = builder.channel(
            ReactionChannel::new(
                format!("photo_{suffix}"),
                [StoichiometricTerm::new(product, 1)],
                activation_delta,
            )
            .with_mode(ReactionChannelMode::Photochemical)
            .with_condition_effect(ChannelConditionEffect::Light {
                band: LightBand::Ultraviolet,
                minimum_power: 1.0e-9,
                multiplier: 1.0,
            })
            .with_selectivity_profile(
                SelectivityProfile::new(
                    ReactionType::PhotochemicalIsomerization,
                    alkene_desc.clone(),
                )
                .never_suppress(),
            ),
        );
    }
    Ok(Some(builder.build()))
}

fn retro_diels_alder_fragments(
    structure: &MolecularStructure,
    diene_atoms: [usize; 4],
    dienophile_atoms: [usize; 2],
) -> Option<(MolecularStructure, MolecularStructure)> {
    let mut edited = structure.clone();
    set_bond_order_in_structure(&mut edited, diene_atoms[0], diene_atoms[1], 2.0)?;
    set_bond_order_in_structure(&mut edited, diene_atoms[1], diene_atoms[2], 1.0)?;
    set_bond_order_in_structure(&mut edited, diene_atoms[2], diene_atoms[3], 2.0)?;
    set_bond_order_in_structure(&mut edited, dienophile_atoms[0], dienophile_atoms[1], 2.0)?;
    remove_bond_in_structure(&mut edited, diene_atoms[0], dienophile_atoms[0])?;
    remove_bond_in_structure(&mut edited, diene_atoms[3], dienophile_atoms[1])?;
    let components = connected_components(&edited);
    if components.len() != 2 {
        return None;
    }
    let first = substructure_for_atoms(&edited, &components[0])?;
    let second = substructure_for_atoms(&edited, &components[1])?;
    Some((first, second))
}

fn set_bond_order_in_structure(
    structure: &mut MolecularStructure,
    first: usize,
    second: usize,
    order: f64,
) -> Option<()> {
    structure
        .bonds
        .iter_mut()
        .find(|bond| {
            (bond.from == first && bond.to == second) || (bond.from == second && bond.to == first)
        })?
        .order = order;
    structure.stereochemistry.clear();
    Some(())
}

fn remove_bond_in_structure(
    structure: &mut MolecularStructure,
    first: usize,
    second: usize,
) -> Option<()> {
    let position = structure.bonds.iter().position(|bond| {
        (bond.from == first && bond.to == second) || (bond.from == second && bond.to == first)
    })?;
    structure.bonds.remove(position);
    structure.stereochemistry.clear();
    Some(())
}

fn connected_components(structure: &MolecularStructure) -> Vec<Vec<usize>> {
    let mut seen = vec![false; structure.atoms.len()];
    let mut components = Vec::new();
    for start in 0..structure.atoms.len() {
        if seen[start] {
            continue;
        }
        let mut stack = vec![start];
        seen[start] = true;
        let mut component = Vec::new();
        while let Some(atom) = stack.pop() {
            component.push(atom);
            for (neighbor, _) in structure.neighbors(atom) {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        component.sort_unstable();
        components.push(component);
    }
    components
}

fn substructure_for_atoms(
    structure: &MolecularStructure,
    atoms: &[usize],
) -> Option<MolecularStructure> {
    let mut mapping = vec![None; structure.atoms.len()];
    for (new_index, old_index) in atoms.iter().copied().enumerate() {
        mapping[old_index] = Some(new_index);
    }
    let bonds = structure
        .bonds
        .iter()
        .filter_map(|bond| {
            Some(MolecularBond {
                from: mapping[bond.from]?,
                to: mapping[bond.to]?,
                order: bond.order,
            })
        })
        .collect::<Vec<_>>();
    let fragment = MolecularStructure {
        source_code: "generated".to_string(),
        atoms: atoms
            .iter()
            .map(|atom| structure.atoms[*atom].clone())
            .collect(),
        bonds,
        stereochemistry: Vec::new(),
    };
    fragment.validate().ok()?;
    Some(fragment)
}

fn stereo_substituent(
    structure: &MolecularStructure,
    atom: usize,
    other_double_bond_atom: usize,
) -> Option<usize> {
    structure
        .neighbors(atom)
        .into_iter()
        .filter_map(|(neighbor, order)| {
            (neighbor != other_double_bond_atom && bond_order_matches(order, 1.0))
                .then_some(neighbor)
        })
        .max_by_key(|neighbor| stereo_priority(structure.atoms[*neighbor].element.as_str()))
}

fn stereo_side_is_symmetric(
    structure: &MolecularStructure,
    atom: usize,
    other_double_bond_atom: usize,
) -> bool {
    let substituents = structure
        .neighbors(atom)
        .into_iter()
        .filter_map(|(neighbor, order)| {
            (neighbor != other_double_bond_atom && bond_order_matches(order, 1.0))
                .then_some(neighbor)
        })
        .collect::<Vec<_>>();
    substituents.len() < 2
        || substituents
            .iter()
            .all(|atom| structure.atoms[*atom].element == structure.atoms[substituents[0]].element)
}

fn stereo_priority(element: &str) -> u8 {
    match element {
        "H" => 1,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Cl" => 17,
        "Br" => 35,
        "I" => 53,
        _ => 0,
    }
}
