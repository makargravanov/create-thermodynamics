use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityContext, SelectivityProfile},
};

pub(crate) fn generate_halide_hydroxide_substitution(
    site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    let base_ea = 25.0;
    let halide_desc = SiteDescriptorBuilder::from_halide_site(site);

    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let halogen = site.halogen;
    // Delegate to the shared halide_ion mapping so every halogen the perception
    // layer can tag (incl. Br) yields its ion here — no local Cl/F/I subset that
    // silently diverges from centers.rs / common.rs.
    let halide_ion = halide_ion(
        structure,
        halogen,
        "halide_hydroxide_substitution",
        &site.participant,
    )?;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[halogen])?;
    let carbon = mapped_atom(&mapping, carbon, "halide substitution carbon")?;
    let oxygen = editor.add_atom(carbon, "O", 0.0, 1.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "halide_hydroxide_substitution",
            &site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant("destroy:hydroxide", 1, 1)
        .product(product, 1)
        .product(halide_ion, 1)
        .activation_energy_kj_per_mol(base_ea)
        .selectivity_profile(SelectivityProfile::new(ReactionType::SN2, halide_desc))
        .build(),
    ))
}

pub(crate) fn generate_halide_dehydrohalogenation(
    site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // E2 elimination: a base removes a beta hydrogen anti-periplanar to the
    // leaving halide, the C-X bond breaks, and the resulting electrons form the
    // new pi bond. Mechanistically the mirror of acid-catalyzed alcohol
    // dehydration, with halide (not water) as the leaving group. Every beta
    // carbon bearing an abstractable hydrogen is a distinct regiochemical
    // outcome (Zaitsev vs Hofmann), so we collect all elimination products into
    // a single Reaction — the selectivity layer then ranks them by the
    // mixture conditions.
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let halogen = site.halogen;
    let halide_ion = halide_ion(structure, halogen, "dehydrohalogenation", &site.participant)?;
    let mut products = Vec::new();
    for (beta_carbon, order) in structure.neighbors(carbon) {
        if structure.atoms[beta_carbon].element != "C" || !bond_order_matches(order, 1.0) {
            continue;
        }
        let Some(beta_hydrogen) = first_bonded_hydrogen(structure, beta_carbon) else {
            continue;
        };
        let mut editor = MolecularEditor::new(structure);
        let mapping = editor.remove_atoms(&[halogen, beta_hydrogen])?;
        let alpha = mapped_atom(&mapping, carbon, "dehydrohalogenation alpha carbon")?;
        let beta = mapped_atom(&mapping, beta_carbon, "dehydrohalogenation beta carbon")?;
        editor.set_bond_order(alpha, beta, 2.0)?;
        products.push(resolver.resolve(editor.finish()?)?);
    }
    if products.is_empty() {
        return Ok(None);
    }
    let mut builder = Reaction::builder(generated_site_reaction_id(
        "dehydrohalogenation",
        &site.participant,
    ))
    .reactant(
        site.participant.substance.id.clone(),
        products.len() as u32,
        1,
    )
    .reactant("destroy:hydroxide", products.len() as u32, 1)
    .product(halide_ion, products.len() as u32)
    .product("destroy:water", products.len() as u32)
    .condition(
        ReactionCondition::new("E2 elimination requires a base to abstract a beta proton")
            .acidity(AcidityCondition::Basic),
    )
    .activation_energy_kj_per_mol(30.0)
    .selectivity_profile(SelectivityProfile::new(
        ReactionType::E2,
        SiteDescriptorBuilder::from_halide_site(site),
    ));
    for product in products {
        builder = builder.product(product, 1);
    }
    Ok(Some(builder.build()))
}

pub(crate) fn generate_alkoxide_protonation(
    site: &AlkoxideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let oxygen = site.oxygen;
    let mut editor = MolecularEditor::new(structure);
    editor.replace_atom(oxygen, "O", 0.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "alkoxide_protonation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:proton", 1, 1)
    .product(product, 1)
    .build())
}

pub(crate) fn generate_thionyl_chloride_substitution(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let oxygen = site.oxygen;
    let proton = site.hydrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[oxygen, proton])?;
    let carbon = mapped_atom(&mapping, carbon, "alcohol carbon")?;
    editor.add_atom(carbon, "Cl", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "thionyl_chloride_substitution",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:thionyl_chloride", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .product("destroy:sulfur_dioxide", 1)
    .build())
}

#[derive(Clone, Copy)]
pub(crate) enum Hydrohalogen {
    Chloride,
    Bromide,
    Iodide,
}

impl Hydrohalogen {
    fn acid_id(self) -> &'static str {
        match self {
            Self::Chloride => "destroy:hydrochloric_acid",
            Self::Bromide => "destroy:hydrobromic_acid",
            Self::Iodide => "destroy:hydroiodic_acid",
        }
    }

    fn element(self) -> &'static str {
        match self {
            Self::Chloride => "Cl",
            Self::Bromide => "Br",
            Self::Iodide => "I",
        }
    }

    fn id_suffix(self) -> &'static str {
        match self {
            Self::Chloride => "chloride",
            Self::Bromide => "bromide",
            Self::Iodide => "iodide",
        }
    }

    fn activation_energy(self) -> f64 {
        match self {
            Self::Chloride => 34.0,
            Self::Bromide => 28.0,
            Self::Iodide => 24.0,
        }
    }
}

pub(crate) fn generate_alcohol_hydrohalogenation(
    site: &AlcoholSite<'_>,
    halogen: Hydrohalogen,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[site.oxygen, site.hydrogen])?;
    let carbon = mapped_atom(&mapping, site.carbon, "alcohol carbon")?;
    editor.add_atom(carbon, halogen.element(), 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(format!(
        "{}/{}",
        generated_site_reaction_id("alcohol_hydrohalogenation", &site.participant),
        halogen.id_suffix()
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant(halogen.acid_id(), 1, 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .condition(
        ReactionCondition::new(
            "alcohol hydrohalogenation requires acidic conditions and available halide",
        )
        .acidity(AcidityCondition::Acidic),
    )
    .activation_energy_kj_per_mol(halogen.activation_energy())
    .selectivity_profile(SelectivityProfile::new(
        if site.degree >= 3 {
            ReactionType::SN1
        } else {
            ReactionType::SN2
        },
        SiteDescriptorBuilder::from_alcohol_site(site),
    ))
    .build())
}

pub(crate) fn generate_halide_ammonia_substitution(
    site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    let base_ea = 25.0;
    let halide_desc = SiteDescriptorBuilder::from_halide_site(site);

    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let halogen = site.halogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[halogen])?;
    let carbon = mapped_atom(&mapping, carbon, "halide carbon")?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "halide_ammonia_substitution",
            &site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant("destroy:ammonia", 2, 2)
        .product(product, 1)
        .product(
            halide_ion(
                structure,
                halogen,
                "halide_ammonia_substitution",
                &site.participant,
            )?,
            1,
        )
        .product("destroy:ammonium", 1)
        .activation_energy_kj_per_mol(base_ea)
        .selectivity_profile(SelectivityProfile::new(ReactionType::SN2, halide_desc))
        .build(),
    ))
}

pub(crate) fn generate_halide_cyanide_substitution(
    site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    let base_ea = 25.0;
    let halide_desc = SiteDescriptorBuilder::from_halide_site(site);

    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let halogen = site.halogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[halogen])?;
    let carbon = mapped_atom(&mapping, carbon, "halide carbon")?;
    let nitrile_carbon = editor.add_atom(carbon, "C", 0.0, 1.0)?;
    editor.add_atom(nitrile_carbon, "N", 0.0, 3.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "halide_cyanide_substitution",
            &site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant("destroy:cyanide", 1, 1)
        .product(product, 1)
        .product(
            halide_ion(
                structure,
                halogen,
                "halide_cyanide_substitution",
                &site.participant,
            )?,
            1,
        )
        .activation_energy_kj_per_mol(base_ea)
        .selectivity_profile(SelectivityProfile::new(ReactionType::SN2, halide_desc))
        .build(),
    ))
}

pub(crate) fn generate_halide_amine_substitution(
    halide_site: &HalideSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    let base_ea = 25.0;
    let halide_desc = SiteDescriptorBuilder::from_halide_site(halide_site);

    // Reject sp2 (vinyl/aryl) halide carbons: an SN2 needs back-side attack on a
    // tetrahedral sp3 carbon, geometrically impossible at a planar sp2 carbon. The
    // functional-group detector tags vinyl halides as plain Halide sites, so this
    // generator must screen them out itself.
    let halide_carbon_structure = halide_site.participant.structure;
    if halide_carbon_structure
        .neighbors(halide_site.carbon)
        .into_iter()
        .any(|(n, order)| {
            bond_order_matches(order, 2.0) && halide_carbon_structure.atoms[n].element == "C"
        })
    {
        return Ok(None);
    }

    let halide = halide_site.participant.substance;
    let halide_structure = halide_site.participant.structure;
    let amine = amine_site.participant.substance;
    let amine_structure = amine_site.participant.structure;
    let halide_carbon = halide_site.carbon;
    let halogen = halide_site.halogen;
    let amine_nitrogen = amine_site.nitrogen;
    let amine_hydrogen =
        *amine_site
            .hydrogens
            .first()
            .ok_or_else(|| ChemistryError::InvalidReaction {
                reaction_id: generated_pair_site_reaction_id(
                    "halide_amine_substitution",
                    &halide_site.participant,
                    &amine_site.participant,
                ),
                reason: "amine has no explicit hydrogen".to_string(),
            })?;
    let mut halide_editor = MolecularEditor::new(halide_structure);
    let halide_mapping = halide_editor.remove_atoms(&[halogen])?;
    let halide_carbon = mapped_atom(&halide_mapping, halide_carbon, "halide carbon")?;
    let halide_fragment = halide_editor.finish()?;

    let mut amine_editor = MolecularEditor::new(amine_structure);
    let amine_mapping = amine_editor.remove_atoms(&[amine_hydrogen])?;
    let amine_nitrogen = mapped_atom(&amine_mapping, amine_nitrogen, "amine nitrogen")?;
    let amine_fragment = amine_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &halide_fragment,
        halide_carbon,
        &amine_fragment,
        amine_nitrogen,
        1.0,
    )?)?;
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "halide_amine_substitution",
            &halide_site.participant,
            &amine_site.participant,
        ))
        .reactant(halide.id.clone(), 1, 1)
        .reactant(amine.id.clone(), 1, 2)
        .product(product, 1)
        .product(
            halide_ion(
                halide_structure,
                halogen,
                "halide_amine_substitution",
                &halide_site.participant,
            )?,
            1,
        )
        .product("destroy:proton", 1)
        .activation_energy_kj_per_mol(base_ea)
        .selectivity_profile(SelectivityProfile::new(ReactionType::SN2, halide_desc))
        .build(),
    ))
}

/// N-alkylation of an amide/imide/lactam N-H by an alkyl halide. Mechanically the
/// same SN2 as `generate_halide_amine_substitution`, but the nucleophile is a weak
/// amide nitrogen rather than a basic amine: it must be deprotonated first, so the
/// reaction needs a basic medium and carries a higher barrier. Reuses the AmineSite
/// shape (N + hydrogens) via the widened `amine_site()` accessor. This is the step
/// that methylates xanthine's ring N-H groups to caffeine — generic over any cyclic
/// or acyclic amide N-H, with no caffeine-specific code.
pub(crate) fn generate_amide_n_alkylation(
    halide_site: &HalideSite<'_>,
    amide_nitrogen_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Reject sp2 (vinyl/aryl) halide carbons: an SN2 needs back-side attack on a
    // tetrahedral sp3 carbon, impossible at a planar sp2 carbon.
    let halide_structure = halide_site.participant.structure;
    if halide_structure
        .neighbors(halide_site.carbon)
        .into_iter()
        .any(|(n, order)| {
            bond_order_matches(order, 2.0) && halide_structure.atoms[n].element == "C"
        })
    {
        return Ok(None);
    }
    let Some(&nitrogen_hydrogen) = amide_nitrogen_site.hydrogens.first() else {
        return Ok(None);
    };
    let halide_desc = SiteDescriptorBuilder::from_halide_site(halide_site);
    let nitrogen_desc = SiteDescriptorBuilder::from_amine_site(amide_nitrogen_site);

    let halide = halide_site.participant.substance;
    let amide = amide_nitrogen_site.participant.substance;
    let amide_structure = amide_nitrogen_site.participant.structure;

    let mut halide_editor = MolecularEditor::new(halide_structure);
    let halide_mapping = halide_editor.remove_atoms(&[halide_site.halogen])?;
    let halide_carbon = mapped_atom(&halide_mapping, halide_site.carbon, "halide carbon")?;
    let halide_fragment = halide_editor.finish()?;

    let mut amide_editor = MolecularEditor::new(amide_structure);
    let amide_mapping = amide_editor.remove_atoms(&[nitrogen_hydrogen])?;
    let amide_nitrogen = mapped_atom(
        &amide_mapping,
        amide_nitrogen_site.nitrogen,
        "amide nitrogen",
    )?;
    let amide_fragment = amide_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &halide_fragment,
        halide_carbon,
        &amide_fragment,
        amide_nitrogen,
        1.0,
    )?)?;
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "amide_n_alkylation",
            &halide_site.participant,
            &amide_nitrogen_site.participant,
        ))
        .reactant(halide.id.clone(), 1, 1)
        .reactant(amide.id.clone(), 1, 1)
        .product(product, 1)
        .product(
            halide_ion(
                halide_structure,
                halide_site.halogen,
                "amide_n_alkylation",
                &halide_site.participant,
            )?,
            1,
        )
        .product("destroy:proton", 1)
        // Needs base to deprotonate the amide N-H before it can act as a nucleophile.
        .condition(
            ReactionCondition::new("amide N-alkylation requires a basic medium")
                .acidity(AcidityCondition::Basic),
        )
        // Higher barrier than a basic-amine SN2 (base 25): amide N is a far weaker
        // nucleophile even once deprotonated.
        .activation_energy_kj_per_mol(45.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::NAlkylation, halide_desc)
                .with_secondary_site(nitrogen_desc),
        )
        .build(),
    ))
}
