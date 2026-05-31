use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile, SubstitutionDegree},
    NucleophileStrength,
};

pub(crate) fn generate_aldol_addition(
    enol: SiteParticipant<'_>,
    acceptor: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let center = enol.alpha_carbon_center()?;
    let alpha_carbon = center.alpha_carbon;
    let alpha_hydrogen = center.alpha_hydrogens[0];
    let (acceptor_carbon, acceptor_oxygen) =
        carbonyl_atoms_from_site(acceptor.structure, &acceptor.site, "aldol addition")?;

    let mut donor_editor = MolecularEditor::new(center.participant.structure);
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
        &center.participant,
        &acceptor,
    ))
    .reactant(center.participant.substance.id.clone(), 1, 1)
    .reactant(acceptor.substance.id.clone(), 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .condition(
        ReactionCondition::new("aldol addition requires basic carbonyl enolization")
            .acidity(AcidityCondition::Basic)
            .max_temperature_kelvin(323.15),
    )
    .activation_energy_kj_per_mol(28.0)
    .selectivity_profile(alpha_selectivity_profile(
        ReactionType::AldolAddition,
        &center,
    ))
    .build())
}

pub(crate) fn generate_alpha_halogenation(
    center: &AlphaCarbonCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    let mut reactions = Vec::new();
    for spec in [
        AlphaHalogenationSpec {
            prefix: "alpha_chlorination",
            halogen_reactant: "destroy:chlorine",
            halogen_element: "Cl",
            acid_product: "destroy:hydrochloric_acid",
        },
        AlphaHalogenationSpec {
            prefix: "alpha_iodination",
            halogen_reactant: "destroy:iodine",
            halogen_element: "I",
            acid_product: "destroy:hydrogen_iodide",
        },
    ] {
        let alpha_hydrogen = center.alpha_hydrogens[0];
        let mut editor = MolecularEditor::new(center.participant.structure);
        let mapping = editor.remove_atoms(&[alpha_hydrogen])?;
        let alpha_carbon = mapped_atom(&mapping, center.alpha_carbon, "alpha carbon")?;
        editor.add_atom(alpha_carbon, spec.halogen_element, 0.0, 1.0)?;
        let product = resolver.resolve(editor.finish()?)?;
        reactions.push(
            Reaction::builder(generated_site_reaction_id(spec.prefix, &center.participant))
                .reactant(center.participant.substance.id.clone(), 1, 1)
                .reactant(spec.halogen_reactant, 1, 1)
                .product(product, 1)
                .product(spec.acid_product, 1)
                .catalyst_order("destroy:proton", 1)
                .condition(
                    ReactionCondition::new("alpha halogenation requires enolization")
                        .acidity(AcidityCondition::Acidic),
                )
                .activation_energy_kj_per_mol(30.0)
                .selectivity_profile(alpha_selectivity_profile(
                    ReactionType::AlphaHalogenation,
                    center,
                ))
                .build(),
        );
    }
    Ok(reactions)
}

pub(crate) fn generate_aldol_dehydration(
    center: &AlphaCarbonCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let Some(beta) = beta_hydroxy_neighbor(center)? else {
        return Ok(None);
    };
    let alpha_hydrogen = center.alpha_hydrogens[0];
    let mut editor = MolecularEditor::new(center.participant.structure);
    let mapping =
        editor.remove_atoms(&[alpha_hydrogen, beta.hydroxyl_oxygen, beta.hydroxyl_hydrogen])?;
    let alpha_carbon = mapped_atom(
        &mapping,
        center.alpha_carbon,
        "aldol dehydration alpha carbon",
    )?;
    let beta_carbon = mapped_atom(&mapping, beta.beta_carbon, "aldol dehydration beta carbon")?;
    editor.set_bond_order(alpha_carbon, beta_carbon, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "aldol_dehydration",
            &center.participant,
        ))
        .reactant(center.participant.substance.id.clone(), 1, 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .catalyst_order("destroy:proton", 1)
        .condition(
            ReactionCondition::new("aldol dehydration requires acid-catalyzed elimination")
                .acidity(AcidityCondition::Acidic)
                .rate_multiplier(0.5),
        )
        .activation_energy_kj_per_mol(32.0)
        .selectivity_profile(alpha_selectivity_profile(
            ReactionType::AldolDehydration,
            center,
        ))
        .build(),
    ))
}

pub(crate) fn generate_enamine_formation(
    carbonyl_site: &CarbonylSite<'_>,
    amine_site: &AmineSite<'_>,
    alpha_center: &AlphaCarbonCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if amine_site.hydrogens.len() != 1 {
        return Ok(None);
    }
    let alpha_hydrogen = alpha_center.alpha_hydrogens[0];
    let amine_hydrogen = amine_site.hydrogens[0];

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_site.participant.structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_site.oxygen, alpha_hydrogen])?;
    let carbonyl_carbon = mapped_atom(
        &carbonyl_mapping,
        carbonyl_site.carbon,
        "enamine carbonyl carbon",
    )?;
    let alpha_carbon = mapped_atom(
        &carbonyl_mapping,
        alpha_center.alpha_carbon,
        "enamine alpha carbon",
    )?;
    carbonyl_editor.set_bond_order(carbonyl_carbon, alpha_carbon, 2.0)?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let mut amine_editor = MolecularEditor::new(amine_site.participant.structure);
    let amine_mapping = amine_editor.remove_atoms(&[amine_hydrogen])?;
    let nitrogen = mapped_atom(&amine_mapping, amine_site.nitrogen, "enamine nitrogen")?;
    let amine_fragment = amine_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &amine_fragment,
        nitrogen,
        1.0,
    )?)?;
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "enamine_formation",
            &carbonyl_site.participant,
            &amine_site.participant,
        ))
        .reactant(carbonyl_site.participant.substance.id.clone(), 1, 1)
        .reactant(amine_site.participant.substance.id.clone(), 1, 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .catalyst_order("destroy:proton", 1)
        .condition(
            ReactionCondition::new(
                "enamine formation requires a secondary amine and water-poor acid",
            )
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.5),
        )
        .activation_energy_kj_per_mol(30.0)
        .selectivity_profile(alpha_selectivity_profile(
            ReactionType::EnamineFormation,
            alpha_center,
        ))
        .build(),
    ))
}

pub(crate) fn generate_enolate_alkylation(
    center: &AlphaCarbonCenter<'_>,
    halide_site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let halide_ion = halide_ion(
        halide_site.participant.structure,
        halide_site.halogen,
        "enolate_alkylation",
        &halide_site.participant,
    )?;
    let alpha_hydrogen = center.alpha_hydrogens[0];
    let mut donor_editor = MolecularEditor::new(center.participant.structure);
    let donor_mapping = donor_editor.remove_atoms(&[alpha_hydrogen])?;
    let alpha_carbon = mapped_atom(&donor_mapping, center.alpha_carbon, "enolate alpha carbon")?;
    let donor_fragment = donor_editor.finish()?;

    let mut halide_editor = MolecularEditor::new(halide_site.participant.structure);
    let halide_mapping = halide_editor.remove_atoms(&[halide_site.halogen])?;
    let alkyl_carbon = mapped_atom(&halide_mapping, halide_site.carbon, "alkyl halide carbon")?;
    let halide_fragment = halide_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &donor_fragment,
        alpha_carbon,
        &halide_fragment,
        alkyl_carbon,
        1.0,
    )?)?;
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "enolate_alkylation",
            &center.participant,
            &halide_site.participant,
        ))
        .reactant(center.participant.substance.id.clone(), 1, 1)
        .reactant(halide_site.participant.substance.id.clone(), 1, 1)
        .reactant("destroy:hydroxide", 1, 1)
        .product(product, 1)
        .product(halide_ion, 1)
        .product("destroy:water", 1)
        .condition(
            ReactionCondition::new("enolate alkylation requires basic enolate formation")
                .acidity(AcidityCondition::Basic),
        )
        .activation_energy_kj_per_mol(26.0)
        .selectivity_profile(
            alpha_selectivity_profile(ReactionType::EnolateAlkylation, center)
                .with_secondary_site(SiteDescriptorBuilder::from_halide_site(halide_site))
                .with_nucleophile_strength(NucleophileStrength::Strong),
        )
        .build(),
    ))
}

pub(crate) fn generate_michael_addition(
    center: &AlphaCarbonCenter<'_>,
    acceptor_site: &UnsaturatedBondSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if acceptor_site.is_alkyne {
        return Ok(None);
    }
    let Some(acceptor) = michael_acceptor_atoms(acceptor_site) else {
        return Ok(None);
    };
    let alpha_hydrogen = center.alpha_hydrogens[0];

    let mut donor_editor = MolecularEditor::new(center.participant.structure);
    let donor_mapping = donor_editor.remove_atoms(&[alpha_hydrogen])?;
    let donor_alpha = mapped_atom(
        &donor_mapping,
        center.alpha_carbon,
        "Michael donor alpha carbon",
    )?;
    let donor_fragment = donor_editor.finish()?;

    let mut acceptor_editor = MolecularEditor::new(acceptor_site.participant.structure);
    acceptor_editor.set_bond_order(acceptor.alpha_carbon, acceptor.beta_carbon, 1.0)?;
    acceptor_editor.add_atom(acceptor.alpha_carbon, "H", 0.0, 1.0)?;
    let acceptor_fragment = acceptor_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &donor_fragment,
        donor_alpha,
        &acceptor_fragment,
        acceptor.beta_carbon,
        1.0,
    )?)?;
    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "michael_addition",
            &center.participant,
            &acceptor_site.participant,
        ))
        .reactant(center.participant.substance.id.clone(), 1, 1)
        .reactant(acceptor_site.participant.substance.id.clone(), 1, 1)
        .catalyst_order("destroy:hydroxide", 1)
        .condition(
            ReactionCondition::new("Michael addition requires basic enolate formation")
                .acidity(AcidityCondition::Basic),
        )
        .product(product, 1)
        .activation_energy_kj_per_mol(28.0)
        .selectivity_profile(alpha_selectivity_profile(
            ReactionType::MichaelAddition,
            center,
        ))
        .build(),
    ))
}

pub(crate) fn generate_claisen_condensation(
    center: &AlphaCarbonCenter<'_>,
    ester_site: &EsterSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if !matches!(center.carbonyl_kind, AlphaCarbonylKind::Ester) {
        return Ok(None);
    }
    let alpha_hydrogen = center.alpha_hydrogens[0];
    let mut donor_editor = MolecularEditor::new(center.participant.structure);
    let donor_mapping = donor_editor.remove_atoms(&[alpha_hydrogen])?;
    let donor_alpha = mapped_atom(
        &donor_mapping,
        center.alpha_carbon,
        "Claisen donor alpha carbon",
    )?;
    let donor_fragment = donor_editor.finish()?;

    let (acyl_fragment, acyl_mapping, alkoxy_fragment, alkoxy_mapping) =
        MolecularEditor::split_at_bond(
            ester_site.participant.structure,
            ester_site.carbon,
            ester_site.alkoxy_oxygen,
        )?;
    let acyl_carbon = mapped_atom(&acyl_mapping, ester_site.carbon, "Claisen acyl carbon")?;
    let _acyl_oxygen = mapped_atom(
        &acyl_mapping,
        ester_site.carbonyl_oxygen,
        "Claisen carbonyl oxygen",
    )?;
    let alkoxy_oxygen = mapped_atom(
        &alkoxy_mapping,
        ester_site.alkoxy_oxygen,
        "Claisen alkoxy oxygen",
    )?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &donor_fragment,
        donor_alpha,
        &acyl_fragment,
        acyl_carbon,
        1.0,
    )?)?;
    let mut alcohol_editor = MolecularEditor::new(&alkoxy_fragment);
    alcohol_editor.add_atom(alkoxy_oxygen, "H", 0.0, 1.0)?;
    let alcohol = resolver.resolve(alcohol_editor.finish()?)?;

    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "claisen_condensation",
            &center.participant,
            &ester_site.participant,
        ))
        .reactant(center.participant.substance.id.clone(), 1, 1)
        .reactant(ester_site.participant.substance.id.clone(), 1, 1)
        .catalyst_order("destroy:hydroxide", 1)
        .condition(
            ReactionCondition::new("Claisen condensation requires basic ester enolate formation")
                .acidity(AcidityCondition::Basic),
        )
        .product(product, 1)
        .product(alcohol, 1)
        .activation_energy_kj_per_mol(31.0)
        .selectivity_profile(alpha_selectivity_profile(
            ReactionType::ClaisenCondensation,
            center,
        ))
        .build(),
    ))
}

struct AlphaHalogenationSpec {
    prefix: &'static str,
    halogen_reactant: &'static str,
    halogen_element: &'static str,
    acid_product: &'static str,
}

struct BetaHydroxyNeighbor {
    beta_carbon: usize,
    hydroxyl_oxygen: usize,
    hydroxyl_hydrogen: usize,
}

struct MichaelAcceptorAtoms {
    alpha_carbon: usize,
    beta_carbon: usize,
}

fn michael_acceptor_atoms(site: &UnsaturatedBondSite<'_>) -> Option<MichaelAcceptorAtoms> {
    let structure = site.participant.structure;
    let first = site.high_degree_carbon;
    let second = site.low_degree_carbon;
    let first_is_alpha = alkene_carbon_is_next_to_carbonyl(structure, first, second);
    let second_is_alpha = alkene_carbon_is_next_to_carbonyl(structure, second, first);
    match (first_is_alpha, second_is_alpha) {
        (true, false) => Some(MichaelAcceptorAtoms {
            alpha_carbon: first,
            beta_carbon: second,
        }),
        (false, true) => Some(MichaelAcceptorAtoms {
            alpha_carbon: second,
            beta_carbon: first,
        }),
        _ => None,
    }
}

fn alkene_carbon_is_next_to_carbonyl(
    structure: &crate::chemistry::molecule::MolecularStructure,
    alkene_carbon: usize,
    other_alkene_carbon: usize,
) -> bool {
    structure
        .neighbors(alkene_carbon)
        .into_iter()
        .any(|(neighbor, order)| {
            neighbor != other_alkene_carbon
                && structure.atoms[neighbor].element == "C"
                && crate::chemistry::molecule::bond_order_matches(order, 1.0)
                && structure
                    .neighbors(neighbor)
                    .into_iter()
                    .any(|(other, bond_order)| {
                        structure.atoms[other].element == "O"
                            && crate::chemistry::molecule::bond_order_matches(bond_order, 2.0)
                    })
        })
}

fn beta_hydroxy_neighbor(
    center: &AlphaCarbonCenter<'_>,
) -> ChemistryResult<Option<BetaHydroxyNeighbor>> {
    for (beta_carbon, alpha_bond_order) in
        center.participant.structure.neighbors(center.alpha_carbon)
    {
        if beta_carbon == center.carbonyl_carbon
            || center.participant.structure.atoms[beta_carbon].element != "C"
            || !crate::chemistry::molecule::bond_order_matches(alpha_bond_order, 1.0)
        {
            continue;
        }
        for (oxygen, order) in center.participant.structure.neighbors(beta_carbon) {
            if center.participant.structure.atoms[oxygen].element != "O"
                || !crate::chemistry::molecule::bond_order_matches(order, 1.0)
            {
                continue;
            }
            if let Some(hydroxyl_hydrogen) =
                first_bonded_hydrogen(center.participant.structure, oxygen)
            {
                return Ok(Some(BetaHydroxyNeighbor {
                    beta_carbon,
                    hydroxyl_oxygen: oxygen,
                    hydroxyl_hydrogen,
                }));
            }
        }
    }
    Ok(None)
}

fn alpha_selectivity_profile(
    reaction_type: ReactionType,
    center: &AlphaCarbonCenter<'_>,
) -> SelectivityProfile {
    let degree = match center.steric_class {
        AlphaStericClass::Primary => SubstitutionDegree::Primary,
        AlphaStericClass::Secondary => SubstitutionDegree::Secondary,
        AlphaStericClass::Tertiary => SubstitutionDegree::Tertiary,
    };
    let carbonyl_oxygen_is_present =
        center.participant.structure.atoms[center.carbonyl_oxygen].element == "O";
    let descriptor = SiteDescriptorBuilder::build(
        crate::chemistry::reactive_site::ReactiveSiteKind::Enol,
        degree,
        u32::from(matches!(center.conjugation, AlphaConjugation::Allylic)),
        match center.acidity {
            AlphaAcidityClass::Ordinary => u32::from(carbonyl_oxygen_is_present),
            AlphaAcidityClass::Activated => 1 + u32::from(carbonyl_oxygen_is_present),
        },
        match center.carbonyl_kind {
            AlphaCarbonylKind::Aldehyde => 0,
            AlphaCarbonylKind::Ketone => 1,
            AlphaCarbonylKind::Ester => 1,
        },
        !center.alpha_hydrogens.is_empty(),
        !matches!(center.conjugation, AlphaConjugation::None),
        matches!(center.conjugation, AlphaConjugation::Benzylic),
    );
    SelectivityProfile::new(reaction_type, descriptor).never_suppress()
}
