use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use super::*;
use crate::chemistry::condition::{AcidityCondition, AtmosphereCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::kinetics::ReactionChannel;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::{Reaction, StoichiometricTerm};
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityContext, SelectivityProfile},
    NucleophileStrength,
};

pub(crate) fn generate_acetal_formation(
    carbonyl_site: &CarbonylSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Reaction> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(carbonyl_site);
    let base_ea = 25.0;

    let carbonyl = carbonyl_site.participant.substance;
    let carbonyl_structure = carbonyl_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let carbonyl_carbon = carbonyl_site.carbon;
    let carbonyl_oxygen = carbonyl_site.oxygen;
    let (alcohol_fragment, alcohol_oxygen) =
        deprotonated_alcohol_fragment(alcohol_site, "acetal formation")?;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_oxygen])?;
    let carbonyl_carbon = mapped_atom(&carbonyl_mapping, carbonyl_carbon, "carbonyl carbon")?;
    let mut product_editor = MolecularEditor::new(&carbonyl_editor.finish()?);
    product_editor.add_group(carbonyl_carbon, &alcohol_fragment, alcohol_oxygen, 1.0)?;
    product_editor.add_group(carbonyl_carbon, &alcohol_fragment, alcohol_oxygen, 1.0)?;
    product_editor.mark_tetrahedral_stereo_mixture_if_valid(carbonyl_carbon)?;
    let product_variants = expand_stereo_product_distribution(product_editor.finish()?)?;

    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "acetal_formation",
        &carbonyl_site.participant,
        &alcohol_site.participant,
    ))
    .reactant(carbonyl.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 2, 1)
    .catalyst_order("destroy:proton", 1)
    .condition(
        ReactionCondition::new("acetal formation requires acidic, water-poor conditions")
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.35),
    );
    if product_variants.len() == 1 {
        let product = resolver.resolve(
            product_variants
                .into_iter()
                .next()
                .expect("length checked")
                .structure,
        )?;
        builder = builder.product(product, 1);
        builder = builder.product("destroy:water", 1);
        builder = builder.activation_energy_kj_per_mol(base_ea);
        builder = builder.selectivity_profile(
            SelectivityProfile::new(ReactionType::AcetalFormation, carbonyl_desc)
                .with_secondary_site(SiteDescriptorBuilder::from_alcohol_site(alcohol_site))
                .with_nucleophile_strength(NucleophileStrength::Weak)
                .never_suppress(),
        );
    } else {
        for variant in product_variants {
            builder = builder.channel(
                ReactionChannel::new(
                    format!("acetal_formation:stereo:{}", variant.channel_suffix),
                    [
                        StoichiometricTerm::new(resolver.resolve(variant.structure)?, 1),
                        StoichiometricTerm::new("destroy:water", 1),
                    ],
                    base_ea + variant.activation_delta_kj_per_mol,
                )
                .with_selectivity_profile(
                    SelectivityProfile::new(ReactionType::AcetalFormation, carbonyl_desc.clone())
                        .with_secondary_site(SiteDescriptorBuilder::from_alcohol_site(alcohol_site))
                        .with_nucleophile_strength(NucleophileStrength::Weak)
                        .never_suppress(),
                ),
            );
        }
        return Ok(builder.build());
    }
    Ok(builder.build())
}

pub(crate) fn generate_imine_formation(
    carbonyl_site: &CarbonylSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Reaction> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(carbonyl_site);
    let base_ea = 25.0;

    let carbonyl = carbonyl_site.participant.substance;
    let carbonyl_structure = carbonyl_site.participant.structure;
    let amine = amine_site.participant.substance;
    let amine_structure = amine_site.participant.structure;
    let carbonyl_carbon = carbonyl_site.carbon;
    let carbonyl_oxygen = carbonyl_site.oxygen;
    let amine_nitrogen = amine_site.nitrogen;
    let hydrogens = &amine_site.hydrogens;
    if hydrogens.len() < 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_pair_site_reaction_id(
                "imine_formation",
                &carbonyl_site.participant,
                &amine_site.participant,
            ),
            reason: "primary amine must have two explicit hydrogens for imine formation"
                .to_string(),
        });
    }

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_oxygen])?;
    let carbonyl_carbon = mapped_atom(&carbonyl_mapping, carbonyl_carbon, "carbonyl carbon")?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let mut amine_editor = MolecularEditor::new(amine_structure);
    let amine_mapping = amine_editor.remove_atoms(&[hydrogens[0], hydrogens[1]])?;
    let amine_nitrogen = mapped_atom(&amine_mapping, amine_nitrogen, "amine nitrogen")?;
    let amine_fragment = amine_editor.finish()?;

    let product_structure = MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &amine_fragment,
        amine_nitrogen,
        2.0,
    )?;
    let product = resolver.resolve(product_structure)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "imine_formation",
        &carbonyl_site.participant,
        &amine_site.participant,
    ))
    .reactant(carbonyl.id.clone(), 1, 1)
    .reactant(amine.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .catalyst_order("destroy:proton", 1)
    .condition(
        ReactionCondition::new("imine formation requires acidic, water-poor conditions")
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.5),
    )
    .activation_energy_kj_per_mol(base_ea)
    .selectivity_profile(
        SelectivityProfile::new(ReactionType::CarbonylAddition, carbonyl_desc)
            .with_nucleophile_strength(NucleophileStrength::Moderate),
    )
    .build())
}

pub(crate) fn generate_hydrazone_formation(
    carbonyl_site: &CarbonylSite<'_>,
    nucleophile_site: &BisNucleophileCenter<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Option<Reaction>> {
    if nucleophile_site.class != BisNucleophileClass::HydrazineLike {
        return Ok(None);
    }

    let nucleophile_structure = nucleophile_site.participant.structure;
    let imine_nitrogen = [
        nucleophile_site.first_nucleophile,
        nucleophile_site.second_nucleophile,
    ]
    .into_iter()
    .find(|nitrogen| bonded_hydrogens(nucleophile_structure, *nitrogen).len() >= 2);
    let Some(imine_nitrogen) = imine_nitrogen else {
        return Ok(None);
    };
    let hydrogens = bonded_hydrogens(nucleophile_structure, imine_nitrogen);

    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(carbonyl_site);
    let carbonyl = carbonyl_site.participant.substance;
    let nucleophile = nucleophile_site.participant.substance;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_site.participant.structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_site.oxygen])?;
    let carbonyl_carbon = mapped_atom(
        &carbonyl_mapping,
        carbonyl_site.carbon,
        "hydrazone carbonyl carbon",
    )?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let mut nucleophile_editor = MolecularEditor::new(nucleophile_structure);
    let nucleophile_mapping = nucleophile_editor.remove_atoms(&[hydrogens[0], hydrogens[1]])?;
    let imine_nitrogen = mapped_atom(
        &nucleophile_mapping,
        imine_nitrogen,
        "hydrazone imine nitrogen",
    )?;
    let nucleophile_fragment = nucleophile_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &nucleophile_fragment,
        imine_nitrogen,
        2.0,
    )?)?;

    Ok(Some(
        Reaction::builder(generated_pair_site_reaction_id(
            "hydrazone_formation",
            &carbonyl_site.participant,
            &nucleophile_site.participant,
        ))
        .reactant(carbonyl.id.clone(), 1, 1)
        .reactant(nucleophile.id.clone(), 1, 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .catalyst_order("destroy:proton", 1)
        .condition(
            ReactionCondition::new("hydrazone formation requires acidic, water-poor conditions")
                .acidity(AcidityCondition::Acidic)
                .max_water_activity(0.5),
        )
        .activation_energy_kj_per_mol(26.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::CarbonylAddition, carbonyl_desc)
                .with_nucleophile_strength(NucleophileStrength::Moderate),
        )
        .build(),
    ))
}

pub(crate) fn generate_organometallic_carbonyl_addition(
    carbonyl: SiteParticipant<'_>,
    organometallic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Reaction> {
    let carbonyl_site = carbonyl.clone().carbonyl_site()?;
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(&carbonyl_site);
    let base_ea = 15.0;

    let (carbonyl_carbon, carbonyl_oxygen) = carbonyl_atoms_from_site(
        carbonyl.structure,
        &carbonyl.site,
        "organometallic addition",
    )?;
    let (organo_carbon, metal, residue_atoms) =
        organometallic_atoms(organometallic.structure, &organometallic.site)?;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl.structure);
    carbonyl_editor.set_bond_order(carbonyl_carbon, carbonyl_oxygen, 1.0)?;
    carbonyl_editor.add_atom(carbonyl_oxygen, "H", 0.0, 1.0)?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let mut organo_editor = MolecularEditor::new(organometallic.structure);
    let mapping = organo_editor.remove_atoms(&residue_atoms)?;
    let organo_carbon = mapped_atom(&mapping, organo_carbon, "organometallic carbon")?;
    let organo_fragment = organo_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &organo_fragment,
        organo_carbon,
        1.0,
    )?)?;
    let residue_mass = atom_mass_sum(organometallic.structure, &residue_atoms)?;
    let residue_charge = atom_charge_sum(organometallic.structure, &residue_atoms)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "organometallic_carbonyl_addition",
        &carbonyl,
        &organometallic,
    ))
    .reactant(carbonyl.substance.id.clone(), 1, 1)
    .reactant(organometallic.substance.id.clone(), 1, 1)
    .chemical_external_reactant("proton donor hydrogen", 1.0, 1.01, 0)
    .chemical_external_product(
        format!(
            "{} salt residue",
            organometallic.structure.atoms[metal].element
        ),
        1.0,
        residue_mass,
        residue_charge,
    )
    .product(product, 1)
    .condition(
        ReactionCondition::new("organometallic carbonyl addition requires dry inert conditions")
            .max_water_activity(0.02)
            .max_oxygen_activity(0.02)
            .atmosphere(AtmosphereCondition::Inert),
    )
    .activation_energy_kj_per_mol(base_ea)
    .selectivity_profile(
        SelectivityProfile::new(ReactionType::CarbonylAddition, carbonyl_desc)
            .with_nucleophile_strength(NucleophileStrength::VeryStrong),
    )
    .build())
}

pub(crate) fn generate_cyanide_nucleophilic_addition(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Reaction> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(site);
    let base_ea = 20.0;

    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let oxygen = site.oxygen;
    let mut editor = MolecularEditor::new(structure);
    editor.set_bond_order(carbon, oxygen, 1.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let nitrile_carbon = editor.add_atom(carbon, "C", 0.0, 1.0)?;
    editor.add_atom(nitrile_carbon, "N", 0.0, 3.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "cyanide_nucleophilic_addition",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen_cyanide", 1, 1)
    .catalyst_order("destroy:cyanide", 1)
    .product(product, 1)
    .activation_energy_kj_per_mol(base_ea)
    .selectivity_profile(
        SelectivityProfile::new(ReactionType::CarbonylAddition, carbonyl_desc)
            .with_nucleophile_strength(NucleophileStrength::Strong),
    )
    .build())
}

pub(crate) fn generate_borohydride_carbonyl_reduction(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Reaction> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(site);
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let oxygen = site.oxygen;

    let mut editor = MolecularEditor::new(structure);
    editor.set_bond_order(carbon, oxygen, 1.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;

    Ok(Reaction::builder(generated_site_reaction_id(
        "borohydride_carbonyl_reduction",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 4, 1)
    .reactant("destroy:borohydride", 1, 1)
    .reactant("destroy:water", 4, 1)
    .product(product, 4)
    .product("destroy:tetrahydroxyborate", 1)
    .condition(
        ReactionCondition::new(
            "borohydride selectively reduces aldehydes and ketones in protic workup",
        )
        .min_water_activity(0.05),
    )
    .activation_energy_kj_per_mol(16.0)
    .selectivity_profile(
        SelectivityProfile::new(ReactionType::CarbonylReduction, carbonyl_desc)
            .with_nucleophile_strength(NucleophileStrength::Strong)
            .never_suppress(),
    )
    .build())
}

pub(crate) fn generate_wolff_kishner_reduction(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
    _context: &SelectivityContext,
) -> ChemistryResult<Reaction> {
    let carbonyl_desc = SiteDescriptorBuilder::from_carbonyl_site(site);
    let base_ea = 30.0;

    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let oxygen = site.oxygen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[oxygen])?;
    let carbon = mapped_atom(&mapping, carbon, "carbonyl carbon")?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "wolff_kishner_reduction",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrazine", 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .product("destroy:nitrogen", 1)
    .activation_energy_kj_per_mol(base_ea)
    .selectivity_profile(
        SelectivityProfile::new(ReactionType::CarbonylAddition, carbonyl_desc)
            .with_nucleophile_strength(NucleophileStrength::Strong),
    )
    .build())
}

pub(crate) fn generate_keto_enol_tautomerization(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // An enol is a hydroxyl borne on a C=C carbon (e.g. the vinyl alcohol that
    // alkyne hydration produces). It tautomerizes to the carbonyl: the O-H proton
    // leaves the oxygen, the C-O bond becomes the carbonyl double bond, the C=C
    // collapses to a single bond, and the beta carbon picks up the proton. The
    // equilibrium lies far toward the keto form, so this is modeled as a one-way
    // conversion. Requiring a genuine order-2.0 C=C neighbor (not aromatic 1.5)
    // means phenols are correctly excluded — their aromatic ring is not an enol.
    let structure = site.participant.structure;
    let enol_carbon = site.carbon;
    let oxygen = site.oxygen;
    let hydroxyl_hydrogen = site.hydrogen;
    let Some(beta_carbon) =
        structure
            .neighbors(enol_carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 2.0))
                .then_some(neighbor)
            })
    else {
        return Ok(None);
    };
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydroxyl_hydrogen])?;
    let enol_carbon = mapped_atom(&mapping, enol_carbon, "tautomerization enol carbon")?;
    let oxygen = mapped_atom(&mapping, oxygen, "tautomerization oxygen")?;
    let beta_carbon = mapped_atom(&mapping, beta_carbon, "tautomerization beta carbon")?;
    editor.set_bond_order(enol_carbon, oxygen, 2.0)?;
    editor.set_bond_order(enol_carbon, beta_carbon, 1.0)?;
    editor.add_atom(beta_carbon, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "keto_enol_tautomerization",
            &site.participant,
        ))
        .reactant(site.participant.substance.id.clone(), 1, 1)
        .product(product, 1)
        .catalyst_order("destroy:proton", 1)
        .condition(ReactionCondition::new(
            "keto-enol tautomerization proceeds under acid or base catalysis",
        ))
        .activation_energy_kj_per_mol(20.0)
        .build(),
    ))
}
