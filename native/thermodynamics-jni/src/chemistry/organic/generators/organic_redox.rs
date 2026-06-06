use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::condition::{AcidityCondition, ReactionCondition};
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::molecule::MolecularEditor;
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};

#[derive(Debug, Clone, Copy)]
enum OxidationReagent {
    AcidicDichromate,
    HydrogenPeroxide,
}

impl OxidationReagent {
    fn alcohol_activation_energy(self) -> f64 {
        match self {
            Self::AcidicDichromate => 25.0,
            Self::HydrogenPeroxide => 34.0,
        }
    }

    fn aldehyde_activation_energy(self) -> f64 {
        match self {
            Self::AcidicDichromate => 25.0,
            Self::HydrogenPeroxide => 28.0,
        }
    }

    fn epoxidation_activation_energy(self) -> f64 {
        match self {
            Self::AcidicDichromate => 0.0,
            Self::HydrogenPeroxide => 32.0,
        }
    }
}

pub(crate) fn generate_alcohol_oxidations(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    let Some(carbonyl_product) = oxidize_alcohol_to_carbonyl(site, resolver)? else {
        return Ok(Vec::new());
    };
    let mut reactions = Vec::new();
    reactions.push(alcohol_dichromate_oxidation(site, carbonyl_product.clone()));
    reactions.push(alcohol_peroxide_oxidation(
        site,
        carbonyl_product.clone(),
        false,
    ));
    if site.degree == 1 {
        if let Some(acid_product) = oxidize_primary_alcohol_to_acid(site, resolver)? {
            reactions.push(alcohol_peroxide_oxidation(site, acid_product, true));
        }
    }
    Ok(reactions)
}

pub(crate) fn generate_aldehyde_oxidations(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    let Some(product) = oxidize_aldehyde_to_acid(site, resolver)? else {
        return Ok(Vec::new());
    };
    Ok(vec![
        aldehyde_dichromate_oxidation(site, product.clone()),
        aldehyde_peroxide_oxidation(site, product),
    ])
}

pub(crate) fn generate_alkene_epoxidation(
    site: &UnsaturatedBondSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if site.is_alkyne {
        return Ok(None);
    }
    let substance = site.participant.substance;
    let mut editor = MolecularEditor::new(site.participant.structure);
    editor.set_bond_order(site.high_degree_carbon, site.low_degree_carbon, 1.0)?;
    let oxygen = editor.add_atom(site.high_degree_carbon, "O", 0.0, 1.0)?;
    editor.add_bond(oxygen, site.low_degree_carbon, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    let reagent = OxidationReagent::HydrogenPeroxide;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "alkene_epoxidation",
            &site.participant,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant("destroy:hydrogen_peroxide", 1, 1)
        .product(product, 1)
        .product("destroy:water", 1)
        .condition(ReactionCondition::new(
            "alkene epoxidation requires an oxygen-transfer oxidant",
        ))
        .activation_energy_kj_per_mol(reagent.epoxidation_activation_energy())
        .selectivity_profile(SelectivityProfile::new(
            ReactionType::OrganicOxidation,
            SiteDescriptorBuilder::from_unsaturated_bond_site(site),
        ))
        .build(),
    ))
}

fn oxidize_alcohol_to_carbonyl(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<crate::chemistry::substance::SubstanceId>> {
    if site.degree >= 3 {
        return Ok(None);
    }
    let structure = site.participant.structure;
    let Some(carbon_hydrogen) = first_bonded_hydrogen(structure, site.carbon) else {
        return Ok(None);
    };
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[carbon_hydrogen, site.hydrogen])?;
    let carbon = mapped_atom(&mapping, site.carbon, "alcohol carbon")?;
    let oxygen = mapped_atom(&mapping, site.oxygen, "alcohol oxygen")?;
    editor.set_bond_order(carbon, oxygen, 2.0)?;
    Ok(Some(resolver.resolve(editor.finish()?)?))
}

fn oxidize_primary_alcohol_to_acid(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<crate::chemistry::substance::SubstanceId>> {
    let structure = site.participant.structure;
    let carbon_hydrogens = bonded_hydrogens(structure, site.carbon);
    if carbon_hydrogens.len() < 2 {
        return Ok(None);
    }
    let second_carbon_hydrogen = carbon_hydrogens[1];
    let first_carbon_hydrogen = carbon_hydrogens[0];
    let mut editor = MolecularEditor::new(structure);
    let mapping =
        editor.remove_atoms(&[first_carbon_hydrogen, second_carbon_hydrogen, site.hydrogen])?;
    let carbon = mapped_atom(&mapping, site.carbon, "alcohol carbon")?;
    let oxygen = mapped_atom(&mapping, site.oxygen, "alcohol oxygen")?;
    editor.set_bond_order(carbon, oxygen, 2.0)?;
    add_hydroxyl(&mut editor, carbon)?;
    Ok(Some(resolver.resolve(editor.finish()?)?))
}

fn oxidize_aldehyde_to_acid(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<crate::chemistry::substance::SubstanceId>> {
    if site.is_ketone {
        return Ok(None);
    }
    let structure = site.participant.structure;
    let Some(hydrogen) = first_bonded_hydrogen(structure, site.carbon) else {
        return Ok(None);
    };
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let carbon = mapped_atom(&mapping, site.carbon, "aldehyde carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    Ok(Some(resolver.resolve(editor.finish()?)?))
}

fn alcohol_dichromate_oxidation(
    site: &AlcoholSite<'_>,
    product: crate::chemistry::substance::SubstanceId,
) -> Reaction {
    let reagent = OxidationReagent::AcidicDichromate;
    Reaction::builder(generated_site_reaction_id(
        "alcohol_oxidation",
        &site.participant,
    ))
    .reactant(site.participant.substance.id.clone(), 3, 1)
    .reactant("destroy:dichromate", 1, 1)
    .reactant("destroy:proton", 8, 1)
    .product(product, 3)
    .product("destroy:chromium_iii", 2)
    .product("destroy:water", 7)
    .condition(
        ReactionCondition::new("dichromate alcohol oxidation requires acidic conditions")
            .acidity(AcidityCondition::Acidic),
    )
    .activation_energy_kj_per_mol(reagent.alcohol_activation_energy())
    .selectivity_profile(SelectivityProfile::new(
        ReactionType::OrganicOxidation,
        SiteDescriptorBuilder::from_alcohol_site(site),
    ))
    .build()
}

fn alcohol_peroxide_oxidation(
    site: &AlcoholSite<'_>,
    product: crate::chemistry::substance::SubstanceId,
    overoxidized: bool,
) -> Reaction {
    let reagent = OxidationReagent::HydrogenPeroxide;
    Reaction::builder(generated_site_reaction_id(
        if overoxidized {
            "alcohol_peroxide_overoxidation"
        } else {
            "alcohol_peroxide_oxidation"
        },
        &site.participant,
    ))
    .reactant(site.participant.substance.id.clone(), 1, 1)
    .reactant(
        "destroy:hydrogen_peroxide",
        if overoxidized { 2 } else { 1 },
        1,
    )
    .product(product, 1)
    .product("destroy:water", if overoxidized { 3 } else { 2 })
    .condition(ReactionCondition::new(
        "peroxide alcohol oxidation requires an oxidizing medium",
    ))
    .activation_energy_kj_per_mol(if overoxidized {
        reagent.alcohol_activation_energy() + 4.0
    } else {
        reagent.alcohol_activation_energy()
    })
    .selectivity_profile(SelectivityProfile::new(
        ReactionType::OrganicOxidation,
        SiteDescriptorBuilder::from_alcohol_site(site),
    ))
    .build()
}

fn aldehyde_dichromate_oxidation(
    site: &CarbonylSite<'_>,
    product: crate::chemistry::substance::SubstanceId,
) -> Reaction {
    let reagent = OxidationReagent::AcidicDichromate;
    Reaction::builder(generated_site_reaction_id(
        "aldehyde_oxidation",
        &site.participant,
    ))
    .reactant(site.participant.substance.id.clone(), 3, 1)
    .reactant("destroy:dichromate", 1, 1)
    .reactant("destroy:proton", 8, 1)
    .product(product, 3)
    .product("destroy:chromium_iii", 2)
    .product("destroy:water", 4)
    .condition(
        ReactionCondition::new("dichromate aldehyde oxidation requires acidic conditions")
            .acidity(AcidityCondition::Acidic),
    )
    .activation_energy_kj_per_mol(reagent.aldehyde_activation_energy())
    .selectivity_profile(SelectivityProfile::new(
        ReactionType::OrganicOxidation,
        SiteDescriptorBuilder::from_carbonyl_site(site),
    ))
    .build()
}

fn aldehyde_peroxide_oxidation(
    site: &CarbonylSite<'_>,
    product: crate::chemistry::substance::SubstanceId,
) -> Reaction {
    let reagent = OxidationReagent::HydrogenPeroxide;
    Reaction::builder(generated_site_reaction_id(
        "aldehyde_peroxide_oxidation",
        &site.participant,
    ))
    .reactant(site.participant.substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen_peroxide", 1, 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .condition(ReactionCondition::new(
        "peroxide aldehyde oxidation requires an oxidizing medium",
    ))
    .activation_energy_kj_per_mol(reagent.aldehyde_activation_energy())
    .selectivity_profile(SelectivityProfile::new(
        ReactionType::OrganicOxidation,
        SiteDescriptorBuilder::from_carbonyl_site(site),
    ))
    .build()
}
