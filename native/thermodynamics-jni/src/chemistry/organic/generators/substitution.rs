use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::MolecularEditor;
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
    let halide_ion = match structure.atoms[halogen].element.as_str() {
        "Cl" => "destroy:chloride",
        "F" => "destroy:fluoride",
        "I" => "destroy:iodide",
        _ => {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: generated_site_reaction_id(
                    "halide_hydroxide_substitution",
                    &site.participant,
                ),
                reason: "halide group does not contain a supported halogen".to_string(),
            })
        }
    };
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
