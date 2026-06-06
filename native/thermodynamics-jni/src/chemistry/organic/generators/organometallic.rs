use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use crate::chemistry::condition::{AtmosphereCondition, ReactionCondition};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::molecule::{MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
    NucleophileStrength,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum OrganometallicFormationMetal {
    Magnesium,
    Lithium,
}

impl OrganometallicFormationMetal {
    fn element(self) -> &'static str {
        match self {
            Self::Magnesium => "Mg",
            Self::Lithium => "Li",
        }
    }

    fn external_metal_name(self) -> &'static str {
        match self {
            Self::Magnesium => "external:magnesium_metal",
            Self::Lithium => "external:lithium_metal",
        }
    }

    fn external_metal_mass(self) -> f64 {
        match self {
            Self::Magnesium => 24.31,
            Self::Lithium => 6.94,
        }
    }
}

pub(crate) fn generate_organometallic_formation(
    site: &HalideSite<'_>,
    metal: OrganometallicFormationMetal,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if site.degree >= 3 {
        return Ok(None);
    }
    let halogen_element = site.participant.structure.atoms[site.halogen].element.as_str();
    if !matches!(halogen_element, "Cl" | "Br" | "I") {
        return Ok(None);
    }

    let mut editor = MolecularEditor::new(site.participant.structure);
    let mapping = editor.remove_atoms(&[site.halogen])?;
    let carbon = mapped_atom(&mapping, site.carbon, "organometallic formation carbon")?;
    let metal_atom = editor.add_atom(carbon, metal.element(), 0.0, 1.0)?;
    if metal == OrganometallicFormationMetal::Magnesium {
        editor.add_atom(metal_atom, halogen_element, 0.0, 1.0)?;
    }
    let product = resolver.resolve(editor.finish()?)?;

    let mut builder = Reaction::builder(generated_site_reaction_id(
        match metal {
            OrganometallicFormationMetal::Magnesium => "organomagnesium_formation",
            OrganometallicFormationMetal::Lithium => "organolithium_formation",
        },
        &site.participant,
    ))
    .reactant(site.participant.substance.id.clone(), 1, 1)
    .chemical_external_reactant(
        metal.external_metal_name(),
        if metal == OrganometallicFormationMetal::Lithium {
            2.0
        } else {
            1.0
        },
        metal.external_metal_mass(),
        0,
    )
    .product(product, 1)
    .condition(
        ReactionCondition::new("organometallic formation requires dry inert conditions")
            .max_water_activity(0.01)
            .max_oxygen_activity(0.01)
            .atmosphere(AtmosphereCondition::Inert),
    )
    .activation_energy_kj_per_mol(match metal {
        OrganometallicFormationMetal::Magnesium => 28.0,
        OrganometallicFormationMetal::Lithium => 18.0,
    });

    if metal == OrganometallicFormationMetal::Lithium {
        builder = builder.chemical_external_product(
            format!("external:lithium_{halogen_element}_salt"),
            1.0,
            metal.external_metal_mass()
                + crate::chemistry::molecule::element_mass(halogen_element)?,
            0,
        );
    }
    Ok(Some(builder.build()))
}

pub(crate) fn generate_organometallic_nitrile_addition(
    nitrile: &NitrileSite<'_>,
    organometallic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let (organo_fragment, organo_carbon, residue) =
        organometallic_fragment(organometallic.structure, &organometallic.site)?;

    let mut nitrile_editor = MolecularEditor::new(nitrile.participant.structure);
    let mapping = nitrile_editor.remove_atoms(&[nitrile.nitrogen])?;
    let nitrile_carbon = mapped_atom(&mapping, nitrile.carbon, "nitrile carbon")?;
    nitrile_editor.add_atom(nitrile_carbon, "O", 0.0, 2.0)?;
    let carbonyl_fragment = nitrile_editor.finish()?;
    let product = resolver.resolve(MolecularEditor::join_structures(
        &carbonyl_fragment,
        nitrile_carbon,
        &organo_fragment,
        organo_carbon,
        1.0,
    )?)?;

    let nitrile_desc = SiteDescriptorBuilder::build(
        crate::chemistry::reactive_site::ReactiveSiteKind::Nitrile,
        crate::chemistry::selectivity::types::SubstitutionDegree::Primary,
        0,
        1,
        0,
        false,
        false,
        false,
    );
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "organometallic_nitrile_addition",
        &nitrile.participant,
        &organometallic,
    ))
    .reactant(nitrile.participant.substance.id.clone(), 1, 1)
    .reactant(organometallic.substance.id.clone(), 1, 1)
    .reactant("destroy:water", 2, 1)
    .product(product, 1)
    .product("destroy:ammonia", 1)
    .chemical_external_product(
        "organometallic hydrolysis salt residue",
        1.0,
        residue.mass + 17.01,
        residue.charge,
    )
    .condition(
        ReactionCondition::new("organometallic nitrile addition requires dry inert addition and wet workup")
            .atmosphere(AtmosphereCondition::Inert),
    )
    .activation_energy_kj_per_mol(20.0)
    .selectivity_profile(
        SelectivityProfile::new(ReactionType::CarbonylAddition, nitrile_desc)
            .with_nucleophile_strength(NucleophileStrength::VeryStrong),
    )
    .build())
}

pub(crate) fn generate_organometallic_epoxide_opening(
    epoxide: SiteParticipant<'_>,
    organometallic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let oxygen = epoxide
        .site
        .atoms
        .iter()
        .copied()
        .find(|atom| epoxide.structure.atoms[*atom].element == "O")
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_pair_site_reaction_id(
                "organometallic_epoxide_opening",
                &epoxide,
                &organometallic,
            ),
            reason: "epoxide site has no oxygen atom".to_string(),
        })?;
    let mut carbons = epoxide
        .site
        .atoms
        .iter()
        .copied()
        .filter(|atom| epoxide.structure.atoms[*atom].element == "C")
        .collect::<Vec<_>>();
    if carbons.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_pair_site_reaction_id(
                "organometallic_epoxide_opening",
                &epoxide,
                &organometallic,
            ),
            reason: "epoxide site must contain exactly two carbon atoms".to_string(),
        });
    }
    carbons.sort_by_key(|carbon| epoxide.structure.carbon_degree(*carbon));
    let attack_carbon = carbons[0];
    let (organo_fragment, organo_carbon, residue) =
        organometallic_fragment(organometallic.structure, &organometallic.site)?;

    let mut epoxide_editor = MolecularEditor::new(epoxide.structure);
    epoxide_editor.remove_bond(oxygen, attack_carbon)?;
    epoxide_editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let epoxide_fragment = epoxide_editor.finish()?;
    let product = resolver.resolve(MolecularEditor::join_structures(
        &epoxide_fragment,
        attack_carbon,
        &organo_fragment,
        organo_carbon,
        1.0,
    )?)?;

    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "organometallic_epoxide_opening",
        &epoxide,
        &organometallic,
    ))
    .reactant(epoxide.substance.id.clone(), 1, 1)
    .reactant(organometallic.substance.id.clone(), 1, 1)
    .chemical_external_reactant("proton donor hydrogen", 1.0, 1.01, 0)
    .product(product, 1)
    .chemical_external_product(
        "organometallic salt residue",
        1.0,
        residue.mass,
        residue.charge,
    )
    .condition(
        ReactionCondition::new("organometallic epoxide opening requires dry inert conditions")
            .max_water_activity(0.02)
            .max_oxygen_activity(0.02)
            .atmosphere(AtmosphereCondition::Inert),
    )
    .activation_energy_kj_per_mol(14.0)
    .build())
}

struct OrganometallicResidue {
    mass: f64,
    charge: i32,
}

fn organometallic_fragment(
    structure: &MolecularStructure,
    site: &crate::chemistry::reactive_site::ReactiveSite,
) -> ChemistryResult<(MolecularStructure, usize, OrganometallicResidue)> {
    let (organo_carbon, _, residue_atoms) = organometallic_atoms(structure, site)?;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&residue_atoms)?;
    let organo_carbon = mapped_atom(&mapping, organo_carbon, "organometallic carbon")?;
    let fragment = editor.finish()?;
    Ok((
        fragment,
        organo_carbon,
        OrganometallicResidue {
            mass: atom_mass_sum(structure, &residue_atoms)?,
            charge: atom_charge_sum(structure, &residue_atoms)?,
        },
    ))
}
