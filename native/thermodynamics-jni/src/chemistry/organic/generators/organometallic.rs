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
    types::{ReactionType, SelectivityProfile, SubstitutionDegree},
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
    if site.degree >= 3
        && !is_bridgehead_organometallic_carbon(site.participant.structure, site.carbon)
    {
        return Ok(None);
    }
    let halogen_element = site.participant.structure.atoms[site.halogen]
        .element
        .as_str();
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

fn is_bridgehead_organometallic_carbon(structure: &MolecularStructure, carbon: usize) -> bool {
    if structure.atoms[carbon].element != "C" || structure.carbon_degree(carbon) < 3 {
        return false;
    }
    let carbon_neighbors = structure
        .neighbors(carbon)
        .into_iter()
        .filter_map(|(neighbor, order)| {
            (structure.atoms[neighbor].element == "C"
                && crate::chemistry::molecule::bond_order_matches(order, 1.0))
            .then_some(neighbor)
        })
        .collect::<Vec<_>>();
    if carbon_neighbors.len() < 3 {
        return false;
    }
    for first_index in 0..carbon_neighbors.len() {
        for second_index in (first_index + 1)..carbon_neighbors.len() {
            if !path_exists_avoiding_atom(
                structure,
                carbon_neighbors[first_index],
                carbon_neighbors[second_index],
                carbon,
            ) {
                return false;
            }
        }
    }
    true
}

fn path_exists_avoiding_atom(
    structure: &MolecularStructure,
    start: usize,
    target: usize,
    excluded: usize,
) -> bool {
    let mut seen = vec![false; structure.atoms.len()];
    let mut stack = vec![start];
    seen[excluded] = true;
    while let Some(atom) = stack.pop() {
        if atom == target {
            return true;
        }
        if seen[atom] {
            continue;
        }
        seen[atom] = true;
        for (neighbor, _) in structure.neighbors(atom) {
            if !seen[neighbor] && structure.atoms[neighbor].element != "H" {
                stack.push(neighbor);
            }
        }
    }
    false
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
        ReactionCondition::new(
            "organometallic nitrile addition requires dry inert addition and wet workup",
        )
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

pub(crate) fn generate_organometallic_carboxylation(
    organometallic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if organometallic.site.kind == crate::chemistry::reactive_site::ReactiveSiteKind::Organocopper {
        return Ok(None);
    }
    let (organo_fragment, organo_carbon, residue) =
        organometallic_fragment(organometallic.structure, &organometallic.site)?;
    let carbon_dioxide = crate::chemistry::frowns::parse_frowns("O=C=O")?;
    let carbonyl_carbon = carbon_dioxide
        .atoms
        .iter()
        .enumerate()
        .find_map(|(index, atom)| (atom.element == "C").then_some(index))
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id(
                "organometallic_carboxylation",
                &organometallic,
            ),
            reason: "carbon dioxide graph has no carbon atom".to_string(),
        })?;
    let mut editor = MolecularEditor::new(&organo_fragment);
    let carbon_dioxide_mapping =
        editor.add_group_with_mapping(organo_carbon, &carbon_dioxide, carbonyl_carbon, 1.0)?;
    let carboxyl_carbon = carbon_dioxide_mapping[carbonyl_carbon];
    let hydroxyl_oxygen = carbon_dioxide
        .neighbors(carbonyl_carbon)
        .into_iter()
        .find_map(|(oxygen, order)| {
            (carbon_dioxide.atoms[oxygen].element == "O"
                && crate::chemistry::molecule::bond_order_matches(order, 2.0))
            .then_some(carbon_dioxide_mapping[oxygen])
        })
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id(
                "organometallic_carboxylation",
                &organometallic,
            ),
            reason: "carbon dioxide graph has no oxygen to protonate".to_string(),
        })?;
    editor.set_bond_order(carboxyl_carbon, hydroxyl_oxygen, 1.0)?;
    editor.add_atom(hydroxyl_oxygen, "H", 0.0, 1.0)?;
    let acid = editor.finish()?;
    let product = resolver.resolve(acid)?;

    let organometallic_desc =
        organometallic_site_descriptor(organometallic.structure, &organometallic.site)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "organometallic_carboxylation",
            &organometallic,
        ))
        .reactant(organometallic.substance.id.clone(), 1, 1)
        .reactant("destroy:carbon_dioxide", 1, 1)
        .reactant("destroy:water", 1, 1)
        .product(product, 1)
        .chemical_external_product(
            "organometallic carboxylation metal hydroxide salt",
            1.0,
            residue.mass + 17.01,
            residue.charge,
        )
        .condition(
            ReactionCondition::new(
                "organometallic carboxylation requires dry CO2 addition and wet workup",
            )
            .max_oxygen_activity(0.02)
            .atmosphere(AtmosphereCondition::Inert),
        )
        .activation_energy_kj_per_mol(16.0)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::CarbonylAddition, organometallic_desc)
                .with_nucleophile_strength(NucleophileStrength::VeryStrong),
        )
        .build(),
    ))
}

fn organometallic_site_descriptor(
    structure: &MolecularStructure,
    site: &crate::chemistry::reactive_site::ReactiveSite,
) -> ChemistryResult<crate::chemistry::selectivity::SiteDescriptor> {
    let (organo_carbon, _, _) = organometallic_atoms(structure, site)?;
    let degree = match structure.carbon_degree(organo_carbon) {
        0 | 1 => SubstitutionDegree::Primary,
        2 => SubstitutionDegree::Secondary,
        _ => SubstitutionDegree::Tertiary,
    };
    let has_beta_hydrogen =
        structure
            .neighbors(organo_carbon)
            .into_iter()
            .any(|(neighbor, order)| {
                structure.atoms[neighbor].element == "C"
                    && crate::chemistry::molecule::bond_order_matches(order, 1.0)
                    && structure.hydrogen_count(neighbor) > 0
            });
    let electron_withdrawing = structure
        .neighbors(organo_carbon)
        .into_iter()
        .filter(|(neighbor, _)| {
            matches!(
                structure.atoms[*neighbor].element.as_str(),
                "O" | "N" | "S" | "F" | "Cl" | "Br" | "I"
            )
        })
        .count() as u32;
    Ok(SiteDescriptorBuilder::build(
        site.kind.clone(),
        degree,
        structure.carbon_degree(organo_carbon).saturating_sub(1) as u32,
        electron_withdrawing,
        0,
        has_beta_hydrogen,
        false,
        false,
    ))
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
