use super::super::centers::*;
use super::super::resolver::DerivedSubstanceResolver;
use super::common::*;
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::kinetics::ReactionChannel;
use crate::chemistry::molecule::{
    MolecularEditor, MolecularStructure, StereoDescriptor, StereoMixtureKind, Stereochemistry,
    TetrahedralStereo,
};
use crate::chemistry::reaction::{Reaction, StoichiometricTerm};
use crate::chemistry::selectivity::{
    engine::SiteDescriptorBuilder,
    types::{ReactionType, SelectivityProfile},
};
use std::collections::BTreeSet;

#[derive(Debug, Clone, Copy)]
pub(crate) struct ElectrophilicAdditionSpec {
    prefix: &'static str,
    electrophile: &'static str,
    high_degree_group: AdditionGroup,
    low_degree_group: AdditionGroup,
    alkyne_stereo_rule: Option<AlkyneStereoRule>,
    nucleophile_ratio: u32,
    activation_energy: f64,
    catalyst: Option<(&'static str, u32)>,
    external_catalyst: Option<&'static str>,
    display_as_reversible: bool,
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum AdditionGroup {
    Atom(&'static str),
    Hydroxyl,
    Borane,
}

#[derive(Debug, Clone, Copy)]
pub(crate) enum AlkyneStereoRule {
    Anti,
}

#[derive(Debug, Clone)]
pub(crate) struct StereoProductVariant {
    pub(crate) structure: MolecularStructure,
    pub(crate) channel_suffix: String,
    pub(crate) activation_delta_kj_per_mol: f64,
    pub(crate) pre_exponential_factor_multiplier: f64,
}

pub(crate) fn electrophilic_addition_specs(alkyne: bool) -> Vec<ElectrophilicAdditionSpec> {
    let activation_energy = if alkyne { 10.0 } else { 25.0 };
    vec![
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_chlorination"
            } else {
                "alkene_chlorination"
            },
            electrophile: "destroy:chlorine",
            high_degree_group: AdditionGroup::Atom("Cl"),
            low_degree_group: AdditionGroup::Atom("Cl"),
            alkyne_stereo_rule: alkyne.then_some(AlkyneStereoRule::Anti),
            nucleophile_ratio: 1,
            activation_energy,
            catalyst: None,
            external_catalyst: None,
            display_as_reversible: false,
        },
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_chlorohydrination"
            } else {
                "alkene_chlorohydrination"
            },
            electrophile: "destroy:hypochlorous_acid",
            high_degree_group: AdditionGroup::Hydroxyl,
            low_degree_group: AdditionGroup::Atom("Cl"),
            alkyne_stereo_rule: alkyne.then_some(AlkyneStereoRule::Anti),
            nucleophile_ratio: 1,
            activation_energy,
            catalyst: None,
            external_catalyst: None,
            display_as_reversible: false,
        },
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_hydrolysis"
            } else {
                "alkene_hydrolysis"
            },
            electrophile: "destroy:water",
            high_degree_group: AdditionGroup::Hydroxyl,
            low_degree_group: AdditionGroup::Atom("H"),
            alkyne_stereo_rule: None,
            nucleophile_ratio: 1,
            activation_energy: 20.0,
            catalyst: Some(("destroy:proton", 2)),
            external_catalyst: None,
            display_as_reversible: true,
        },
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_borane_hydroboration"
            } else {
                "alkene_borane_hydroboration"
            },
            electrophile: "destroy:diborane",
            high_degree_group: AdditionGroup::Atom("H"),
            low_degree_group: AdditionGroup::Borane,
            alkyne_stereo_rule: None,
            nucleophile_ratio: 2,
            activation_energy,
            catalyst: None,
            external_catalyst: None,
            display_as_reversible: false,
        },
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_hydrochlorination"
            } else {
                "alkene_hydrochlorination"
            },
            electrophile: "destroy:hydrochloric_acid",
            high_degree_group: AdditionGroup::Atom("Cl"),
            low_degree_group: AdditionGroup::Atom("H"),
            alkyne_stereo_rule: None,
            nucleophile_ratio: 1,
            activation_energy,
            catalyst: None,
            external_catalyst: None,
            display_as_reversible: false,
        },
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_hydrogenation"
            } else {
                "alkene_hydrogenation"
            },
            electrophile: "destroy:hydrogen",
            high_degree_group: AdditionGroup::Atom("H"),
            low_degree_group: AdditionGroup::Atom("H"),
            alkyne_stereo_rule: None,
            nucleophile_ratio: 1,
            activation_energy,
            catalyst: None,
            external_catalyst: Some("forge:dusts/nickel"),
            display_as_reversible: false,
        },
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_hydroiodination"
            } else {
                "alkene_hydroiodination"
            },
            electrophile: "destroy:hydrogen_iodide",
            high_degree_group: AdditionGroup::Atom("I"),
            low_degree_group: AdditionGroup::Atom("H"),
            alkyne_stereo_rule: None,
            nucleophile_ratio: 1,
            activation_energy,
            catalyst: None,
            external_catalyst: None,
            display_as_reversible: false,
        },
        ElectrophilicAdditionSpec {
            prefix: if alkyne {
                "alkyne_iodination"
            } else {
                "alkene_iodination"
            },
            electrophile: "destroy:iodine",
            high_degree_group: AdditionGroup::Atom("I"),
            low_degree_group: AdditionGroup::Atom("I"),
            alkyne_stereo_rule: alkyne.then_some(AlkyneStereoRule::Anti),
            nucleophile_ratio: 1,
            activation_energy,
            catalyst: None,
            external_catalyst: None,
            display_as_reversible: false,
        },
    ]
}

pub(crate) fn generate_electrophilic_addition(
    site: &UnsaturatedBondSite<'_>,
    spec: ElectrophilicAdditionSpec,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let high_degree_carbon = site.high_degree_carbon;
    let low_degree_carbon = site.low_degree_carbon;
    let is_alkyne = site.is_alkyne;
    let mut editor = MolecularEditor::new(structure);
    editor.set_bond_order(
        high_degree_carbon,
        low_degree_carbon,
        if is_alkyne { 2.0 } else { 1.0 },
    )?;
    add_addition_group(&mut editor, high_degree_carbon, spec.high_degree_group)?;
    add_addition_group(&mut editor, low_degree_carbon, spec.low_degree_group)?;
    if is_alkyne {
        if let Some(rule) = spec.alkyne_stereo_rule {
            apply_alkyne_stereo_rule(&mut editor, high_degree_carbon, low_degree_carbon, rule)?;
        } else {
            editor
                .mark_double_bond_stereo_mixture_if_valid(high_degree_carbon, low_degree_carbon)?;
        }
    } else {
        editor.mark_tetrahedral_stereo_mixture_if_valid(high_degree_carbon)?;
        editor.mark_tetrahedral_stereo_mixture_if_valid(low_degree_carbon)?;
    }
    let product_variants = expand_stereo_product_distribution(editor.finish()?)?;
    let mut products = Vec::new();
    for variant in product_variants {
        products.push((
            resolver.resolve(variant.structure)?,
            variant.channel_suffix,
            variant.activation_delta_kj_per_mol,
            variant.pre_exponential_factor_multiplier,
        ));
    }
    let mut builder = Reaction::builder(generated_site_reaction_id(spec.prefix, &site.participant))
        .reactant(substance.id.clone(), spec.nucleophile_ratio, 1)
        .reactant(spec.electrophile, 1, 1)
        .activation_energy_kj_per_mol(spec.activation_energy);
    if products.len() == 1 {
        builder = builder.product(products.remove(0).0, spec.nucleophile_ratio);
    } else {
        for (product, suffix, activation_delta, pre_exponential_multiplier) in products {
            builder = builder.channel(
                ReactionChannel::new(
                    format!("{}:stereo:{}", spec.prefix, suffix),
                    [StoichiometricTerm::new(product, spec.nucleophile_ratio)],
                    spec.activation_energy + activation_delta,
                )
                .with_pre_exponential_factor(10_000.0 * pre_exponential_multiplier),
            );
        }
    }
    if let Some((catalyst, order)) = spec.catalyst {
        builder = builder.catalyst_order(catalyst, order);
    }
    if let Some(catalyst) = spec.external_catalyst {
        builder = builder.external_catalyst(catalyst, 1.0);
    }
    if spec.display_as_reversible {
        builder = builder.display_as_reversible();
    }
    Ok(builder.build())
}

pub(crate) fn generate_chain_growth_polymerization(
    site: &UnsaturatedBondSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    // Addition polymerization opens the C=C and links monomers into a chain. The
    // product is a polymer material carrying its repeat unit and representative
    // chain length, not the repeat unit pretending to be a standalone molecule.
    // Alkynes and dienes are out of scope (the perception returns None for >1 C=C).
    if site.is_alkyne {
        return Ok(None);
    }
    let substance = site.participant.substance;
    let Some((polymer, monomer_count)) = crate::chemistry::polymer::chain_growth_polymer_substance(
        substance.id.clone(),
        site.participant.structure,
    )?
    else {
        return Ok(None);
    };
    let product = resolver.resolve_substance(polymer)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "chain_growth_polymerization",
            &site.participant,
        ))
        .reactant(substance.id.clone(), monomer_count, 1)
        .product(product, 1)
        // Initiated radical/ionic chain growth is thermally promoted but has no
        // sharp cutoff, so the Arrhenius barrier alone gates it (no binary
        // temperature ban). Sized for a typical vinyl propagation.
        .activation_energy_kj_per_mol(40.0)
        .selectivity_profile(SelectivityProfile::new(
            ReactionType::ChainGrowthPolymerization,
            SiteDescriptorBuilder::from_unsaturated_bond_site(site),
        ))
        .build(),
    ))
}

pub(crate) fn apply_alkyne_stereo_rule(
    editor: &mut MolecularEditor,
    first: usize,
    second: usize,
    rule: AlkyneStereoRule,
) -> ChemistryResult<bool> {
    let structure = editor.structure();
    let Some((first_substituent, second_substituent)) =
        double_bond_stereo_substituents(&structure, first, second)
    else {
        return Ok(false);
    };
    match rule {
        AlkyneStereoRule::Anti => editor.set_double_bond_stereo(
            first,
            second,
            first_substituent,
            second_substituent,
            StereoDescriptor::Trans,
        )?,
    }
    Ok(true)
}

pub(crate) fn double_bond_stereo_substituents(
    structure: &MolecularStructure,
    first: usize,
    second: usize,
) -> Option<(usize, usize)> {
    let first_substituents = bonded_substituents_except(structure, first, second);
    let second_substituents = bonded_substituents_except(structure, second, first);
    if first_substituents.len() != 2 || second_substituents.len() != 2 {
        return None;
    }
    let first_substituent = preferred_stereo_substituent(structure, &first_substituents)?;
    let second_substituent = preferred_stereo_substituent(structure, &second_substituents)?;
    Some((first_substituent, second_substituent))
}

pub(crate) fn bonded_substituents_except(
    structure: &MolecularStructure,
    atom: usize,
    excluded: usize,
) -> Vec<usize> {
    structure
        .bonds
        .iter()
        .filter_map(|bond| {
            if bond.from == atom && bond.to != excluded {
                Some(bond.to)
            } else if bond.to == atom && bond.from != excluded {
                Some(bond.from)
            } else {
                None
            }
        })
        .collect()
}

pub(crate) fn preferred_stereo_substituent(
    structure: &MolecularStructure,
    substituents: &[usize],
) -> Option<usize> {
    substituents.iter().copied().max_by_key(|index| {
        let atom = &structure.atoms[*index];
        (atomic_stereo_priority(&atom.element), atom.r_group_number)
    })
}

pub(crate) fn atomic_stereo_priority(element: &str) -> u16 {
    match element {
        "H" => 1,
        "B" => 5,
        "C" => 6,
        "N" => 7,
        "O" => 8,
        "F" => 9,
        "Cl" => 17,
        "Br" => 35,
        "I" => 53,
        "R" => 200,
        _ => 0,
    }
}

pub(crate) fn expand_stereo_product_distribution(
    structure: MolecularStructure,
) -> ChemistryResult<Vec<StereoProductVariant>> {
    expand_stereo_product_distribution_with_parameters(structure, "single".to_string(), 0.0, 1.0)
}

pub(crate) fn expand_stereo_product_distribution_with_parameters(
    structure: MolecularStructure,
    suffix: String,
    activation_delta_kj_per_mol: f64,
    pre_exponential_factor_multiplier: f64,
) -> ChemistryResult<Vec<StereoProductVariant>> {
    let Some(position) = structure
        .stereochemistry
        .iter()
        .position(|stereo| matches!(stereo, Stereochemistry::Mixture { .. }))
    else {
        return Ok(vec![StereoProductVariant {
            structure,
            channel_suffix: suffix,
            activation_delta_kj_per_mol,
            pre_exponential_factor_multiplier,
        }]);
    };
    let mut base = structure;
    let mixture = base.stereochemistry.remove(position);
    let variants = match mixture {
        Stereochemistry::Mixture {
            atoms,
            kind: StereoMixtureKind::Tetrahedral,
        } => {
            if atoms.len() != 5 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: "<generated-organic>".to_string(),
                    reason:
                        "tetrahedral stereo mixture must contain one center and four substituents"
                            .to_string(),
                });
            }
            let substituents = [atoms[1], atoms[2], atoms[3], atoms[4]];
            vec![
                (
                    Stereochemistry::Tetrahedral(TetrahedralStereo {
                        center: atoms[0],
                        substituents,
                        descriptor: StereoDescriptor::Clockwise,
                    }),
                    "tetra_cw".to_string(),
                    0.0,
                    1.0,
                ),
                (
                    Stereochemistry::Tetrahedral(TetrahedralStereo {
                        center: atoms[0],
                        substituents,
                        descriptor: StereoDescriptor::CounterClockwise,
                    }),
                    "tetra_ccw".to_string(),
                    0.0,
                    1.0,
                ),
            ]
        }
        Stereochemistry::Mixture {
            atoms,
            kind: StereoMixtureKind::DoubleBond,
        } => {
            if atoms.len() != 4 {
                return Err(ChemistryError::InvalidReaction {
                    reaction_id: "<generated-organic>".to_string(),
                    reason: "double-bond stereo mixture must contain bond atoms and substituents"
                        .to_string(),
                });
            }
            let steric_penalty = geometric_z_steric_penalty_kj_per_mol(
                &base, atoms[0], atoms[1], atoms[2], atoms[3],
            );
            vec![
                (
                    Stereochemistry::DoubleBond(crate::chemistry::molecule::DoubleBondStereo {
                        first: atoms[0],
                        second: atoms[1],
                        first_substituent: atoms[2],
                        second_substituent: atoms[3],
                        descriptor: StereoDescriptor::E,
                    }),
                    "db_e".to_string(),
                    0.0,
                    1.0,
                ),
                (
                    Stereochemistry::DoubleBond(crate::chemistry::molecule::DoubleBondStereo {
                        first: atoms[0],
                        second: atoms[1],
                        first_substituent: atoms[2],
                        second_substituent: atoms[3],
                        descriptor: StereoDescriptor::Z,
                    }),
                    "db_z".to_string(),
                    steric_penalty,
                    z_pre_exponential_multiplier(steric_penalty),
                ),
            ]
        }
        Stereochemistry::Mixture {
            kind: StereoMixtureKind::General,
            ..
        } => {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<generated-organic>".to_string(),
                reason: "general stereo mixture has no quantitative distribution rule".to_string(),
            });
        }
        _ => unreachable!("position selected a stereo mixture"),
    };
    let mut result = Vec::new();
    for (stereo, variant_suffix, variant_activation_delta, variant_pre_exponential_multiplier) in
        variants
    {
        let mut variant = base.clone();
        variant.stereochemistry.push(stereo);
        variant.validate()?;
        result.extend(expand_stereo_product_distribution_with_parameters(
            variant,
            format!("{suffix}_{variant_suffix}"),
            activation_delta_kj_per_mol + variant_activation_delta,
            pre_exponential_factor_multiplier * variant_pre_exponential_multiplier,
        )?);
    }
    Ok(result)
}

pub(crate) fn geometric_z_steric_penalty_kj_per_mol(
    structure: &MolecularStructure,
    first: usize,
    second: usize,
    first_substituent: usize,
    second_substituent: usize,
) -> f64 {
    let first_bulk = substituent_steric_bulk(structure, first_substituent, first);
    let second_bulk = substituent_steric_bulk(structure, second_substituent, second);
    (1.5 + 0.35 * (first_bulk + second_bulk)).clamp(1.5, 8.0)
}

pub(crate) fn z_pre_exponential_multiplier(steric_penalty_kj_per_mol: f64) -> f64 {
    (1.0 - steric_penalty_kj_per_mol / 20.0).clamp(0.55, 0.95)
}

pub(crate) fn substituent_steric_bulk(
    structure: &MolecularStructure,
    substituent: usize,
    blocked_atom: usize,
) -> f64 {
    let mut visited = BTreeSet::new();
    substituent_steric_bulk_inner(structure, substituent, blocked_atom, &mut visited)
}

pub(crate) fn substituent_steric_bulk_inner(
    structure: &MolecularStructure,
    atom_index: usize,
    blocked_atom: usize,
    visited: &mut BTreeSet<usize>,
) -> f64 {
    if atom_index == blocked_atom || !visited.insert(atom_index) {
        return 0.0;
    }
    let atom = &structure.atoms[atom_index];
    let mut bulk = match atom.element.as_str() {
        "H" => 0.2,
        "B" | "C" | "N" | "O" | "F" => 1.0,
        "Cl" => 1.8,
        "Br" => 2.1,
        "I" => 2.4,
        "R" => 1.5,
        _ => 1.0,
    };
    for neighbor in bonded_substituents_except(structure, atom_index, blocked_atom) {
        bulk += 0.35 * substituent_steric_bulk_inner(structure, neighbor, atom_index, visited);
    }
    bulk
}
