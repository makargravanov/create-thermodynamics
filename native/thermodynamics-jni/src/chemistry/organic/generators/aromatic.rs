use super::super::resolver::DerivedSubstanceResolver;
use super::super::space::SiteParticipant;
use super::common::*;
use crate::chemistry::condition::ReactionCondition;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::kinetics::ReactionChannel;
use crate::chemistry::molecule::{MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::{Reaction, StoichiometricTerm};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SubstituentClass {
    StronglyActivating,   // -OH, -NH2
    ModeratelyActivating, // -OR, -NHCOR, -OCOR
    WeaklyActivating,     // -R, -Ar
    WeaklyDeactivating,   // -F, -Cl, -Br, -I
    ModeratelyDeactivating,// -CHO, -COR, -COOH, -COOR, -SO3H
    StronglyDeactivating, // -NO2, -CF3, -CN, -NH3+
}

impl SubstituentClass {
    fn effect_for_distance(self, distance: usize) -> f64 {
        match self {
            Self::StronglyActivating => match distance {
                1 => -15.0,
                2 => -3.0,
                _ => -16.0,
            },
            Self::ModeratelyActivating => match distance {
                1 => -10.0,
                2 => -2.0,
                _ => -11.0,
            },
            Self::WeaklyActivating => match distance {
                1 => -4.0,
                2 => -1.0,
                _ => -5.0,
            },
            Self::WeaklyDeactivating => match distance {
                1 => 2.0,
                2 => 6.0,
                _ => 1.5,
            },
            Self::ModeratelyDeactivating => match distance {
                1 => 8.0,
                2 => 4.0,
                _ => 8.5,
            },
            Self::StronglyDeactivating => match distance {
                1 => 15.0,
                2 => 9.0,
                _ => 16.0,
            },
        }
    }

    fn is_deactivating(self) -> bool {
        matches!(self, Self::ModeratelyDeactivating | Self::StronglyDeactivating)
    }
}

fn classify_substituent(
    structure: &MolecularStructure,
    _ring_carbon: usize,
    sub_atom: usize,
) -> SubstituentClass {
    let atom_element = structure.atoms[sub_atom].element.as_str();
    match atom_element {
        "O" => {
            let has_h = structure.neighbors(sub_atom).iter().any(|(n, _)| {
                structure.atoms[*n].element == "H"
            });
            if has_h {
                SubstituentClass::StronglyActivating
            } else {
                let is_ester = structure.neighbors(sub_atom).iter().any(|(n, _)| {
                    structure.atoms[*n].element == "C" && has_double_bonded_oxygen(structure, *n)
                });
                if is_ester {
                    SubstituentClass::ModeratelyActivating
                } else {
                    SubstituentClass::ModeratelyActivating
                }
            }
        }
        "N" => {
            if structure.atoms[sub_atom].charge > 0.1 {
                SubstituentClass::StronglyDeactivating
            } else {
                let is_amide = structure.neighbors(sub_atom).iter().any(|(n, _)| {
                    structure.atoms[*n].element == "C" && has_double_bonded_oxygen(structure, *n)
                });
                if is_amide {
                    SubstituentClass::ModeratelyActivating
                } else {
                    SubstituentClass::StronglyActivating
                }
            }
        }
        "C" => {
            if has_double_bonded_oxygen(structure, sub_atom) {
                SubstituentClass::ModeratelyDeactivating
            } else if has_triple_bonded_nitrogen(structure, sub_atom) {
                SubstituentClass::StronglyDeactivating
            } else if halogen_count(structure, sub_atom) >= 3 {
                SubstituentClass::StronglyDeactivating
            } else {
                SubstituentClass::WeaklyActivating
            }
        }
        "F" | "Cl" | "Br" | "I" => SubstituentClass::WeaklyDeactivating,
        "S" => {
            let double_bonded_oxygens = structure.neighbors(sub_atom).iter().filter(|(n, order)| {
                structure.atoms[*n].element == "O" && crate::chemistry::molecule::bond_order_matches(*order, 2.0)
            }).count();
            if double_bonded_oxygens >= 2 {
                SubstituentClass::ModeratelyDeactivating
            } else {
                SubstituentClass::WeaklyActivating
            }
        }
        _ => SubstituentClass::WeaklyActivating,
    }
}

fn has_double_bonded_oxygen(structure: &MolecularStructure, carbon: usize) -> bool {
    structure.neighbors(carbon).iter().any(|(n, order)| {
        structure.atoms[*n].element == "O" && crate::chemistry::molecule::bond_order_matches(*order, 2.0)
    })
}

fn has_triple_bonded_nitrogen(structure: &MolecularStructure, carbon: usize) -> bool {
    structure.neighbors(carbon).iter().any(|(n, order)| {
        structure.atoms[*n].element == "N" && crate::chemistry::molecule::bond_order_matches(*order, 3.0)
    })
}

fn halogen_count(structure: &MolecularStructure, carbon: usize) -> usize {
    structure.neighbors(carbon).iter().filter(|(n, _)| {
        matches!(structure.atoms[*n].element.as_str(), "F" | "Cl" | "Br" | "I")
    }).count()
}

fn ring_distance(
    structure: &MolecularStructure,
    ring_atoms: &[usize],
    start: usize,
    end: usize,
) -> Option<usize> {
    if start == end {
        return Some(0);
    }
    let mut queue = std::collections::VecDeque::new();
    let mut visited = std::collections::BTreeMap::new();
    queue.push_back((start, 0));
    visited.insert(start, 0);

    while let Some((curr, dist)) = queue.pop_front() {
        if curr == end {
            return Some(dist);
        }
        for (neighbor, order) in structure.neighbors(curr) {
            if ring_atoms.contains(&neighbor) && crate::chemistry::molecule::bond_order_matches(order, 1.5) {
                if !visited.contains_key(&neighbor) {
                    visited.insert(neighbor, dist + 1);
                    queue.push_back((neighbor, dist + 1));
                }
            }
        }
    }
    None
}

fn is_substituent_bulky(structure: &MolecularStructure, sub_atom: usize) -> bool {
    let element = structure.atoms[sub_atom].element.as_str();
    if element == "Br" || element == "I" {
        return true;
    }
    if element == "C" {
        let non_h_neighbors = structure.neighbors(sub_atom).iter().filter(|(n, _)| {
            structure.atoms[*n].element != "H"
        }).count();
        return non_h_neighbors >= 3;
    }
    false
}

pub(crate) fn compute_eas_activation_delta(
    structure: &MolecularStructure,
    ring_atoms: &[usize],
    target_carbon: usize,
) -> f64 {
    let mut total_delta = 0.0;
    for &ring_car in ring_atoms {
        for (neighbor, order) in structure.neighbors(ring_car) {
            if !ring_atoms.contains(&neighbor) && structure.atoms[neighbor].element != "H" && !crate::chemistry::molecule::bond_order_matches(order, 1.5) {
                let sub_class = classify_substituent(structure, ring_car, neighbor);
                if let Some(dist) = ring_distance(structure, ring_atoms, ring_car, target_carbon) {
                    let mut effect = sub_class.effect_for_distance(dist);
                    if dist == 1 && is_substituent_bulky(structure, neighbor) {
                        effect += 3.0;
                    }
                    total_delta += effect;
                }
            }
        }
    }
    total_delta
}

pub(crate) fn is_ring_deactivated_for_fc(
    structure: &MolecularStructure,
    ring_atoms: &[usize],
) -> bool {
    for &ring_car in ring_atoms {
        for (neighbor, order) in structure.neighbors(ring_car) {
            if !ring_atoms.contains(&neighbor) && structure.atoms[neighbor].element != "H" && !crate::chemistry::molecule::bond_order_matches(order, 1.5) {
                let sub_class = classify_substituent(structure, ring_car, neighbor);
                if sub_class.is_deactivating() || matches!(sub_class, SubstituentClass::WeaklyDeactivating) {
                    return true;
                }
            }
        }
    }
    false
}

pub(crate) fn generate_eas_reaction<F>(
    aromatic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
    reaction_id_prefix: &str,
    reagents: &[(&str, u32)],
    catalysts: &[&str],
    base_ea: f64,
    is_friedel_crafts: bool,
    condition: Option<ReactionCondition>,
    editor_transform: F,
    co_products: &[(&str, u32)],
) -> ChemistryResult<Option<Reaction>>
where
    F: Fn(&mut MolecularEditor, usize) -> ChemistryResult<()>,
{
    let ring_atoms = &aromatic.site.atoms;
    if is_friedel_crafts && is_ring_deactivated_for_fc(aromatic.structure, ring_atoms) {
        return Ok(None);
    }

    let mut variants = Vec::new();
    for &carbon in ring_atoms {
        if aromatic.structure.atoms[carbon].element != "C" {
            continue;
        }
        let Some(hydrogen) = first_bonded_hydrogen(aromatic.structure, carbon) else {
            continue;
        };
        let mut editor = MolecularEditor::new(aromatic.structure);
        let mapping = editor.remove_atoms(&[hydrogen])?;
        let carbon_mapped = mapped_atom(&mapping, carbon, "aromatic substitution carbon")?;
        editor_transform(&mut editor, carbon_mapped)?;
        let product = resolver.resolve(editor.finish()?)?;
        let activation_delta = compute_eas_activation_delta(aromatic.structure, ring_atoms, carbon);
        variants.push((product, activation_delta));
    }

    if variants.is_empty() {
        return Ok(None);
    }

    let mut builder = Reaction::builder(generated_site_reaction_id(reaction_id_prefix, &aromatic))
        .reactant(aromatic.substance.id.clone(), 1, 1);
    for &(reagent_id, coef) in reagents {
        builder = builder.reactant(reagent_id, coef, coef);
    }
    for &cat_id in catalysts {
        builder = builder.catalyst_order(cat_id, 1);
    }
    if let Some(cond) = condition {
        builder = builder.condition(cond);
    }

    if variants.len() == 1 {
        builder = builder.product(variants[0].0.clone(), 1);
        for &(co_prod, coef) in co_products {
            builder = builder.product(co_prod, coef);
        }
        builder = builder.activation_energy_kj_per_mol((base_ea + variants[0].1).max(5.0));
    } else {
        for (index, (product, activation_delta)) in variants.into_iter().enumerate() {
            let mut terms = vec![StoichiometricTerm::new(product, 1)];
            for &(co_prod, coef) in co_products {
                terms.push(StoichiometricTerm::new(co_prod, coef));
            }
            builder = builder.channel(ReactionChannel::new(
                format!("{reaction_id_prefix}:position_{index}"),
                terms,
                (base_ea + activation_delta).max(5.0),
            ));
        }
    }
    Ok(Some(builder.build()))
}

pub(crate) fn generate_aromatic_nitration(
    aromatic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    use crate::chemistry::condition::AcidityCondition;
    let condition = Some(
        ReactionCondition::new("aromatic nitration requires strongly acidic conditions")
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.65),
    );
    generate_eas_reaction(
        aromatic,
        resolver,
        "aromatic_nitration",
        &[("destroy:nitric_acid", 1)],
        &["destroy:sulfuric_acid"],
        30.0,
        false,
        condition,
        |editor, carbon| {
            let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
            editor.add_atom(nitrogen, "O", 0.0, 1.5)?;
            editor.add_atom(nitrogen, "O", 0.0, 1.5)?;
            Ok(())
        },
        &[("destroy:water", 1)],
    )
}

pub(crate) fn generate_aromatic_chlorination(
    aromatic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    generate_eas_reaction(
        aromatic,
        resolver,
        "aromatic_chlorination",
        &[("destroy:chlorine", 1)],
        &["destroy:ferric_chloride"],
        28.0,
        false,
        None,
        |editor, carbon| {
            editor.add_atom(carbon, "Cl", 0.0, 1.0)?;
            Ok(())
        },
        &[("destroy:hydrochloric_acid", 1)],
    )
}

pub(crate) fn generate_aromatic_bromination(
    aromatic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    generate_eas_reaction(
        aromatic,
        resolver,
        "aromatic_bromination",
        &[("destroy:bromine", 1)],
        &["destroy:ferric_bromide"],
        28.0,
        false,
        None,
        |editor, carbon| {
            editor.add_atom(carbon, "Br", 0.0, 1.0)?;
            Ok(())
        },
        &[("destroy:hydrobromic_acid", 1)],
    )
}

pub(crate) fn generate_aromatic_sulfonation(
    aromatic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    generate_eas_reaction(
        aromatic,
        resolver,
        "aromatic_sulfonation",
        &[("destroy:sulfuric_acid", 1)],
        &[],
        32.0,
        false,
        None,
        |editor, carbon| {
            let sulfur = editor.add_atom(carbon, "S", 0.0, 1.0)?;
            editor.add_atom(sulfur, "O", 0.0, 2.0)?;
            editor.add_atom(sulfur, "O", 0.0, 2.0)?;
            let oxygen = editor.add_atom(sulfur, "O", 0.0, 1.0)?;
            editor.add_atom(oxygen, "H", 0.0, 1.0)?;
            Ok(())
        },
        &[("destroy:water", 1)],
    )
}

pub(crate) fn generate_fc_alkylation(
    aromatic: SiteParticipant<'_>,
    halide: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let ring_atoms = &aromatic.site.atoms;
    if is_ring_deactivated_for_fc(aromatic.structure, ring_atoms) {
        return Ok(None);
    }

    let halide_site = halide.clone().halide_site()?;
    let alkyl_carbon = halide_site.carbon;
    let halogen = halide_site.halogen;

    let halogen_element = halide.structure.atoms[halogen].element.as_str();
    let acid_co_product = match halogen_element {
        "Cl" => "destroy:hydrochloric_acid",
        "Br" => "destroy:hydrobromic_acid",
        "I" => "destroy:hydroiodic_acid",
        _ => "destroy:hydrochloric_acid",
    };

    let mut halide_editor = MolecularEditor::new(halide.structure);
    let halide_mapping = halide_editor.remove_atoms(&[halogen])?;
    let alkyl_carbon_mapped = mapped_atom(&halide_mapping, alkyl_carbon, "halide alkyl carbon")?;
    let halide_fragment = halide_editor.finish()?;

    let mut variants = Vec::new();
    for &carbon in ring_atoms {
        if aromatic.structure.atoms[carbon].element != "C" {
            continue;
        }
        let Some(hydrogen) = first_bonded_hydrogen(aromatic.structure, carbon) else {
            continue;
        };

        let mut aromatic_editor = MolecularEditor::new(aromatic.structure);
        let aromatic_mapping = aromatic_editor.remove_atoms(&[hydrogen])?;
        let aromatic_carbon_mapped = mapped_atom(&aromatic_mapping, carbon, "aromatic substitution carbon")?;
        let aromatic_fragment = aromatic_editor.finish()?;

        let product_structure = MolecularEditor::join_structures(
            &aromatic_fragment,
            aromatic_carbon_mapped,
            &halide_fragment,
            alkyl_carbon_mapped,
            1.0,
        )?;
        let product = resolver.resolve(product_structure)?;
        let activation_delta = compute_eas_activation_delta(aromatic.structure, ring_atoms, carbon);
        variants.push((product, activation_delta));
    }

    if variants.is_empty() {
        return Ok(None);
    }

    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "fc_alkylation",
        &aromatic,
        &halide,
    ))
    .reactant(aromatic.substance.id.clone(), 1, 1)
    .reactant(halide.substance.id.clone(), 1, 1)
    .catalyst_order("destroy:aluminum_trichloride", 1);

    let base_ea = 35.0;
    if variants.len() == 1 {
        builder = builder
            .product(variants[0].0.clone(), 1)
            .product(acid_co_product, 1)
            .activation_energy_kj_per_mol((base_ea + variants[0].1).max(5.0));
    } else {
        for (index, (product, activation_delta)) in variants.into_iter().enumerate() {
            builder = builder.channel(ReactionChannel::new(
                format!("fc_alkylation:position_{index}"),
                [
                    StoichiometricTerm::new(product, 1),
                    StoichiometricTerm::new(acid_co_product, 1),
                ],
                (base_ea + activation_delta).max(5.0),
            ));
        }
    }
    Ok(Some(builder.build()))
}

pub(crate) fn generate_fc_acylation(
    aromatic: SiteParticipant<'_>,
    acyl_chloride: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let ring_atoms = &aromatic.site.atoms;
    if is_ring_deactivated_for_fc(aromatic.structure, ring_atoms) {
        return Ok(None);
    }

    let acyl_site = acyl_chloride.clone().acyl_chloride_site()?;
    let acyl_carbon = acyl_site.carbon;
    let chlorine = acyl_site.chlorine;

    let mut acyl_editor = MolecularEditor::new(acyl_chloride.structure);
    let acyl_mapping = acyl_editor.remove_atoms(&[chlorine])?;
    let acyl_carbon_mapped = mapped_atom(&acyl_mapping, acyl_carbon, "acyl chloride carbon")?;
    let acyl_fragment = acyl_editor.finish()?;

    let mut variants = Vec::new();
    for &carbon in ring_atoms {
        if aromatic.structure.atoms[carbon].element != "C" {
            continue;
        }
        let Some(hydrogen) = first_bonded_hydrogen(aromatic.structure, carbon) else {
            continue;
        };

        let mut aromatic_editor = MolecularEditor::new(aromatic.structure);
        let aromatic_mapping = aromatic_editor.remove_atoms(&[hydrogen])?;
        let aromatic_carbon_mapped = mapped_atom(&aromatic_mapping, carbon, "aromatic substitution carbon")?;
        let aromatic_fragment = aromatic_editor.finish()?;

        let product_structure = MolecularEditor::join_structures(
            &aromatic_fragment,
            aromatic_carbon_mapped,
            &acyl_fragment,
            acyl_carbon_mapped,
            1.0,
        )?;
        let product = resolver.resolve(product_structure)?;
        let activation_delta = compute_eas_activation_delta(aromatic.structure, ring_atoms, carbon);
        variants.push((product, activation_delta));
    }

    if variants.is_empty() {
        return Ok(None);
    }

    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "fc_acylation",
        &aromatic,
        &acyl_chloride,
    ))
    .reactant(aromatic.substance.id.clone(), 1, 1)
    .reactant(acyl_chloride.substance.id.clone(), 1, 1)
    .catalyst_order("destroy:aluminum_trichloride", 1);

    let base_ea = 30.0;
    if variants.len() == 1 {
        builder = builder
            .product(variants[0].0.clone(), 1)
            .product("destroy:hydrochloric_acid", 1)
            .activation_energy_kj_per_mol((base_ea + variants[0].1).max(5.0));
    } else {
        for (index, (product, activation_delta)) in variants.into_iter().enumerate() {
            builder = builder.channel(ReactionChannel::new(
                format!("fc_acylation:position_{index}"),
                [
                    StoichiometricTerm::new(product, 1),
                    StoichiometricTerm::new("destroy:hydrochloric_acid", 1),
                ],
                (base_ea + activation_delta).max(5.0),
            ));
        }
    }
    Ok(Some(builder.build()))
}
