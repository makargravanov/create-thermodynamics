use super::super::resolver::DerivedSubstanceResolver;
use super::common::{first_bonded_hydrogen, sanitize_id};
use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::reactive_site::ReactiveSiteKind;
use crate::chemistry::selectivity::types::{
    ElectronicEnvironment, ReactionType, SelectivityProfile, SiteDescriptor, SubstitutionDegree,
};
use crate::chemistry::substance::{Substance, SubstanceTagId};

const HYDROGEN_ID: &str = "destroy:hydrogen";
/// Thermal cracking is endothermic and only runs hot — gate it at a typical
/// pyrolysis onset (~700 C). Below this the Arrhenius rate is negligible anyway,
/// but the hard gate keeps the reaction out of room-temperature mixtures.
/// Breaking a C-C sigma bond homolytically is steep; sized so cracking is slow
/// even once the temperature gate opens, matching its industrial severity.
const CRACKING_ACTIVATION_ENERGY_KJ_PER_MOL: f64 = 260.0;
const CRACKING_PRE_EXPONENTIAL_FACTOR: f64 = 1.0e12;
const PYROLYSIS_ACTIVATION_ENERGY_KJ_PER_MOL: f64 = 220.0;
const PYROLYSIS_PRE_EXPONENTIAL_FACTOR: f64 = 5.0e11;
/// An alkane needs at least three carbons to split into a smaller alkane plus an
/// alkene (ethane would have to make a carbene, which this mechanism cannot).
const MIN_CRACKABLE_CARBONS: usize = 3;

/// Generates thermal cracking reactions for a saturated hydrocarbon: every
/// acyclic C-C single bond, in both orientations, undergoes beta-scission into a
/// smaller alkane plus a terminal alkene. The bond breaks, one carbon is capped
/// with H (the alkane fragment), and the other forms a C=C with a neighbor that
/// sheds an H (the alkene fragment) — so no atoms are gained or lost overall.
/// Returns an empty vec for anything that is not a concrete neutral alkane.
pub(crate) fn generate_cracking(
    feedstock: &Substance,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    if feedstock.charge != 0 || is_hypothetical(feedstock) {
        return Ok(Vec::new());
    }
    let Some(structure) = feedstock.molecular_structure.as_ref() else {
        return Ok(Vec::new());
    };
    if !is_concrete_alkane(structure) {
        return Ok(Vec::new());
    }
    let mut reactions = Vec::new();
    for bond in &structure.bonds {
        if structure.atoms[bond.from].element != "C" || structure.atoms[bond.to].element != "C" {
            continue;
        }
        for (alkane_carbon, alkene_carbon) in [(bond.from, bond.to), (bond.to, bond.from)] {
            if let Some(reaction) =
                crack_bond(feedstock, structure, alkane_carbon, alkene_carbon, resolver)?
            {
                reactions.push(reaction);
            }
        }
    }
    Ok(reactions)
}

/// General thermal dehydrogenation for concrete hydrocarbons. Any C-C single or
/// double bond whose two carbon atoms both carry explicit hydrogen can lose H2
/// and increase the bond order by one:
///
/// - `C-C -> C=C + H2`
/// - `C=C -> C#C + H2`
///
/// This uses graph structure, not substance names, so repeated dynamic
/// generation can discover alkenes, dienes, alkynes and cyclic products from the
/// same small rule.
pub(crate) fn generate_pyrolysis(
    feedstock: &Substance,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Vec<Reaction>> {
    if feedstock.charge != 0 || is_hypothetical(feedstock) {
        return Ok(Vec::new());
    }
    let Some(structure) = feedstock.molecular_structure.as_ref() else {
        return Ok(Vec::new());
    };
    if !is_concrete_hydrocarbon(structure) {
        return Ok(Vec::new());
    }

    let mut best_by_product = std::collections::BTreeMap::new();
    for bond in &structure.bonds {
        if structure.atoms[bond.from].element != "C" || structure.atoms[bond.to].element != "C" {
            continue;
        }
        if bond_order_matches(bond.order, 1.5) || bond.order >= 3.0 {
            continue;
        }
        let next_order = if bond_order_matches(bond.order, 1.0) {
            2.0
        } else if bond_order_matches(bond.order, 2.0) {
            3.0
        } else {
            continue;
        };
        let Some(reaction) =
            dehydrogenate_bond(feedstock, structure, bond.from, bond.to, next_order, resolver)?
        else {
            continue;
        };
        let product_id = reaction.products[0].substance_id.clone();
        best_by_product.entry(product_id).or_insert(reaction);
    }
    Ok(best_by_product.into_values().collect())
}

/// Beta-scission of the `alkane_carbon`-`alkene_carbon` bond: the alkane fragment
/// caps `alkane_carbon` with H, the alkene fragment forms a C=C between
/// `alkene_carbon` and a carbon neighbor that sheds one H. Returns `None` when the
/// bond is in a ring (cracking does not open rings) or when the alkene side has no
/// suitable neighbor to receive the double bond (e.g. a terminal methyl).
fn crack_bond(
    feedstock: &Substance,
    structure: &MolecularStructure,
    alkane_carbon: usize,
    alkene_carbon: usize,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if bond_is_in_ring(structure, alkane_carbon, alkene_carbon) {
        return Ok(None);
    }
    // The new double bond forms between the alkene carbon and one of its OTHER
    // carbon neighbors, which must give up a hydrogen. No such neighbor (terminal
    // carbon) means this orientation cannot make an alkene — skip it.
    let Some(double_bond_partner) = structure
        .neighbors(alkene_carbon)
        .into_iter()
        .find(|(neighbor, order)| {
            *neighbor != alkane_carbon
                && structure.atoms[*neighbor].element == "C"
                && bond_order_matches(*order, 1.0)
                && first_bonded_hydrogen(structure, *neighbor).is_some()
        })
        .map(|(neighbor, _)| neighbor)
    else {
        return Ok(None);
    };

    let (alkane_fragment, alkane_map, alkene_fragment, alkene_map) =
        MolecularEditor::split_at_bond(structure, alkane_carbon, alkene_carbon)?;

    // Alkane side: cap the now-dangling carbon with a hydrogen.
    let mut alkane_editor = MolecularEditor::new(&alkane_fragment);
    alkane_editor.add_atom(mapped(&alkane_map, alkane_carbon)?, "H", 0.0, 1.0)?;
    let alkane_id = resolver.resolve(alkane_editor.finish()?)?;

    // Alkene side: promote the carbon's bond to its partner to a double bond and
    // remove one hydrogen from that partner to keep the valence (and mass) right.
    let mut alkene_editor = MolecularEditor::new(&alkene_fragment);
    let partner_new = mapped(&alkene_map, double_bond_partner)?;
    alkene_editor.set_bond_order(mapped(&alkene_map, alkene_carbon)?, partner_new, 2.0)?;
    let partner_hydrogen = first_bonded_hydrogen(&alkene_editor.structure(), partner_new)
        .ok_or_else(|| ChemistryError::GenerationInvariantViolation {
            generator: "cracking".to_string(),
            substance_id: feedstock.id.to_string(),
            reason: "alkene partner lost its explicit hydrogen after graph split".to_string(),
        })?;
    alkene_editor.remove_atoms(&[partner_hydrogen])?;
    let alkene_id = resolver.resolve(alkene_editor.finish()?)?;

    Ok(Some(
        Reaction::builder(format!(
            "cracking/{}/{alkane_carbon}_{alkene_carbon}",
            sanitize_id(feedstock.id.as_str())
        ))
        .reactant(feedstock.id.clone(), 1, 1)
        .product(alkane_id, 1)
        .product(alkene_id, 1)
        .reactant_phase_access(
            feedstock.id.clone(),
            [
                MixturePhase::Organic,
                MixturePhase::Gas,
                MixturePhase::SupercriticalFluid,
            ],
        )
        .pre_exponential_factor(CRACKING_PRE_EXPONENTIAL_FACTOR)
        .activation_energy_kj_per_mol(CRACKING_ACTIVATION_ENERGY_KJ_PER_MOL)
        .selectivity_profile(
            SelectivityProfile::new(
                ReactionType::HydrocarbonCracking,
                cracking_site_descriptor(structure, alkane_carbon, alkene_carbon),
            )
            .never_suppress(),
        )
        .build(),
    ))
}

fn dehydrogenate_bond(
    feedstock: &Substance,
    structure: &MolecularStructure,
    first: usize,
    second: usize,
    next_order: f64,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let Some(first_hydrogen) = first_bonded_hydrogen(structure, first) else {
        return Ok(None);
    };
    let Some(second_hydrogen) = first_bonded_hydrogen(structure, second) else {
        return Ok(None);
    };

    let descriptor = cracking_site_descriptor(structure, first, second);
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[first_hydrogen, second_hydrogen])?;
    let first = mapped_for_generator("pyrolysis", &mapping, first)?;
    let second = mapped_for_generator("pyrolysis", &mapping, second)?;
    editor.set_bond_order(first, second, next_order)?;
    let product_id = resolver.resolve(editor.finish()?)?;
    if product_id == feedstock.id {
        return Ok(None);
    }

    Ok(Some(
        Reaction::builder(format!(
            "pyrolysis/{}/{first}_{second}_{next_order:.0}",
            sanitize_id(feedstock.id.as_str())
        ))
        .reactant(feedstock.id.clone(), 1, 1)
        .product(product_id.clone(), 1)
        .product(HYDROGEN_ID, 1)
        .reactant_phase_access(
            feedstock.id.clone(),
            [
                MixturePhase::Organic,
                MixturePhase::Gas,
                MixturePhase::SupercriticalFluid,
            ],
        )
        .product_phase(product_id, MixturePhase::Organic)
        .product_phase(HYDROGEN_ID, MixturePhase::Gas)
        .pre_exponential_factor(PYROLYSIS_PRE_EXPONENTIAL_FACTOR)
        .activation_energy_kj_per_mol(PYROLYSIS_ACTIVATION_ENERGY_KJ_PER_MOL)
        .selectivity_profile(
            SelectivityProfile::new(ReactionType::HydrocarbonPyrolysis, descriptor)
                .never_suppress(),
        )
        .build(),
    ))
}

/// A concrete alkane: neutral, no R-group placeholders, only C and H, every bond
/// a single bond, and at least [`MIN_CRACKABLE_CARBONS`] carbons.
fn is_concrete_alkane(structure: &MolecularStructure) -> bool {
    let mut carbons = 0;
    for atom in &structure.atoms {
        match atom.element.as_str() {
            "C" => carbons += 1,
            "H" => {}
            _ => return false,
        }
    }
    if carbons < MIN_CRACKABLE_CARBONS {
        return false;
    }
    structure
        .bonds
        .iter()
        .all(|bond| bond_order_matches(bond.order, 1.0))
}

fn is_concrete_hydrocarbon(structure: &MolecularStructure) -> bool {
    structure
        .atoms
        .iter()
        .all(|atom| matches!(atom.element.as_str(), "C" | "H"))
}

/// Whether the `first`-`second` bond lies on a ring: true when `second` is still
/// reachable from `first` after that bond is ignored. Cracking only splits chain
/// bonds, so ring bonds are skipped.
fn bond_is_in_ring(structure: &MolecularStructure, first: usize, second: usize) -> bool {
    let mut seen = vec![false; structure.atoms.len()];
    let mut stack = vec![first];
    seen[first] = true;
    while let Some(atom) = stack.pop() {
        for (neighbor, _) in structure.neighbors(atom) {
            let skip_split_bond =
                (atom == first && neighbor == second) || (atom == second && neighbor == first);
            if skip_split_bond || seen[neighbor] {
                continue;
            }
            if neighbor == second {
                return true;
            }
            seen[neighbor] = true;
            stack.push(neighbor);
        }
    }
    false
}

fn cracking_site_descriptor(
    structure: &MolecularStructure,
    first: usize,
    second: usize,
) -> SiteDescriptor {
    let first_degree = structure.carbon_degree(first);
    let second_degree = structure.carbon_degree(second);
    let max_degree = first_degree.max(second_degree);
    let degree = match max_degree {
        0 | 1 => SubstitutionDegree::Primary,
        2 => SubstitutionDegree::Secondary,
        _ => SubstitutionDegree::Tertiary,
    };
    let bulky = bulky_substituent_count(structure, first)
        .saturating_add(bulky_substituent_count(structure, second));
    SiteDescriptor::new(
        ReactiveSiteKind::AlkylHydrogen,
        degree,
        ElectronicEnvironment {
            electron_donating_groups: first_degree.saturating_add(second_degree).saturating_sub(2)
                as u32,
            electron_withdrawing_groups: heteroatom_neighbor_count(structure, first)
                + heteroatom_neighbor_count(structure, second),
            resonance_stabilization: is_allylic(structure, first) || is_allylic(structure, second),
            aromatic: false,
        },
        bulky,
        false,
    )
}

fn is_allylic(structure: &MolecularStructure, atom: usize) -> bool {
    structure
        .neighbors(atom)
        .into_iter()
        .any(|(neighbor, order)| {
            bond_order_matches(order, 1.0)
                && structure
                    .neighbors(neighbor)
                    .into_iter()
                    .any(|(other, other_order)| {
                        other != atom
                            && structure.atoms[other].element == "C"
                            && bond_order_matches(other_order, 2.0)
                    })
        })
}

fn heteroatom_neighbor_count(structure: &MolecularStructure, atom: usize) -> u32 {
    structure
        .neighbors(atom)
        .into_iter()
        .filter(|(neighbor, _)| !matches!(structure.atoms[*neighbor].element.as_str(), "C" | "H"))
        .count() as u32
}

fn bulky_substituent_count(structure: &MolecularStructure, atom: usize) -> u32 {
    structure
        .neighbors(atom)
        .into_iter()
        .filter(|(neighbor, _)| {
            structure.atoms[*neighbor].element != "H"
                && structure
                    .neighbors(*neighbor)
                    .into_iter()
                    .filter(|(second, _)| {
                        *second != atom && structure.atoms[*second].element != "H"
                    })
                    .count()
                    >= 2
        })
        .count() as u32
}

fn mapped(mapping: &[Option<usize>], old_index: usize) -> ChemistryResult<usize> {
    mapped_for_generator("cracking", mapping, old_index)
}

fn mapped_for_generator(
    generator: &str,
    mapping: &[Option<usize>],
    old_index: usize,
) -> ChemistryResult<usize> {
    mapping[old_index].ok_or_else(|| {
        crate::chemistry::error::ChemistryError::GenerationInvariantViolation {
            generator: generator.to_string(),
            substance_id: String::new(),
            reason: format!("atom {old_index} missing from generated fragment"),
        }
    })
}

fn is_hypothetical(substance: &Substance) -> bool {
    substance
        .tags
        .iter()
        .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;
    use crate::chemistry::substance::SubstanceId;

    fn alkane(id: &str, code: &str) -> Substance {
        let structure = parse_frowns(code).unwrap();
        let summary = structure.summary().unwrap();
        Substance::new(
            id,
            summary.charge,
            summary.molar_mass_grams,
            700.0,
            250.0,
            80.0,
            20_000.0,
        )
        .with_molecular_structure(structure)
    }

    fn run(feedstock: &Substance) -> Vec<Reaction> {
        let mut resolver =
            DerivedSubstanceResolver::new_from_canonical_to_id(std::collections::BTreeMap::new());
        generate_cracking(feedstock, &mut resolver).unwrap()
    }

    fn run_pyrolysis(feedstock: &Substance) -> (Vec<Reaction>, Vec<Substance>) {
        let mut resolver =
            DerivedSubstanceResolver::new_from_canonical_to_id(std::collections::BTreeMap::new());
        let reactions = generate_pyrolysis(feedstock, &mut resolver).unwrap();
        (reactions, resolver.substances)
    }

    #[test]
    fn propane_cracks_into_methane_and_ethene() {
        let reactions = run(&alkane("test:propane", "CCC"));
        // Propane's two equivalent C-C bonds each crack to methane + ethene; the
        // edge always has exactly one feedstock reactant and two fragment products.
        assert!(!reactions.is_empty());
        for reaction in &reactions {
            assert_eq!(reaction.reactants.len(), 1);
            assert_eq!(reaction.products.len(), 2);
            assert!(reaction.conditions.is_empty());
            assert!(matches!(
                reaction
                    .selectivity_profile
                    .as_ref()
                    .map(|profile| profile.mechanism),
                Some(ReactionType::HydrocarbonCracking)
            ));
            assert!(reaction
                .phase_access
                .get(&SubstanceId::from("test:propane"))
                .is_some_and(|access| access.phases.contains(&MixturePhase::SupercriticalFluid)));
        }
    }

    #[test]
    fn ethane_is_too_short_to_crack() {
        assert!(run(&alkane("test:ethane", "CC")).is_empty());
    }

    #[test]
    fn cracking_conserves_mass() {
        let feedstock = alkane("test:pentane", "CCCCC");
        let mut resolver =
            DerivedSubstanceResolver::new_from_canonical_to_id(std::collections::BTreeMap::new());
        let reactions = generate_cracking(&feedstock, &mut resolver).unwrap();
        assert!(!reactions.is_empty());
        let mass_of = |id: &SubstanceId| -> f64 {
            resolver
                .substances
                .iter()
                .find(|s| &s.id == id)
                .map(|s| s.molar_mass_grams)
                .unwrap_or_else(|| feedstock.molar_mass_grams)
        };
        for reaction in &reactions {
            let products: f64 = reaction
                .products
                .iter()
                .map(|p| mass_of(&p.substance_id))
                .sum();
            assert!(
                (products - feedstock.molar_mass_grams).abs() < 1.0e-6,
                "cracking must conserve mass"
            );
        }
    }

    #[test]
    fn unsaturated_or_heteroatom_feedstock_does_not_crack() {
        assert!(run(&alkane("test:propene", "C=CC")).is_empty());
        assert!(run(&alkane("test:ethanol", "CCO")).is_empty());
    }

    #[test]
    fn ethane_pyrolysis_dehydrogenates_to_ethene_and_hydrogen() {
        let (reactions, products) = run_pyrolysis(&alkane("test:ethane", "CC"));
        assert_eq!(reactions.len(), 1);
        let reaction = &reactions[0];
        assert_eq!(
            reaction.reactants[0].substance_id,
            SubstanceId::from("test:ethane")
        );
        assert!(reaction
            .products
            .iter()
            .any(|term| term.substance_id == SubstanceId::from(HYDROGEN_ID)));
        assert!(matches!(
            reaction
                .selectivity_profile
                .as_ref()
                .map(|profile| profile.mechanism),
            Some(ReactionType::HydrocarbonPyrolysis)
        ));
        let product_id = reaction
            .products
            .iter()
            .find(|term| term.substance_id != SubstanceId::from(HYDROGEN_ID))
            .unwrap()
            .substance_id
            .clone();
        let product = products
            .iter()
            .find(|substance| substance.id == product_id)
            .expect("pyrolysis must create an organic product");
        assert!(product
            .molecular_structure
            .as_ref()
            .unwrap()
            .bonds
            .iter()
            .any(|bond| bond_order_matches(bond.order, 2.0)));
    }

    #[test]
    fn propane_pyrolysis_deduplicates_equivalent_dehydrogenation_sites() {
        let (reactions, _) = run_pyrolysis(&alkane("test:propane", "CCC"));
        assert_eq!(reactions.len(), 1);
    }

    #[test]
    fn alkene_pyrolysis_can_make_an_alkyne() {
        let (reactions, products) = run_pyrolysis(&alkane("test:ethene", "C=C"));
        assert_eq!(reactions.len(), 1);
        let product_id = reactions[0]
            .products
            .iter()
            .find(|term| term.substance_id != SubstanceId::from(HYDROGEN_ID))
            .unwrap()
            .substance_id
            .clone();
        let product = products
            .iter()
            .find(|substance| substance.id == product_id)
            .unwrap();
        assert!(product
            .molecular_structure
            .as_ref()
            .unwrap()
            .bonds
            .iter()
            .any(|bond| bond_order_matches(bond.order, 3.0)));
    }

    #[test]
    fn heteroatom_feedstock_does_not_get_hydrocarbon_pyrolysis() {
        assert!(run_pyrolysis(&alkane("test:ethanol", "CCO")).0.is_empty());
    }
}
