use std::collections::BTreeMap;

use super::common::sanitize_id;
use super::super::resolver::DerivedSubstanceResolver;
use crate::chemistry::error::ChemistryResult;
use crate::chemistry::selectivity::types::{
    ElectronicEnvironment, ReactionType, SelectivityProfile, SiteDescriptor, SubstitutionDegree,
};
use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::molecule::{bond_order_matches, MolecularEditor, MolecularStructure};
use crate::chemistry::reaction::Reaction;
use crate::chemistry::reactive_site::ReactiveSiteKind;
use crate::chemistry::substance::{Substance, SubstanceId, SubstanceTagId};

const CHLORINE_ID: &str = "destroy:chlorine";
const BROMINE_ID: &str = "destroy:bromine";
const HYDROCHLORIC_ACID_ID: &str = "destroy:hydrochloric_acid";
const HYDROBROMIC_ACID_ID: &str = "destroy:hydrobromic_acid";

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum RadicalHalogen {
    Chlorine,
    Bromine,
}

#[derive(Debug, Clone)]
struct RadicalSite {
    hydrogen: usize,
    descriptor: SiteDescriptor,
    stability: f64,
}

#[derive(Debug, Clone)]
struct Candidate {
    product_id: SubstanceId,
    descriptor: SiteDescriptor,
    activation_energy_kj_per_mol: f64,
    pre_exponential_factor: f64,
}

pub(crate) fn generate_radical_halogenations(
    fuel: &Substance,
    resolver: &mut DerivedSubstanceResolver,
    available_halogens: &[RadicalHalogen],
) -> ChemistryResult<Vec<Reaction>> {
    if available_halogens.is_empty() {
        return Ok(Vec::new());
    }
    if fuel.charge != 0 || is_hypothetical(fuel) {
        return Ok(Vec::new());
    }
    let Some(structure) = fuel.molecular_structure.as_ref() else {
        return Ok(Vec::new());
    };
    if structure.atoms.iter().any(|atom| atom.element == "R") {
        return Ok(Vec::new());
    }
    let sites = radical_abstraction_sites(structure);
    if sites.is_empty() {
        return Ok(Vec::new());
    }

    let mut reactions = Vec::new();
    for halogen in available_halogens {
        let mut best_by_product = BTreeMap::<SubstanceId, Candidate>::new();
        for site in &sites {
            let mut editor = MolecularEditor::new(structure);
            editor.replace_atom(site.hydrogen, halogen.symbol(), 0.0)?;
            let product_id = resolver.resolve(editor.finish()?)?;
            let candidate = Candidate {
                product_id: product_id.clone(),
                descriptor: site.descriptor.clone(),
                activation_energy_kj_per_mol: halogen.activation_energy(site.stability),
                pre_exponential_factor: halogen.pre_exponential_factor(site.stability),
            };
            best_by_product
                .entry(product_id)
                .and_modify(|current| {
                    if candidate.activation_energy_kj_per_mol
                        < current.activation_energy_kj_per_mol
                    {
                        *current = candidate.clone();
                    }
                })
                .or_insert(candidate);
        }

        for candidate in best_by_product.into_values() {
            reactions.push(radical_halogenation_reaction(fuel, *halogen, candidate));
        }
    }
    Ok(reactions)
}

fn radical_halogenation_reaction(
    fuel: &Substance,
    halogen: RadicalHalogen,
    candidate: Candidate,
) -> Reaction {
    Reaction::builder(format!(
        "radical_halogenation/{}/{}/{}",
        sanitize_id(fuel.id.as_str()),
        halogen.id_token(),
        sanitize_id(candidate.product_id.as_str())
    ))
    .reactant(fuel.id.clone(), 1, 1)
    .reactant(halogen.reagent_id(), 1, 1)
    .product(candidate.product_id.clone(), 1)
    .product(halogen.acid_product_id(), 1)
    .reactant_phase_access(fuel.id.clone(), [MixturePhase::Gas, MixturePhase::Organic])
    .reactant_phase_access(
        halogen.reagent_id(),
        [MixturePhase::Gas, MixturePhase::Organic],
    )
    .product_phase(candidate.product_id, MixturePhase::Organic)
    .product_phase(halogen.acid_product_id(), MixturePhase::Gas)
    .pre_exponential_factor(candidate.pre_exponential_factor)
    .activation_energy_kj_per_mol(candidate.activation_energy_kj_per_mol)
    .selectivity_profile(
        SelectivityProfile::new(ReactionType::RadicalHalogenation, candidate.descriptor)
            .never_suppress(),
    )
    .build()
}

fn radical_abstraction_sites(structure: &MolecularStructure) -> Vec<RadicalSite> {
    let mut sites = Vec::new();
    for carbon in 0..structure.atoms.len() {
        if structure.atoms[carbon].element != "C" || !is_radical_substitution_carbon(structure, carbon)
        {
            continue;
        }
        let Some(hydrogen) = structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (bond_order_matches(order, 1.0) && structure.atoms[neighbor].element == "H")
                    .then_some(neighbor)
            })
        else {
            continue;
        };
        let stability = radical_stability_score(structure, carbon);
        sites.push(RadicalSite {
            hydrogen,
            descriptor: radical_site_descriptor(structure, carbon, stability),
            stability,
        });
    }
    sites
}

fn is_radical_substitution_carbon(structure: &MolecularStructure, carbon: usize) -> bool {
    if is_aromatic_carbon(structure, carbon) {
        return false;
    }
    structure.neighbors(carbon).into_iter().all(|(neighbor, order)| {
        structure.atoms[neighbor].element == "H" || bond_order_matches(order, 1.0)
    })
}

fn radical_site_descriptor(
    structure: &MolecularStructure,
    carbon: usize,
    stability: f64,
) -> SiteDescriptor {
    let degree = if is_benzylic(structure, carbon) {
        SubstitutionDegree::Benzylic
    } else if is_allylic(structure, carbon) {
        SubstitutionDegree::Allylic
    } else {
        match structure.carbon_degree(carbon) {
            0 | 1 => SubstitutionDegree::Primary,
            2 => SubstitutionDegree::Secondary,
            _ => SubstitutionDegree::Tertiary,
        }
    };
    let electronics = ElectronicEnvironment {
        electron_donating_groups: structure.carbon_degree(carbon).saturating_sub(1) as u32,
        electron_withdrawing_groups: heteroatom_neighbor_count(structure, carbon),
        resonance_stabilization: stability >= 2.0,
        aromatic: false,
    };
    SiteDescriptor::new(
        ReactiveSiteKind::AlkylHydrogen,
        degree,
        electronics,
        bulky_substituent_count(structure, carbon),
        false,
    )
}

fn radical_stability_score(structure: &MolecularStructure, carbon: usize) -> f64 {
    let mut score = structure.carbon_degree(carbon) as f64;
    if is_benzylic(structure, carbon) {
        score += 2.5;
    }
    if is_allylic(structure, carbon) {
        score += 1.8;
    }
    score -= heteroatom_neighbor_count(structure, carbon) as f64 * 0.4;
    score.max(0.0)
}

fn heteroatom_neighbor_count(structure: &MolecularStructure, atom: usize) -> u32 {
    structure
        .neighbors(atom)
        .into_iter()
        .filter(|(neighbor, _)| {
            matches!(
                structure.atoms[*neighbor].element.as_str(),
                "O" | "N" | "S" | "F" | "Cl" | "Br" | "I"
            )
        })
        .count() as u32
}

fn bulky_substituent_count(structure: &MolecularStructure, carbon: usize) -> u32 {
    structure
        .neighbors(carbon)
        .into_iter()
        .filter(|(neighbor, _)| {
            structure.atoms[*neighbor].element != "H"
                && structure
                    .neighbors(*neighbor)
                    .into_iter()
                    .filter(|(second, _)| *second != carbon && structure.atoms[*second].element != "H")
                    .count()
                    >= 2
        })
        .count() as u32
}

fn is_benzylic(structure: &MolecularStructure, atom: usize) -> bool {
    structure.neighbors(atom).into_iter().any(|(neighbor, order)| {
        bond_order_matches(order, 1.0) && is_aromatic_carbon(structure, neighbor)
    })
}

fn is_allylic(structure: &MolecularStructure, atom: usize) -> bool {
    structure.neighbors(atom).into_iter().any(|(neighbor, order)| {
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

fn is_aromatic_carbon(structure: &MolecularStructure, atom: usize) -> bool {
    structure.atoms.get(atom).is_some_and(|atom| atom.element == "C")
        && structure
            .neighbors(atom)
            .iter()
            .filter(|(_, order)| bond_order_matches(*order, 1.5))
            .count()
            >= 2
}

fn is_hypothetical(substance: &Substance) -> bool {
    substance
        .tags
        .iter()
        .any(|tag| tag == &SubstanceTagId::from("destroy:hypothetical"))
}

impl RadicalHalogen {
    fn symbol(self) -> &'static str {
        match self {
            Self::Chlorine => "Cl",
            Self::Bromine => "Br",
        }
    }

    fn id_token(self) -> &'static str {
        match self {
            Self::Chlorine => "chlorination",
            Self::Bromine => "bromination",
        }
    }

    fn reagent_id(self) -> &'static str {
        match self {
            Self::Chlorine => CHLORINE_ID,
            Self::Bromine => BROMINE_ID,
        }
    }

    fn acid_product_id(self) -> &'static str {
        match self {
            Self::Chlorine => HYDROCHLORIC_ACID_ID,
            Self::Bromine => HYDROBROMIC_ACID_ID,
        }
    }

    fn activation_energy(self, stability: f64) -> f64 {
        match self {
            Self::Chlorine => (48.0 - stability * 4.0).max(26.0),
            Self::Bromine => (58.0 - stability * 7.0).max(24.0),
        }
    }

    fn pre_exponential_factor(self, stability: f64) -> f64 {
        match self {
            Self::Chlorine => 8.0e7 * (1.0 + stability * 0.15),
            Self::Bromine => 3.0e7 * (1.0 + stability * 0.35),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;

    fn substance(id: &str, code: &str) -> Substance {
        let structure = parse_frowns(code).unwrap();
        let summary = structure.summary().unwrap();
        Substance::new(
            id,
            summary.charge,
            summary.molar_mass_grams,
            700.0,
            350.0,
            100.0,
            20_000.0,
        )
        .with_molecular_structure(structure)
    }

    #[test]
    fn radical_halogenation_replaces_explicit_hydrogen() {
        let methane = substance("test:methane", "C");
        let mut resolver = DerivedSubstanceResolver::new_from_canonical_to_id(BTreeMap::new());
        let reactions = generate_radical_halogenations(
            &methane,
            &mut resolver,
            &[RadicalHalogen::Chlorine, RadicalHalogen::Bromine],
        )
        .unwrap();

        assert_eq!(reactions.len(), 2);
        assert!(reactions.iter().any(|reaction| {
            reaction.reactants.iter().any(|term| term.substance_id.as_str() == CHLORINE_ID)
                && reaction
                    .products
                    .iter()
                    .any(|term| term.substance_id.as_str() == HYDROCHLORIC_ACID_ID)
        }));
        assert!(resolver
            .substances
            .iter()
            .any(|substance| substance.molecular_structure.as_ref().is_some_and(|structure| {
                structure.atoms.iter().any(|atom| atom.element == "Cl")
            })));
    }

    #[test]
    fn radical_halogenation_ignores_aromatic_ring_hydrogen() {
        let benzene = substance("test:benzene", "destroy:benzene:,,,,,");
        let mut resolver = DerivedSubstanceResolver::new_from_canonical_to_id(BTreeMap::new());
        let reactions = generate_radical_halogenations(
            &benzene,
            &mut resolver,
            &[RadicalHalogen::Chlorine, RadicalHalogen::Bromine],
        )
        .unwrap();

        assert!(reactions.is_empty());
    }
}
