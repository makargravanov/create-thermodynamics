use std::collections::{BTreeMap, BTreeSet};

use super::error::{ChemistryError, ChemistryResult};
use super::functional_group::{FunctionalGroup, FunctionalGroupType};
use super::molecule::{MolecularEditor, MolecularStructure};
use super::reaction::Reaction;
use super::registry::{ChemistryRegistry, ChemistryRegistryBuilder};
use super::substance::{Substance, SubstanceId};

const DEFAULT_DERIVED_DENSITY: f64 = 1000.0;
const DEFAULT_DERIVED_HEAT_CAPACITY: f64 = 100.0;
const DEFAULT_DERIVED_LATENT_HEAT: f64 = 20_000.0;

pub fn destroy_registry_with_generated_reactions_builder(
) -> ChemistryResult<ChemistryRegistryBuilder> {
    let base_registry = super::destroy_registry_builder()?.build()?;
    let generated = generate_organic_reactions(&base_registry)?;
    let mut builder = ChemistryRegistryBuilder::from_registry(&base_registry);
    for substance in generated.substances {
        builder = builder.substance(substance);
    }
    for reaction in generated.reactions {
        builder = builder.reaction(reaction);
    }
    Ok(builder)
}

#[derive(Debug, Default)]
struct GeneratedOrganicCatalog {
    substances: Vec<Substance>,
    reactions: Vec<Reaction>,
}

struct DerivedSubstanceResolver {
    canonical_to_id: BTreeMap<String, SubstanceId>,
    generated_id_to_canonical: BTreeMap<SubstanceId, String>,
    substances: Vec<Substance>,
}

impl DerivedSubstanceResolver {
    fn new(registry: &ChemistryRegistry) -> Self {
        let mut canonical_to_id = BTreeMap::new();
        for substance in registry.substances() {
            if let Some(structure) = &substance.molecular_structure {
                canonical_to_id
                    .entry(structure.canonical_code())
                    .or_insert_with(|| substance.id.clone());
            }
        }
        Self {
            canonical_to_id,
            generated_id_to_canonical: BTreeMap::new(),
            substances: Vec::new(),
        }
    }

    fn resolve(&mut self, structure: MolecularStructure) -> ChemistryResult<SubstanceId> {
        let canonical = structure.canonical_code();
        if let Some(id) = self.canonical_to_id.get(&canonical) {
            return Ok(id.clone());
        }
        let id = SubstanceId::new(format!("destroy:derived_{:016x}", stable_hash(&canonical)))?;
        if let Some(existing) = self.generated_id_to_canonical.get(&id) {
            if existing != &canonical {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: id.to_string(),
                    reason: "derived substance id collision".to_string(),
                });
            }
        }
        let summary = structure.summary()?;
        let substance = Substance::new(
            id.clone(),
            summary.charge,
            summary.molar_mass_grams,
            DEFAULT_DERIVED_DENSITY,
            if summary.charge == 0 {
                1000.0
            } else {
                f64::MAX
            },
            DEFAULT_DERIVED_HEAT_CAPACITY,
            DEFAULT_DERIVED_LATENT_HEAT,
        )
        .with_catalog_metadata(
            Some(format!("generated:{}", structure.canonical_code())),
            None,
            0x20FF_FFFF,
            Vec::new(),
        )
        .with_molecular_structure(structure);
        self.canonical_to_id.insert(canonical.clone(), id.clone());
        self.generated_id_to_canonical.insert(id.clone(), canonical);
        self.substances.push(substance);
        Ok(id)
    }
}

fn generate_organic_reactions(
    registry: &ChemistryRegistry,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let mut resolver = DerivedSubstanceResolver::new(registry);
    let mut reactions = Vec::new();
    let mut reaction_ids = BTreeSet::new();
    let reactants = registry.substances().cloned().collect::<Vec<_>>();

    for substance in &reactants {
        let Some(structure) = &substance.molecular_structure else {
            continue;
        };
        for group in &substance.functional_groups {
            match group.group_type {
                FunctionalGroupType::Halide => {
                    if let Some(reaction) = generate_halide_hydroxide_substitution(
                        substance,
                        structure,
                        group,
                        &mut resolver,
                    )? {
                        push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                    }
                }
                FunctionalGroupType::Alcohol => {
                    if let Some(reaction) =
                        generate_alcohol_oxidation(substance, structure, group, &mut resolver)?
                    {
                        push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                    }
                }
                FunctionalGroupType::Nitrile => {
                    let reaction =
                        generate_nitrile_hydrolysis(substance, structure, group, &mut resolver)?;
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                FunctionalGroupType::Nitro => {
                    let reaction =
                        generate_nitro_hydrogenation(substance, structure, group, &mut resolver)?;
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                _ => {}
            }
        }
    }

    for acid in &reactants {
        let Some(acid_structure) = &acid.molecular_structure else {
            continue;
        };
        for acid_group in acid
            .functional_groups
            .iter()
            .filter(|group| group.group_type == FunctionalGroupType::CarboxylicAcid)
        {
            for alcohol in &reactants {
                let Some(alcohol_structure) = &alcohol.molecular_structure else {
                    continue;
                };
                for alcohol_group in alcohol
                    .functional_groups
                    .iter()
                    .filter(|group| group.group_type == FunctionalGroupType::Alcohol)
                {
                    let reaction = generate_carboxylic_acid_esterification(
                        acid,
                        acid_structure,
                        acid_group,
                        alcohol,
                        alcohol_structure,
                        alcohol_group,
                        &mut resolver,
                    )?;
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
        }
    }

    Ok(GeneratedOrganicCatalog {
        substances: resolver.substances,
        reactions,
    })
}

fn push_unique_reaction(
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
    reaction: Reaction,
) -> ChemistryResult<()> {
    let id = reaction.id.to_string();
    if reaction_ids.insert(id.clone()) {
        reactions.push(reaction);
        Ok(())
    } else {
        Err(ChemistryError::DuplicateReaction(id))
    }
}

fn generate_halide_hydroxide_substitution(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let carbon = group.atoms[0];
    let halogen = group.atoms[1];
    let halide_ion = match structure.atoms[halogen].element.as_str() {
        "Cl" => "destroy:chloride",
        "F" => "destroy:fluoride",
        "I" => "destroy:iodide",
        _ => {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: generated_group_reaction_id(
                    "halide_hydroxide_substitution",
                    substance,
                    group,
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
        Reaction::builder(generated_group_reaction_id(
            "halide_hydroxide_substitution",
            substance,
            group,
        ))
        .reactant(substance.id.clone(), 1, 1)
        .reactant(
            "destroy:hydroxide",
            1,
            if group.degree.unwrap_or_default() == 3 {
                0
            } else {
                1
            },
        )
        .product(product, 1)
        .product(halide_ion, 1)
        .build(),
    ))
}

fn generate_alcohol_oxidation(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if group.degree.unwrap_or_default() >= 3 {
        return Ok(None);
    }
    let carbon = group.atoms[0];
    let oxygen = group.atoms[1];
    let Some(carbon_hydrogen) = first_bonded_hydrogen(structure, carbon) else {
        return Ok(None);
    };
    let oxygen_hydrogen = first_bonded_hydrogen(structure, oxygen).ok_or_else(|| {
        ChemistryError::InvalidReaction {
            reaction_id: generated_group_reaction_id("alcohol_oxidation", substance, group),
            reason: "alcohol oxygen has no explicit hydrogen".to_string(),
        }
    })?;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[carbon_hydrogen, oxygen_hydrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "alcohol carbon")?;
    let oxygen = mapped_atom(&mapping, oxygen, "alcohol oxygen")?;
    editor.set_bond_order(carbon, oxygen, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_group_reaction_id(
            "alcohol_oxidation",
            substance,
            group,
        ))
        .reactant(substance.id.clone(), 3, 1)
        .reactant("destroy:dichromate", 1, 1)
        .reactant("destroy:proton", 8, 1)
        .product(product, 3)
        .product("destroy:chromium_iii", 2)
        .product("destroy:water", 7)
        .activation_energy_kj_per_mol(25.0)
        .build(),
    ))
}

fn generate_carboxylic_acid_esterification(
    acid: &Substance,
    acid_structure: &MolecularStructure,
    acid_group: &FunctionalGroup,
    alcohol: &Substance,
    alcohol_structure: &MolecularStructure,
    alcohol_group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let acid_carbon = acid_group.atoms[0];
    let acid_hydroxyl_oxygen = acid_group.atoms[2];
    let acid_proton =
        first_bonded_hydrogen(acid_structure, acid_hydroxyl_oxygen).ok_or_else(|| {
            ChemistryError::InvalidReaction {
                reaction_id: generated_pair_group_reaction_id(
                    "carboxylic_acid_esterification",
                    acid,
                    acid_group,
                    alcohol,
                    alcohol_group,
                ),
                reason: "carboxylic acid has no explicit hydroxyl hydrogen".to_string(),
            }
        })?;
    let alcohol_oxygen = alcohol_group.atoms[1];
    let alcohol_proton =
        first_bonded_hydrogen(alcohol_structure, alcohol_oxygen).ok_or_else(|| {
            ChemistryError::InvalidReaction {
                reaction_id: generated_pair_group_reaction_id(
                    "carboxylic_acid_esterification",
                    acid,
                    acid_group,
                    alcohol,
                    alcohol_group,
                ),
                reason: "alcohol has no explicit hydroxyl hydrogen".to_string(),
            }
        })?;

    let mut acid_editor = MolecularEditor::new(acid_structure);
    let acid_mapping = acid_editor.remove_atoms(&[acid_proton, acid_hydroxyl_oxygen])?;
    let acid_carbon = mapped_atom(&acid_mapping, acid_carbon, "acid carbon")?;
    let acid_fragment = acid_editor.finish()?;

    let mut alcohol_editor = MolecularEditor::new(alcohol_structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol_proton])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol_oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let product_structure = MolecularEditor::join_structures(
        &acid_fragment,
        acid_carbon,
        &alcohol_fragment,
        alcohol_oxygen,
        1.0,
    )?;
    let product = resolver.resolve(product_structure)?;
    Ok(Reaction::builder(generated_pair_group_reaction_id(
        "carboxylic_acid_esterification",
        acid,
        acid_group,
        alcohol,
        alcohol_group,
    ))
    .reactant(acid.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 1, 0)
    .catalyst_order("destroy:sulfuric_acid", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .build())
}

fn generate_nitrile_hydrolysis(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let nitrogen = group.atoms[1];
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "nitrile carbon")?;
    editor.add_atom(carbon, "O", 0.0, 2.0)?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "nitrile_hydrolysis",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(product, 1)
    .build())
}

fn generate_nitro_hydrogenation(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let nitrogen = group.atoms[1];
    let first_oxygen = group.atoms[2];
    let second_oxygen = group.atoms[3];
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[first_oxygen, second_oxygen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "nitro nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "nitro_hydrogenation",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen", 3, 1)
    .product("destroy:water", 2)
    .product(product, 1)
    .external_catalyst("forge:dusts/palladium", 1.0)
    .build())
}

fn first_bonded_hydrogen(structure: &MolecularStructure, atom: usize) -> Option<usize> {
    structure
        .neighbors(atom)
        .into_iter()
        .map(|(neighbor, _)| neighbor)
        .find(|neighbor| structure.atoms[*neighbor].element == "H")
}

fn mapped_atom(mapping: &[Option<usize>], old_index: usize, role: &str) -> ChemistryResult<usize> {
    mapping
        .get(old_index)
        .and_then(|value| *value)
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: "<generated-organic>".to_string(),
            reason: format!("{role} was removed during graph transformation"),
        })
}

fn generated_reaction_id(prefix: &str, substance: &Substance) -> String {
    format!("{prefix}/{}", sanitize_id(substance.id.as_str()))
}

fn generated_group_reaction_id(
    prefix: &str,
    substance: &Substance,
    group: &FunctionalGroup,
) -> String {
    format!(
        "{}/{}",
        generated_reaction_id(prefix, substance),
        atom_suffix(group)
    )
}

fn generated_pair_reaction_id(prefix: &str, first: &Substance, second: &Substance) -> String {
    format!(
        "{prefix}/{}/{}",
        sanitize_id(first.id.as_str()),
        sanitize_id(second.id.as_str())
    )
}

fn generated_pair_group_reaction_id(
    prefix: &str,
    first: &Substance,
    first_group: &FunctionalGroup,
    second: &Substance,
    second_group: &FunctionalGroup,
) -> String {
    format!(
        "{}/{}/{}",
        generated_pair_reaction_id(prefix, first, second),
        atom_suffix(first_group),
        atom_suffix(second_group)
    )
}

fn atom_suffix(group: &FunctionalGroup) -> String {
    group
        .atoms
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join("_")
}

fn sanitize_id(value: &str) -> String {
    value
        .chars()
        .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
        .collect()
}

fn stable_hash(value: &str) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    for byte in value.bytes() {
        hash ^= byte as u64;
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::OnceLock;

    fn generated_registry() -> ChemistryRegistry {
        static REGISTRY: OnceLock<ChemistryRegistry> = OnceLock::new();
        REGISTRY
            .get_or_init(|| {
                destroy_registry_with_generated_reactions_builder()
                    .unwrap()
                    .build()
                    .unwrap()
            })
            .clone()
    }

    #[test]
    fn generated_registry_builds_without_duplicate_derived_substances() {
        let registry = generated_registry();
        let mut canonical_codes = BTreeSet::new();
        for substance in registry.substances() {
            if !substance.id.as_str().starts_with("destroy:derived_") {
                continue;
            }
            if let Some(structure) = &substance.molecular_structure {
                assert!(canonical_codes.insert(structure.canonical_code()));
            }
        }
        assert!(registry.reactions().count() > super::super::DESTROY_REGISTERED_REACTION_COUNT);
    }

    #[test]
    fn halide_hydroxide_substitution_generates_ethanol_from_chloroethane() {
        let registry = generated_registry();
        let reaction = registry
            .reactions()
            .find(|reaction| {
                reaction
                    .id
                    .as_str()
                    .starts_with("halide_hydroxide_substitution/destroy_chloroethane/")
            })
            .unwrap();
        assert!(reaction
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:ethanol"));
        assert!(reaction
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:chloride"));
    }

    #[test]
    fn alcohol_oxidation_generates_acetone_from_isopropanol() {
        let registry = generated_registry();
        let reaction = registry
            .reactions()
            .find(|reaction| {
                reaction
                    .id
                    .as_str()
                    .starts_with("alcohol_oxidation/destroy_isopropanol/")
            })
            .unwrap();
        assert!(reaction
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:acetone"));
    }

    #[test]
    fn esterification_generates_product_from_acetic_acid_and_ethanol() {
        let registry = generated_registry();
        let reaction = registry
            .reactions()
            .find(|reaction| {
                reaction.id.as_str().starts_with(
                    "carboxylic_acid_esterification/destroy_acetic_acid/destroy_ethanol/",
                )
            })
            .unwrap();
        assert!(reaction
            .products
            .iter()
            .any(|term| term.substance_id.as_str() != "destroy:water"));
    }

    #[test]
    fn nitrile_hydrolysis_and_nitro_hydrogenation_are_registered() {
        let registry = generated_registry();
        assert!(registry.reactions().any(|reaction| {
            reaction
                .id
                .as_str()
                .starts_with("nitrile_hydrolysis/destroy_generic_nitrile/")
        }));
        let nitro = registry
            .reactions()
            .find(|reaction| {
                reaction
                    .id
                    .as_str()
                    .starts_with("nitro_hydrogenation/destroy_dinitrotoluene/")
            })
            .unwrap();
        assert!(!nitro.external_catalysts.is_empty());
    }
}
