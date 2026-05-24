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
pub(crate) struct GeneratedOrganicCatalog {
    pub(crate) substances: Vec<Substance>,
    pub(crate) reactions: Vec<Reaction>,
}

struct GenerationScope {
    substances: BTreeSet<SubstanceId>,
}

impl GenerationScope {
    fn all(registry: &ChemistryRegistry) -> Self {
        Self {
            substances: registry
                .substances()
                .map(|substance| substance.id.clone())
                .collect(),
        }
    }

    #[cfg(test)]
    fn from_substances(substances: &BTreeSet<SubstanceId>) -> Self {
        Self {
            substances: substances.clone(),
        }
    }

    fn contains(&self, substance_id: &SubstanceId) -> bool {
        self.substances.contains(substance_id)
    }
}

#[derive(Clone, Copy)]
pub(crate) struct GroupParticipant<'a> {
    pub(crate) substance: &'a Substance,
    pub(crate) structure: &'a MolecularStructure,
    pub(crate) group_index: usize,
}

impl<'a> GroupParticipant<'a> {
    fn group(self) -> ChemistryResult<&'a FunctionalGroup> {
        self.substance
            .functional_groups
            .get(self.group_index)
            .ok_or_else(|| ChemistryError::InvalidSubstance {
                substance_id: self.substance.id.to_string(),
                reason: format!(
                    "functional group index {} is outside the substance group list",
                    self.group_index
                ),
            })
    }

    fn is_seed(self, seeds: Option<&BTreeSet<SubstanceId>>) -> bool {
        seeds.is_none_or(|seeds| seeds.contains(&self.substance.id))
    }
}

pub(crate) struct OrganicGenerationSpace<'a> {
    all_substances: Vec<&'a Substance>,
    participants_by_group: BTreeMap<FunctionalGroupType, Vec<GroupParticipant<'a>>>,
}

impl<'a> OrganicGenerationSpace<'a> {
    fn new(
        substances: impl IntoIterator<Item = &'a Substance>,
        scope: &GenerationScope,
    ) -> ChemistryResult<Self> {
        let mut all_substances = Vec::new();
        let mut participants_by_group: BTreeMap<FunctionalGroupType, Vec<GroupParticipant<'a>>> =
            BTreeMap::new();

        for substance in substances {
            all_substances.push(substance);
            if !scope.contains(&substance.id) {
                continue;
            }
            let Some(structure) = substance.molecular_structure.as_ref() else {
                continue;
            };
            for (group_index, group) in substance.functional_groups.iter().enumerate() {
                participants_by_group
                    .entry(group.group_type.clone())
                    .or_default()
                    .push(GroupParticipant {
                        substance,
                        structure,
                        group_index,
                    });
            }
        }

        Ok(Self {
            all_substances,
            participants_by_group,
        })
    }

    pub(crate) fn from_participants(
        participants_by_group: BTreeMap<FunctionalGroupType, Vec<GroupParticipant<'a>>>,
    ) -> Self {
        Self {
            all_substances: Vec::new(),
            participants_by_group,
        }
    }

    fn participants(&self) -> impl Iterator<Item = GroupParticipant<'a>> + '_ {
        self.participants_by_group
            .values()
            .flat_map(|participants| participants.iter().copied())
    }

    fn participants_of(
        &self,
        group_type: &FunctionalGroupType,
    ) -> impl Iterator<Item = GroupParticipant<'a>> + '_ {
        self.participants_by_group
            .get(group_type)
            .into_iter()
            .flat_map(|participants| participants.iter().copied())
    }
}

struct DerivedSubstanceResolver {
    canonical_to_id: BTreeMap<String, SubstanceId>,
    generated_id_to_canonical: BTreeMap<SubstanceId, String>,
    substances: Vec<Substance>,
}

impl DerivedSubstanceResolver {
    fn new_from_canonical_to_id(canonical_to_id: BTreeMap<String, SubstanceId>) -> Self {
        Self {
            canonical_to_id,
            generated_id_to_canonical: BTreeMap::new(),
            substances: Vec::new(),
        }
    }

    fn resolve(&mut self, structure: MolecularStructure) -> ChemistryResult<SubstanceId> {
        let canonical = super::frowns::write_frowns(&structure)?;
        if let Some(id) = self.canonical_to_id.get(&canonical) {
            return Ok(id.clone());
        }
        let id = SubstanceId::new(canonical.clone())?;
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
        .with_catalog_metadata(Some(canonical.clone()), None, 0x20FF_FFFF, Vec::new())
        .with_molecular_structure(structure);
        self.canonical_to_id.insert(canonical.clone(), id.clone());
        self.generated_id_to_canonical.insert(id.clone(), canonical);
        self.substances.push(substance);
        Ok(id)
    }
}

pub(crate) fn generate_organic_reactions(
    registry: &ChemistryRegistry,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let scope = GenerationScope::all(registry);
    let space = OrganicGenerationSpace::new(registry.substances(), &scope)?;
    generate_organic_reactions_with_space(&space, None)
}

#[cfg(test)]
pub(crate) fn generate_organic_reactions_for_substances(
    substances: &[&Substance],
    seeds: &BTreeSet<SubstanceId>,
    scope: &BTreeSet<SubstanceId>,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let scope = GenerationScope::from_substances(scope);
    let space = OrganicGenerationSpace::new(substances.iter().copied(), &scope)?;
    generate_organic_reactions_with_space(&space, Some(seeds))
}

fn generate_organic_reactions_with_space(
    space: &OrganicGenerationSpace<'_>,
    seeds: Option<&BTreeSet<SubstanceId>>,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let canonical_to_id = canonical_to_id_from_substances(space.all_substances.iter().copied())?;
    let seed_participants = space
        .participants()
        .filter(|participant| participant.is_seed(seeds))
        .collect::<Vec<_>>();
    generate_organic_reactions_for_seed_participants(space, seed_participants, canonical_to_id)
}

pub(crate) fn generate_organic_reactions_for_seed_participants<'a>(
    space: &OrganicGenerationSpace<'a>,
    seed_participants: impl IntoIterator<Item = GroupParticipant<'a>>,
    canonical_to_id: BTreeMap<String, SubstanceId>,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let mut resolver = DerivedSubstanceResolver::new_from_canonical_to_id(canonical_to_id);
    let mut reactions = Vec::new();
    let mut reaction_ids = BTreeSet::new();

    for participant in seed_participants {
        let substance = participant.substance;
        let structure = participant.structure;
        let group = participant.group()?;
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
                let reaction = generate_halide_ammonia_substitution(
                    substance,
                    structure,
                    group,
                    &mut resolver,
                )?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_halide_cyanide_substitution(
                    substance,
                    structure,
                    group,
                    &mut resolver,
                )?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Alcohol => {
                if let Some(reaction) =
                    generate_alcohol_oxidation(substance, structure, group, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) =
                    generate_alcohol_dehydration(substance, structure, group, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                let reaction = generate_thionyl_chloride_substitution(
                    substance,
                    structure,
                    group,
                    &mut resolver,
                )?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Alkoxide => {
                let reaction =
                    generate_alkoxide_protonation(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Nitrile => {
                let reaction =
                    generate_nitrile_hydrolysis(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction =
                    generate_nitrile_hydrogenation(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Nitro => {
                let reaction =
                    generate_nitro_hydrogenation(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::AcylChloride => {
                let reaction =
                    generate_acyl_chloride_hydrolysis(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::CarboxylicAcid => {
                let reaction =
                    generate_acyl_chloride_formation(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Carbonyl => {
                if let Some(reaction) =
                    generate_aldehyde_oxidation(substance, structure, group, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                let reaction = generate_cyanide_nucleophilic_addition(
                    substance,
                    structure,
                    group,
                    &mut resolver,
                )?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction =
                    generate_wolff_kishner_reduction(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::UnsubstitutedAmide => {
                let reaction =
                    generate_amide_hydrolysis(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::PrimaryAmine => {
                let reaction =
                    generate_amine_phosgenation(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::NonTertiaryAmine => {
                let reaction =
                    generate_cyanamide_addition(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Isocyanate => {
                let reaction =
                    generate_isocyanate_hydrolysis(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Borane => {
                let reaction =
                    generate_borane_oxidation(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::BorateEster => {
                let reaction =
                    generate_borate_ester_hydrolysis(substance, structure, group, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            FunctionalGroupType::Alkene => {
                for spec in electrophilic_addition_specs(false) {
                    let reaction = generate_electrophilic_addition(
                        substance,
                        structure,
                        group,
                        spec,
                        &mut resolver,
                    )?;
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            FunctionalGroupType::Alkyne => {
                for spec in electrophilic_addition_specs(true) {
                    let reaction = generate_electrophilic_addition(
                        substance,
                        structure,
                        group,
                        spec,
                        &mut resolver,
                    )?;
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            _ => {}
        }

        generate_pair_reactions_for_seed(
            participant,
            space,
            &mut resolver,
            &mut reactions,
            &mut reaction_ids,
        )?;
    }

    Ok(GeneratedOrganicCatalog {
        substances: resolver.substances,
        reactions,
    })
}

fn canonical_to_id_from_substances<'a>(
    substances: impl IntoIterator<Item = &'a Substance>,
) -> ChemistryResult<BTreeMap<String, SubstanceId>> {
    let mut canonical_to_id = BTreeMap::new();
    for substance in substances {
        if let Some(structure) = &substance.molecular_structure {
            canonical_to_id
                .entry(super::frowns::write_frowns(structure)?)
                .or_insert_with(|| substance.id.clone());
        }
    }
    Ok(canonical_to_id)
}

fn generate_pair_reactions_for_seed(
    seed: GroupParticipant<'_>,
    space: &OrganicGenerationSpace<'_>,
    resolver: &mut DerivedSubstanceResolver,
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
) -> ChemistryResult<()> {
    match seed.group()?.group_type {
        FunctionalGroupType::CarboxylicAcid => {
            for alcohol in space.participants_of(&FunctionalGroupType::Alcohol) {
                let reaction = generate_carboxylic_acid_esterification(
                    seed.substance,
                    seed.structure,
                    seed.group()?,
                    alcohol.substance,
                    alcohol.structure,
                    alcohol.group()?,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        FunctionalGroupType::Alcohol => {
            for acid in space.participants_of(&FunctionalGroupType::CarboxylicAcid) {
                let reaction = generate_carboxylic_acid_esterification(
                    acid.substance,
                    acid.structure,
                    acid.group()?,
                    seed.substance,
                    seed.structure,
                    seed.group()?,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for acyl_chloride in space.participants_of(&FunctionalGroupType::AcylChloride) {
                let reaction = generate_acyl_chloride_esterification(
                    acyl_chloride.substance,
                    acyl_chloride.structure,
                    acyl_chloride.group()?,
                    seed.substance,
                    seed.structure,
                    seed.group()?,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        FunctionalGroupType::AcylChloride => {
            for alcohol in space.participants_of(&FunctionalGroupType::Alcohol) {
                let reaction = generate_acyl_chloride_esterification(
                    seed.substance,
                    seed.structure,
                    seed.group()?,
                    alcohol.substance,
                    alcohol.structure,
                    alcohol.group()?,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        FunctionalGroupType::Halide => {
            for amine in space.participants_of(&FunctionalGroupType::NonTertiaryAmine) {
                let reaction = generate_halide_amine_substitution(
                    seed.substance,
                    seed.structure,
                    seed.group()?,
                    amine.substance,
                    amine.structure,
                    amine.group()?,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        FunctionalGroupType::NonTertiaryAmine => {
            for halide in space.participants_of(&FunctionalGroupType::Halide) {
                let reaction = generate_halide_amine_substitution(
                    halide.substance,
                    halide.structure,
                    halide.group()?,
                    seed.substance,
                    seed.structure,
                    seed.group()?,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        _ => {}
    }
    Ok(())
}

fn push_unique_reaction(
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
    reaction: Reaction,
) -> ChemistryResult<()> {
    let id = reaction.id.to_string();
    if reaction_ids.insert(id.clone()) {
        reactions.push(reaction);
    }
    Ok(())
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

fn generate_acyl_chloride_formation(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let hydroxyl_oxygen = group.atoms[2];
    let proton = explicit_group_hydrogen(structure, group, 3, hydroxyl_oxygen, "carboxylic acid")?;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydroxyl_oxygen, proton])?;
    let carbon = mapped_atom(&mapping, carbon, "carboxylic acid carbon")?;
    editor.add_atom(carbon, "Cl", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "acyl_chloride_formation",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:phosgene", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .product("destroy:carbon_dioxide", 1)
    .build())
}

fn generate_acyl_chloride_hydrolysis(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let chlorine = group.atoms[2];
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[chlorine])?;
    let carbon = mapped_atom(&mapping, carbon, "acyl chloride carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "acyl_chloride_hydrolysis",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
}

fn generate_acyl_chloride_esterification(
    acyl_chloride: &Substance,
    acyl_chloride_structure: &MolecularStructure,
    acyl_chloride_group: &FunctionalGroup,
    alcohol: &Substance,
    alcohol_structure: &MolecularStructure,
    alcohol_group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let acyl_carbon = acyl_chloride_group.atoms[0];
    let chlorine = acyl_chloride_group.atoms[2];
    let alcohol_oxygen = alcohol_group.atoms[1];
    let alcohol_proton = explicit_group_hydrogen(
        alcohol_structure,
        alcohol_group,
        2,
        alcohol_oxygen,
        "alcohol",
    )?;
    let mut acyl_editor = MolecularEditor::new(acyl_chloride_structure);
    let acyl_mapping = acyl_editor.remove_atoms(&[chlorine])?;
    let acyl_carbon = mapped_atom(&acyl_mapping, acyl_carbon, "acyl chloride carbon")?;
    let acyl_fragment = acyl_editor.finish()?;

    let mut alcohol_editor = MolecularEditor::new(alcohol_structure);
    let alcohol_mapping = alcohol_editor.remove_atoms(&[alcohol_proton])?;
    let alcohol_oxygen = mapped_atom(&alcohol_mapping, alcohol_oxygen, "alcohol oxygen")?;
    let alcohol_fragment = alcohol_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &acyl_fragment,
        acyl_carbon,
        &alcohol_fragment,
        alcohol_oxygen,
        1.0,
    )?)?;
    Ok(Reaction::builder(generated_pair_group_reaction_id(
        "acyl_chloride_esterification",
        acyl_chloride,
        acyl_chloride_group,
        alcohol,
        alcohol_group,
    ))
    .reactant(acyl_chloride.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
}

fn generate_alcohol_dehydration(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let alcohol_carbon = group.atoms[0];
    let oxygen = group.atoms[1];
    let proton = explicit_group_hydrogen(structure, group, 2, oxygen, "alcohol")?;
    let mut products = Vec::new();
    for (neighbor, order) in structure.neighbors(alcohol_carbon) {
        if structure.atoms[neighbor].element != "C"
            || !super::molecule::bond_order_matches(order, 1.0)
        {
            continue;
        }
        let Some(beta_hydrogen) = first_bonded_hydrogen(structure, neighbor) else {
            continue;
        };
        let mut editor = MolecularEditor::new(structure);
        let mapping = editor.remove_atoms(&[beta_hydrogen, oxygen, proton])?;
        let carbon = mapped_atom(&mapping, alcohol_carbon, "alcohol carbon")?;
        let neighbor = mapped_atom(&mapping, neighbor, "dehydration neighbor carbon")?;
        editor.set_bond_order(carbon, neighbor, 2.0)?;
        products.push(resolver.resolve(editor.finish()?)?);
    }
    if products.is_empty() {
        return Ok(None);
    }
    let mut builder = Reaction::builder(generated_group_reaction_id(
        "alcohol_dehydration",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), products.len() as u32, 1)
    .reactant("destroy:oleum", products.len() as u32, 1)
    .product("destroy:sulfuric_acid", (products.len() * 2) as u32);
    for product in products {
        builder = builder.product(product, 1);
    }
    Ok(Some(builder.build()))
}

fn generate_alkoxide_protonation(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let oxygen = group.atoms[1];
    let mut editor = MolecularEditor::new(structure);
    editor.replace_atom(oxygen, "O", 0.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "alkoxide_protonation",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:proton", 1, 1)
    .product(product, 1)
    .build())
}

fn generate_thionyl_chloride_substitution(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let oxygen = group.atoms[1];
    let proton = explicit_group_hydrogen(structure, group, 2, oxygen, "alcohol")?;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[oxygen, proton])?;
    let carbon = mapped_atom(&mapping, carbon, "alcohol carbon")?;
    editor.add_atom(carbon, "Cl", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "thionyl_chloride_substitution",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:thionyl_chloride", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .product("destroy:sulfur_dioxide", 1)
    .build())
}

fn generate_aldehyde_oxidation(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    if group.is_ketone.unwrap_or(true) {
        return Ok(None);
    }
    let carbon = group.atoms[0];
    let Some(hydrogen) = first_bonded_hydrogen(structure, carbon) else {
        return Ok(None);
    };
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "aldehyde carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_group_reaction_id(
            "aldehyde_oxidation",
            substance,
            group,
        ))
        .reactant(substance.id.clone(), 3, 1)
        .reactant("destroy:dichromate", 1, 1)
        .reactant("destroy:proton", 8, 1)
        .product(product, 3)
        .product("destroy:chromium_iii", 2)
        .product("destroy:water", 4)
        .activation_energy_kj_per_mol(25.0)
        .build(),
    ))
}

fn generate_cyanide_nucleophilic_addition(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let oxygen = group.atoms[1];
    let mut editor = MolecularEditor::new(structure);
    editor.set_bond_order(carbon, oxygen, 1.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let nitrile_carbon = editor.add_atom(carbon, "C", 0.0, 1.0)?;
    editor.add_atom(nitrile_carbon, "N", 0.0, 3.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "cyanide_nucleophilic_addition",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen_cyanide", 1, 1)
    .catalyst_order("destroy:cyanide", 1)
    .product(product, 1)
    .build())
}

fn generate_wolff_kishner_reduction(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let oxygen = group.atoms[1];
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[oxygen])?;
    let carbon = mapped_atom(&mapping, carbon, "carbonyl carbon")?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "wolff_kishner_reduction",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrazine", 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .product("destroy:nitrogen", 1)
    .build())
}

fn generate_amide_hydrolysis(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let nitrogen = group.atoms[2];
    let hydrogens = explicit_group_hydrogens(structure, group, nitrogen, "amide nitrogen")?;
    if hydrogens.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_group_reaction_id("amide_hydrolysis", substance, group),
            reason: "unsubstituted amide must have exactly two explicit nitrogen hydrogens"
                .to_string(),
        });
    }
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen, hydrogens[0], hydrogens[1]])?;
    let carbon = mapped_atom(&mapping, carbon, "amide carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "amide_hydrolysis",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(product, 1)
    .product("destroy:ammonia", 1)
    .build())
}

fn generate_amine_phosgenation(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let nitrogen = group.atoms[1];
    let hydrogens = explicit_group_hydrogens(structure, group, nitrogen, "primary amine")?;
    if hydrogens.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_group_reaction_id("amine_phosgenation", substance, group),
            reason: "primary amine must have exactly two explicit hydrogens".to_string(),
        });
    }
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogens[0], hydrogens[1]])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "amine nitrogen")?;
    let carbon = editor.add_atom(nitrogen, "C", 0.0, 2.0)?;
    editor.add_atom(carbon, "O", 0.0, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "amine_phosgenation",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:phosgene", 1, 1)
    .product("destroy:hydrochloric_acid", 2)
    .product(product, 1)
    .build())
}

fn generate_cyanamide_addition(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let nitrogen = group.atoms[1];
    let hydrogen = explicit_group_hydrogen(structure, group, 2, nitrogen, "amine")?;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "amine nitrogen")?;
    let carbon = editor.add_atom(nitrogen, "C", 0.0, 1.0)?;
    let imine_nitrogen = editor.add_atom(carbon, "N", 0.0, 2.0)?;
    editor.add_atom(imine_nitrogen, "H", 0.0, 1.0)?;
    let amine_nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(amine_nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(amine_nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "cyanamide_addition",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:cyanamide", 1, 1)
    .product(product, 1)
    .build())
}

fn generate_halide_ammonia_substitution(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let halogen = group.atoms[1];
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[halogen])?;
    let carbon = mapped_atom(&mapping, carbon, "halide carbon")?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "halide_ammonia_substitution",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant(
        "destroy:ammonia",
        2,
        if group.degree.unwrap_or_default() == 3 {
            1
        } else {
            2
        },
    )
    .product(product, 1)
    .product(
        halide_ion(
            structure,
            halogen,
            "halide_ammonia_substitution",
            substance,
            group,
        )?,
        1,
    )
    .product("destroy:ammonium", 1)
    .build())
}

fn generate_halide_cyanide_substitution(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let halogen = group.atoms[1];
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[halogen])?;
    let carbon = mapped_atom(&mapping, carbon, "halide carbon")?;
    let nitrile_carbon = editor.add_atom(carbon, "C", 0.0, 1.0)?;
    editor.add_atom(nitrile_carbon, "N", 0.0, 3.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "halide_cyanide_substitution",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant(
        "destroy:cyanide",
        1,
        if group.degree.unwrap_or_default() == 3 {
            0
        } else {
            1
        },
    )
    .product(product, 1)
    .product(
        halide_ion(
            structure,
            halogen,
            "halide_cyanide_substitution",
            substance,
            group,
        )?,
        1,
    )
    .build())
}

fn generate_halide_amine_substitution(
    halide: &Substance,
    halide_structure: &MolecularStructure,
    halide_group: &FunctionalGroup,
    amine: &Substance,
    amine_structure: &MolecularStructure,
    amine_group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let halide_carbon = halide_group.atoms[0];
    let halogen = halide_group.atoms[1];
    let amine_nitrogen = amine_group.atoms[1];
    let amine_hydrogen =
        explicit_group_hydrogen(amine_structure, amine_group, 2, amine_nitrogen, "amine")?;
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
    Ok(Reaction::builder(generated_pair_group_reaction_id(
        "halide_amine_substitution",
        halide,
        halide_group,
        amine,
        amine_group,
    ))
    .reactant(halide.id.clone(), 1, 1)
    .reactant(amine.id.clone(), 1, 2)
    .product(product, 1)
    .product(
        halide_ion(
            halide_structure,
            halogen,
            "halide_amine_substitution",
            halide,
            halide_group,
        )?,
        1,
    )
    .product("destroy:proton", 1)
    .build())
}

fn generate_isocyanate_hydrolysis(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let nitrogen = group.atoms[1];
    let functional_carbon = group.atoms[2];
    let oxygen = group.atoms[3];
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[functional_carbon, oxygen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "isocyanate nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "isocyanate_hydrolysis",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product("destroy:carbon_dioxide", 1)
    .product(product, 1)
    .build())
}

fn generate_nitrile_hydrogenation(
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
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "nitrile_hydrogenation",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen", 2, 1)
    .product(product, 1)
    .external_catalyst("forge:dusts/nickel", 1.0)
    .build())
}

fn generate_borane_oxidation(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbon = group.atoms[0];
    let boron = group.atoms[1];
    let mut editor = MolecularEditor::new(structure);
    editor.insert_bridging_atom(carbon, boron, "O", 0.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_group_reaction_id(
        "borane_oxidation",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen_peroxide", 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .build())
}

fn generate_borate_ester_hydrolysis(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let oxygen = group.atoms[1];
    let boron = group.atoms[2];
    let (first, first_mapping, second, second_mapping) =
        MolecularEditor::split_at_bond(structure, oxygen, boron)?;
    let (boron_fragment, boron_mapping, alcohol_fragment, oxygen_mapping) =
        if first_mapping[boron].is_some() {
            (first, first_mapping, second, second_mapping)
        } else {
            (second, second_mapping, first, first_mapping)
        };

    let mut boron_editor = MolecularEditor::new(&boron_fragment);
    let boron = mapped_atom(&boron_mapping, boron, "borate boron")?;
    add_hydroxyl(&mut boron_editor, boron)?;
    let boron_product = resolver.resolve(boron_editor.finish()?)?;

    let mut alcohol_editor = MolecularEditor::new(&alcohol_fragment);
    let oxygen = mapped_atom(&oxygen_mapping, oxygen, "borate ester oxygen")?;
    alcohol_editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let alcohol_product = resolver.resolve(alcohol_editor.finish()?)?;

    Ok(Reaction::builder(generated_group_reaction_id(
        "borate_ester_hydrolysis",
        substance,
        group,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product(boron_product, 1)
    .product(alcohol_product, 1)
    .build())
}

#[derive(Debug, Clone, Copy)]
struct ElectrophilicAdditionSpec {
    prefix: &'static str,
    electrophile: &'static str,
    high_degree_group: AdditionGroup,
    low_degree_group: AdditionGroup,
    nucleophile_ratio: u32,
    activation_energy: f64,
    catalyst: Option<(&'static str, u32)>,
    external_catalyst: Option<&'static str>,
    display_as_reversible: bool,
}

#[derive(Debug, Clone, Copy)]
enum AdditionGroup {
    Atom(&'static str),
    Hydroxyl,
    Borane,
}

fn electrophilic_addition_specs(alkyne: bool) -> Vec<ElectrophilicAdditionSpec> {
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
            nucleophile_ratio: 1,
            activation_energy,
            catalyst: None,
            external_catalyst: None,
            display_as_reversible: false,
        },
    ]
}

fn generate_electrophilic_addition(
    substance: &Substance,
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    spec: ElectrophilicAdditionSpec,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let high_degree_carbon = group.atoms[0];
    let low_degree_carbon = group.atoms[1];
    let is_alkyne = group.group_type == FunctionalGroupType::Alkyne;
    let mut editor = MolecularEditor::new(structure);
    editor.set_bond_order(
        high_degree_carbon,
        low_degree_carbon,
        if is_alkyne { 2.0 } else { 1.0 },
    )?;
    add_addition_group(&mut editor, high_degree_carbon, spec.high_degree_group)?;
    add_addition_group(&mut editor, low_degree_carbon, spec.low_degree_group)?;
    let product = resolver.resolve(editor.finish()?)?;
    let mut builder = Reaction::builder(generated_group_reaction_id(spec.prefix, substance, group))
        .reactant(substance.id.clone(), spec.nucleophile_ratio, 1)
        .reactant(spec.electrophile, 1, 1)
        .product(product, spec.nucleophile_ratio)
        .activation_energy_kj_per_mol(spec.activation_energy);
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

fn add_hydroxyl(editor: &mut MolecularEditor, parent: usize) -> ChemistryResult<usize> {
    let oxygen = editor.add_atom(parent, "O", 0.0, 1.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    Ok(oxygen)
}

fn add_addition_group(
    editor: &mut MolecularEditor,
    parent: usize,
    group: AdditionGroup,
) -> ChemistryResult<()> {
    match group {
        AdditionGroup::Atom(element) => {
            editor.add_atom(parent, element, 0.0, 1.0)?;
        }
        AdditionGroup::Hydroxyl => {
            add_hydroxyl(editor, parent)?;
        }
        AdditionGroup::Borane => {
            let boron = editor.add_atom(parent, "B", 0.0, 1.0)?;
            editor.add_atom(boron, "H", 0.0, 1.0)?;
            editor.add_atom(boron, "H", 0.0, 1.0)?;
        }
    }
    Ok(())
}

fn explicit_group_hydrogen(
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    preferred_atom_position: usize,
    parent: usize,
    role: &str,
) -> ChemistryResult<usize> {
    if let Some(index) = group.atoms.get(preferred_atom_position) {
        if structure.atoms[*index].element == "H" {
            return Ok(*index);
        }
    }
    first_bonded_hydrogen(structure, parent).ok_or_else(|| ChemistryError::InvalidReaction {
        reaction_id: "<generated-organic>".to_string(),
        reason: format!("{role} has no explicit hydrogen"),
    })
}

fn explicit_group_hydrogens(
    structure: &MolecularStructure,
    group: &FunctionalGroup,
    parent: usize,
    role: &str,
) -> ChemistryResult<Vec<usize>> {
    let hydrogens = group
        .atoms
        .iter()
        .copied()
        .filter(|index| structure.atoms[*index].element == "H")
        .collect::<Vec<_>>();
    if hydrogens.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<generated-organic>".to_string(),
            reason: format!("{role} has no explicit hydrogen"),
        });
    }
    for hydrogen in &hydrogens {
        if !structure
            .neighbors(parent)
            .into_iter()
            .any(|(neighbor, _)| neighbor == *hydrogen)
        {
            return Err(ChemistryError::InvalidReaction {
                reaction_id: "<generated-organic>".to_string(),
                reason: format!("{role} group hydrogen is not bonded to the expected atom"),
            });
        }
    }
    Ok(hydrogens)
}

fn halide_ion(
    structure: &MolecularStructure,
    halogen: usize,
    prefix: &str,
    substance: &Substance,
    group: &FunctionalGroup,
) -> ChemistryResult<&'static str> {
    match structure.atoms[halogen].element.as_str() {
        "Cl" => Ok("destroy:chloride"),
        "F" => Ok("destroy:fluoride"),
        "I" => Ok("destroy:iodide"),
        _ => Err(ChemistryError::InvalidReaction {
            reaction_id: generated_group_reaction_id(prefix, substance, group),
            reason: "halide group does not contain a supported halogen".to_string(),
        }),
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::OnceLock;

    const ACTIVE_DESTROY_GENERIC_REACTIONS: &[&str] = &[
        "acyl_chloride_esterification",
        "acyl_chloride_formation",
        "acyl_chloride_hydrolysis",
        "alcohol_dehydration",
        "alcohol_oxidation",
        "aldehyde_oxidation",
        "alkene_chlorination",
        "alkene_chlorohydrination",
        "alkene_hydrolysis",
        "alkene_borane_hydroboration",
        "alkene_hydrochlorination",
        "alkene_hydrogenation",
        "alkene_hydroiodination",
        "alkene_iodination",
        "alkoxide_protonation",
        "alkyne_chlorination",
        "alkyne_chlorohydrination",
        "alkyne_hydrolysis",
        "alkyne_borane_hydroboration",
        "alkyne_hydrochlorination",
        "alkyne_hydrogenation",
        "alkyne_hydroiodination",
        "alkyne_iodination",
        "amide_hydrolysis",
        "amine_phosgenation",
        "borane_oxidation",
        "borate_ester_hydrolysis",
        "cyanamide_addition",
        "carboxylic_acid_esterification",
        "cyanide_nucleophilic_addition",
        "halide_amine_substitution",
        "halide_ammonia_substitution",
        "halide_cyanide_substitution",
        "halide_hydroxide_substitution",
        "isocyanate_hydrolysis",
        "nitrile_hydrogenation",
        "nitrile_hydrolysis",
        "nitro_hydrogenation",
        "thionyl_chloride_substitution",
        "wolff_kishner_reduction",
    ];

    const EXCLUDED_DESTROY_GENERIC_REACTIONS: &[&str] = &[
        "electrophilic_hydroboration",
        "borate_esterification",
        "borohydride_carbonyl_reduction",
        "carboxylic_acid_reduction",
    ];

    const ACTIVE_GENERATORS_WITHOUT_CATALOG_SUBSTRATE: &[&str] = &["aldehyde_oxidation"];

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

    fn reaction_with_prefix<'a>(registry: &'a ChemistryRegistry, prefix: &str) -> &'a Reaction {
        registry
            .reactions()
            .find(|reaction| reaction.id.as_str().starts_with(prefix))
            .unwrap_or_else(|| panic!("missing generated reaction with prefix {prefix}"))
    }

    #[test]
    fn generation_space_indexes_only_substances_inside_scope() {
        let registry = super::super::destroy_registry_builder()
            .unwrap()
            .build()
            .unwrap();
        let substances = registry.substances().collect::<Vec<_>>();
        let scope = GenerationScope::from_substances(&BTreeSet::from([SubstanceId::from(
            "destroy:acetic_acid",
        )]));
        let space = OrganicGenerationSpace::new(substances.iter().copied(), &scope).unwrap();

        let acids = space
            .participants_of(&FunctionalGroupType::CarboxylicAcid)
            .collect::<Vec<_>>();
        assert_eq!(acids.len(), 1);
        assert_eq!(acids[0].substance.id.as_str(), "destroy:acetic_acid");
        assert_eq!(
            space.participants_of(&FunctionalGroupType::Alcohol).count(),
            0
        );
    }

    #[test]
    fn scoped_generation_matches_full_static_generation() {
        let registry = super::super::destroy_registry_builder()
            .unwrap()
            .build()
            .unwrap();
        let full = generate_organic_reactions(&registry).unwrap();
        let substances = registry.substances().collect::<Vec<_>>();
        let all_ids = substances
            .iter()
            .map(|substance| substance.id.clone())
            .collect::<BTreeSet<_>>();
        let scoped =
            generate_organic_reactions_for_substances(&substances, &all_ids, &all_ids).unwrap();

        let full_substance_ids = full
            .substances
            .iter()
            .map(|substance| substance.id.as_str())
            .collect::<BTreeSet<_>>();
        let scoped_substance_ids = scoped
            .substances
            .iter()
            .map(|substance| substance.id.as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(full_substance_ids, scoped_substance_ids);

        let full_reaction_ids = full
            .reactions
            .iter()
            .map(|reaction| reaction.id.as_str())
            .collect::<BTreeSet<_>>();
        let scoped_reaction_ids = scoped
            .reactions
            .iter()
            .map(|reaction| reaction.id.as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(full_reaction_ids, scoped_reaction_ids);
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
                assert!(canonical_codes.insert(structure.canonical_code().unwrap()));
            }
        }
        assert!(registry.reactions().count() > super::super::DESTROY_REGISTERED_REACTION_COUNT);
    }

    #[test]
    fn active_destroy_generic_reactions_are_accounted_for() {
        assert_eq!(ACTIVE_DESTROY_GENERIC_REACTIONS.len(), 40);
        assert_eq!(EXCLUDED_DESTROY_GENERIC_REACTIONS.len(), 4);

        let registry = generated_registry();
        for prefix in ACTIVE_DESTROY_GENERIC_REACTIONS {
            if ACTIVE_GENERATORS_WITHOUT_CATALOG_SUBSTRATE.contains(prefix) {
                continue;
            }
            assert!(
                registry
                    .reactions()
                    .any(|reaction| reaction.id.as_str().starts_with(prefix)),
                "missing generated reaction for active Destroy generator {prefix}",
            );
        }
        for prefix in EXCLUDED_DESTROY_GENERIC_REACTIONS {
            assert!(
                !registry
                    .reactions()
                    .any(|reaction| reaction.id.as_str().starts_with(prefix)),
                "excluded Destroy generator {prefix} should not be registered",
            );
        }
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

    #[test]
    fn acyl_chloride_generators_are_registered() {
        let registry = generated_registry();
        let hydrolysis = reaction_with_prefix(
            &registry,
            "acyl_chloride_hydrolysis/destroy_generic_acyl_chloride/",
        );
        assert!(hydrolysis
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:hydrochloric_acid"));

        let formation =
            reaction_with_prefix(&registry, "acyl_chloride_formation/destroy_acetic_acid/");
        assert!(formation
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:carbon_dioxide"));

        let esterification = reaction_with_prefix(
            &registry,
            "acyl_chloride_esterification/destroy_generic_acyl_chloride/destroy_ethanol/",
        );
        assert!(esterification
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:hydrochloric_acid"));
    }

    #[test]
    fn halide_substitution_generators_are_registered() {
        let registry = generated_registry();
        let ammonia = reaction_with_prefix(
            &registry,
            "halide_ammonia_substitution/destroy_chloroethane/",
        );
        assert!(ammonia
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:ammonium"));

        let cyanide = reaction_with_prefix(
            &registry,
            "halide_cyanide_substitution/destroy_chloroethane/",
        );
        assert!(cyanide
            .reactants
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:cyanide"));

        let amine = reaction_with_prefix(
            &registry,
            "halide_amine_substitution/destroy_chloroethane/destroy_methylamine/",
        );
        assert!(amine
            .products
            .iter()
            .any(|term| term.substance_id.as_str() == "destroy:proton"));
    }

    #[test]
    fn electrophilic_addition_generators_are_registered() {
        let registry = generated_registry();
        for prefix in [
            "alkene_chlorination/destroy_ethene/",
            "alkene_chlorohydrination/destroy_ethene/",
            "alkene_hydrolysis/destroy_ethene/",
            "alkene_borane_hydroboration/destroy_ethene/",
            "alkene_hydrochlorination/destroy_ethene/",
            "alkene_hydrogenation/destroy_ethene/",
            "alkene_hydroiodination/destroy_ethene/",
            "alkene_iodination/destroy_ethene/",
            "alkyne_chlorination/destroy_acetylene/",
            "alkyne_hydrogenation/destroy_acetylene/",
        ] {
            reaction_with_prefix(&registry, prefix);
        }
        let hydrogenation = reaction_with_prefix(&registry, "alkene_hydrogenation/destroy_ethene/");
        assert!(!hydrogenation.external_catalysts.is_empty());
    }

    #[test]
    fn heteroatom_generators_are_registered() {
        let registry = generated_registry();
        reaction_with_prefix(&registry, "amide_hydrolysis/destroy_acetamide/");
        reaction_with_prefix(&registry, "amine_phosgenation/destroy_methylamine/");
        reaction_with_prefix(&registry, "cyanamide_addition/destroy_methylamine/");
        reaction_with_prefix(
            &registry,
            "isocyanate_hydrolysis/destroy_generic_isocyanate/",
        );
        reaction_with_prefix(&registry, "borane_oxidation/destroy_generic_borane/");
        reaction_with_prefix(
            &registry,
            "borate_ester_hydrolysis/destroy_generic_borate_ester/",
        );
        reaction_with_prefix(&registry, "nitrile_hydrogenation/destroy_generic_nitrile/");
        reaction_with_prefix(&registry, "thionyl_chloride_substitution/destroy_ethanol/");
        reaction_with_prefix(&registry, "wolff_kishner_reduction/destroy_acetone/");
    }
}
