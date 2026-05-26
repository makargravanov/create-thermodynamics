use std::collections::{BTreeMap, BTreeSet};

use super::condition::{AcidityCondition, AtmosphereCondition, ReactionCondition};
use super::error::{ChemistryError, ChemistryResult};
use super::kinetics::ReactionChannel;
use super::molecule::{
    MolecularEditor, MolecularStructure, StereoDescriptor, StereoMixtureKind, Stereochemistry,
    TetrahedralStereo,
};
use super::reaction::{Reaction, StoichiometricTerm};
use super::reactive_site::{try_find_reactive_sites, ReactiveSite, ReactiveSiteKind};
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

#[derive(Clone)]
pub(crate) struct SiteParticipant<'a> {
    pub(crate) substance: &'a Substance,
    pub(crate) structure: &'a MolecularStructure,
    pub(crate) site: ReactiveSite,
}

impl<'a> SiteParticipant<'a> {
    fn is_seed(&self, seeds: Option<&BTreeSet<SubstanceId>>) -> bool {
        seeds.is_none_or(|seeds| seeds.contains(&self.substance.id))
    }
}

pub(crate) struct OrganicGenerationSpace<'a> {
    all_substances: Vec<&'a Substance>,
    participants_by_site: BTreeMap<ReactiveSiteKind, Vec<SiteParticipant<'a>>>,
}

impl<'a> OrganicGenerationSpace<'a> {
    fn new(
        substances: impl IntoIterator<Item = &'a Substance>,
        scope: &GenerationScope,
    ) -> ChemistryResult<Self> {
        let mut all_substances = Vec::new();
        let mut participants_by_site: BTreeMap<ReactiveSiteKind, Vec<SiteParticipant<'a>>> =
            BTreeMap::new();

        for substance in substances {
            all_substances.push(substance);
            if !scope.contains(&substance.id) {
                continue;
            }
            let Some(structure) = substance.molecular_structure.as_ref() else {
                continue;
            };
            for site in try_find_reactive_sites(structure)? {
                participants_by_site
                    .entry(site.kind.clone())
                    .or_default()
                    .push(SiteParticipant {
                        substance,
                        structure,
                        site,
                    });
            }
        }

        Ok(Self {
            all_substances,
            participants_by_site,
        })
    }

    pub(crate) fn from_substances_for_scope(
        substances: impl IntoIterator<Item = &'a Substance>,
        scope: &BTreeSet<SubstanceId>,
    ) -> ChemistryResult<Self> {
        Self::new(
            substances,
            &GenerationScope {
                substances: scope.clone(),
            },
        )
    }

    fn site_participants(&self) -> impl Iterator<Item = SiteParticipant<'a>> + '_ {
        self.participants_by_site
            .values()
            .flat_map(|participants| participants.iter().cloned())
    }

    fn sites_of(&self, kind: &ReactiveSiteKind) -> impl Iterator<Item = SiteParticipant<'a>> + '_ {
        self.participants_by_site
            .get(kind)
            .into_iter()
            .flat_map(|participants| participants.iter().cloned())
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

#[derive(Clone)]
struct HalideSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    halogen: usize,
    degree: usize,
}

#[derive(Clone)]
struct AlcoholSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    oxygen: usize,
    hydrogen: usize,
    degree: usize,
}

#[derive(Clone)]
struct AlkoxideSite<'a> {
    participant: SiteParticipant<'a>,
    oxygen: usize,
}

#[derive(Clone)]
struct CarbonylSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    oxygen: usize,
    is_ketone: bool,
}

#[derive(Clone)]
struct CarboxylicAcidSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    hydroxyl_oxygen: usize,
    hydroxyl_hydrogen: usize,
}

#[derive(Clone)]
struct AcylChlorideSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    chlorine: usize,
}

#[derive(Clone)]
struct AmideSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    nitrogen: usize,
    nitrogen_hydrogens: Vec<usize>,
}

#[derive(Clone)]
struct AmineSite<'a> {
    participant: SiteParticipant<'a>,
    nitrogen: usize,
    hydrogens: Vec<usize>,
}

#[derive(Clone)]
struct NitrileSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    nitrogen: usize,
}

#[derive(Clone)]
struct NitroSite<'a> {
    participant: SiteParticipant<'a>,
    nitrogen: usize,
    oxygens: [usize; 2],
}

#[derive(Clone)]
struct UnsaturatedBondSite<'a> {
    participant: SiteParticipant<'a>,
    high_degree_carbon: usize,
    low_degree_carbon: usize,
    is_alkyne: bool,
}

#[derive(Clone)]
struct BoraneSite<'a> {
    participant: SiteParticipant<'a>,
    carbon: usize,
    boron: usize,
}

#[derive(Clone)]
struct BorateEsterSite<'a> {
    participant: SiteParticipant<'a>,
    oxygen: usize,
    boron: usize,
}

#[derive(Clone)]
struct IsocyanateSite<'a> {
    participant: SiteParticipant<'a>,
    nitrogen: usize,
    functional_carbon: usize,
    oxygen: usize,
}

impl<'a> SiteParticipant<'a> {
    fn require_kind(&self, expected: ReactiveSiteKind) -> ChemistryResult<()> {
        if self.site.kind == expected {
            Ok(())
        } else {
            Err(ChemistryError::InvalidReaction {
                reaction_id: generated_site_reaction_id("typed_site", self),
                reason: format!(
                    "expected {:?} reactive site, got {:?}",
                    expected, self.site.kind
                ),
            })
        }
    }

    fn halide_site(self) -> ChemistryResult<HalideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Halide)?;
        let carbon = self.site_atom_by_element("C", "halide carbon")?;
        let halogen = self
            .site
            .leaving_atom
            .or_else(|| {
                self.site.atoms.iter().copied().find(|atom| {
                    matches!(
                        self.structure.atoms[*atom].element.as_str(),
                        "F" | "Cl" | "I"
                    )
                })
            })
            .ok_or_else(|| self.site_error("halide site has no supported halogen"))?;
        Ok(HalideSite {
            degree: self
                .site
                .substitution_degree
                .unwrap_or_else(|| self.structure.carbon_degree(carbon)),
            participant: self,
            carbon,
            halogen,
        })
    }

    fn alcohol_site(self) -> ChemistryResult<AlcoholSite<'a>> {
        self.require_kind(ReactiveSiteKind::Alcohol)?;
        let oxygen = self.site_atom_by_element("O", "alcohol oxygen")?;
        let carbon = self.bonded_site_atom(oxygen, "C", 1.0, "alcohol carbon")?;
        let hydrogen = first_bonded_hydrogen(self.structure, oxygen)
            .ok_or_else(|| self.site_error("alcohol oxygen has no explicit hydrogen"))?;
        Ok(AlcoholSite {
            degree: self
                .site
                .substitution_degree
                .unwrap_or_else(|| self.structure.carbon_degree(carbon)),
            participant: self,
            carbon,
            oxygen,
            hydrogen,
        })
    }

    fn alkoxide_site(self) -> ChemistryResult<AlkoxideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Alkoxide)?;
        let oxygen = self.site_atom_by_element("O", "alkoxide oxygen")?;
        self.bonded_site_atom(oxygen, "C", 1.0, "alkoxide carbon")?;
        Ok(AlkoxideSite {
            participant: self,
            oxygen,
        })
    }

    fn carbonyl_site(self) -> ChemistryResult<CarbonylSite<'a>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl
        ) {
            return Err(self.site_error("site is not a carbonyl center"));
        }
        let (carbon, oxygen) = carbonyl_atoms_from_site(self.structure, &self.site, "carbonyl")?;
        let carbon_neighbors = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .filter(|(neighbor, order)| {
                *neighbor != oxygen
                    && self.structure.atoms[*neighbor].element == "C"
                    && super::molecule::bond_order_matches(*order, 1.0)
            })
            .count();
        Ok(CarbonylSite {
            is_ketone: self.site.kind == ReactiveSiteKind::Ketone
                || (self.site.kind == ReactiveSiteKind::Carbonyl && carbon_neighbors >= 2),
            participant: self,
            carbon,
            oxygen,
        })
    }

    fn carboxylic_acid_site(self) -> ChemistryResult<CarboxylicAcidSite<'a>> {
        self.require_kind(ReactiveSiteKind::CarboxylicAcid)?;
        let (carbon, carbonyl_oxygen) =
            carbonyl_atoms_from_site(self.structure, &self.site, "carboxylic acid")?;
        let hydroxyl_oxygen = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != carbonyl_oxygen
                    && self.structure.atoms[neighbor].element == "O"
                    && super::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("carboxylic acid has no hydroxyl oxygen"))?;
        let hydroxyl_hydrogen = first_bonded_hydrogen(self.structure, hydroxyl_oxygen)
            .ok_or_else(|| self.site_error("carboxylic acid has no explicit hydroxyl hydrogen"))?;
        Ok(CarboxylicAcidSite {
            participant: self,
            carbon,
            hydroxyl_oxygen,
            hydroxyl_hydrogen,
        })
    }

    fn acyl_chloride_site(self) -> ChemistryResult<AcylChlorideSite<'a>> {
        self.require_kind(ReactiveSiteKind::AcylChloride)?;
        let (carbon, _) = carbonyl_atoms_from_site(self.structure, &self.site, "acyl chloride")?;
        let chlorine = self.bonded_site_atom(carbon, "Cl", 1.0, "acyl chloride chlorine")?;
        Ok(AcylChlorideSite {
            participant: self,
            carbon,
            chlorine,
        })
    }

    fn amide_site(self) -> ChemistryResult<AmideSite<'a>> {
        self.require_kind(ReactiveSiteKind::Amide)?;
        let (carbon, oxygen) = carbonyl_atoms_from_site(self.structure, &self.site, "amide")?;
        let nitrogen = self
            .structure
            .neighbors(carbon)
            .into_iter()
            .find_map(|(neighbor, order)| {
                (neighbor != oxygen
                    && self.structure.atoms[neighbor].element == "N"
                    && super::molecule::bond_order_matches(order, 1.0))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error("amide has no nitrogen atom"))?;
        let nitrogen_hydrogens = bonded_hydrogens(self.structure, nitrogen);
        Ok(AmideSite {
            participant: self,
            carbon,
            nitrogen,
            nitrogen_hydrogens,
        })
    }

    fn amine_site(self) -> ChemistryResult<AmineSite<'a>> {
        if !matches!(
            self.site.kind,
            ReactiveSiteKind::PrimaryAmine | ReactiveSiteKind::NonTertiaryAmine
        ) {
            return Err(self.site_error("site is not an amine center"));
        }
        let nitrogen = self.site_atom_by_element("N", "amine nitrogen")?;
        let hydrogens = bonded_hydrogens(self.structure, nitrogen);
        if hydrogens.is_empty() {
            return Err(self.site_error("amine has no explicit nitrogen hydrogen"));
        }
        Ok(AmineSite {
            participant: self,
            nitrogen,
            hydrogens,
        })
    }

    fn nitrile_site(self) -> ChemistryResult<NitrileSite<'a>> {
        self.require_kind(ReactiveSiteKind::Nitrile)?;
        let carbon = self.site_atom_by_element("C", "nitrile carbon")?;
        let nitrogen = self.bonded_site_atom(carbon, "N", 3.0, "nitrile nitrogen")?;
        Ok(NitrileSite {
            participant: self,
            carbon,
            nitrogen,
        })
    }

    fn nitro_site(self) -> ChemistryResult<NitroSite<'a>> {
        self.require_kind(ReactiveSiteKind::Nitro)?;
        let nitrogen = self.site_atom_by_element("N", "nitro nitrogen")?;
        let oxygens = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| self.structure.atoms[*atom].element == "O")
            .collect::<Vec<_>>();
        let oxygens: [usize; 2] = oxygens
            .try_into()
            .map_err(|_| self.site_error("nitro center must contain exactly two oxygens"))?;
        Ok(NitroSite {
            participant: self,
            nitrogen,
            oxygens,
        })
    }

    fn unsaturated_bond_site(self) -> ChemistryResult<UnsaturatedBondSite<'a>> {
        let is_alkyne = match self.site.kind {
            ReactiveSiteKind::Alkene => false,
            ReactiveSiteKind::Alkyne => true,
            _ => return Err(self.site_error("site is not an unsaturated bond")),
        };
        let carbons = self
            .site
            .atoms
            .iter()
            .copied()
            .filter(|atom| self.structure.atoms[*atom].element == "C")
            .collect::<Vec<_>>();
        if carbons.len() != 2 {
            return Err(self.site_error("unsaturated bond must contain exactly two carbons"));
        }
        let first_degree = self.structure.carbon_degree(carbons[0]).saturating_sub(1);
        let second_degree = self.structure.carbon_degree(carbons[1]).saturating_sub(1);
        let (high_degree_carbon, low_degree_carbon) = if second_degree > first_degree {
            (carbons[1], carbons[0])
        } else {
            (carbons[0], carbons[1])
        };
        Ok(UnsaturatedBondSite {
            participant: self,
            high_degree_carbon,
            low_degree_carbon,
            is_alkyne,
        })
    }

    fn borane_site(self) -> ChemistryResult<BoraneSite<'a>> {
        self.require_kind(ReactiveSiteKind::Borane)?;
        let carbon = self.site_atom_by_element("C", "borane carbon")?;
        let boron = self.bonded_site_atom(carbon, "B", 1.0, "borane boron")?;
        Ok(BoraneSite {
            participant: self,
            carbon,
            boron,
        })
    }

    fn borate_ester_site(self) -> ChemistryResult<BorateEsterSite<'a>> {
        self.require_kind(ReactiveSiteKind::BorateEster)?;
        let oxygen = self.site_atom_by_element("O", "borate ester oxygen")?;
        let boron = self.bonded_site_atom(oxygen, "B", 1.0, "borate ester boron")?;
        Ok(BorateEsterSite {
            participant: self,
            oxygen,
            boron,
        })
    }

    fn isocyanate_site(self) -> ChemistryResult<IsocyanateSite<'a>> {
        self.require_kind(ReactiveSiteKind::Isocyanate)?;
        let nitrogen = self.site_atom_by_element("N", "isocyanate nitrogen")?;
        let functional_carbon = self.bonded_site_atom(nitrogen, "C", 2.0, "isocyanate carbon")?;
        let oxygen = self.bonded_site_atom(functional_carbon, "O", 2.0, "isocyanate oxygen")?;
        Ok(IsocyanateSite {
            participant: self,
            nitrogen,
            functional_carbon,
            oxygen,
        })
    }

    fn site_atom_by_element(&self, element: &str, label: &str) -> ChemistryResult<usize> {
        self.site
            .atoms
            .iter()
            .copied()
            .find(|atom| self.structure.atoms[*atom].element == element)
            .ok_or_else(|| self.site_error(&format!("reactive site is missing {label}")))
    }

    fn bonded_site_atom(
        &self,
        parent: usize,
        element: &str,
        order: f64,
        label: &str,
    ) -> ChemistryResult<usize> {
        self.structure
            .neighbors(parent)
            .into_iter()
            .find_map(|(neighbor, bond_order)| {
                (self.site.atoms.contains(&neighbor)
                    && self.structure.atoms[neighbor].element == element
                    && super::molecule::bond_order_matches(bond_order, order))
                .then_some(neighbor)
            })
            .ok_or_else(|| self.site_error(&format!("reactive site is missing {label}")))
    }

    fn site_error(&self, reason: &str) -> ChemistryError {
        ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("typed_site", self),
            reason: reason.to_string(),
        }
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
    let seed_ids = seeds.cloned().unwrap_or_else(|| {
        space
            .all_substances
            .iter()
            .map(|substance| substance.id.clone())
            .collect()
    });
    generate_organic_reactions_for_seed_substances(space, &seed_ids, canonical_to_id)
}

pub(crate) fn generate_organic_reactions_for_seed_participants<'a>(
    space: &OrganicGenerationSpace<'a>,
    seed_participants: impl IntoIterator<Item = SiteParticipant<'a>>,
    canonical_to_id: BTreeMap<String, SubstanceId>,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let mut resolver = DerivedSubstanceResolver::new_from_canonical_to_id(canonical_to_id);
    let mut reactions = Vec::new();
    let mut reaction_ids = BTreeSet::new();

    for participant in seed_participants {
        match participant.site.kind {
            ReactiveSiteKind::Halide => {
                let site = participant.clone().halide_site()?;
                if let Some(reaction) =
                    generate_halide_hydroxide_substitution(&site, &mut resolver)?
                {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                let reaction = generate_halide_ammonia_substitution(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_halide_cyanide_substitution(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Alcohol => {
                let site = participant.clone().alcohol_site()?;
                if let Some(reaction) = generate_alcohol_oxidation(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                if let Some(reaction) = generate_alcohol_dehydration(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                let reaction = generate_thionyl_chloride_substitution(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Alkoxide => {
                let site = participant.clone().alkoxide_site()?;
                let reaction = generate_alkoxide_protonation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Nitrile => {
                let site = participant.clone().nitrile_site()?;
                let reaction = generate_nitrile_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_nitrile_hydrogenation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Nitro => {
                let site = participant.clone().nitro_site()?;
                let reaction = generate_nitro_hydrogenation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::AcylChloride => {
                let site = participant.clone().acyl_chloride_site()?;
                let reaction = generate_acyl_chloride_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::CarboxylicAcid => {
                let site = participant.clone().carboxylic_acid_site()?;
                let reaction = generate_acyl_chloride_formation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => {
                let site = participant.clone().carbonyl_site()?;
                if let Some(reaction) = generate_aldehyde_oxidation(&site, &mut resolver)? {
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
                let reaction = generate_cyanide_nucleophilic_addition(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                let reaction = generate_wolff_kishner_reduction(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Amide => {
                let site = participant.clone().amide_site()?;
                let reaction = generate_amide_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::PrimaryAmine => {
                let site = participant.clone().amine_site()?;
                let reaction = generate_amine_phosgenation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::NonTertiaryAmine => {
                let site = participant.clone().amine_site()?;
                let reaction = generate_cyanamide_addition(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Isocyanate => {
                let site = participant.clone().isocyanate_site()?;
                let reaction = generate_isocyanate_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Borane => {
                let site = participant.clone().borane_site()?;
                let reaction = generate_borane_oxidation(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::BorateEster => {
                let site = participant.clone().borate_ester_site()?;
                let reaction = generate_borate_ester_hydrolysis(&site, &mut resolver)?;
                push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Alkene => {
                let site = participant.clone().unsaturated_bond_site()?;
                for spec in electrophilic_addition_specs(false) {
                    let reaction = match generate_electrophilic_addition(&site, spec, &mut resolver)
                    {
                        Ok(reaction) => reaction,
                        Err(error) if is_unknown_stereo_distribution(&error) => continue,
                        Err(error) => return Err(error),
                    };
                    push_unique_reaction(&mut reactions, &mut reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Alkyne => {
                let site = participant.clone().unsaturated_bond_site()?;
                for spec in electrophilic_addition_specs(true) {
                    let reaction = match generate_electrophilic_addition(&site, spec, &mut resolver)
                    {
                        Ok(reaction) => reaction,
                        Err(error) if is_unknown_stereo_distribution(&error) => continue,
                        Err(error) => return Err(error),
                    };
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

pub(crate) fn generate_organic_reactions_for_seed_substances<'a>(
    space: &OrganicGenerationSpace<'a>,
    seeds: &BTreeSet<SubstanceId>,
    canonical_to_id: BTreeMap<String, SubstanceId>,
) -> ChemistryResult<GeneratedOrganicCatalog> {
    let seed_participants = space
        .site_participants()
        .filter(|participant| participant.is_seed(Some(seeds)))
        .filter(|participant| is_generator_seed_site(&participant.site.kind))
        .collect::<Vec<_>>();
    let mut generated = generate_organic_reactions_for_seed_participants(
        space,
        seed_participants,
        canonical_to_id,
    )?;
    let canonical_to_id = canonical_to_id_from_generated(space, &generated)?;
    let mut resolver = DerivedSubstanceResolver::new_from_canonical_to_id(canonical_to_id);
    let mut reaction_ids = generated
        .reactions
        .iter()
        .map(|reaction| reaction.id.to_string())
        .collect::<BTreeSet<_>>();
    generate_site_reactions_for_seed_participants(
        space,
        space
            .site_participants()
            .filter(|participant| participant.is_seed(Some(seeds)))
            .collect::<Vec<_>>(),
        &mut resolver,
        &mut generated.reactions,
        &mut reaction_ids,
    )?;
    generated.substances.extend(resolver.substances);
    Ok(generated)
}

fn is_generator_seed_site(kind: &ReactiveSiteKind) -> bool {
    matches!(
        kind,
        ReactiveSiteKind::Halide
            | ReactiveSiteKind::Alcohol
            | ReactiveSiteKind::Alkoxide
            | ReactiveSiteKind::Nitrile
            | ReactiveSiteKind::Nitro
            | ReactiveSiteKind::AcylChloride
            | ReactiveSiteKind::CarboxylicAcid
            | ReactiveSiteKind::Aldehyde
            | ReactiveSiteKind::Ketone
            | ReactiveSiteKind::Carbonyl
            | ReactiveSiteKind::Amide
            | ReactiveSiteKind::PrimaryAmine
            | ReactiveSiteKind::NonTertiaryAmine
            | ReactiveSiteKind::Isocyanate
            | ReactiveSiteKind::Borane
            | ReactiveSiteKind::BorateEster
            | ReactiveSiteKind::Alkene
            | ReactiveSiteKind::Alkyne
    )
}

fn canonical_to_id_from_generated(
    space: &OrganicGenerationSpace<'_>,
    generated: &GeneratedOrganicCatalog,
) -> ChemistryResult<BTreeMap<String, SubstanceId>> {
    let mut canonical_to_id =
        canonical_to_id_from_substances(space.all_substances.iter().copied())?;
    for substance in &generated.substances {
        if let Some(structure) = &substance.molecular_structure {
            canonical_to_id
                .entry(super::frowns::write_frowns(structure)?)
                .or_insert_with(|| substance.id.clone());
        }
    }
    Ok(canonical_to_id)
}

fn is_unknown_stereo_distribution(error: &ChemistryError) -> bool {
    matches!(
        error,
        ChemistryError::InvalidReaction { reason, .. }
            if reason.contains("stereo distribution")
                || reason.contains("stereo mixture has no quantitative distribution")
    )
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
    seed: SiteParticipant<'_>,
    space: &OrganicGenerationSpace<'_>,
    resolver: &mut DerivedSubstanceResolver,
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
) -> ChemistryResult<()> {
    match seed.site.kind {
        ReactiveSiteKind::CarboxylicAcid => {
            let acid_site = seed.clone().carboxylic_acid_site()?;
            for alcohol in space.sites_of(&ReactiveSiteKind::Alcohol) {
                let alcohol_site = alcohol.alcohol_site()?;
                let reaction =
                    generate_carboxylic_acid_esterification(&acid_site, &alcohol_site, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        ReactiveSiteKind::Alcohol => {
            let alcohol_site = seed.clone().alcohol_site()?;
            for acid in space.sites_of(&ReactiveSiteKind::CarboxylicAcid) {
                let acid_site = acid.carboxylic_acid_site()?;
                let reaction =
                    generate_carboxylic_acid_esterification(&acid_site, &alcohol_site, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for acyl_chloride in space.sites_of(&ReactiveSiteKind::AcylChloride) {
                let acyl_chloride_site = acyl_chloride.acyl_chloride_site()?;
                let reaction = generate_acyl_chloride_esterification(
                    &acyl_chloride_site,
                    &alcohol_site,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for carbonyl_kind in carbonyl_site_kinds() {
                for carbonyl in space.sites_of(&carbonyl_kind) {
                    let carbonyl_site = carbonyl.carbonyl_site()?;
                    let reaction =
                        generate_acetal_formation(&carbonyl_site, &alcohol_site, resolver)?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        ReactiveSiteKind::AcylChloride => {
            let acyl_chloride_site = seed.clone().acyl_chloride_site()?;
            for alcohol in space.sites_of(&ReactiveSiteKind::Alcohol) {
                let alcohol_site = alcohol.alcohol_site()?;
                let reaction = generate_acyl_chloride_esterification(
                    &acyl_chloride_site,
                    &alcohol_site,
                    resolver,
                )?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        ReactiveSiteKind::Halide => {
            let halide_site = seed.clone().halide_site()?;
            for amine in space.sites_of(&ReactiveSiteKind::NonTertiaryAmine) {
                let amine_site = amine.amine_site()?;
                let reaction =
                    generate_halide_amine_substitution(&halide_site, &amine_site, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        ReactiveSiteKind::NonTertiaryAmine => {
            let amine_site = seed.clone().amine_site()?;
            for halide in space.sites_of(&ReactiveSiteKind::Halide) {
                let halide_site = halide.halide_site()?;
                let reaction =
                    generate_halide_amine_substitution(&halide_site, &amine_site, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => {
            let carbonyl_site = seed.clone().carbonyl_site()?;
            for alcohol in space.sites_of(&ReactiveSiteKind::Alcohol) {
                let alcohol_site = alcohol.alcohol_site()?;
                let reaction = generate_acetal_formation(&carbonyl_site, &alcohol_site, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            for amine in space.sites_of(&ReactiveSiteKind::PrimaryAmine) {
                let amine_site = amine.amine_site()?;
                let reaction = generate_imine_formation(&carbonyl_site, &amine_site, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
        }
        ReactiveSiteKind::PrimaryAmine => {
            let amine_site = seed.clone().amine_site()?;
            for carbonyl_kind in carbonyl_site_kinds() {
                for carbonyl in space.sites_of(&carbonyl_kind) {
                    let carbonyl_site = carbonyl.carbonyl_site()?;
                    let reaction = generate_imine_formation(&carbonyl_site, &amine_site, resolver)?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
        }
        _ => {}
    }
    Ok(())
}

fn carbonyl_site_kinds() -> [ReactiveSiteKind; 3] {
    [
        ReactiveSiteKind::Aldehyde,
        ReactiveSiteKind::Ketone,
        ReactiveSiteKind::Carbonyl,
    ]
}

fn generate_site_reactions_for_seed_participants<'a>(
    space: &OrganicGenerationSpace<'a>,
    seed_sites: impl IntoIterator<Item = SiteParticipant<'a>>,
    resolver: &mut DerivedSubstanceResolver,
    reactions: &mut Vec<Reaction>,
    reaction_ids: &mut BTreeSet<String>,
) -> ChemistryResult<()> {
    for seed in seed_sites {
        match seed.site.kind {
            ReactiveSiteKind::Aldehyde | ReactiveSiteKind::Ketone | ReactiveSiteKind::Carbonyl => {
                for organometallic_kind in [
                    ReactiveSiteKind::Organomagnesium,
                    ReactiveSiteKind::Organolithium,
                    ReactiveSiteKind::Organocopper,
                ] {
                    for organometallic in space.sites_of(&organometallic_kind) {
                        let reaction = generate_organometallic_carbonyl_addition(
                            seed.clone(),
                            organometallic,
                            resolver,
                        )?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
                for enol in space.sites_of(&ReactiveSiteKind::Enol) {
                    let reaction = generate_aldol_addition(enol, seed.clone(), resolver)?;
                    push_unique_reaction(reactions, reaction_ids, reaction)?;
                }
            }
            ReactiveSiteKind::Organomagnesium
            | ReactiveSiteKind::Organolithium
            | ReactiveSiteKind::Organocopper => {
                for carbonyl_kind in [
                    ReactiveSiteKind::Aldehyde,
                    ReactiveSiteKind::Ketone,
                    ReactiveSiteKind::Carbonyl,
                ] {
                    for carbonyl in space.sites_of(&carbonyl_kind) {
                        let reaction = generate_organometallic_carbonyl_addition(
                            carbonyl,
                            seed.clone(),
                            resolver,
                        )?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::Enol => {
                for carbonyl_kind in [
                    ReactiveSiteKind::Aldehyde,
                    ReactiveSiteKind::Ketone,
                    ReactiveSiteKind::Carbonyl,
                ] {
                    for carbonyl in space.sites_of(&carbonyl_kind) {
                        let reaction = generate_aldol_addition(seed.clone(), carbonyl, resolver)?;
                        push_unique_reaction(reactions, reaction_ids, reaction)?;
                    }
                }
            }
            ReactiveSiteKind::AromaticRing => {
                let reaction = generate_aromatic_nitration(seed, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            ReactiveSiteKind::Epoxide => {
                let reaction = generate_epoxide_hydrolysis(seed, resolver)?;
                push_unique_reaction(reactions, reaction_ids, reaction)?;
            }
            _ => {}
        }
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
    site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
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
        .reactant("destroy:hydroxide", 1, if site.degree == 3 { 0 } else { 1 })
        .product(product, 1)
        .product(halide_ion, 1)
        .build(),
    ))
}

fn generate_alcohol_oxidation(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    if site.degree >= 3 {
        return Ok(None);
    }
    let carbon = site.carbon;
    let oxygen = site.oxygen;
    let Some(carbon_hydrogen) = first_bonded_hydrogen(structure, carbon) else {
        return Ok(None);
    };
    let oxygen_hydrogen = site.hydrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[carbon_hydrogen, oxygen_hydrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "alcohol carbon")?;
    let oxygen = mapped_atom(&mapping, oxygen, "alcohol oxygen")?;
    editor.set_bond_order(carbon, oxygen, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "alcohol_oxidation",
            &site.participant,
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
    acid_site: &CarboxylicAcidSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let acid = acid_site.participant.substance;
    let acid_structure = acid_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let alcohol_structure = alcohol_site.participant.structure;
    let acid_carbon = acid_site.carbon;
    let acid_hydroxyl_oxygen = acid_site.hydroxyl_oxygen;
    let acid_proton = acid_site.hydroxyl_hydrogen;
    let alcohol_oxygen = alcohol_site.oxygen;
    let alcohol_proton = alcohol_site.hydrogen;

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
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "carboxylic_acid_esterification",
        &acid_site.participant,
        &alcohol_site.participant,
    ))
    .reactant(acid.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 1, 0)
    .catalyst_order("destroy:sulfuric_acid", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .build())
}

fn generate_acetal_formation(
    carbonyl_site: &CarbonylSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbonyl = carbonyl_site.participant.substance;
    let carbonyl_structure = carbonyl_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let carbonyl_carbon = carbonyl_site.carbon;
    let carbonyl_oxygen = carbonyl_site.oxygen;
    let (alcohol_fragment, alcohol_oxygen) =
        deprotonated_alcohol_fragment(alcohol_site, "acetal formation")?;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_oxygen])?;
    let carbonyl_carbon = mapped_atom(&carbonyl_mapping, carbonyl_carbon, "carbonyl carbon")?;
    let mut product_editor = MolecularEditor::new(&carbonyl_editor.finish()?);
    product_editor.add_group(carbonyl_carbon, &alcohol_fragment, alcohol_oxygen, 1.0)?;
    product_editor.add_group(carbonyl_carbon, &alcohol_fragment, alcohol_oxygen, 1.0)?;
    product_editor.mark_tetrahedral_stereo_mixture_if_valid(carbonyl_carbon)?;
    let product_variants = expand_stereo_product_distribution(product_editor.finish()?)?;

    let mut builder = Reaction::builder(generated_pair_site_reaction_id(
        "acetal_formation",
        &carbonyl_site.participant,
        &alcohol_site.participant,
    ))
    .reactant(carbonyl.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 2, 1)
    .catalyst_order("destroy:proton", 1)
    .condition(
        ReactionCondition::new("acetal formation requires acidic, water-poor conditions")
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.35),
    );
    if product_variants.len() == 1 {
        let product = resolver.resolve(
            product_variants
                .into_iter()
                .next()
                .expect("length checked")
                .structure,
        )?;
        builder = builder.product(product, 1);
        builder = builder.product("destroy:water", 1);
    } else {
        for variant in product_variants {
            builder = builder.channel(ReactionChannel::new(
                format!("acetal_formation:stereo:{}", variant.channel_suffix),
                [
                    StoichiometricTerm::new(resolver.resolve(variant.structure)?, 1),
                    StoichiometricTerm::new("destroy:water", 1),
                ],
                25.0 + variant.activation_delta_kj_per_mol,
            ));
        }
        return Ok(builder.build());
    }
    Ok(builder.build())
}

fn generate_imine_formation(
    carbonyl_site: &CarbonylSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let carbonyl = carbonyl_site.participant.substance;
    let carbonyl_structure = carbonyl_site.participant.structure;
    let amine = amine_site.participant.substance;
    let amine_structure = amine_site.participant.structure;
    let carbonyl_carbon = carbonyl_site.carbon;
    let carbonyl_oxygen = carbonyl_site.oxygen;
    let amine_nitrogen = amine_site.nitrogen;
    let hydrogens = &amine_site.hydrogens;
    if hydrogens.len() < 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_pair_site_reaction_id(
                "imine_formation",
                &carbonyl_site.participant,
                &amine_site.participant,
            ),
            reason: "primary amine must have two explicit hydrogens for imine formation"
                .to_string(),
        });
    }

    let mut carbonyl_editor = MolecularEditor::new(carbonyl_structure);
    let carbonyl_mapping = carbonyl_editor.remove_atoms(&[carbonyl_oxygen])?;
    let carbonyl_carbon = mapped_atom(&carbonyl_mapping, carbonyl_carbon, "carbonyl carbon")?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let mut amine_editor = MolecularEditor::new(amine_structure);
    let amine_mapping = amine_editor.remove_atoms(&[hydrogens[0], hydrogens[1]])?;
    let amine_nitrogen = mapped_atom(&amine_mapping, amine_nitrogen, "amine nitrogen")?;
    let amine_fragment = amine_editor.finish()?;

    let product_structure = MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &amine_fragment,
        amine_nitrogen,
        2.0,
    )?;
    let product = resolver.resolve(product_structure)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "imine_formation",
        &carbonyl_site.participant,
        &amine_site.participant,
    ))
    .reactant(carbonyl.id.clone(), 1, 1)
    .reactant(amine.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .catalyst_order("destroy:proton", 1)
    .condition(
        ReactionCondition::new("imine formation requires acidic, water-poor conditions")
            .acidity(AcidityCondition::Acidic)
            .max_water_activity(0.5),
    )
    .build())
}

fn generate_organometallic_carbonyl_addition(
    carbonyl: SiteParticipant<'_>,
    organometallic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let (carbonyl_carbon, carbonyl_oxygen) = carbonyl_atoms_from_site(
        carbonyl.structure,
        &carbonyl.site,
        "organometallic addition",
    )?;
    let (organo_carbon, metal, residue_atoms) =
        organometallic_atoms(organometallic.structure, &organometallic.site)?;

    let mut carbonyl_editor = MolecularEditor::new(carbonyl.structure);
    carbonyl_editor.set_bond_order(carbonyl_carbon, carbonyl_oxygen, 1.0)?;
    carbonyl_editor.add_atom(carbonyl_oxygen, "H", 0.0, 1.0)?;
    let carbonyl_fragment = carbonyl_editor.finish()?;

    let mut organo_editor = MolecularEditor::new(organometallic.structure);
    let mapping = organo_editor.remove_atoms(&residue_atoms)?;
    let organo_carbon = mapped_atom(&mapping, organo_carbon, "organometallic carbon")?;
    let organo_fragment = organo_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &carbonyl_fragment,
        carbonyl_carbon,
        &organo_fragment,
        organo_carbon,
        1.0,
    )?)?;
    let residue_mass = atom_mass_sum(organometallic.structure, &residue_atoms)?;
    let residue_charge = atom_charge_sum(organometallic.structure, &residue_atoms)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "organometallic_carbonyl_addition",
        &carbonyl,
        &organometallic,
    ))
    .reactant(carbonyl.substance.id.clone(), 1, 1)
    .reactant(organometallic.substance.id.clone(), 1, 1)
    .chemical_external_reactant("proton donor hydrogen", 1.0, 1.01, 0)
    .chemical_external_product(
        format!(
            "{} salt residue",
            organometallic.structure.atoms[metal].element
        ),
        1.0,
        residue_mass,
        residue_charge,
    )
    .product(product, 1)
    .condition(
        ReactionCondition::new("organometallic carbonyl addition requires dry inert conditions")
            .max_water_activity(0.02)
            .max_oxygen_activity(0.02)
            .atmosphere(AtmosphereCondition::Inert),
    )
    .build())
}

fn generate_aldol_addition(
    enol: SiteParticipant<'_>,
    acceptor: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let (_, alpha_carbon) = enol_atoms(enol.structure, &enol.site)?;
    let alpha_hydrogen = first_bonded_hydrogen(enol.structure, alpha_carbon).ok_or_else(|| {
        ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("aldol_addition", &enol),
            reason: "aldol donor alpha carbon has no explicit hydrogen".to_string(),
        }
    })?;
    let (acceptor_carbon, acceptor_oxygen) =
        carbonyl_atoms_from_site(acceptor.structure, &acceptor.site, "aldol addition")?;

    let mut donor_editor = MolecularEditor::new(enol.structure);
    let donor_mapping = donor_editor.remove_atoms(&[alpha_hydrogen])?;
    let alpha_carbon = mapped_atom(&donor_mapping, alpha_carbon, "aldol alpha carbon")?;
    let donor_fragment = donor_editor.finish()?;

    let mut acceptor_editor = MolecularEditor::new(acceptor.structure);
    acceptor_editor.set_bond_order(acceptor_carbon, acceptor_oxygen, 1.0)?;
    acceptor_editor.add_atom(acceptor_oxygen, "H", 0.0, 1.0)?;
    let acceptor_fragment = acceptor_editor.finish()?;

    let product = resolver.resolve(MolecularEditor::join_structures(
        &donor_fragment,
        alpha_carbon,
        &acceptor_fragment,
        acceptor_carbon,
        1.0,
    )?)?;
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "aldol_addition",
        &enol,
        &acceptor,
    ))
    .reactant(enol.substance.id.clone(), 1, 1)
    .reactant(acceptor.substance.id.clone(), 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .condition(
        ReactionCondition::new("aldol addition requires basic carbonyl enolization")
            .acidity(AcidityCondition::Basic)
            .max_temperature_kelvin(323.15),
    )
    .build())
}

fn generate_aromatic_nitration(
    aromatic: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let mut variants = Vec::new();
    for carbon in aromatic_substitution_carbons(aromatic.structure, &aromatic.site) {
        let Some(hydrogen) = first_bonded_hydrogen(aromatic.structure, carbon) else {
            continue;
        };
        let mut editor = MolecularEditor::new(aromatic.structure);
        let mapping = editor.remove_atoms(&[hydrogen])?;
        let carbon = mapped_atom(&mapping, carbon, "aromatic nitration carbon")?;
        let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
        editor.add_atom(nitrogen, "O", 0.0, 1.5)?;
        editor.add_atom(nitrogen, "O", 0.0, 1.5)?;
        let product = resolver.resolve(editor.finish()?)?;
        variants.push((
            product,
            aromatic_activation_delta(aromatic.structure, carbon),
        ));
    }
    if variants.is_empty() {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("aromatic_nitration", &aromatic),
            reason: "aromatic nitration found no aromatic carbon with explicit hydrogen"
                .to_string(),
        });
    }
    let mut builder =
        Reaction::builder(generated_site_reaction_id("aromatic_nitration", &aromatic))
            .reactant(aromatic.substance.id.clone(), 1, 1)
            .reactant("destroy:nitric_acid", 1, 1)
            .catalyst_order("destroy:sulfuric_acid", 1)
            .condition(
                ReactionCondition::new("aromatic nitration requires strongly acidic conditions")
                    .acidity(AcidityCondition::Acidic)
                    .max_water_activity(0.65),
            );
    if variants.len() == 1 {
        builder = builder
            .product(variants[0].0.clone(), 1)
            .product("destroy:water", 1);
    } else {
        for (index, (product, activation_delta)) in variants.into_iter().enumerate() {
            builder = builder.channel(ReactionChannel::new(
                format!("aromatic_nitration:position_{index}"),
                [
                    StoichiometricTerm::new(product, 1),
                    StoichiometricTerm::new("destroy:water", 1),
                ],
                30.0 + activation_delta,
            ));
        }
    }
    Ok(builder.build())
}

fn generate_epoxide_hydrolysis(
    epoxide: SiteParticipant<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let oxygen = epoxide
        .site
        .atoms
        .iter()
        .copied()
        .find(|atom| epoxide.structure.atoms[*atom].element == "O")
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("epoxide_hydrolysis", &epoxide),
            reason: "epoxide site has no oxygen atom".to_string(),
        })?;
    let carbons = epoxide
        .site
        .atoms
        .iter()
        .copied()
        .filter(|atom| epoxide.structure.atoms[*atom].element == "C")
        .collect::<Vec<_>>();
    if carbons.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("epoxide_hydrolysis", &epoxide),
            reason: "epoxide site must contain exactly two carbon atoms".to_string(),
        });
    }
    let mut editor = MolecularEditor::new(epoxide.structure);
    editor.remove_bond(oxygen, carbons[0])?;
    add_hydroxyl(&mut editor, carbons[0])?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(
        Reaction::builder(generated_site_reaction_id("epoxide_hydrolysis", &epoxide))
            .reactant(epoxide.substance.id.clone(), 1, 1)
            .reactant("destroy:water", 1, 1)
            .catalyst_order("destroy:proton", 1)
            .product(product, 1)
            .condition(
                ReactionCondition::new("epoxide hydrolysis requires aqueous acid")
                    .acidity(AcidityCondition::Acidic)
                    .min_water_activity(0.1),
            )
            .build(),
    )
}

fn deprotonated_alcohol_fragment(
    site: &AlcoholSite<'_>,
    _role: &str,
) -> ChemistryResult<(MolecularStructure, usize)> {
    let structure = site.participant.structure;
    let oxygen = site.oxygen;
    let hydrogen = site.hydrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let oxygen = mapped_atom(&mapping, oxygen, "alcohol oxygen")?;
    Ok((editor.finish()?, oxygen))
}

fn generate_nitrile_hydrolysis(
    site: &NitrileSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let nitrogen = site.nitrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "nitrile carbon")?;
    editor.add_atom(carbon, "O", 0.0, 2.0)?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "nitrile_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(product, 1)
    .build())
}

fn generate_nitro_hydrogenation(
    site: &NitroSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let [first_oxygen, second_oxygen] = site.oxygens;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[first_oxygen, second_oxygen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "nitro nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "nitro_hydrogenation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen", 3, 1)
    .product("destroy:water", 2)
    .product(product, 1)
    .external_catalyst("forge:dusts/palladium", 1.0)
    .build())
}

fn generate_acyl_chloride_formation(
    site: &CarboxylicAcidSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let hydroxyl_oxygen = site.hydroxyl_oxygen;
    let proton = site.hydroxyl_hydrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydroxyl_oxygen, proton])?;
    let carbon = mapped_atom(&mapping, carbon, "carboxylic acid carbon")?;
    editor.add_atom(carbon, "Cl", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "acyl_chloride_formation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:phosgene", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .product("destroy:carbon_dioxide", 1)
    .build())
}

fn generate_acyl_chloride_hydrolysis(
    site: &AcylChlorideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let chlorine = site.chlorine;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[chlorine])?;
    let carbon = mapped_atom(&mapping, carbon, "acyl chloride carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "acyl_chloride_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
}

fn generate_acyl_chloride_esterification(
    acyl_chloride_site: &AcylChlorideSite<'_>,
    alcohol_site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let acyl_chloride = acyl_chloride_site.participant.substance;
    let acyl_chloride_structure = acyl_chloride_site.participant.structure;
    let alcohol = alcohol_site.participant.substance;
    let alcohol_structure = alcohol_site.participant.structure;
    let acyl_carbon = acyl_chloride_site.carbon;
    let chlorine = acyl_chloride_site.chlorine;
    let alcohol_oxygen = alcohol_site.oxygen;
    let alcohol_proton = alcohol_site.hydrogen;
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
    Ok(Reaction::builder(generated_pair_site_reaction_id(
        "acyl_chloride_esterification",
        &acyl_chloride_site.participant,
        &alcohol_site.participant,
    ))
    .reactant(acyl_chloride.id.clone(), 1, 1)
    .reactant(alcohol.id.clone(), 1, 1)
    .product(product, 1)
    .product("destroy:hydrochloric_acid", 1)
    .build())
}

fn generate_alcohol_dehydration(
    site: &AlcoholSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let alcohol_carbon = site.carbon;
    let oxygen = site.oxygen;
    let proton = site.hydrogen;
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
    let mut builder = Reaction::builder(generated_site_reaction_id(
        "alcohol_dehydration",
        &site.participant,
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

fn generate_thionyl_chloride_substitution(
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

fn generate_aldehyde_oxidation(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Option<Reaction>> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    if site.is_ketone {
        return Ok(None);
    }
    let carbon = site.carbon;
    let Some(hydrogen) = first_bonded_hydrogen(structure, carbon) else {
        return Ok(None);
    };
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "aldehyde carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Some(
        Reaction::builder(generated_site_reaction_id(
            "aldehyde_oxidation",
            &site.participant,
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
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let oxygen = site.oxygen;
    let mut editor = MolecularEditor::new(structure);
    editor.set_bond_order(carbon, oxygen, 1.0)?;
    editor.add_atom(oxygen, "H", 0.0, 1.0)?;
    let nitrile_carbon = editor.add_atom(carbon, "C", 0.0, 1.0)?;
    editor.add_atom(nitrile_carbon, "N", 0.0, 3.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "cyanide_nucleophilic_addition",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen_cyanide", 1, 1)
    .catalyst_order("destroy:cyanide", 1)
    .product(product, 1)
    .build())
}

fn generate_wolff_kishner_reduction(
    site: &CarbonylSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let oxygen = site.oxygen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[oxygen])?;
    let carbon = mapped_atom(&mapping, carbon, "carbonyl carbon")?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "wolff_kishner_reduction",
        &site.participant,
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
    site: &AmideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let nitrogen = site.nitrogen;
    let hydrogens = &site.nitrogen_hydrogens;
    if hydrogens.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("amide_hydrolysis", &site.participant),
            reason: "unsubstituted amide must have exactly two explicit nitrogen hydrogens"
                .to_string(),
        });
    }
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen, hydrogens[0], hydrogens[1]])?;
    let carbon = mapped_atom(&mapping, carbon, "amide carbon")?;
    add_hydroxyl(&mut editor, carbon)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "amide_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .catalyst_order("destroy:proton", 1)
    .product(product, 1)
    .product("destroy:ammonia", 1)
    .build())
}

fn generate_amine_phosgenation(
    site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let hydrogens = &site.hydrogens;
    if hydrogens.len() != 2 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("amine_phosgenation", &site.participant),
            reason: "primary amine must have exactly two explicit hydrogens".to_string(),
        });
    }
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[hydrogens[0], hydrogens[1]])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "amine nitrogen")?;
    let carbon = editor.add_atom(nitrogen, "C", 0.0, 2.0)?;
    editor.add_atom(carbon, "O", 0.0, 2.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "amine_phosgenation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:phosgene", 1, 1)
    .product("destroy:hydrochloric_acid", 2)
    .product(product, 1)
    .build())
}

fn generate_cyanamide_addition(
    site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let hydrogen = *site
        .hydrogens
        .first()
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id("cyanamide_addition", &site.participant),
            reason: "amine has no explicit hydrogen".to_string(),
        })?;
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
    Ok(Reaction::builder(generated_site_reaction_id(
        "cyanamide_addition",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:cyanamide", 1, 1)
    .product(product, 1)
    .build())
}

fn generate_halide_ammonia_substitution(
    site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
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
    Ok(Reaction::builder(generated_site_reaction_id(
        "halide_ammonia_substitution",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:ammonia", 2, if site.degree == 3 { 1 } else { 2 })
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
    .build())
}

fn generate_halide_cyanide_substitution(
    site: &HalideSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
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
    Ok(Reaction::builder(generated_site_reaction_id(
        "halide_cyanide_substitution",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:cyanide", 1, if site.degree == 3 { 0 } else { 1 })
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
    .build())
}

fn generate_halide_amine_substitution(
    halide_site: &HalideSite<'_>,
    amine_site: &AmineSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
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
    Ok(Reaction::builder(generated_pair_site_reaction_id(
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
    .build())
}

fn generate_isocyanate_hydrolysis(
    site: &IsocyanateSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let nitrogen = site.nitrogen;
    let functional_carbon = site.functional_carbon;
    let oxygen = site.oxygen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[functional_carbon, oxygen])?;
    let nitrogen = mapped_atom(&mapping, nitrogen, "isocyanate nitrogen")?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "isocyanate_hydrolysis",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:water", 1, 1)
    .product("destroy:carbon_dioxide", 1)
    .product(product, 1)
    .build())
}

fn generate_nitrile_hydrogenation(
    site: &NitrileSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let nitrogen = site.nitrogen;
    let mut editor = MolecularEditor::new(structure);
    let mapping = editor.remove_atoms(&[nitrogen])?;
    let carbon = mapped_atom(&mapping, carbon, "nitrile carbon")?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    editor.add_atom(carbon, "H", 0.0, 1.0)?;
    let nitrogen = editor.add_atom(carbon, "N", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    editor.add_atom(nitrogen, "H", 0.0, 1.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "nitrile_hydrogenation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen", 2, 1)
    .product(product, 1)
    .external_catalyst("forge:dusts/nickel", 1.0)
    .build())
}

fn generate_borane_oxidation(
    site: &BoraneSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let carbon = site.carbon;
    let boron = site.boron;
    let mut editor = MolecularEditor::new(structure);
    editor.insert_bridging_atom(carbon, boron, "O", 0.0)?;
    let product = resolver.resolve(editor.finish()?)?;
    Ok(Reaction::builder(generated_site_reaction_id(
        "borane_oxidation",
        &site.participant,
    ))
    .reactant(substance.id.clone(), 1, 1)
    .reactant("destroy:hydrogen_peroxide", 1, 1)
    .catalyst_order("destroy:hydroxide", 1)
    .product(product, 1)
    .product("destroy:water", 1)
    .build())
}

fn generate_borate_ester_hydrolysis(
    site: &BorateEsterSite<'_>,
    resolver: &mut DerivedSubstanceResolver,
) -> ChemistryResult<Reaction> {
    let substance = site.participant.substance;
    let structure = site.participant.structure;
    let oxygen = site.oxygen;
    let boron = site.boron;
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

    Ok(Reaction::builder(generated_site_reaction_id(
        "borate_ester_hydrolysis",
        &site.participant,
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
    alkyne_stereo_rule: Option<AlkyneStereoRule>,
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

#[derive(Debug, Clone, Copy)]
enum AlkyneStereoRule {
    Anti,
}

#[derive(Debug, Clone)]
struct StereoProductVariant {
    structure: MolecularStructure,
    channel_suffix: String,
    activation_delta_kj_per_mol: f64,
    pre_exponential_factor_multiplier: f64,
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

fn generate_electrophilic_addition(
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

fn apply_alkyne_stereo_rule(
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

fn double_bond_stereo_substituents(
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

fn bonded_substituents_except(
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

fn preferred_stereo_substituent(
    structure: &MolecularStructure,
    substituents: &[usize],
) -> Option<usize> {
    substituents.iter().copied().max_by_key(|index| {
        let atom = &structure.atoms[*index];
        (atomic_stereo_priority(&atom.element), atom.r_group_number)
    })
}

fn atomic_stereo_priority(element: &str) -> u16 {
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

fn expand_stereo_product_distribution(
    structure: MolecularStructure,
) -> ChemistryResult<Vec<StereoProductVariant>> {
    expand_stereo_product_distribution_with_parameters(structure, "single".to_string(), 0.0, 1.0)
}

fn expand_stereo_product_distribution_with_parameters(
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
                    Stereochemistry::DoubleBond(super::molecule::DoubleBondStereo {
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
                    Stereochemistry::DoubleBond(super::molecule::DoubleBondStereo {
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

fn geometric_z_steric_penalty_kj_per_mol(
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

fn z_pre_exponential_multiplier(steric_penalty_kj_per_mol: f64) -> f64 {
    (1.0 - steric_penalty_kj_per_mol / 20.0).clamp(0.55, 0.95)
}

fn substituent_steric_bulk(
    structure: &MolecularStructure,
    substituent: usize,
    blocked_atom: usize,
) -> f64 {
    let mut visited = BTreeSet::new();
    substituent_steric_bulk_inner(structure, substituent, blocked_atom, &mut visited)
}

fn substituent_steric_bulk_inner(
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

fn bonded_hydrogens(structure: &MolecularStructure, parent: usize) -> Vec<usize> {
    structure
        .neighbors(parent)
        .into_iter()
        .map(|(neighbor, _)| neighbor)
        .filter(|neighbor| structure.atoms[*neighbor].element == "H")
        .collect()
}

fn halide_ion(
    structure: &MolecularStructure,
    halogen: usize,
    prefix: &str,
    participant: &SiteParticipant<'_>,
) -> ChemistryResult<&'static str> {
    match structure.atoms[halogen].element.as_str() {
        "Cl" => Ok("destroy:chloride"),
        "F" => Ok("destroy:fluoride"),
        "I" => Ok("destroy:iodide"),
        _ => Err(ChemistryError::InvalidReaction {
            reaction_id: generated_site_reaction_id(prefix, participant),
            reason: "halide group does not contain a supported halogen".to_string(),
        }),
    }
}

fn carbonyl_atoms_from_site(
    structure: &MolecularStructure,
    site: &ReactiveSite,
    role: &str,
) -> ChemistryResult<(usize, usize)> {
    for carbon in site
        .atoms
        .iter()
        .copied()
        .filter(|atom| structure.atoms[*atom].element == "C")
    {
        if let Some((oxygen, _)) =
            structure
                .neighbors(carbon)
                .into_iter()
                .find(|(neighbor, order)| {
                    structure.atoms[*neighbor].element == "O"
                        && super::molecule::bond_order_matches(*order, 2.0)
                })
        {
            return Ok((carbon, oxygen));
        }
    }
    Err(ChemistryError::InvalidReaction {
        reaction_id: "<generated-organic>".to_string(),
        reason: format!("{role} site does not contain a carbonyl bond"),
    })
}

fn enol_atoms(
    structure: &MolecularStructure,
    site: &ReactiveSite,
) -> ChemistryResult<(usize, usize)> {
    let (carbonyl, _) = carbonyl_atoms_from_site(structure, site, "enol")?;
    let alpha = site
        .atoms
        .iter()
        .copied()
        .find(|atom| {
            *atom != carbonyl
                && structure.atoms[*atom].element == "C"
                && structure
                    .neighbors(carbonyl)
                    .iter()
                    .any(|(neighbor, order)| {
                        *neighbor == *atom && super::molecule::bond_order_matches(*order, 1.0)
                    })
        })
        .ok_or_else(|| ChemistryError::InvalidReaction {
            reaction_id: "<generated-organic>".to_string(),
            reason: "enol site does not contain an alpha carbon".to_string(),
        })?;
    Ok((carbonyl, alpha))
}

fn organometallic_atoms(
    structure: &MolecularStructure,
    site: &ReactiveSite,
) -> ChemistryResult<(usize, usize, Vec<usize>)> {
    let mut organo_carbon = None;
    let mut metal = None;
    for atom in &site.atoms {
        match structure.atoms[*atom].element.as_str() {
            "C" => organo_carbon = Some(*atom),
            "Mg" | "Li" | "Cu" => metal = Some(*atom),
            _ => {}
        }
    }
    let organo_carbon = organo_carbon.ok_or_else(|| ChemistryError::InvalidReaction {
        reaction_id: "<generated-organic>".to_string(),
        reason: "organometallic site has no carbon atom".to_string(),
    })?;
    let metal = metal.ok_or_else(|| ChemistryError::InvalidReaction {
        reaction_id: "<generated-organic>".to_string(),
        reason: "organometallic site has no metal atom".to_string(),
    })?;
    let mut residue_atoms = vec![metal];
    for (neighbor, order) in structure.neighbors(metal) {
        if neighbor != organo_carbon && super::molecule::bond_order_matches(order, 1.0) {
            residue_atoms.push(neighbor);
        }
    }
    residue_atoms.sort_unstable();
    residue_atoms.dedup();
    Ok((organo_carbon, metal, residue_atoms))
}

fn atom_mass_sum(structure: &MolecularStructure, atoms: &[usize]) -> ChemistryResult<f64> {
    atoms.iter().try_fold(0.0, |sum, atom| {
        Ok(sum + super::molecule::element_mass(&structure.atoms[*atom].element)?)
    })
}

fn atom_charge_sum(structure: &MolecularStructure, atoms: &[usize]) -> ChemistryResult<i32> {
    let charge = atoms
        .iter()
        .map(|atom| structure.atoms[*atom].charge)
        .sum::<f64>();
    if (charge - charge.round()).abs() > 1.0e-9 {
        return Err(ChemistryError::InvalidReaction {
            reaction_id: "<generated-organic>".to_string(),
            reason: "external residue has non-integral charge".to_string(),
        });
    }
    Ok(charge.round() as i32)
}

fn aromatic_substitution_carbons(
    structure: &MolecularStructure,
    site: &ReactiveSite,
) -> Vec<usize> {
    site.atoms
        .iter()
        .copied()
        .filter(|atom| {
            structure.atoms[*atom].element == "C"
                && first_bonded_hydrogen(structure, *atom).is_some()
                && structure
                    .neighbors(*atom)
                    .iter()
                    .filter(|(_, order)| super::molecule::bond_order_matches(*order, 1.5))
                    .count()
                    >= 2
        })
        .collect()
}

fn aromatic_activation_delta(structure: &MolecularStructure, carbon: usize) -> f64 {
    let mut delta: f64 = 10.0;
    for (neighbor, order) in structure.neighbors(carbon) {
        if !super::molecule::bond_order_matches(order, 1.5) {
            continue;
        }
        for (substituent, substituent_order) in structure.neighbors(neighbor) {
            if substituent == carbon || super::molecule::bond_order_matches(substituent_order, 1.5)
            {
                continue;
            }
            delta = delta.min(match structure.atoms[substituent].element.as_str() {
                "O" | "N" => 0.0,
                "C" => 2.5,
                "Cl" | "Br" | "I" | "F" => 4.0,
                _ => 8.0,
            });
        }
    }
    delta
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

fn generated_pair_reaction_id(prefix: &str, first: &Substance, second: &Substance) -> String {
    format!(
        "{prefix}/{}/{}",
        sanitize_id(first.id.as_str()),
        sanitize_id(second.id.as_str())
    )
}

fn generated_site_reaction_id(prefix: &str, participant: &SiteParticipant<'_>) -> String {
    format!(
        "{}/{}/{}",
        prefix,
        sanitize_id(participant.substance.id.as_str()),
        participant
            .site
            .atoms
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("_")
    )
}

fn generated_pair_site_reaction_id(
    prefix: &str,
    first: &SiteParticipant<'_>,
    second: &SiteParticipant<'_>,
) -> String {
    format!(
        "{}/{}/{}/{}",
        generated_pair_reaction_id(prefix, first.substance, second.substance),
        first
            .site
            .atoms
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("_"),
        second
            .site
            .atoms
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("_"),
        site_kind_suffix(&first.site.kind)
    )
}

fn site_kind_suffix(kind: &ReactiveSiteKind) -> &'static str {
    match kind {
        ReactiveSiteKind::AcidAnhydride => "acid_anhydride",
        ReactiveSiteKind::AcylChloride => "acyl_chloride",
        ReactiveSiteKind::Alcohol => "alcohol",
        ReactiveSiteKind::Alkene => "alkene",
        ReactiveSiteKind::Alkoxide => "alkoxide",
        ReactiveSiteKind::Alkyne => "alkyne",
        ReactiveSiteKind::Aldehyde => "aldehyde",
        ReactiveSiteKind::Amide => "amide",
        ReactiveSiteKind::AromaticCarbon => "aromatic_carbon",
        ReactiveSiteKind::AromaticRing => "aromatic_ring",
        ReactiveSiteKind::ArylHalide => "aryl_halide",
        ReactiveSiteKind::Azide => "azide",
        ReactiveSiteKind::Borane => "borane",
        ReactiveSiteKind::BoricAcid => "boric_acid",
        ReactiveSiteKind::BorateEster => "borate_ester",
        ReactiveSiteKind::Carbonyl => "carbonyl",
        ReactiveSiteKind::CarboxylicAcid => "carboxylic_acid",
        ReactiveSiteKind::Diazonium => "diazonium",
        ReactiveSiteKind::Enol => "enol",
        ReactiveSiteKind::Enolate => "enolate",
        ReactiveSiteKind::Epoxide => "epoxide",
        ReactiveSiteKind::Ester => "ester",
        ReactiveSiteKind::Ether => "ether",
        ReactiveSiteKind::Halide => "halide",
        ReactiveSiteKind::Imine => "imine",
        ReactiveSiteKind::Isocyanate => "isocyanate",
        ReactiveSiteKind::Ketone => "ketone",
        ReactiveSiteKind::Nitrile => "nitrile",
        ReactiveSiteKind::Nitro => "nitro",
        ReactiveSiteKind::NonTertiaryAmine => "non_tertiary_amine",
        ReactiveSiteKind::Organocopper => "organocopper",
        ReactiveSiteKind::Organolithium => "organolithium",
        ReactiveSiteKind::Organomagnesium => "organomagnesium",
        ReactiveSiteKind::Phenol => "phenol",
        ReactiveSiteKind::PrimaryAmine => "primary_amine",
        ReactiveSiteKind::Sulfide => "sulfide",
        ReactiveSiteKind::SulfonylChloride => "sulfonyl_chloride",
        ReactiveSiteKind::Thiol => "thiol",
    }
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
    const ACTIVE_GENERATORS_WITH_UNKNOWN_STEREO_DISTRIBUTION: &[&str] = &[];

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
            .sites_of(&ReactiveSiteKind::CarboxylicAcid)
            .collect::<Vec<_>>();
        assert_eq!(acids.len(), 1);
        assert_eq!(acids[0].substance.id.as_str(), "destroy:acetic_acid");
        assert_eq!(space.sites_of(&ReactiveSiteKind::Alcohol).count(), 0);
    }

    #[test]
    fn organic_generation_has_no_functional_group_transition_layer() {
        let source = include_str!("mod.rs");
        assert!(!source.contains(concat!("legacy", "_group", "_from", "_site")));
        assert!(!source.contains(concat!("sites", "_of", "_legacy", "_group")));
        assert!(!source.contains(concat!("Functional", "Group")));
    }

    #[test]
    fn acetal_and_imine_generators_create_concrete_products_with_conditions() {
        let registry = generated_registry();
        let acetal = reaction_with_prefix(&registry, "acetal_formation/");
        assert!(acetal
            .conditions
            .iter()
            .any(|condition| condition.acidity == Some(AcidityCondition::Acidic)));
        assert!(acetal
            .products
            .iter()
            .chain(
                acetal
                    .channels
                    .iter()
                    .flat_map(|channel| channel.products.iter())
            )
            .any(|term| term.substance_id.as_str() == "destroy:water"));

        let imine = reaction_with_prefix(&registry, "imine_formation/");
        assert!(imine
            .conditions
            .iter()
            .any(|condition| condition.max_water_activity.is_some()));
        assert!(imine.products.len() >= 2);
    }

    #[test]
    fn reactive_site_generators_add_aromatic_nitration_and_epoxide_hydrolysis() {
        let registry = generated_registry();
        let nitration = reaction_with_prefix(&registry, "aromatic_nitration/destroy_benzene/");
        assert!(nitration
            .conditions
            .iter()
            .any(|condition| condition.acidity == Some(AcidityCondition::Acidic)));
        assert!(!nitration.channels.is_empty() || !nitration.products.is_empty());

        let mut dynamic =
            super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let epoxide = dynamic
            .resolve_frowns(
                "destroy:graph:atoms=C.C.O.H.H.H.H;bonds=0-s-1,0-s-2,1-s-2,0-s-3,0-s-4,1-s-5,1-s-6",
            )
            .unwrap();
        let report = dynamic.generate_reactions_for(&epoxide, 1).unwrap();
        assert!(report.added_reactions > 0);
        assert!(dynamic
            .reactions()
            .any(|reaction| reaction.id.as_str().starts_with("epoxide_hydrolysis/")));
    }

    #[test]
    fn organometallic_and_aldol_generators_create_carbon_carbon_bonds() {
        let mut dynamic =
            super::super::dynamic::DynamicChemistryRegistry::from_destroy_catalog().unwrap();
        let methyl_magnesium_chloride = dynamic.resolve_frowns("CMgCl").unwrap();
        dynamic
            .generate_reactions_for_substances(
                [
                    SubstanceId::from("destroy:acetone"),
                    methyl_magnesium_chloride,
                ],
                1,
            )
            .unwrap();
        let organometallic = dynamic
            .reactions()
            .find(|reaction| {
                reaction
                    .id
                    .as_str()
                    .starts_with("organometallic_carbonyl_addition/")
            })
            .unwrap();
        assert!(organometallic
            .conditions
            .iter()
            .any(|condition| condition.atmosphere == Some(AtmosphereCondition::Inert)));
        assert!(!organometallic.external_products.is_empty());

        let acetaldehyde = dynamic.resolve_frowns("CC=O").unwrap();
        dynamic
            .generate_reactions_for_substances(
                [SubstanceId::from("destroy:acetone"), acetaldehyde],
                1,
            )
            .unwrap();
        assert!(dynamic
            .reactions()
            .any(|reaction| reaction.id.as_str().starts_with("aldol_addition/")));
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
            if ACTIVE_GENERATORS_WITH_UNKNOWN_STEREO_DISTRIBUTION.contains(prefix) {
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
            "alkyne_hydrogenation/destroy_acetylene/",
        ] {
            reaction_with_prefix(&registry, prefix);
        }
        let hydrogenation = reaction_with_prefix(&registry, "alkene_hydrogenation/destroy_ethene/");
        assert!(!hydrogenation.external_catalysts.is_empty());
    }

    #[test]
    fn geometric_stereo_products_use_kinetic_channels() {
        let structure = super::super::frowns::parse_frowns(
            "destroy:graph:atoms=C.C.H.Cl.H.I;bonds=0-2-1,0-s-2,0-s-3,1-s-4,1-s-5;stereo=mix:db:0.1.3.5",
        )
        .unwrap();
        let variants = expand_stereo_product_distribution(structure).unwrap();
        let e_variant = variants
            .iter()
            .find(|variant| variant.channel_suffix.contains("db_e"))
            .unwrap();
        let z_variant = variants
            .iter()
            .find(|variant| variant.channel_suffix.contains("db_z"))
            .unwrap();

        assert!(z_variant.activation_delta_kj_per_mol > e_variant.activation_delta_kj_per_mol);
        assert!(
            z_variant.pre_exponential_factor_multiplier
                < e_variant.pre_exponential_factor_multiplier
        );
        assert!(e_variant
            .structure
            .stereochemistry
            .iter()
            .any(|stereo| matches!(stereo, Stereochemistry::DoubleBond(double_bond) if double_bond.descriptor == StereoDescriptor::E)));
        assert!(z_variant
            .structure
            .stereochemistry
            .iter()
            .any(|stereo| matches!(stereo, Stereochemistry::DoubleBond(double_bond) if double_bond.descriptor == StereoDescriptor::Z)));
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
