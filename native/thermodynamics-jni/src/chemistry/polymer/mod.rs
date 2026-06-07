use std::fmt::{Display, Formatter};

use super::error::{ChemistryError, ChemistryResult};
use super::functional_group::{find_functional_groups, FunctionalGroupType};
use super::molecule::{
    bond_order_matches, MolecularAtom, MolecularBond, MolecularEditor, MolecularStructure,
    ValenceSaturation,
};
use super::substance::{
    LiquidPhasePreference, SolventRole, Substance, SubstanceId, SubstancePhaseBehavior,
    SubstanceRepresentation,
};

const MIN_POLYMERIZATION_CONVERSION: f64 = 1.0e-9;
const MAX_REASONABLE_DISPERSITY: f64 = 100.0;
const REPRESENTATIVE_CHAIN_GROWTH_DEGREE: u32 = 100;
const REPRESENTATIVE_STEP_GROWTH_DEGREE: u32 = 20;

/// Bulk macromolecules do not vaporize — they thermally degrade first. An
/// effectively unreachable boiling point marks them as non-volatile, matching
/// the registry's established `f64::MAX` "never boils" convention (used for
/// ions like silver/chloride) while staying finite so thermodynamic math is
/// well-defined.
const POLYMER_BOILING_POINT_KELVIN: f64 = f64::MAX;
const DEFAULT_POLYMER_DENSITY_GRAMS_PER_BUCKET: f64 = 1000.0;
const DEFAULT_POLYMER_HEAT_CAPACITY_J_PER_MOL_KELVIN: f64 = 100.0;
const DEFAULT_POLYMER_LATENT_HEAT_J_PER_MOL: f64 = 20_000.0;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PolymerId(String);

impl PolymerId {
    pub fn new(value: impl Into<String>) -> ChemistryResult<Self> {
        let value = value.into();
        if value.trim().is_empty() {
            return Err(polymer_error("polymer id must not be empty"));
        }
        Ok(Self(value))
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl Display for PolymerId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl From<&str> for PolymerId {
    fn from(value: &str) -> Self {
        Self(value.to_string())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolymerizationMechanism {
    ChainGrowthRadical,
    StepGrowthCondensation,
    RingOpening,
    CoordinationInsertion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PolymerArchitecture {
    Linear,
    Branched {
        branch_points_per_chain: u32,
    },
    Network {
        crosslink_mol_fraction: OrderedFinite,
    },
}

impl PolymerArchitecture {
    fn id_token(&self) -> String {
        match self {
            Self::Linear => "linear".to_string(),
            Self::Branched {
                branch_points_per_chain,
            } => format!("branched:{branch_points_per_chain}"),
            Self::Network {
                crosslink_mol_fraction,
            } => format!("network:{:.6}", crosslink_mol_fraction.value()),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct OrderedFinite(u64);

impl OrderedFinite {
    pub fn new(value: f64) -> ChemistryResult<Self> {
        if !value.is_finite() || value < 0.0 {
            return Err(polymer_error("finite non-negative value required"));
        }
        Ok(Self(value.to_bits()))
    }

    pub fn value(self) -> f64 {
        f64::from_bits(self.0)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PolymerEndGroup {
    pub label: String,
    pub molar_mass_grams: OrderedFinite,
}

impl PolymerEndGroup {
    pub fn new(label: impl Into<String>, molar_mass_grams: f64) -> ChemistryResult<Self> {
        let label = label.into();
        if label.trim().is_empty() {
            return Err(polymer_error("polymer end group label must not be empty"));
        }
        Ok(Self {
            label,
            molar_mass_grams: OrderedFinite::new(molar_mass_grams)?,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolymerRepeatUnit {
    pub source_monomer: SubstanceId,
    pub structure: MolecularStructure,
    pub connection_atoms: [usize; 2],
    pub molar_mass_grams: f64,
}

impl PolymerRepeatUnit {
    pub fn new(
        source_monomer: impl Into<SubstanceId>,
        structure: MolecularStructure,
        connection_atoms: [usize; 2],
    ) -> ChemistryResult<Self> {
        validate_repeat_unit_structure(&structure, connection_atoms)?;
        let summary = structure.summary()?;
        if summary.charge != 0 {
            return Err(polymer_error("polymer repeat unit must be neutral"));
        }
        Ok(Self {
            source_monomer: source_monomer.into(),
            structure,
            connection_atoms,
            molar_mass_grams: summary.molar_mass_grams,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChainLengthDistribution {
    pub number_average_degree: f64,
    pub weight_average_degree: f64,
    pub dispersity: f64,
}

impl ChainLengthDistribution {
    pub fn new(number_average_degree: f64, weight_average_degree: f64) -> ChemistryResult<Self> {
        validate_positive_finite(number_average_degree, "number-average polymerization degree")?;
        validate_positive_finite(weight_average_degree, "weight-average polymerization degree")?;
        if weight_average_degree < number_average_degree {
            return Err(polymer_error(
                "weight-average polymerization degree must not be below number-average degree",
            ));
        }
        let dispersity = weight_average_degree / number_average_degree;
        validate_positive_finite(dispersity, "polymer dispersity")?;
        if dispersity > MAX_REASONABLE_DISPERSITY {
            return Err(polymer_error("polymer dispersity is outside supported range"));
        }
        Ok(Self {
            number_average_degree,
            weight_average_degree,
            dispersity,
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolymerMaterial {
    pub id: PolymerId,
    pub mechanism: PolymerizationMechanism,
    pub repeat_unit: PolymerRepeatUnit,
    pub distribution: ChainLengthDistribution,
    pub architecture: PolymerArchitecture,
    pub end_groups: Vec<PolymerEndGroup>,
}

impl PolymerMaterial {
    pub fn validate(&self) -> ChemistryResult<()> {
        if self.id.as_str().trim().is_empty() {
            return Err(polymer_error("polymer id must not be empty"));
        }
        validate_repeat_unit_structure(&self.repeat_unit.structure, self.repeat_unit.connection_atoms)?;
        ChainLengthDistribution::new(
            self.distribution.number_average_degree,
            self.distribution.weight_average_degree,
        )?;
        match &self.architecture {
            PolymerArchitecture::Linear => {}
            PolymerArchitecture::Branched {
                branch_points_per_chain,
            } => {
                if *branch_points_per_chain == 0 {
                    return Err(polymer_error(
                        "branched polymer architecture requires at least one branch point",
                    ));
                }
            }
            PolymerArchitecture::Network {
                crosslink_mol_fraction,
            } => {
                let fraction = crosslink_mol_fraction.value();
                if !(0.0..=1.0).contains(&fraction) {
                    return Err(polymer_error(
                        "network crosslink fraction must be between 0 and 1",
                    ));
                }
            }
        }
        Ok(())
    }

    pub fn number_average_molar_mass_grams(&self) -> f64 {
        self.repeat_unit.molar_mass_grams * self.distribution.number_average_degree
            + self.end_group_molar_mass_grams()
    }

    pub fn weight_average_molar_mass_grams(&self) -> f64 {
        self.repeat_unit.molar_mass_grams * self.distribution.weight_average_degree
            + self.end_group_molar_mass_grams()
    }

    fn end_group_molar_mass_grams(&self) -> f64 {
        self.end_groups
            .iter()
            .map(|group| group.molar_mass_grams.value())
            .sum()
    }
}

/// Bridges a derived [`PolymerMaterial`] into a registry [`Substance`] so the
/// macromolecule can appear in a mixture. The repeat unit is stored as catalog
/// metadata, not as `molecular_structure`: a polymer is not a small molecule and
/// ordinary organic generators must not react with its repeat unit as if it were
/// an isolated substance.
pub fn polymer_material_to_substance(material: &PolymerMaterial) -> ChemistryResult<Substance> {
    material.validate()?;
    let id = SubstanceId::new(material.id.as_str())?;
    let repeat_frowns = super::frowns::write_frowns(&material.repeat_unit.structure)?;
    let substance = Substance::new(
        id,
        0,
        material.number_average_molar_mass_grams(),
        DEFAULT_POLYMER_DENSITY_GRAMS_PER_BUCKET,
        POLYMER_BOILING_POINT_KELVIN,
        DEFAULT_POLYMER_HEAT_CAPACITY_J_PER_MOL_KELVIN,
        DEFAULT_POLYMER_LATENT_HEAT_J_PER_MOL,
    )
    .with_catalog_metadata(Some(repeat_frowns), None, 0x20FF_FFFF, Vec::new())
    .with_phase_properties(SubstancePhaseBehavior {
        preferred_liquid_phase: LiquidPhasePreference::Organic,
        aqueous_solubility_mol_per_bucket: Some(0.0),
        organic_solubility_mol_per_bucket: Some(0.0),
        can_precipitate: true,
        can_form_liquid_phase: false,
        solvent_role: SolventRole::NotSolvent,
    })
    .with_representation(SubstanceRepresentation::Polymer {
        repeat_unit_id: material.repeat_unit.source_monomer.clone(),
        number_average_degree: material.distribution.number_average_degree,
        weight_average_degree: material.distribution.weight_average_degree,
        architecture: material.architecture.id_token(),
    });
    Ok(substance)
}

#[derive(Debug, Clone, PartialEq)]
pub struct ChainGrowthInput {
    pub monomer_id: SubstanceId,
    pub monomer_structure: MolecularStructure,
    pub conversion_fraction: f64,
    pub initiator_to_monomer_mol_fraction: f64,
    pub chain_transfer_fraction: f64,
}

impl ChainGrowthInput {
    pub fn new(
        monomer_id: impl Into<SubstanceId>,
        monomer_structure: MolecularStructure,
        conversion_fraction: f64,
        initiator_to_monomer_mol_fraction: f64,
        chain_transfer_fraction: f64,
    ) -> Self {
        Self {
            monomer_id: monomer_id.into(),
            monomer_structure,
            conversion_fraction,
            initiator_to_monomer_mol_fraction,
            chain_transfer_fraction,
        }
    }
}

pub fn derive_chain_growth_polymer(input: ChainGrowthInput) -> ChemistryResult<Option<PolymerMaterial>> {
    validate_fraction(input.conversion_fraction, "polymerization conversion")?;
    validate_fraction(
        input.initiator_to_monomer_mol_fraction,
        "initiator-to-monomer fraction",
    )?;
    validate_fraction(input.chain_transfer_fraction, "chain transfer fraction")?;
    if input.conversion_fraction < MIN_POLYMERIZATION_CONVERSION {
        return Ok(None);
    }
    let Some((first, second)) = polymerizable_carbon_carbon_double_bond(&input.monomer_structure)
    else {
        return Ok(None);
    };
    let repeat_unit = alkene_repeat_unit(
        input.monomer_id.clone(),
        &input.monomer_structure,
        first,
        second,
    )?;
    let active_chain_fraction =
        input.initiator_to_monomer_mol_fraction + input.chain_transfer_fraction;
    if active_chain_fraction <= 0.0 {
        return Err(polymer_error(
            "chain-growth polymerization requires initiator or chain transfer source",
        ));
    }
    let number_average_degree = input.conversion_fraction / active_chain_fraction;
    let dispersity = 1.0 + input.conversion_fraction;
    let distribution =
        ChainLengthDistribution::new(number_average_degree, number_average_degree * dispersity)?;
    let material = PolymerMaterial {
        id: PolymerId::new(format!(
            "polymer:chain_growth:{}:{}-{}",
            sanitize_polymer_id(input.monomer_id.as_str()),
            first,
            second
        ))?,
        mechanism: PolymerizationMechanism::ChainGrowthRadical,
        repeat_unit,
        distribution,
        architecture: PolymerArchitecture::Linear,
        end_groups: vec![
            PolymerEndGroup::new("initiator-derived", 0.0)?,
            PolymerEndGroup::new("terminated-chain-end", 1.01)?,
        ],
    };
    material.validate()?;
    Ok(Some(material))
}

/// Builds the *repeat-unit structure* for a chain-growth (addition) polymer
/// directly from a monomer, or `None` if the monomer cannot polymerize this way.
///
/// This is the simulator-facing counterpart to [`derive_chain_growth_polymer`]:
/// where that returns the fractional [`PolymerMaterial`] distribution for the
/// planner/analytics, this returns the single discrete repeat unit (the alkene
/// C=C collapsed to C–C with both carbons marked `UnsaturatedAllowed` to carry
/// the chain bonds). Resolving that structure to a substance yields a mass-exact
/// `1 monomer → 1 repeat-unit` reaction edge: no atoms are gained or lost, so it
/// balances to the registry tolerance and dedupes by canonical structure for free.
pub fn chain_growth_repeat_unit_structure(
    monomer_id: SubstanceId,
    monomer: &MolecularStructure,
) -> ChemistryResult<Option<MolecularStructure>> {
    // Abstract catalog templates carry R-group placeholders and cannot be a
    // concrete repeat unit — skip them rather than erroring, the same way an
    // unrecognized monomer is skipped. Only fully specified monomers polymerize.
    if monomer.atoms.iter().any(|atom| atom.element == "R") {
        return Ok(None);
    }
    let Some((first, second)) = polymerizable_carbon_carbon_double_bond(monomer) else {
        return Ok(None);
    };
    let repeat_unit = alkene_repeat_unit(monomer_id, monomer, first, second)?;
    Ok(Some(repeat_unit.structure))
}

pub fn chain_growth_polymer_substance(
    monomer_id: SubstanceId,
    monomer: &MolecularStructure,
) -> ChemistryResult<Option<(Substance, u32)>> {
    if has_r_group(monomer) {
        return Ok(None);
    }
    let Some((first, second)) = polymerizable_carbon_carbon_double_bond(monomer) else {
        return Ok(None);
    };
    let repeat_unit = alkene_repeat_unit(monomer_id.clone(), monomer, first, second)?;
    let degree = REPRESENTATIVE_CHAIN_GROWTH_DEGREE;
    let material = PolymerMaterial {
        id: PolymerId::new(format!(
            "polymer:chain_growth:{}:{}-{}:n{}",
            sanitize_polymer_id(monomer_id.as_str()),
            first,
            second,
            degree
        ))?,
        mechanism: PolymerizationMechanism::ChainGrowthRadical,
        repeat_unit,
        distribution: ChainLengthDistribution::new(degree as f64, degree as f64)?,
        architecture: PolymerArchitecture::Linear,
        end_groups: Vec::new(),
    };
    Ok(Some((polymer_material_to_substance(&material)?, degree)))
}

/// The nucleophilic partner in an AA+BB polycondensation: a diol gives a
/// polyester (O-acyl linkage), a diamine gives a polyamide (N-acyl linkage).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepGrowthLinkage {
    Polyester,
    Polyamide,
}

/// One condensable end of a difunctional monomer: the atom that forms the new
/// linkage and the leaving hydrogen that departs with it (as part of water).
#[derive(Debug, Clone, Copy)]
struct CondensableEnd {
    /// Index of the atom that bonds into the chain (acyl C, diol O, or amine N).
    link_atom: usize,
    /// A hydrogen on (or on the hydroxyl of) this end that leaves as water.
    leaving_hydrogen: usize,
    /// For a carboxylic acid, the hydroxyl oxygen that also leaves; `None` for the
    /// nucleophile, which only sheds a hydrogen.
    leaving_oxygen: Option<usize>,
}

/// Builds the repeat-unit structure for an AA+BB step-growth polycondensation of
/// a diacid with a diol (polyester) or diamine (polyamide), or `None` if the two
/// monomers are not a clean difunctional pair. One internal linkage forms and two
/// waters leave, so the reaction `1 diacid + 1 comonomer → 1 repeat-unit + 2 H2O`
/// is mass-exact; the unit's two far ends are left open (`UnsaturatedAllowed`) to
/// carry the chain. Mirrors the analytics-vs-simulator split of
/// [`chain_growth_repeat_unit_structure`].
pub fn step_growth_repeat_unit_structure(
    diacid: &MolecularStructure,
    comonomer: &MolecularStructure,
) -> ChemistryResult<Option<(MolecularStructure, StepGrowthLinkage)>> {
    if has_r_group(diacid) || has_r_group(comonomer) {
        return Ok(None);
    }
    let Some(acid_ends) = difunctional_acid_ends(diacid) else {
        return Ok(None);
    };
    let Some((linkage, nucleophile_ends)) = difunctional_nucleophile_ends(comonomer) else {
        return Ok(None);
    };
    let structure = assemble_polycondensation_unit(diacid, &acid_ends, comonomer, &nucleophile_ends)?;
    Ok(Some((structure, linkage)))
}

pub fn step_growth_polymer_substance(
    diacid_id: &SubstanceId,
    diacid: &MolecularStructure,
    comonomer_id: &SubstanceId,
    comonomer: &MolecularStructure,
) -> ChemistryResult<Option<(Substance, StepGrowthLinkage, u32)>> {
    let Some((repeat_unit_structure, linkage)) =
        step_growth_repeat_unit_structure(diacid, comonomer)?
    else {
        return Ok(None);
    };
    let connection_atoms = open_connection_atoms(&repeat_unit_structure)?;
    let repeat_unit = PolymerRepeatUnit::new(
        SubstanceId::new(format!(
            "polymer_repeat:{}+{}",
            sanitize_polymer_id(diacid_id.as_str()),
            sanitize_polymer_id(comonomer_id.as_str())
        ))?,
        repeat_unit_structure,
        connection_atoms,
    )?;
    let degree = REPRESENTATIVE_STEP_GROWTH_DEGREE;
    let material = PolymerMaterial {
        id: PolymerId::new(format!(
            "polymer:step_growth:{}:{}:{}:n{}",
            sanitize_polymer_id(diacid_id.as_str()),
            sanitize_polymer_id(comonomer_id.as_str()),
            match linkage {
                StepGrowthLinkage::Polyester => "polyester",
                StepGrowthLinkage::Polyamide => "polyamide",
            },
            degree
        ))?,
        mechanism: PolymerizationMechanism::StepGrowthCondensation,
        repeat_unit,
        distribution: ChainLengthDistribution::new(degree as f64, degree as f64)?,
        architecture: PolymerArchitecture::Linear,
        end_groups: Vec::new(),
    };
    Ok(Some((polymer_material_to_substance(&material)?, linkage, degree)))
}

fn open_connection_atoms(structure: &MolecularStructure) -> ChemistryResult<[usize; 2]> {
    let atoms = structure
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| {
            (atom.valence_saturation == ValenceSaturation::UnsaturatedAllowed).then_some(index)
        })
        .collect::<Vec<_>>();
    match atoms.as_slice() {
        [first, second] => Ok([*first, *second]),
        _ => Err(polymer_error("polymer repeat unit must have exactly two open connection atoms")),
    }
}

fn has_r_group(structure: &MolecularStructure) -> bool {
    structure.atoms.iter().any(|atom| atom.element == "R")
}

/// Finds the two carboxylic-acid ends of a diacid. Returns `None` unless the
/// molecule carries exactly two `-COOH` groups (a clean difunctional diacid).
fn difunctional_acid_ends(structure: &MolecularStructure) -> Option<[CondensableEnd; 2]> {
    let mut ends = Vec::new();
    for group in find_functional_groups(structure) {
        if group.group_type == FunctionalGroupType::CarboxylicAcid {
            // CarboxylicAcid atoms = [carbon, carbonyl_oxygen, hydroxyl_oxygen, hydrogen].
            ends.push(CondensableEnd {
                link_atom: group.atoms[0],
                leaving_hydrogen: group.atoms[3],
                leaving_oxygen: Some(group.atoms[2]),
            });
        }
    }
    match ends.len() {
        2 => Some([ends[0], ends[1]]),
        _ => None,
    }
}

/// Finds the two nucleophilic ends of a comonomer and classifies the linkage:
/// two `-OH` alcohols → polyester, two primary-amine `-NH2` → polyamide. Returns
/// `None` for anything that is not cleanly one difunctional nucleophile type.
fn difunctional_nucleophile_ends(
    structure: &MolecularStructure,
) -> Option<(StepGrowthLinkage, [CondensableEnd; 2])> {
    let mut alcohols = Vec::new();
    let mut amines = Vec::new();
    for group in find_functional_groups(structure) {
        match group.group_type {
            // Alcohol atoms = [carbon, oxygen, hydrogen].
            FunctionalGroupType::Alcohol => alcohols.push(CondensableEnd {
                link_atom: group.atoms[1],
                leaving_hydrogen: group.atoms[2],
                leaving_oxygen: None,
            }),
            // PrimaryAmine atoms = [carbon, nitrogen, hydrogen, hydrogen].
            FunctionalGroupType::PrimaryAmine => amines.push(CondensableEnd {
                link_atom: group.atoms[1],
                leaving_hydrogen: group.atoms[2],
                leaving_oxygen: None,
            }),
            _ => {}
        }
    }
    match (alcohols.len(), amines.len()) {
        (2, 0) => Some((StepGrowthLinkage::Polyester, [alcohols[0], alcohols[1]])),
        (0, 2) => Some((StepGrowthLinkage::Polyamide, [amines[0], amines[1]])),
        _ => None,
    }
}

/// Joins a diacid and its difunctional comonomer into a single repeat-unit
/// structure. Exactly ONE linkage forms (acid end 0 to nucleophile end 0); both
/// acid ends shed their -OH and both nucleophile ends shed one H, so two waters
/// leave overall. The remaining acid and nucleophile ends are marked open
/// (`UnsaturatedAllowed`) to carry the chain — the same convention alkene repeat
/// units use for their backbone carbons.
fn assemble_polycondensation_unit(
    diacid: &MolecularStructure,
    acid_ends: &[CondensableEnd; 2],
    comonomer: &MolecularStructure,
    nucleophile_ends: &[CondensableEnd; 2],
) -> ChemistryResult<MolecularStructure> {
    let offset = diacid.atoms.len();
    let mut editor = MolecularEditor::new(diacid);
    // Append every comonomer atom (shifted by `offset`) and bond the two link
    // atoms. The molecule is transiently over-valent until the leaving groups go.
    editor.add_group(
        acid_ends[0].link_atom,
        comonomer,
        nucleophile_ends[0].link_atom,
        1.0,
    )?;
    let mut leaving = vec![
        acid_ends[0].leaving_hydrogen,
        acid_ends[1].leaving_hydrogen,
        nucleophile_ends[0].leaving_hydrogen + offset,
        nucleophile_ends[1].leaving_hydrogen + offset,
    ];
    for end in acid_ends {
        if let Some(oxygen) = end.leaving_oxygen {
            leaving.push(oxygen);
        }
    }
    let mapping = editor.remove_atoms(&leaving)?;
    let mut structure = editor.finish()?;
    let acid_open = mapping[acid_ends[1].link_atom]
        .ok_or_else(|| polymer_error("polycondensation lost its open acid end"))?;
    let nucleophile_open = mapping[nucleophile_ends[1].link_atom + offset]
        .ok_or_else(|| polymer_error("polycondensation lost its open nucleophile end"))?;
    mark_repeat_connection(&mut structure.atoms[acid_open]);
    mark_repeat_connection(&mut structure.atoms[nucleophile_open]);
    structure.source_code = "polymer-repeat-unit".to_string();
    Ok(structure)
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepGrowthInput {
    pub repeat_unit: PolymerRepeatUnit,
    pub functional_group_conversion: f64,
    pub stoichiometric_balance: f64,
    pub architecture: PolymerArchitecture,
}

pub fn derive_step_growth_polymer(input: StepGrowthInput) -> ChemistryResult<Option<PolymerMaterial>> {
    validate_fraction(
        input.functional_group_conversion,
        "step-growth functional group conversion",
    )?;
    validate_positive_finite(input.stoichiometric_balance, "step-growth stoichiometric balance")?;
    if input.stoichiometric_balance > 1.0 {
        return Err(polymer_error(
            "step-growth stoichiometric balance must be normalized to <= 1",
        ));
    }
    if input.functional_group_conversion < MIN_POLYMERIZATION_CONVERSION {
        return Ok(None);
    }
    let r = input.stoichiometric_balance;
    let p = input.functional_group_conversion;
    let denominator = 1.0 + r - 2.0 * r * p;
    if denominator <= 0.0 {
        return Err(polymer_error(
            "step-growth conversion and stoichiometry produce divergent chain length",
        ));
    }
    let number_average_degree = (1.0 + r) / denominator;
    let weight_average_degree = number_average_degree * (1.0 + p);
    let source = input.repeat_unit.source_monomer.clone();
    let material = PolymerMaterial {
        id: PolymerId::new(format!(
            "polymer:step_growth:{}",
            sanitize_polymer_id(source.as_str())
        ))?,
        mechanism: PolymerizationMechanism::StepGrowthCondensation,
        repeat_unit: input.repeat_unit,
        distribution: ChainLengthDistribution::new(number_average_degree, weight_average_degree)?,
        architecture: input.architecture,
        end_groups: Vec::new(),
    };
    material.validate()?;
    Ok(Some(material))
}

fn alkene_repeat_unit(
    monomer_id: SubstanceId,
    monomer: &MolecularStructure,
    first: usize,
    second: usize,
) -> ChemistryResult<PolymerRepeatUnit> {
    let mut structure = monomer.clone();
    for bond in &mut structure.bonds {
        if (bond.from == first && bond.to == second) || (bond.from == second && bond.to == first) {
            bond.order = 1.0;
            break;
        }
    }
    structure.source_code = "polymer-repeat-unit".to_string();
    mark_repeat_connection(&mut structure.atoms[first]);
    mark_repeat_connection(&mut structure.atoms[second]);
    PolymerRepeatUnit::new(monomer_id, structure, [first, second])
}

fn mark_repeat_connection(atom: &mut MolecularAtom) {
    atom.valence_saturation = ValenceSaturation::UnsaturatedAllowed;
}

fn polymerizable_carbon_carbon_double_bond(structure: &MolecularStructure) -> Option<(usize, usize)> {
    let mut result = None;
    for MolecularBond { from, to, order } in &structure.bonds {
        if !bond_order_matches(*order, 2.0) {
            continue;
        }
        if structure.atoms[*from].element != "C" || structure.atoms[*to].element != "C" {
            continue;
        }
        if is_aromatic_bond(structure, *from, *to) {
            continue;
        }
        if result.is_some() {
            return None;
        }
        result = Some((*from, *to));
    }
    result
}

fn is_aromatic_bond(structure: &MolecularStructure, first: usize, second: usize) -> bool {
    structure
        .bonds
        .iter()
        .any(|bond| {
            ((bond.from == first && bond.to == second) || (bond.from == second && bond.to == first))
                && bond_order_matches(bond.order, 1.5)
        })
}

fn validate_repeat_unit_structure(
    structure: &MolecularStructure,
    connection_atoms: [usize; 2],
) -> ChemistryResult<()> {
    structure.validate()?;
    if connection_atoms[0] == connection_atoms[1] {
        return Err(polymer_error("repeat unit connection atoms must be distinct"));
    }
    for atom in connection_atoms {
        let Some(connection) = structure.atoms.get(atom) else {
            return Err(polymer_error("repeat unit connection atom does not exist"));
        };
        if connection.element == "H" {
            return Err(polymer_error("hydrogen cannot be a polymer repeat connection atom"));
        }
    }
    if structure.atoms.iter().any(|atom| atom.element == "R") {
        return Err(polymer_error(
            "polymer repeat units must be concrete structures without R groups",
        ));
    }
    Ok(())
}

fn validate_fraction(value: f64, name: &str) -> ChemistryResult<()> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(polymer_error(&format!("{name} must be between 0 and 1")));
    }
    Ok(())
}

fn validate_positive_finite(value: f64, name: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(polymer_error(&format!("{name} must be positive and finite")));
    }
    Ok(())
}

fn sanitize_polymer_id(value: &str) -> String {
    value
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || ch == '_' || ch == '-' {
                ch
            } else {
                '_'
            }
        })
        .collect()
}

fn polymer_error(reason: &str) -> ChemistryError {
    ChemistryError::InvalidReaction {
        reaction_id: "polymer".to_string(),
        reason: reason.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::frowns::parse_frowns;

    #[test]
    fn polymer_material_becomes_a_nonvolatile_substance_carrying_its_repeat_unit() {
        let ethene = parse_frowns("C=C").unwrap();
        let material = derive_chain_growth_polymer(ChainGrowthInput::new(
            "test:ethene",
            ethene,
            0.80,
            0.01,
            0.0,
        ))
        .unwrap()
        .expect("ethene should produce a chain-growth polymer");

        let substance = polymer_material_to_substance(&material).unwrap();

        assert_eq!(substance.id.as_str(), material.id.as_str());
        assert_eq!(substance.charge, 0);
        assert!(
            (substance.molar_mass_grams - material.number_average_molar_mass_grams()).abs()
                < 1.0e-9
        );
        assert!(
            substance.boiling_point_kelvin >= POLYMER_BOILING_POINT_KELVIN,
            "a bulk polymer must be effectively non-volatile"
        );
        assert!(
            substance.molecular_structure.is_none(),
            "a polymer must not expose its repeat unit as an ordinary molecular graph"
        );
        assert!(substance.structure_code.is_some());
        assert!(matches!(
            substance.representation,
            SubstanceRepresentation::Polymer { .. }
        ));
        assert!(substance.validate().is_ok());
    }

    #[test]
    fn chain_growth_alkene_polymer_keeps_repeat_unit_instead_of_large_graph() {
        let ethene = parse_frowns("C=C").unwrap();
        let polymer = derive_chain_growth_polymer(ChainGrowthInput::new(
            "test:ethene",
            ethene.clone(),
            0.80,
            0.01,
            0.0,
        ))
        .unwrap()
        .expect("ethene should produce a chain-growth polymer");

        assert_eq!(polymer.repeat_unit.structure.atom_count(), ethene.atom_count());
        assert_eq!(polymer.repeat_unit.connection_atoms, [0, 1]);
        assert!(polymer.distribution.number_average_degree > 70.0);
        assert!(polymer.number_average_molar_mass_grams() > 1_900.0);
    }

    #[test]
    fn non_alkene_does_not_create_chain_growth_polymer() {
        let methane = parse_frowns("C").unwrap();
        let polymer = derive_chain_growth_polymer(ChainGrowthInput::new(
            "test:methane",
            methane,
            0.80,
            0.01,
            0.0,
        ))
        .unwrap();

        assert!(polymer.is_none());
    }

    #[test]
    fn step_growth_uses_carothers_distribution() {
        let repeat = PolymerRepeatUnit::new(
            "test:repeat",
            parse_frowns("CCO").unwrap(),
            [0, 2],
        )
        .unwrap();
        let polymer = derive_step_growth_polymer(StepGrowthInput {
            repeat_unit: repeat,
            functional_group_conversion: 0.95,
            stoichiometric_balance: 1.0,
            architecture: PolymerArchitecture::Linear,
        })
        .unwrap()
        .expect("high conversion should create step-growth polymer");

        assert!((polymer.distribution.number_average_degree - 20.0).abs() < 1.0e-9);
        assert!((polymer.distribution.dispersity - 1.95).abs() < 1.0e-9);
    }

    #[test]
    fn polyamide_repeat_unit_condenses_a_diacid_and_diamine_minus_two_waters() {
        // Nylon-6,6: adipic acid (C6H10O4) + hexanediamine (C6H16N2) lose two
        // waters to form the repeat unit, leaving one acid end and one amine end
        // open for the chain.
        let adipic = parse_frowns("O=C(OH)CCCCC(=O)OH").unwrap();
        let hexanediamine = parse_frowns("NCCCCCCN").unwrap();
        let (unit, linkage) = step_growth_repeat_unit_structure(&adipic, &hexanediamine)
            .unwrap()
            .expect("a diacid and a diamine should condense to a polyamide unit");

        assert_eq!(linkage, StepGrowthLinkage::Polyamide);
        // Water built with explicit hydrogens — parse_frowns("O") is a bare O atom.
        let water = parse_frowns("destroy:graph:atoms=H.O.H;bonds=0-s-1,1-s-2")
            .unwrap()
            .summary()
            .unwrap()
            .molar_mass_grams;
        let expected = adipic.summary().unwrap().molar_mass_grams
            + hexanediamine.summary().unwrap().molar_mass_grams
            - 2.0 * water;
        assert!(
            (unit.summary().unwrap().molar_mass_grams - expected).abs() < 1.0e-6,
            "repeat unit mass must equal the two monomers minus two waters"
        );
        // Exactly two atoms are left open to extend the chain.
        let open = unit
            .atoms
            .iter()
            .filter(|atom| atom.valence_saturation == ValenceSaturation::UnsaturatedAllowed)
            .count();
        assert_eq!(open, 2, "a linear repeat unit has two open connection atoms");
    }

    #[test]
    fn step_growth_unit_rejects_monofunctional_or_mismatched_monomers() {
        let acetic = parse_frowns("CC(=O)OH").unwrap();
        let hexanediamine = parse_frowns("NCCCCCCN").unwrap();
        // A monoacid is not difunctional, so no repeat unit forms.
        assert!(step_growth_repeat_unit_structure(&acetic, &hexanediamine)
            .unwrap()
            .is_none());
    }

}
