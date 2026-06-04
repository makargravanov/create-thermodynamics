include!("destroy_reactions.generated.rs");

use std::collections::BTreeMap;

use super::condition::{AtmosphereCondition, ReactionCondition};
use super::mixture::MixturePhase;
use super::registry::ChemistryRegistry;
use super::substance::{MaterialFormulaUnit, Substance, SubstanceId, SubstanceRepresentation};

pub const DESTROY_METALLURGY_REACTION_COUNT: usize = 28;

#[derive(Debug, Clone, PartialEq, Eq)]
struct MineralFormula {
    cations: BTreeMap<SubstanceId, u32>,
    oxide: u32,
    sulfide: u32,
    carbonate: u32,
    silicate: u32,
}

impl MineralFormula {
    fn from_units(units: &[MaterialFormulaUnit]) -> Self {
        let mut formula = Self {
            cations: BTreeMap::new(),
            oxide: 0,
            sulfide: 0,
            carbonate: 0,
            silicate: 0,
        };
        for unit in units {
            match unit.substance_id.as_str() {
                "destroy:oxide" => formula.oxide += unit.coefficient,
                "destroy:sulfide" => formula.sulfide += unit.coefficient,
                "destroy:carbonate" => formula.carbonate += unit.coefficient,
                "destroy:silicate" => formula.silicate += unit.coefficient,
                _ => {
                    *formula
                        .cations
                        .entry(unit.substance_id.clone())
                        .or_insert(0) += unit.coefficient;
                }
            }
        }
        formula
    }

    fn is_simple_oxide(&self) -> bool {
        self.oxide > 0 && self.sulfide == 0 && self.carbonate == 0 && self.silicate == 0
    }

    fn is_simple_sulfide(&self) -> bool {
        self.sulfide > 0 && self.oxide == 0 && self.carbonate == 0 && self.silicate == 0
    }

    fn is_simple_carbonate(&self) -> bool {
        self.carbonate > 0 && self.oxide == 0 && self.sulfide == 0 && self.silicate == 0
    }

    fn is_simple_silicate(&self) -> bool {
        self.silicate > 0 && self.oxide == 0 && self.sulfide == 0 && self.carbonate == 0
    }
}

pub fn destroy_metallurgy_reactions_registry_builder(
    builder: ChemistryRegistryBuilder,
) -> ChemistryResult<ChemistryRegistryBuilder> {
    let registry = builder.build()?;
    let mut builder = ChemistryRegistryBuilder::from_registry(&registry);
    for reaction in generated_metallurgy_reactions(&registry)? {
        if registry.reaction(&reaction.id).is_err() {
            builder = builder.reaction(reaction);
        }
    }
    Ok(builder)
}

fn generated_metallurgy_reactions(registry: &ChemistryRegistry) -> ChemistryResult<Vec<Reaction>> {
    let oxides = indexed_formulas(registry, FormulaKind::Oxide)?;
    let sulfides = indexed_formulas(registry, FormulaKind::Sulfide)?;
    let carbonates = indexed_formulas(registry, FormulaKind::Carbonate)?;
    let silicates = indexed_formulas(registry, FormulaKind::Silicate)?;
    let metals = metal_ids_by_element(registry);

    let mut reactions = Vec::new();
    for (carbonate, carbonate_formula) in &carbonates {
        if let Some((oxide, _)) = oxides.iter().find(|(_, oxide_formula)| {
            carbonate_decomposition_matches(carbonate_formula, oxide_formula)
        }) {
            reactions.push(carbonate_calcination_reaction(
                carbonate,
                oxide,
                carbonate_formula.carbonate,
                carbonate_calcination_temperature_kelvin(carbonate),
            ));
        }
    }
    for (sulfide, sulfide_formula) in &sulfides {
        if let Some((oxide, _)) = oxides
            .iter()
            .find(|(_, oxide_formula)| sulfide_roasting_matches(sulfide_formula, oxide_formula))
        {
            reactions.push(sulfide_roasting_reaction(
                sulfide,
                oxide,
                sulfide_formula.sulfide,
            ));
        }
    }
    for (oxide, oxide_formula) in &oxides {
        if let Some(metal_products) = metal_products_for_formula(oxide_formula, &metals) {
            if carbon_monoxide_reduction_allowed(oxide) {
                reactions.push(carbon_monoxide_reduction_reaction(
                    oxide,
                    oxide_formula.oxide,
                    &metal_products,
                ));
            }
            if hydrogen_reduction_allowed(oxide) {
                reactions.push(hydrogen_reduction_reaction(
                    oxide,
                    oxide_formula.oxide,
                    &metal_products,
                ));
            }
        }
    }
    for (silicate, silicate_formula) in &silicates {
        if let Some((base_oxide, base_oxide_count)) =
            oxides.iter().find_map(|(oxide, oxide_formula)| {
                slag_base_oxide_count(silicate_formula, oxide_formula)
                    .map(|coefficient| (oxide, coefficient))
            })
        {
            reactions.push(slag_formation_reaction(
                base_oxide,
                base_oxide_count,
                silicate,
                silicate_formula.silicate,
            ));
        }
    }
    Ok(reactions)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FormulaKind {
    Oxide,
    Sulfide,
    Carbonate,
    Silicate,
}

fn indexed_formulas(
    registry: &ChemistryRegistry,
    kind: FormulaKind,
) -> ChemistryResult<Vec<(&Substance, MineralFormula)>> {
    let mut formulas = Vec::new();
    for substance in registry.substances() {
        let units = match &substance.representation {
            SubstanceRepresentation::Oxide { formula_units } => Some(formula_units.as_slice()),
            SubstanceRepresentation::IonicSolid { formula_units } => Some(formula_units.as_slice()),
            _ => None,
        };
        let Some(units) = units else {
            continue;
        };
        let formula = MineralFormula::from_units(units);
        let matches = match kind {
            FormulaKind::Oxide => formula.is_simple_oxide(),
            FormulaKind::Sulfide => formula.is_simple_sulfide(),
            FormulaKind::Carbonate => formula.is_simple_carbonate(),
            FormulaKind::Silicate => formula.is_simple_silicate(),
        };
        if matches {
            formulas.push((substance, formula));
        }
    }
    Ok(formulas)
}

fn metal_ids_by_element(registry: &ChemistryRegistry) -> BTreeMap<String, SubstanceId> {
    registry
        .substances()
        .filter_map(|substance| match &substance.representation {
            SubstanceRepresentation::Metal { element_symbol } => {
                Some((element_symbol.clone(), substance.id.clone()))
            }
            _ => None,
        })
        .collect()
}

fn carbonate_decomposition_matches(carbonate: &MineralFormula, oxide: &MineralFormula) -> bool {
    carbonate.carbonate > 0
        && oxide.oxide == carbonate.carbonate
        && carbonate.cations == oxide.cations
}

fn sulfide_roasting_matches(sulfide: &MineralFormula, oxide: &MineralFormula) -> bool {
    sulfide.sulfide > 0 && oxide.oxide == sulfide.sulfide && sulfide.cations == oxide.cations
}

fn metal_products_for_formula(
    formula: &MineralFormula,
    metals: &BTreeMap<String, SubstanceId>,
) -> Option<Vec<(SubstanceId, u32)>> {
    let mut products = Vec::new();
    for (ion, coefficient) in &formula.cations {
        let metal = metal_id_for_ion(ion, metals)?;
        products.push((metal, *coefficient));
    }
    Some(products)
}

fn metal_id_for_ion(
    ion: &SubstanceId,
    metals: &BTreeMap<String, SubstanceId>,
) -> Option<SubstanceId> {
    let element = match ion.as_str() {
        "destroy:iron_ii" | "destroy:iron_iii" => "Fe",
        "destroy:copper_i" | "destroy:copper_ii" => "Cu",
        "destroy:zinc_ion" => "Zn",
        "destroy:nickel_ion" => "Ni",
        "destroy:lead_ii" => "Pb",
        "destroy:aluminum_iii" => "Al",
        "destroy:calcium_ion" => "Ca",
        "destroy:magnesium_ion" => "Mg",
        _ => return None,
    };
    metals.get(element).cloned()
}

fn hydrogen_reduction_allowed(oxide: &Substance) -> bool {
    !matches!(
        oxide.id.as_str(),
        "destroy:iron_iii_oxide"
            | "destroy:magnetite"
            | "destroy:copper_i_oxide"
            | "destroy:aluminum_oxide"
            | "destroy:calcium_oxide"
            | "destroy:magnesium_oxide"
            | "destroy:silica"
    )
}

fn carbon_monoxide_reduction_allowed(oxide: &Substance) -> bool {
    !matches!(
        oxide.id.as_str(),
        "destroy:aluminum_oxide"
            | "destroy:calcium_oxide"
            | "destroy:magnesium_oxide"
            | "destroy:silica"
    )
}

fn slag_base_oxide_count(silicate: &MineralFormula, oxide: &MineralFormula) -> Option<u32> {
    if silicate.silicate == 0 || oxide.oxide == 0 {
        return None;
    }
    let mut ratio = None;
    for (cation, silicate_count) in &silicate.cations {
        let oxide_count = oxide.cations.get(cation)?;
        if silicate_count % oxide_count != 0 {
            return None;
        }
        let current = silicate_count / oxide_count;
        if ratio.is_some_and(|known| known != current) {
            return None;
        }
        ratio = Some(current);
    }
    if oxide
        .cations
        .keys()
        .any(|cation| !silicate.cations.contains_key(cation))
    {
        return None;
    }
    let ratio = ratio?;
    if oxide.oxide * ratio == silicate.silicate * 2 {
        Some(ratio)
    } else {
        None
    }
}

fn carbonate_calcination_reaction(
    carbonate: &Substance,
    oxide: &Substance,
    carbonate_count: u32,
    min_temperature_kelvin: f64,
) -> Reaction {
    Reaction::builder(format!("{}.calcination", carbonate.id))
        .reactant(carbonate.id.clone(), 1, 1)
        .product(oxide.id.clone(), 1)
        .product("destroy:carbon_dioxide", carbonate_count)
        .reactant_phase_access(
            carbonate.id.clone(),
            [MixturePhase::Solid, MixturePhase::MoltenSlag],
        )
        .product_phase(oxide.id.clone(), MixturePhase::MoltenSlag)
        .product_phase("destroy:carbon_dioxide", MixturePhase::Gas)
        .condition(high_temperature_reaction(
            "carbonate calcination requires the solid or molten carbonate to be hot enough to decompose",
            min_temperature_kelvin,
        ))
        .pre_exponential_factor(5.0e7)
        .activation_energy_kj_per_mol(95.0)
        .enthalpy_change_kj_per_mol(170.0 * carbonate_count as f64)
        .build()
}

fn sulfide_roasting_reaction(
    sulfide: &Substance,
    oxide: &Substance,
    sulfide_count: u32,
) -> Reaction {
    Reaction::builder(format!("{}.roasting", sulfide.id))
        .reactant(sulfide.id.clone(), 2, 1)
        .reactant("destroy:oxygen", 3 * sulfide_count, 1)
        .product(oxide.id.clone(), 2)
        .product("destroy:sulfur_dioxide", 2 * sulfide_count)
        .reactant_phase_access(
            sulfide.id.clone(),
            [MixturePhase::Solid, MixturePhase::MoltenSlag],
        )
        .reactant_phase_access("destroy:oxygen", [MixturePhase::Gas])
        .product_phase(oxide.id.clone(), MixturePhase::MoltenSlag)
        .product_phase("destroy:sulfur_dioxide", MixturePhase::Gas)
        .condition(high_temperature_reaction(
            "sulfide roasting requires hot solid or molten sulfide and gaseous oxygen",
            850.0,
        ))
        .pre_exponential_factor(2.0e8)
        .activation_energy_kj_per_mol(80.0)
        .enthalpy_change_kj_per_mol(-330.0 * sulfide_count as f64)
        .build()
}

fn carbon_monoxide_reduction_reaction(
    oxide: &Substance,
    oxygen_count: u32,
    metal_products: &[(SubstanceId, u32)],
) -> Reaction {
    let mut builder = Reaction::builder(format!("{}.carbon_monoxide_reduction", oxide.id))
        .reactant(oxide.id.clone(), 1, 1)
        .reactant("destroy:carbon_monoxide", oxygen_count, 1)
        .product("destroy:carbon_dioxide", oxygen_count)
        .reactant_phase_access(
            oxide.id.clone(),
            [MixturePhase::Solid, MixturePhase::MoltenSlag],
        )
        .reactant_phase_access("destroy:carbon_monoxide", [MixturePhase::Gas])
        .product_phase("destroy:carbon_dioxide", MixturePhase::Gas)
        .condition(high_temperature_reaction(
            "carbon monoxide reduction requires hot oxide and gaseous carbon monoxide",
            reduction_temperature_kelvin(oxide),
        ))
        .condition(reducing_atmosphere(
            "carbon monoxide reduction is blocked by an oxidizing gas phase",
        ))
        .pre_exponential_factor(1.0e8)
        .activation_energy_kj_per_mol(85.0)
        .enthalpy_change_kj_per_mol(-25.0 * oxygen_count as f64);
    for (metal, coefficient) in metal_products {
        builder = builder
            .product(metal.clone(), *coefficient)
            .product_phase(metal.clone(), MixturePhase::MoltenMetal);
    }
    builder.build()
}

fn hydrogen_reduction_reaction(
    oxide: &Substance,
    oxygen_count: u32,
    metal_products: &[(SubstanceId, u32)],
) -> Reaction {
    let mut builder = Reaction::builder(format!("{}.hydrogen_reduction", oxide.id))
        .reactant(oxide.id.clone(), 1, 1)
        .reactant("destroy:hydrogen", oxygen_count, 1)
        .product("destroy:water", oxygen_count)
        .reactant_phase_access(
            oxide.id.clone(),
            [MixturePhase::Solid, MixturePhase::MoltenSlag],
        )
        .reactant_phase_access("destroy:hydrogen", [MixturePhase::Gas])
        .product_phase("destroy:water", MixturePhase::Gas)
        .condition(high_temperature_reaction(
            "hydrogen reduction requires hot oxide and gaseous hydrogen",
            reduction_temperature_kelvin(oxide),
        ))
        .condition(reducing_atmosphere(
            "hydrogen reduction is blocked by an oxidizing gas phase",
        ))
        .pre_exponential_factor(8.0e7)
        .activation_energy_kj_per_mol(90.0)
        .enthalpy_change_kj_per_mol(35.0 * oxygen_count as f64);
    for (metal, coefficient) in metal_products {
        builder = builder
            .product(metal.clone(), *coefficient)
            .product_phase(metal.clone(), MixturePhase::MoltenMetal);
    }
    builder.build()
}

fn slag_formation_reaction(
    base_oxide: &Substance,
    base_oxide_count: u32,
    silicate: &Substance,
    silica_count: u32,
) -> Reaction {
    Reaction::builder(format!("{}.slag_formation", silicate.id))
        .reactant(base_oxide.id.clone(), base_oxide_count, 1)
        .reactant("destroy:silica", silica_count, 1)
        .product(silicate.id.clone(), 1)
        .reactant_phase_access(
            base_oxide.id.clone(),
            [MixturePhase::Solid, MixturePhase::MoltenSlag],
        )
        .reactant_phase_access(
            "destroy:silica",
            [MixturePhase::Solid, MixturePhase::MoltenSlag],
        )
        .product_phase(silicate.id.clone(), MixturePhase::MoltenSlag)
        .condition(high_temperature_reaction(
            "slag formation requires hot basic oxide and silica in the solid or molten slag phase",
            slag_temperature_kelvin(base_oxide, silicate),
        ))
        .pre_exponential_factor(3.0e7)
        .activation_energy_kj_per_mol(70.0)
        .enthalpy_change_kj_per_mol(-80.0)
        .build()
}

fn high_temperature_reaction(reason: &str, min_temperature_kelvin: f64) -> ReactionCondition {
    ReactionCondition::new(reason)
        .min_temperature_kelvin(min_temperature_kelvin)
        .rate_multiplier(1.0)
}

fn reducing_atmosphere(reason: &str) -> ReactionCondition {
    ReactionCondition::new(reason)
        .atmosphere(AtmosphereCondition::Inert)
        .max_oxygen_activity(1.0e-8)
}

fn reduction_temperature_kelvin(oxide: &Substance) -> f64 {
    match oxide.id.as_str() {
        "destroy:copper_ii_oxide" => 650.0,
        "destroy:copper_i_oxide" | "destroy:lead_ii_oxide" => 750.0,
        "destroy:nickel_ii_oxide" => 900.0,
        "destroy:iron_ii_oxide" => 1_000.0,
        "destroy:iron_iii_oxide" | "destroy:magnetite" => 1_050.0,
        "destroy:zinc_oxide" => 1_200.0,
        _ => oxide.melting_point_kelvin.min(1_500.0).max(900.0),
    }
}

fn carbonate_calcination_temperature_kelvin(carbonate: &Substance) -> f64 {
    carbonate.melting_point_kelvin.min(1_173.0)
}

fn slag_temperature_kelvin(base_oxide: &Substance, silicate: &Substance) -> f64 {
    match silicate.id.as_str() {
        "destroy:iron_ii_silicate" => 1_200.0,
        _ => base_oxide
            .melting_point_kelvin
            .min(silicate.melting_point_kelvin)
            .min(1_450.0)
            .max(1_200.0),
    }
}
