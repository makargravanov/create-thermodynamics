use super::catalysis::CatalystSurfaceSpec;
use super::complex::{ComplexGeometry, ComplexLigand, ComplexSpec, LigandExchangeLability};
use super::error::{ChemistryError, ChemistryResult};
use super::mixture::MixturePhase;
use super::molecule::{
    element_mass, parse_java_structure, parse_legacy_structure, MolecularStructure,
    MolecularSummary,
};
use super::redox::{RedoxEnvironment, RedoxHalfReaction};
use super::registry::{ChemistryRegistryBuilder, GasSolubilityModel, SolventMiscibility};
use super::solution::{AcidBaseSpec, EquilibriumSpec};
use super::substance::{
    LiquidPhasePreference, MaterialFormulaUnit, SolventRole, Substance, SubstanceId,
    SubstancePhaseProperties, SubstanceRepresentation, SubstanceTagId,
};

const DEFAULT_DENSITY_GRAMS_PER_BUCKET: f64 = 1000.0;
const DEFAULT_MOLAR_HEAT_CAPACITY: f64 = 100.0;
const DEFAULT_LATENT_HEAT: f64 = 20_000.0;
const DEFAULT_COLOR_ARGB: u32 = 0x20FF_FFFF;

pub const DESTROY_SUBSTANCE_COUNT: usize = DESTROY_SUBSTANCES.len();

#[derive(Debug, Clone, Copy)]
struct RawSubstance {
    id: &'static str,
    structure_code: Option<&'static str>,
    java_structure_code: Option<&'static str>,
    translation_key: Option<&'static str>,
    boiling_point_celsius: Option<f64>,
    boiling_point_kelvin: Option<f64>,
    density: Option<f64>,
    molar_heat_capacity: Option<f64>,
    specific_heat_capacity: Option<f64>,
    color_argb: u32,
    tags: &'static [&'static str],
}

#[derive(Debug, Clone, Copy)]
struct RawMetallurgyIon {
    id: &'static str,
    element: Option<&'static str>,
    charge: i32,
    color_argb: u32,
}

#[derive(Debug, Clone, Copy)]
struct RawMetallurgyMetal {
    id: &'static str,
    element: &'static str,
    solid_density: f64,
    melting_point_kelvin: f64,
    boiling_point_kelvin: f64,
    heat_capacity: f64,
    fusion_heat: f64,
    vaporization_heat: f64,
    color_argb: u32,
}

#[derive(Debug, Clone, Copy)]
struct RawMetallurgyMaterial {
    id: &'static str,
    representation: RawMetallurgyMaterialRepresentation,
    formula_units: &'static [(&'static str, u32)],
    solid_density: f64,
    melting_point_kelvin: f64,
    boiling_point_kelvin: f64,
    heat_capacity: f64,
    fusion_heat: f64,
    color_argb: u32,
}

#[derive(Debug, Clone, Copy)]
enum RawMetallurgyMaterialRepresentation {
    IonicSolid,
    Oxide,
}

pub fn destroy_substances_registry_builder() -> ChemistryResult<ChemistryRegistryBuilder> {
    let mut builder = ChemistryRegistryBuilder::new();
    for tag in DESTROY_TAGS {
        builder = builder.substance_tag(*tag);
    }
    for raw in DESTROY_SUBSTANCES {
        builder = builder.substance(raw.to_substance()?);
    }
    builder = register_metallurgy_substances(builder)?;
    builder = register_phase_tables(builder);
    Ok(builder)
}

impl RawMetallurgyIon {
    fn to_substance(self) -> ChemistryResult<Substance> {
        let molar_mass = match self.element {
            Some(element) => element_mass(element)?,
            None => material_formula_mass(self.id, self.formula_units())?,
        };
        let mut substance = Substance::new(
            SubstanceId::new(format!("destroy:{}", self.id))?,
            self.charge,
            molar_mass,
            DEFAULT_DENSITY_GRAMS_PER_BUCKET,
            f64::MAX,
            DEFAULT_MOLAR_HEAT_CAPACITY,
            DEFAULT_LATENT_HEAT,
        )
        .with_phase_properties(SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: Some(10.0),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: false,
            can_form_liquid_phase: false,
            solvent_role: SolventRole::NotSolvent,
        })
        .with_catalog_metadata(None, None, self.color_argb, Vec::new());
        if let Some(element) = self.element {
            substance = substance.with_representation(SubstanceRepresentation::Ion {
                parent_element: Some(element.to_string()),
            });
        }
        Ok(substance)
    }

    fn formula_units(self) -> &'static [(&'static str, u32)] {
        match self.id {
            "carbonate" => &[("carbon", 1), ("oxide_atom", 3)],
            "silicate" => &[("silicon_atom", 1), ("oxide_atom", 4)],
            _ => &[],
        }
    }
}

impl RawMetallurgyMetal {
    fn to_substance(self) -> ChemistryResult<Substance> {
        Ok(Substance::new(
            SubstanceId::new(format!("destroy:{}", self.id))?,
            0,
            element_mass(self.element)?,
            self.solid_density,
            self.boiling_point_kelvin,
            self.heat_capacity,
            self.vaporization_heat,
        )
        .with_solid_density_grams_per_bucket(self.solid_density)
        .with_melting_point_kelvin(self.melting_point_kelvin)
        .with_fusion_heat_j_per_mol(self.fusion_heat)
        .with_phase_properties(solid_material_phase_properties())
        .with_representation(SubstanceRepresentation::Metal {
            element_symbol: self.element.to_string(),
        })
        .with_catalog_metadata(None, None, self.color_argb, Vec::new()))
    }
}

impl RawMetallurgyMaterial {
    fn to_substance(self) -> ChemistryResult<Substance> {
        let formula_units = self
            .formula_units
            .iter()
            .map(|(id, coefficient)| MaterialFormulaUnit::new(SubstanceId::from(*id), *coefficient))
            .collect::<Vec<_>>();
        let representation = match self.representation {
            RawMetallurgyMaterialRepresentation::IonicSolid => {
                SubstanceRepresentation::IonicSolid {
                    formula_units: formula_units.clone(),
                }
            }
            RawMetallurgyMaterialRepresentation::Oxide => SubstanceRepresentation::Oxide {
                formula_units: formula_units.clone(),
            },
        };
        Ok(Substance::new(
            SubstanceId::new(format!("destroy:{}", self.id))?,
            0,
            material_formula_mass(self.id, self.formula_units)?,
            self.solid_density,
            self.boiling_point_kelvin,
            self.heat_capacity,
            DEFAULT_LATENT_HEAT,
        )
        .with_solid_density_grams_per_bucket(self.solid_density)
        .with_melting_point_kelvin(self.melting_point_kelvin)
        .with_fusion_heat_j_per_mol(self.fusion_heat)
        .with_phase_properties(solid_material_phase_properties())
        .with_representation(representation)
        .with_catalog_metadata(None, None, self.color_argb, Vec::new()))
    }
}

fn material_formula_mass(
    material_id: &str,
    formula_units: &[(&'static str, u32)],
) -> ChemistryResult<f64> {
    let mut mass = 0.0;
    for (id, coefficient) in formula_units {
        let component_mass = match *id {
            "destroy:aluminum_iii" => element_mass("Al")?,
            "destroy:magnesium_ion" => element_mass("Mg")?,
            "destroy:silicon_iv" => element_mass("Si")?,
            "destroy:carbonate" => element_mass("C")? + element_mass("O")? * 3.0,
            "destroy:silicate" => element_mass("Si")? + element_mass("O")? * 4.0,
            "destroy:calcium_ion" => element_mass("Ca")?,
            "destroy:copper_i" | "destroy:copper_ii" => element_mass("Cu")?,
            "destroy:iron_ii" | "destroy:iron_iii" => element_mass("Fe")?,
            "destroy:lead_ii" => element_mass("Pb")?,
            "destroy:nickel_ion" => element_mass("Ni")?,
            "destroy:oxide" => element_mass("O")?,
            "destroy:sulfide" => element_mass("S")?,
            "destroy:zinc_ion" => element_mass("Zn")?,
            "carbon" => element_mass("C")?,
            "oxide_atom" => element_mass("O")?,
            "silicon_atom" => element_mass("Si")?,
            _ => {
                return Err(ChemistryError::InvalidSubstance {
                    substance_id: format!("destroy:{material_id}"),
                    reason: format!("unknown metallurgy formula component '{id}'"),
                });
            }
        };
        mass += component_mass * *coefficient as f64;
    }
    if !mass.is_finite() || mass <= 0.0 {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: format!("destroy:{material_id}"),
            reason: "metallurgy formula produced invalid mass".to_string(),
        });
    }
    Ok(mass)
}

fn solid_material_phase_properties() -> SubstancePhaseProperties {
    SubstancePhaseProperties {
        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
        aqueous_solubility_mol_per_bucket: Some(0.0),
        organic_solubility_mol_per_bucket: Some(0.0),
        can_precipitate: true,
        can_form_liquid_phase: false,
        solvent_role: SolventRole::NotSolvent,
    }
}

fn register_metallurgy_substances(
    mut builder: ChemistryRegistryBuilder,
) -> ChemistryResult<ChemistryRegistryBuilder> {
    for ion in METALLURGY_IONS {
        builder = builder.substance(ion.to_substance()?);
    }
    for metal in METALLURGY_METALS {
        builder = builder.substance(metal.to_substance()?);
    }
    for material in METALLURGY_MATERIALS {
        builder = builder.substance(material.to_substance()?);
    }
    Ok(builder)
}

fn register_phase_tables(builder: ChemistryRegistryBuilder) -> ChemistryRegistryBuilder {
    builder
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "forge:dusts/nickel",
            58.6934,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "forge:dusts/palladium",
            106.42,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "forge:dusts/platinum",
            195.084,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "forge:dusts/rhodium",
            102.9055,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical("forge:dusts/iron", 55.845, 0))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "forge:dusts/copper",
            63.546,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical("forge:dusts/zinc", 65.38, 0))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/nickel\"), 1f)",
            58.6934,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/palladium\"), 1f)",
            106.42,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/platinum\"), 1f)",
            195.084,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/platinum\"), 3f)",
            195.084,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/iron\"), 1f)",
            55.845,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/copper\"), 1f)",
            63.546,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/zinc\"), 1f)",
            65.38,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::chemical(
            "addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/rhodium\"), 1f)",
            102.9055,
            0,
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::unchecked(
            "addSimpleItemCatalyst(DestroyItems.ZEOLITE::get, 1f)",
            "zeolite item surface has no explicit chemical formula in the Rust model",
        ))
        .catalyst_surface_spec(CatalystSurfaceSpec::unchecked(
            "addSimpleItemCatalyst(DestroyItems.SILICA::get, 1f)",
            "silica item surface has no explicit chemical formula in the Rust model",
        ))
        .complex_spec(
            ComplexSpec::new(
                "destroy:copper_ii_tetraammine",
                "destroy:copper_ii",
                [ComplexLigand::new("destroy:ammonia", 4)],
                2,
                1.0e13,
            )
            .with_coordination_number(4)
            .with_geometry(ComplexGeometry::SquarePlanar)
            .with_ligand_exchange_lability(LigandExchangeLability::Labile)
            .with_translation_key("copper_ii_tetraammine")
            .with_color_argb(0x804A90E2),
        )
        .complex_spec(
            ComplexSpec::new(
                "destroy:nickel_tetraammine",
                "destroy:nickel_ion",
                [ComplexLigand::new("destroy:ammonia", 4)],
                2,
                1.0e8,
            )
            .with_coordination_number(4)
            .with_geometry(ComplexGeometry::Tetrahedral)
            .with_ligand_exchange_lability(LigandExchangeLability::Labile)
            .with_translation_key("nickel_tetraammine")
            .with_color_argb(0x8062B25D),
        )
        .complex_spec(
            ComplexSpec::new(
                "destroy:ferric_hexacyanide",
                "destroy:iron_iii",
                [ComplexLigand::new("destroy:cyanide", 6)],
                -3,
                1.0e31,
            )
            .with_coordination_number(6)
            .with_geometry(ComplexGeometry::Octahedral)
            .with_ligand_exchange_lability(LigandExchangeLability::Inert)
            .with_translation_key("ferric_hexacyanide")
            .with_color_argb(0x80364D9B),
        )
        .complex_spec(
            ComplexSpec::new(
                "destroy:cuprous_dicyanide",
                "destroy:copper_i",
                [ComplexLigand::new("destroy:cyanide", 2)],
                -1,
                1.0e24,
            )
            .with_coordination_number(2)
            .with_geometry(ComplexGeometry::Linear)
            .with_ligand_exchange_lability(LigandExchangeLability::Intermediate)
            .with_translation_key("cuprous_dicyanide")
            .with_color_argb(0x8077A5A1),
        )
        .complex_spec(
            ComplexSpec::new(
                "destroy:zinc_tetraammine",
                "destroy:zinc_ion",
                [ComplexLigand::new("destroy:ammonia", 4)],
                2,
                3.0e9,
            )
            .with_coordination_number(4)
            .with_geometry(ComplexGeometry::Tetrahedral)
            .with_ligand_exchange_lability(LigandExchangeLability::Labile)
            .with_translation_key("zinc_tetraammine")
            .with_color_argb(0x8090A8C8),
        )
        .complex_spec(
            ComplexSpec::new(
                "destroy:ferric_tetrachloride",
                "destroy:iron_iii",
                [ComplexLigand::new("destroy:chloride", 4)],
                -1,
                1.0e2,
            )
            .with_coordination_number(4)
            .with_geometry(ComplexGeometry::Tetrahedral)
            .with_ligand_exchange_lability(LigandExchangeLability::Labile)
            .with_translation_key("ferric_tetrachloride")
            .with_color_argb(0x80B78238),
        )
        .gas_solubility(
            "destroy:oxygen",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 1.3e-8,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.12,
                transfer_coefficient_per_tick: 0.15,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:nitrogen",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 6.4e-9,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.08,
                transfer_coefficient_per_tick: 0.12,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:carbon_monoxide",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 9.5e-9,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.08,
                transfer_coefficient_per_tick: 0.14,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:carbon_dioxide",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 3.4e-7,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.20,
                transfer_coefficient_per_tick: 0.18,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:chlorine",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 6.0e-5,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.25,
                transfer_coefficient_per_tick: 0.20,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:ammonia",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 5.8e-4,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.15,
                transfer_coefficient_per_tick: 0.25,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:hydrochloric_acid",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 1.0e-2,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.30,
                transfer_coefficient_per_tick: 0.30,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:sulfur_dioxide",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 1.2e-5,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.20,
                transfer_coefficient_per_tick: 0.22,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:nitrogen_dioxide",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 1.0e-5,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.20,
                transfer_coefficient_per_tick: 0.22,
                estimated: true,
            },
        )
        .gas_solubility(
            "destroy:hydrogen",
            GasSolubilityModel::Henry {
                henry_mol_per_bucket_pascal: 7.8e-9,
                temperature_kelvin: 298.0,
                salting_out_coefficient: 0.08,
                transfer_coefficient_per_tick: 0.15,
                estimated: true,
            },
        )
        .solvent_miscibility(
            "destroy:water",
            "destroy:ethanol",
            SolventMiscibility::FullyMiscible,
        )
        .solvent_miscibility(
            "destroy:water",
            "destroy:acetone",
            SolventMiscibility::FullyMiscible,
        )
        .solvent_miscibility(
            "destroy:water",
            "destroy:chloroform",
            SolventMiscibility::PartiallyMiscible {
                limit_mol_per_bucket: 0.1,
            },
        )
        .equilibrium(EquilibriumSpec::new(
            "destroy:water.autoionization",
            [(SubstanceId::from("destroy:water"), 1, MixturePhase::Aqueous)],
            [
                (
                    SubstanceId::from("destroy:proton"),
                    1,
                    MixturePhase::Aqueous,
                ),
                (
                    SubstanceId::from("destroy:hydroxide"),
                    1,
                    MixturePhase::Aqueous,
                ),
            ],
            1.0e-14,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:acetic_acid",
            "destroy:acetic_acid",
            "destroy:acetate",
            4.76,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:ammonium",
            "destroy:ammonium",
            "destroy:ammonia",
            9.25,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:hydrochloric_acid",
            "destroy:hydrochloric_acid",
            "destroy:chloride",
            -6.3,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:hydrofluoric_acid",
            "destroy:hydrofluoric_acid",
            "destroy:fluoride",
            3.17,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:hydrogen_cyanide",
            "destroy:hydrogen_cyanide",
            "destroy:cyanide",
            9.21,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:hydrogen_iodide",
            "destroy:hydrogen_iodide",
            "destroy:iodide",
            -10.0,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:hydrogensulfate",
            "destroy:hydrogensulfate",
            "destroy:sulfate",
            1.99,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:hypochlorous_acid",
            "destroy:hypochlorous_acid",
            "destroy:hypochlorite",
            7.53,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:nitric_acid",
            "destroy:nitric_acid",
            "destroy:nitrate",
            -1.4,
        ))
        .acid_base_pair(AcidBaseSpec::new(
            "destroy:sulfuric_acid",
            "destroy:sulfuric_acid",
            "destroy:hydrogensulfate",
            -3.0,
        ))
        .equilibrium(EquilibriumSpec::new(
            "destroy:boric_acid.hydrolysis",
            [
                (
                    SubstanceId::from("destroy:boric_acid"),
                    1,
                    MixturePhase::Aqueous,
                ),
                (SubstanceId::from("destroy:water"), 1, MixturePhase::Aqueous),
            ],
            [
                (
                    SubstanceId::from("destroy:tetrahydroxyborate"),
                    1,
                    MixturePhase::Aqueous,
                ),
                (
                    SubstanceId::from("destroy:proton"),
                    1,
                    MixturePhase::Aqueous,
                ),
            ],
            10.0_f64.powf(-9.24),
        ))
        .redox_half_reaction(
            RedoxHalfReaction::oxidation(
                "destroy:iron_ii_to_iron_iii",
                [(SubstanceId::from("destroy:iron_ii"), 1)],
                [(SubstanceId::from("destroy:iron_iii"), 1)],
                1,
                RedoxEnvironment::Any,
            )
            .with_standard_potential_volts(-0.771),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:iron_iii_to_iron_ii",
                [(SubstanceId::from("destroy:iron_iii"), 1)],
                [(SubstanceId::from("destroy:iron_ii"), 1)],
                1,
                RedoxEnvironment::Any,
            )
            .with_standard_potential_volts(0.771),
        )
        .redox_half_reaction(
            RedoxHalfReaction::oxidation(
                "destroy:copper_i_to_copper_ii",
                [(SubstanceId::from("destroy:copper_i"), 1)],
                [(SubstanceId::from("destroy:copper_ii"), 1)],
                1,
                RedoxEnvironment::Any,
            )
            .with_standard_potential_volts(-0.153),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:copper_ii_to_copper_i",
                [(SubstanceId::from("destroy:copper_ii"), 1)],
                [(SubstanceId::from("destroy:copper_i"), 1)],
                1,
                RedoxEnvironment::Any,
            )
            .with_standard_potential_volts(0.153),
        )
        .redox_half_reaction(
            RedoxHalfReaction::oxidation(
                "destroy:iodide_to_iodine",
                [(SubstanceId::from("destroy:iodide"), 2)],
                [(SubstanceId::from("destroy:iodine"), 1)],
                2,
                RedoxEnvironment::Any,
            )
            .with_standard_potential_volts(-0.535),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:iodine_to_iodide",
                [(SubstanceId::from("destroy:iodine"), 1)],
                [(SubstanceId::from("destroy:iodide"), 2)],
                2,
                RedoxEnvironment::Any,
            )
            .with_standard_potential_volts(0.535),
        )
        .redox_half_reaction(
            RedoxHalfReaction::oxidation(
                "destroy:hydrogen_to_proton",
                [(SubstanceId::from("destroy:hydrogen"), 1)],
                [(SubstanceId::from("destroy:proton"), 2)],
                2,
                RedoxEnvironment::Acidic,
            )
            .with_standard_potential_volts(0.0),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:proton_to_hydrogen",
                [(SubstanceId::from("destroy:proton"), 2)],
                [(SubstanceId::from("destroy:hydrogen"), 1)],
                2,
                RedoxEnvironment::Acidic,
            )
            .with_standard_potential_volts(0.0),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:oxygen_to_water",
                [
                    (SubstanceId::from("destroy:oxygen"), 1),
                    (SubstanceId::from("destroy:proton"), 4),
                ],
                [(SubstanceId::from("destroy:water"), 2)],
                4,
                RedoxEnvironment::Acidic,
            )
            .with_standard_potential_volts(1.229),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:hydrogen_peroxide_to_water",
                [
                    (SubstanceId::from("destroy:hydrogen_peroxide"), 1),
                    (SubstanceId::from("destroy:proton"), 2),
                ],
                [(SubstanceId::from("destroy:water"), 2)],
                2,
                RedoxEnvironment::Acidic,
            )
            .with_standard_potential_volts(1.776),
        )
        .redox_half_reaction(
            RedoxHalfReaction::oxidation(
                "destroy:hydrogen_peroxide_to_oxygen",
                [(SubstanceId::from("destroy:hydrogen_peroxide"), 1)],
                [
                    (SubstanceId::from("destroy:oxygen"), 1),
                    (SubstanceId::from("destroy:proton"), 2),
                ],
                2,
                RedoxEnvironment::Acidic,
            )
            .with_standard_potential_volts(-0.682),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:chlorine_to_chloride",
                [(SubstanceId::from("destroy:chlorine"), 1)],
                [(SubstanceId::from("destroy:chloride"), 2)],
                2,
                RedoxEnvironment::Any,
            )
            .with_standard_potential_volts(1.358),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:hypochlorous_acid_to_chloride",
                [
                    (SubstanceId::from("destroy:hypochlorous_acid"), 1),
                    (SubstanceId::from("destroy:proton"), 1),
                ],
                [
                    (SubstanceId::from("destroy:chloride"), 1),
                    (SubstanceId::from("destroy:water"), 1),
                ],
                2,
                RedoxEnvironment::Acidic,
            )
            .with_standard_potential_volts(1.49),
        )
        .redox_half_reaction(
            RedoxHalfReaction::reduction(
                "destroy:dichromate_to_chromium_iii",
                [
                    (SubstanceId::from("destroy:dichromate"), 1),
                    (SubstanceId::from("destroy:proton"), 14),
                ],
                [
                    (SubstanceId::from("destroy:chromium_iii"), 2),
                    (SubstanceId::from("destroy:water"), 7),
                ],
                6,
                RedoxEnvironment::Acidic,
            )
            .with_standard_potential_volts(1.33),
        )
}

impl RawSubstance {
    fn to_substance(self) -> ChemistryResult<Substance> {
        let structure =
            parse_raw_structure(self.id, self.structure_code, self.java_structure_code)?;
        let summary = summarize_structure(self.id, &structure)?;
        let boiling_point_kelvin = if summary.charge != 0 {
            f64::MAX
        } else if let Some(value) = self.boiling_point_kelvin {
            value
        } else if let Some(value) = self.boiling_point_celsius {
            value + 273.0
        } else {
            estimate_boiling_point(summary.molar_mass_grams)
        };
        let molar_heat_capacity = if let Some(value) = self.molar_heat_capacity {
            value
        } else if let Some(value) = self.specific_heat_capacity {
            value / summary.molar_mass_grams
        } else {
            DEFAULT_MOLAR_HEAT_CAPACITY
        };
        let color = if self.color_argb == 0 {
            DEFAULT_COLOR_ARGB
        } else {
            self.color_argb
        };
        let tags = self
            .tags
            .iter()
            .map(|tag| SubstanceTagId::new(format!("destroy:{tag}")))
            .collect::<ChemistryResult<Vec<_>>>()?;
        Ok(Substance::new(
            SubstanceId::new(format!("destroy:{}", self.id))?,
            summary.charge,
            summary.molar_mass_grams,
            self.density.unwrap_or(DEFAULT_DENSITY_GRAMS_PER_BUCKET),
            boiling_point_kelvin,
            molar_heat_capacity,
            DEFAULT_LATENT_HEAT,
        )
        .with_phase_properties(estimate_phase_properties(self.id, &summary, self.tags))
        .with_catalog_metadata(
            self.structure_code
                .or(self.java_structure_code)
                .map(str::to_string),
            self.translation_key.map(str::to_string),
            color,
            tags,
        )
        .with_molecular_structure(structure))
    }
}

fn estimate_phase_properties(
    id: &str,
    summary: &MolecularSummary,
    tags: &[&str],
) -> SubstancePhaseProperties {
    if id == "water" {
        return SubstancePhaseProperties::aqueous_solvent();
    }
    if matches!(id, "proton" | "hydroxide") {
        return SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: None,
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: false,
            can_form_liquid_phase: true,
            solvent_role: SolventRole::NotSolvent,
        };
    }
    if summary.charge != 0 {
        return SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: Some(10.0),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: false,
            solvent_role: SolventRole::NotSolvent,
        };
    }
    if tags.contains(&"solvent") {
        return SubstancePhaseProperties::organic_unlimited(0.1);
    }
    if id == "trimethylphosphine" {
        return SubstancePhaseProperties::organic_solute(0.0);
    }
    if id.ends_with("_acid") || id.contains("acid") || id == "ammonia" {
        return SubstancePhaseProperties {
            preferred_liquid_phase: LiquidPhasePreference::Aqueous,
            aqueous_solubility_mol_per_bucket: None,
            organic_solubility_mol_per_bucket: Some(0.25),
            can_precipitate: false,
            can_form_liquid_phase: true,
            solvent_role: SolventRole::KnownSolvent,
        };
    }
    SubstancePhaseProperties::organic_unlimited(0.05)
}

fn estimate_boiling_point(molar_mass_grams: f64) -> f64 {
    2.042_598_921_281_41 * molar_mass_grams + 178.176_866_128_713
}

pub fn summarize_legacy_structure(structure_code: &str) -> ChemistryResult<MolecularSummary> {
    parse_legacy_structure(structure_code)?.summary()
}

fn parse_raw_structure(
    id: &str,
    structure_code: Option<&str>,
    java_structure_code: Option<&str>,
) -> ChemistryResult<MolecularStructure> {
    match (structure_code, java_structure_code) {
        (Some(code), _) => parse_legacy_structure(code),
        (None, Some(code)) => parse_java_structure(code),
        (None, None) => Err(ChemistryError::InvalidSubstance {
            substance_id: format!("destroy:{id}"),
            reason: "substance has no structure".to_string(),
        }),
    }
}

fn summarize_structure(
    id: &str,
    structure: &MolecularStructure,
) -> ChemistryResult<MolecularSummary> {
    let summary = structure.summary()?;
    if !summary.molar_mass_grams.is_finite() || summary.molar_mass_grams <= 0.0 {
        return Err(ChemistryError::InvalidSubstance {
            substance_id: format!("destroy:{id}"),
            reason: "structure produced invalid molar mass".to_string(),
        });
    }
    Ok(summary)
}

pub const DESTROY_TAGS: &[&str] = &[
    "destroy:abundant_in_air",
    "destroy:acutely_toxic",
    "destroy:acid_rain",
    "destroy:adhesive",
    "destroy:bleach",
    "destroy:carcinogen",
    "destroy:explosive",
    "destroy:flame_retardant",
    "destroy:fragrant",
    "destroy:fuel_additive",
    "destroy:greenhouse",
    "destroy:hypothetical",
    "destroy:ozone_depleter",
    "destroy:plasticizer",
    "destroy:refrigerant",
    "destroy:smelly",
    "destroy:smog",
    "destroy:solvent",
];

const METALLURGY_IONS: &[RawMetallurgyIon] = &[
    RawMetallurgyIon {
        id: "aluminum_iii",
        element: Some("Al"),
        charge: 3,
        color_argb: 0x80C7D0D8,
    },
    RawMetallurgyIon {
        id: "magnesium_ion",
        element: Some("Mg"),
        charge: 2,
        color_argb: 0x80D6D6D6,
    },
    RawMetallurgyIon {
        id: "silicon_iv",
        element: Some("Si"),
        charge: 4,
        color_argb: 0x809098A0,
    },
    RawMetallurgyIon {
        id: "carbonate",
        element: None,
        charge: -2,
        color_argb: 0x80D8D0B8,
    },
    RawMetallurgyIon {
        id: "silicate",
        element: None,
        charge: -4,
        color_argb: 0x80B8B8A8,
    },
];

const METALLURGY_METALS: &[RawMetallurgyMetal] = &[
    RawMetallurgyMetal {
        id: "iron_metal",
        element: "Fe",
        solid_density: 7_874.0,
        melting_point_kelvin: 1811.0,
        boiling_point_kelvin: 3134.0,
        heat_capacity: 25.10,
        fusion_heat: 13_800.0,
        vaporization_heat: 340_000.0,
        color_argb: 0xFFB0A8A0,
    },
    RawMetallurgyMetal {
        id: "copper_metal",
        element: "Cu",
        solid_density: 8_960.0,
        melting_point_kelvin: 1357.77,
        boiling_point_kelvin: 2835.0,
        heat_capacity: 24.44,
        fusion_heat: 13_100.0,
        vaporization_heat: 300_000.0,
        color_argb: 0xFFC87533,
    },
    RawMetallurgyMetal {
        id: "zinc_metal",
        element: "Zn",
        solid_density: 7_140.0,
        melting_point_kelvin: 692.68,
        boiling_point_kelvin: 1180.0,
        heat_capacity: 25.47,
        fusion_heat: 7_320.0,
        vaporization_heat: 115_000.0,
        color_argb: 0xFFBFC4C7,
    },
    RawMetallurgyMetal {
        id: "nickel_metal",
        element: "Ni",
        solid_density: 8_908.0,
        melting_point_kelvin: 1728.0,
        boiling_point_kelvin: 3186.0,
        heat_capacity: 26.07,
        fusion_heat: 17_500.0,
        vaporization_heat: 378_000.0,
        color_argb: 0xFFC4BFA8,
    },
    RawMetallurgyMetal {
        id: "lead_metal",
        element: "Pb",
        solid_density: 11_340.0,
        melting_point_kelvin: 600.61,
        boiling_point_kelvin: 2022.0,
        heat_capacity: 26.65,
        fusion_heat: 4_770.0,
        vaporization_heat: 179_500.0,
        color_argb: 0xFF77777F,
    },
    RawMetallurgyMetal {
        id: "aluminum_metal",
        element: "Al",
        solid_density: 2_700.0,
        melting_point_kelvin: 933.47,
        boiling_point_kelvin: 2792.0,
        heat_capacity: 24.20,
        fusion_heat: 10_700.0,
        vaporization_heat: 294_000.0,
        color_argb: 0xFFD0D4D8,
    },
    RawMetallurgyMetal {
        id: "calcium_metal",
        element: "Ca",
        solid_density: 1_550.0,
        melting_point_kelvin: 1115.0,
        boiling_point_kelvin: 1757.0,
        heat_capacity: 25.93,
        fusion_heat: 8_540.0,
        vaporization_heat: 155_000.0,
        color_argb: 0xFFD8D8C8,
    },
    RawMetallurgyMetal {
        id: "magnesium_metal",
        element: "Mg",
        solid_density: 1_738.0,
        melting_point_kelvin: 923.0,
        boiling_point_kelvin: 1363.0,
        heat_capacity: 24.87,
        fusion_heat: 8_480.0,
        vaporization_heat: 128_000.0,
        color_argb: 0xFFD6D6D0,
    },
];

const METALLURGY_MATERIALS: &[RawMetallurgyMaterial] = &[
    RawMetallurgyMaterial {
        id: "iron_ii_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:iron_ii", 1), ("destroy:oxide", 1)],
        solid_density: 5_745.0,
        melting_point_kelvin: 1650.0,
        boiling_point_kelvin: 3_000.0,
        heat_capacity: 49.5,
        fusion_heat: 30_000.0,
        color_argb: 0xFF2E2E28,
    },
    RawMetallurgyMaterial {
        id: "iron_iii_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:iron_iii", 2), ("destroy:oxide", 3)],
        solid_density: 5_260.0,
        melting_point_kelvin: 1838.0,
        boiling_point_kelvin: 3_000.0,
        heat_capacity: 103.9,
        fusion_heat: 70_000.0,
        color_argb: 0xFF8E3329,
    },
    RawMetallurgyMaterial {
        id: "magnetite",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[
            ("destroy:iron_ii", 1),
            ("destroy:iron_iii", 2),
            ("destroy:oxide", 4),
        ],
        solid_density: 5_180.0,
        melting_point_kelvin: 1870.0,
        boiling_point_kelvin: 3_000.0,
        heat_capacity: 150.0,
        fusion_heat: 95_000.0,
        color_argb: 0xFF1C1C1F,
    },
    RawMetallurgyMaterial {
        id: "copper_i_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:copper_i", 2), ("destroy:oxide", 1)],
        solid_density: 6_000.0,
        melting_point_kelvin: 1508.0,
        boiling_point_kelvin: 2_070.0,
        heat_capacity: 64.0,
        fusion_heat: 32_000.0,
        color_argb: 0xFFB23B23,
    },
    RawMetallurgyMaterial {
        id: "copper_ii_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:copper_ii", 1), ("destroy:oxide", 1)],
        solid_density: 6_315.0,
        melting_point_kelvin: 1599.0,
        boiling_point_kelvin: 2_000.0,
        heat_capacity: 42.3,
        fusion_heat: 22_000.0,
        color_argb: 0xFF171717,
    },
    RawMetallurgyMaterial {
        id: "zinc_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:zinc_ion", 1), ("destroy:oxide", 1)],
        solid_density: 5_606.0,
        melting_point_kelvin: 2248.0,
        boiling_point_kelvin: 2_500.0,
        heat_capacity: 40.3,
        fusion_heat: 45_000.0,
        color_argb: 0xFFF0E8D8,
    },
    RawMetallurgyMaterial {
        id: "nickel_ii_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:nickel_ion", 1), ("destroy:oxide", 1)],
        solid_density: 6_670.0,
        melting_point_kelvin: 2220.0,
        boiling_point_kelvin: 3_000.0,
        heat_capacity: 44.0,
        fusion_heat: 44_000.0,
        color_argb: 0xFF48543A,
    },
    RawMetallurgyMaterial {
        id: "lead_ii_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:lead_ii", 1), ("destroy:oxide", 1)],
        solid_density: 9_530.0,
        melting_point_kelvin: 1_161.0,
        boiling_point_kelvin: 1_750.0,
        heat_capacity: 45.8,
        fusion_heat: 19_000.0,
        color_argb: 0xFFC28A2C,
    },
    RawMetallurgyMaterial {
        id: "aluminum_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:aluminum_iii", 2), ("destroy:oxide", 3)],
        solid_density: 3_950.0,
        melting_point_kelvin: 2_345.0,
        boiling_point_kelvin: 3_250.0,
        heat_capacity: 79.0,
        fusion_heat: 109_000.0,
        color_argb: 0xFFE6E3DC,
    },
    RawMetallurgyMaterial {
        id: "calcium_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:calcium_ion", 1), ("destroy:oxide", 1)],
        solid_density: 3_340.0,
        melting_point_kelvin: 2_886.0,
        boiling_point_kelvin: 3_120.0,
        heat_capacity: 42.0,
        fusion_heat: 63_000.0,
        color_argb: 0xFFE8E3D5,
    },
    RawMetallurgyMaterial {
        id: "magnesium_oxide",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:magnesium_ion", 1), ("destroy:oxide", 1)],
        solid_density: 3_580.0,
        melting_point_kelvin: 3_125.0,
        boiling_point_kelvin: 3_873.0,
        heat_capacity: 37.0,
        fusion_heat: 77_000.0,
        color_argb: 0xFFECEAE0,
    },
    RawMetallurgyMaterial {
        id: "silica",
        representation: RawMetallurgyMaterialRepresentation::Oxide,
        formula_units: &[("destroy:silicon_iv", 1), ("destroy:oxide", 2)],
        solid_density: 2_650.0,
        melting_point_kelvin: 1_986.0,
        boiling_point_kelvin: 3_223.0,
        heat_capacity: 44.0,
        fusion_heat: 9_000.0,
        color_argb: 0xFFE5E1D2,
    },
    RawMetallurgyMaterial {
        id: "iron_ii_sulfide",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:iron_ii", 1), ("destroy:sulfide", 1)],
        solid_density: 4_840.0,
        melting_point_kelvin: 1_461.0,
        boiling_point_kelvin: 2_500.0,
        heat_capacity: 50.0,
        fusion_heat: 30_000.0,
        color_argb: 0xFF222220,
    },
    RawMetallurgyMaterial {
        id: "copper_i_sulfide",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:copper_i", 2), ("destroy:sulfide", 1)],
        solid_density: 5_600.0,
        melting_point_kelvin: 1_403.0,
        boiling_point_kelvin: 2_200.0,
        heat_capacity: 76.0,
        fusion_heat: 35_000.0,
        color_argb: 0xFF2B2A24,
    },
    RawMetallurgyMaterial {
        id: "copper_ii_sulfide",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:copper_ii", 1), ("destroy:sulfide", 1)],
        solid_density: 4_760.0,
        melting_point_kelvin: 1_473.0,
        boiling_point_kelvin: 2_200.0,
        heat_capacity: 49.0,
        fusion_heat: 27_000.0,
        color_argb: 0xFF151512,
    },
    RawMetallurgyMaterial {
        id: "zinc_sulfide",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:zinc_ion", 1), ("destroy:sulfide", 1)],
        solid_density: 4_090.0,
        melting_point_kelvin: 2_123.0,
        boiling_point_kelvin: 2_800.0,
        heat_capacity: 48.0,
        fusion_heat: 42_000.0,
        color_argb: 0xFFE5DDB8,
    },
    RawMetallurgyMaterial {
        id: "lead_ii_sulfide",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:lead_ii", 1), ("destroy:sulfide", 1)],
        solid_density: 7_600.0,
        melting_point_kelvin: 1_387.0,
        boiling_point_kelvin: 2_500.0,
        heat_capacity: 49.5,
        fusion_heat: 23_000.0,
        color_argb: 0xFF202020,
    },
    RawMetallurgyMaterial {
        id: "nickel_ii_sulfide",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:nickel_ion", 1), ("destroy:sulfide", 1)],
        solid_density: 5_300.0,
        melting_point_kelvin: 1_573.0,
        boiling_point_kelvin: 2_500.0,
        heat_capacity: 50.0,
        fusion_heat: 33_000.0,
        color_argb: 0xFF2F3023,
    },
    RawMetallurgyMaterial {
        id: "calcium_carbonate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:calcium_ion", 1), ("destroy:carbonate", 1)],
        solid_density: 2_710.0,
        melting_point_kelvin: 1_612.0,
        boiling_point_kelvin: 2_200.0,
        heat_capacity: 81.9,
        fusion_heat: 37_000.0,
        color_argb: 0xFFE5E1D2,
    },
    RawMetallurgyMaterial {
        id: "magnesium_carbonate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:magnesium_ion", 1), ("destroy:carbonate", 1)],
        solid_density: 2_960.0,
        melting_point_kelvin: 813.0,
        boiling_point_kelvin: 1_600.0,
        heat_capacity: 75.0,
        fusion_heat: 30_000.0,
        color_argb: 0xFFE4E0D8,
    },
    RawMetallurgyMaterial {
        id: "zinc_carbonate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:zinc_ion", 1), ("destroy:carbonate", 1)],
        solid_density: 4_440.0,
        melting_point_kelvin: 573.0,
        boiling_point_kelvin: 1_600.0,
        heat_capacity: 82.0,
        fusion_heat: 28_000.0,
        color_argb: 0xFFE5E0C8,
    },
    RawMetallurgyMaterial {
        id: "iron_ii_carbonate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:iron_ii", 1), ("destroy:carbonate", 1)],
        solid_density: 3_960.0,
        melting_point_kelvin: 763.0,
        boiling_point_kelvin: 1_600.0,
        heat_capacity: 82.0,
        fusion_heat: 28_000.0,
        color_argb: 0xFFB8B08A,
    },
    RawMetallurgyMaterial {
        id: "copper_ii_carbonate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:copper_ii", 1), ("destroy:carbonate", 1)],
        solid_density: 4_000.0,
        melting_point_kelvin: 473.0,
        boiling_point_kelvin: 1_600.0,
        heat_capacity: 85.0,
        fusion_heat: 27_000.0,
        color_argb: 0xFF4F9C78,
    },
    RawMetallurgyMaterial {
        id: "lead_ii_carbonate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:lead_ii", 1), ("destroy:carbonate", 1)],
        solid_density: 6_600.0,
        melting_point_kelvin: 673.0,
        boiling_point_kelvin: 1_600.0,
        heat_capacity: 81.0,
        fusion_heat: 26_000.0,
        color_argb: 0xFFE6E0C4,
    },
    RawMetallurgyMaterial {
        id: "calcium_silicate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:calcium_ion", 2), ("destroy:silicate", 1)],
        solid_density: 2_900.0,
        melting_point_kelvin: 1_817.0,
        boiling_point_kelvin: 2_800.0,
        heat_capacity: 120.0,
        fusion_heat: 70_000.0,
        color_argb: 0xFFD8D5C0,
    },
    RawMetallurgyMaterial {
        id: "magnesium_silicate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:magnesium_ion", 2), ("destroy:silicate", 1)],
        solid_density: 3_200.0,
        melting_point_kelvin: 1_850.0,
        boiling_point_kelvin: 2_800.0,
        heat_capacity: 120.0,
        fusion_heat: 75_000.0,
        color_argb: 0xFFD8D8C8,
    },
    RawMetallurgyMaterial {
        id: "iron_ii_silicate",
        representation: RawMetallurgyMaterialRepresentation::IonicSolid,
        formula_units: &[("destroy:iron_ii", 2), ("destroy:silicate", 1)],
        solid_density: 4_300.0,
        melting_point_kelvin: 1_478.0,
        boiling_point_kelvin: 2_800.0,
        heat_capacity: 130.0,
        fusion_heat: 70_000.0,
        color_argb: 0xFF5D5844,
    },
];

const DESTROY_SUBSTANCES: &[RawSubstance] = &[
    RawSubstance {
        id: "acetamide",
        structure_code: Some(r#"destroy:linear:CC(=O)N"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(221.2),
        boiling_point_kelvin: None,
        density: Some(1159.0),
        molar_heat_capacity: Some(91.3),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "carcinogen", "plasticizer", "smog"],
    },
    RawSubstance {
        id: "acetate",
        structure_code: Some(r#"destroy:linear:CC~(~O^-0.5)O^-0.5"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "acetic_acid",
        structure_code: Some(r#"destroy:linear:CC(=O)OH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(118.5),
        boiling_point_kelvin: None,
        density: Some(1049.0),
        molar_heat_capacity: Some(123.1),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "acetic_anhydride",
        structure_code: Some(r#"destroy:linear:CC(=O)OC(=O)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(140.0),
        boiling_point_kelvin: None,
        density: Some(1082.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smelly", "smog"],
    },
    RawSubstance {
        id: "acetone",
        structure_code: Some(r#"destroy:linear:CC(=O)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(56.08),
        boiling_point_kelvin: None,
        density: Some(784.5),
        molar_heat_capacity: Some(126.3),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["solvent", "smelly", "smog"],
    },
    RawSubstance {
        id: "acetone_cyanohydrin",
        structure_code: Some(r#"destroy:linear:CC(OH)(C#N)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(95.0),
        boiling_point_kelvin: None,
        density: Some(932.0),
        molar_heat_capacity: Some(160.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smog"],
    },
    RawSubstance {
        id: "acetylene",
        structure_code: Some(r#"destroy:linear:C#C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-75.0),
        boiling_point_kelvin: None,
        density: Some(613.0),
        molar_heat_capacity: Some(44.036),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "acrylonitrile",
        structure_code: Some(r#"destroy:linear:C=CC#N"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(77.0),
        boiling_point_kelvin: None,
        density: Some(810.0),
        molar_heat_capacity: Some(113.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smelly", "smog"],
    },
    RawSubstance {
        id: "adipic_acid",
        structure_code: Some(r#"destroy:linear:O=C(OH)CCCCC(=O)OH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(337.5),
        boiling_point_kelvin: None,
        density: Some(1360.0),
        molar_heat_capacity: Some(196.5),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "adiponitrile",
        structure_code: Some(r#"destroy:linear:N#CCCCCC#N"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(295.1),
        boiling_point_kelvin: None,
        density: Some(951.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smog"],
    },
    RawSubstance {
        id: "aibn",
        structure_code: Some(r#"destroy:linear:CC(C)(C#N)N=NC(C)(C#N)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(10000.0),
        boiling_point_kelvin: None,
        density: Some(1100.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "ammonia",
        structure_code: Some(r#"destroy:linear:N"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-33.34),
        boiling_point_kelvin: None,
        density: Some(900.0),
        molar_heat_capacity: Some(80.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["refrigerant", "smelly"],
    },
    RawSubstance {
        id: "ammonium",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.NITROGEN, 1).addAtom(LegacyElement.HYDROGEN).addAtom(LegacyElement.HYDROGEN).addAtom(LegacyElement.HYDROGEN).addAtom(LegacyElement.HYDROGEN)"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "argon",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.ARGON)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: Some(87.302),
        density: Some(1395.4),
        molar_heat_capacity: Some(20.85),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "aspirin",
        structure_code: Some(r#"destroy:benzene:OC(=O)C,C(=O)O,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(140.0),
        boiling_point_kelvin: None,
        density: Some(1400.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "benzene",
        structure_code: Some(r#"destroy:benzene:,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(80.1),
        boiling_point_kelvin: None,
        density: Some(876.5),
        molar_heat_capacity: Some(134.8),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["carcinogen", "smog", "solvent"],
    },
    RawSubstance {
        id: "benzyl_chloride",
        structure_code: Some(r#"destroy:benzene:CCl,,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(179.0),
        boiling_point_kelvin: None,
        density: Some(1100.0),
        molar_heat_capacity: Some(182.4),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "ozone_depleter", "smog"],
    },
    RawSubstance {
        id: "boric_acid",
        structure_code: Some(r#"destroy:linear:OB(O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(300.0),
        boiling_point_kelvin: None,
        density: Some(1435.0),
        molar_heat_capacity: Some(81.3),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "flame_retardant"],
    },
    RawSubstance {
        id: "borohydride",
        structure_code: Some(r#"destroy:linear:HB^-1(H)(H)H"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "butadiene",
        structure_code: Some(r#"destroy:linear:C=CC=C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-4.41),
        boiling_point_kelvin: None,
        density: Some(614.9),
        molar_heat_capacity: Some(123.65),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "calcium_ion",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.CALCIUM, 2)"#),
        translation_key: Some(r#"calcium"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "carbon_dioxide",
        structure_code: Some(r#"destroy:linear:O=C=O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-78.4645),
        boiling_point_kelvin: None,
        density: Some(827.3),
        molar_heat_capacity: Some(37.135),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["greenhouse"],
    },
    RawSubstance {
        id: "carbon_monoxide",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.CARBON) .addAtom(LegacyElement.OXYGEN, BondType.TRIPLE)"#,
        ),
        translation_key: None,
        boiling_point_celsius: Some(-191.5),
        boiling_point_kelvin: None,
        density: Some(789.0),
        molar_heat_capacity: Some(29.1),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "greenhouse"],
    },
    RawSubstance {
        id: "carbon_tetrachloride",
        structure_code: Some(r#"destroy:linear:ClC(Cl)(Cl)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(76.72),
        boiling_point_kelvin: None,
        density: Some(1586.7),
        molar_heat_capacity: Some(132.6),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "carcinogen", "ozone_depleter"],
    },
    RawSubstance {
        id: "chloride",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.CHLORINE, -1)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "chlorine",
        structure_code: Some(r#"destroy:linear:ClCl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-34.04),
        boiling_point_kelvin: None,
        density: Some(1562.5),
        molar_heat_capacity: Some(33.949),
        specific_heat_capacity: None,
        color_argb: 0x20F9FCC2,
        tags: &["acutely_toxic", "ozone_depleter", "smelly"],
    },
    RawSubstance {
        id: "chloroaurate",
        structure_code: Some(r#"destroy:linear:ClAu^-1(Cl)(Cl)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0x7FEDCA4A,
        tags: &[],
    },
    RawSubstance {
        id: "chlorodifluoromethane",
        structure_code: Some(r#"destroy:linear:ClC(F)F"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-40.7),
        boiling_point_kelvin: None,
        density: Some(1186.8),
        molar_heat_capacity: Some(112.6),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["ozone_depleter", "refrigerant"],
    },
    RawSubstance {
        id: "chloroethane",
        structure_code: Some(r#"destroy:linear:CCCl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(12.27),
        boiling_point_kelvin: None,
        density: Some(889.8),
        molar_heat_capacity: Some(40.7),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["ozone_depleter", "carcinogen"],
    },
    RawSubstance {
        id: "chloroethene",
        structure_code: Some(r#"destroy:linear:C=CCl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-13.4),
        boiling_point_kelvin: None,
        density: Some(911.0),
        molar_heat_capacity: Some(85.92),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["carcinogen"],
    },
    RawSubstance {
        id: "chloromethane",
        structure_code: Some(r#"destroy:linear:CCl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-23.8),
        boiling_point_kelvin: None,
        density: Some(1003.0),
        molar_heat_capacity: Some(81.2),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "refrigerant"],
    },
    RawSubstance {
        id: "chloroform",
        structure_code: Some(r#"destroy:linear:ClC(Cl)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(61.15),
        boiling_point_kelvin: None,
        density: Some(1489.0),
        molar_heat_capacity: Some(114.25),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "ozone_depleter", "solvent"],
    },
    RawSubstance {
        id: "chromate",
        structure_code: Some(r#"destroy:linear:O=Cr=(-O^-1)(-O^-1)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0xD0F7ED2C,
        tags: &[],
    },
    RawSubstance {
        id: "chromium_iii",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.CHROMIUM, 3)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0xD00D9614,
        tags: &[],
    },
    RawSubstance {
        id: "cisplatin",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom( LegacyElement.PLATINUM) .addAtom(LegacyElement.CHLORINE) .addAtom(LegacyElement.CHLORINE) .addGroup( LegacyMolecularStructure.atom(LegacyElement.NITROGEN) .addAtom(LegacyElement.HYDROGEN) .addAtom(LegacyElement.HYDROGEN) .addAtom(LegacyElement.HYDROGEN), true ).addGroup( LegacyMolecularStructure.atom(LegacyElement.NITROGEN) .addAtom(LegacyElement.HYDROGEN) .addAtom(LegacyElement.HYDROGEN) .addAtom(LegacyElement.HYDROGEN), true )"#,
        ),
        translation_key: None,
        boiling_point_celsius: Some(270.0),
        boiling_point_kelvin: None,
        density: Some(3740.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "copper_i",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.COPPER, 1)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0xE0D30823,
        tags: &[],
    },
    RawSubstance {
        id: "copper_ii",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.COPPER, 2)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0xE00FFCA1,
        tags: &[],
    },
    RawSubstance {
        id: "creatine",
        structure_code: Some(r#"destroy:linear:NC(=N)N(C)CC(=O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: Some(171.1),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "cubane",
        structure_code: Some(r#"destroy:cubane:,,,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(161.6),
        boiling_point_kelvin: None,
        density: Some(1290.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "cubanedicarboxylic_acid",
        structure_code: Some(r#"destroy:cubane:C(=O)OH,,,,,,C(=O)OH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(457.4),
        boiling_point_kelvin: None,
        density: Some(2400.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "cyanamide",
        structure_code: Some(r#"destroy:linear:N#CN"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(260.0),
        boiling_point_kelvin: None,
        density: Some(1280.0),
        molar_heat_capacity: Some(78.2),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog", "carcinogen"],
    },
    RawSubstance {
        id: "cyanamide_ion",
        structure_code: Some(r#"destroy:linear:N^-1=C=N^-1"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "cyanide",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.CARBON, -1) .addAtom(LegacyElement.NITROGEN, BondType.TRIPLE)"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic"],
    },
    RawSubstance {
        id: "cyclohexene",
        structure_code: Some(r#"destroy:cyclohexene:,,,,,,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(82.98),
        boiling_point_kelvin: None,
        density: Some(811.0),
        molar_heat_capacity: Some(152.9),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smelly", "smog"],
    },
    RawSubstance {
        id: "cyclopentadienide",
        structure_code: Some(r#"destroy:cyclopentadienide:,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "diborane",
        structure_code: Some(r#"destroy:diborane:,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-92.49),
        boiling_point_kelvin: None,
        density: Some(1131.0),
        molar_heat_capacity: Some(56.7),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "carcinogen"],
    },
    RawSubstance {
        id: "dichlorodifluoromethane",
        structure_code: Some(r#"destroy:linear:ClC(F)(F)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-28.9),
        boiling_point_kelvin: None,
        density: Some(1486.0),
        molar_heat_capacity: Some(126.8),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["greenhouse", "ozone_depleter", "refrigerant"],
    },
    RawSubstance {
        id: "dichloromethane",
        structure_code: Some(r#"destroy:linear:ClCCl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(39.6),
        boiling_point_kelvin: None,
        density: Some(1326.6),
        molar_heat_capacity: Some(81.2),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "refrigerant", "carcinogen", "solvent"],
    },
    RawSubstance {
        id: "dichromate",
        structure_code: Some(r#"destroy:linear:O=Cr(=O)(O^-1)OCr=(=O)(O^-1)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0xD0DB3D0D,
        tags: &[],
    },
    RawSubstance {
        id: "dinitrotoluene",
        structure_code: Some(r#"destroy:benzene:C,N~(~O)O,,N~(~O)O,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(250.0),
        boiling_point_kelvin: None,
        density: Some(1520.0),
        molar_heat_capacity: Some(243.3),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "carcinogen", "plasticizer"],
    },
    RawSubstance {
        id: "ethanol",
        structure_code: Some(r#"destroy:linear:CCO"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(78.23),
        boiling_point_kelvin: None,
        density: Some(789.45),
        molar_heat_capacity: Some(109.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["solvent", "smog"],
    },
    RawSubstance {
        id: "ethene",
        structure_code: Some(r#"destroy:linear:C=C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-103.7),
        boiling_point_kelvin: None,
        density: Some(567.9),
        molar_heat_capacity: Some(67.3),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "ethylanthraquinone",
        structure_code: Some(r#"destroy:anthraquinone:,,,O,,,CC,,O,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(415.4),
        boiling_point_kelvin: None,
        density: Some(1231.0),
        molar_heat_capacity: Some(286.6),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "ethylanthrahydroquinone",
        structure_code: Some(r#"destroy:anthracene:,,,O,,,CC,,O,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(415.4),
        boiling_point_kelvin: None,
        density: Some(1231.0),
        molar_heat_capacity: Some(286.6),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "ethylbenzene",
        structure_code: Some(r#"destroy:benzene:CC,,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(136.0),
        boiling_point_kelvin: None,
        density: Some(866.5),
        molar_heat_capacity: None,
        specific_heat_capacity: Some(1726.0),
        color_argb: 0,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "ethoxide",
        structure_code: Some(r#"destroy:linear:CCO^-1"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "fluoride",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.FLUORINE, -1)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_acid_anhydride",
        structure_code: Some(r#"destroy:linear:RC(=O)OC(=O)R"#),
        java_structure_code: None,
        translation_key: Some(r#"acid_anhydride"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_acyl_chloride",
        structure_code: Some(r#"destroy:linear:RC(=O)Cl"#),
        java_structure_code: None,
        translation_key: Some(r#"acyl_chloride"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_alcohol",
        structure_code: Some(r#"destroy:linear:RC(R)(R)O"#),
        java_structure_code: None,
        translation_key: Some(r#"alcohol"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_alkene",
        structure_code: Some(r#"destroy:linear:RC=(R)C(R)R"#),
        java_structure_code: None,
        translation_key: Some(r#"alkene"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_alkoxide",
        structure_code: Some(r#"destroy:linear:RC(R)(R)O^-1"#),
        java_structure_code: None,
        translation_key: Some(r#"alkoxide"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_alkyne",
        structure_code: Some(r#"destroy:linear:RC#CR"#),
        java_structure_code: None,
        translation_key: Some(r#"alkyne"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_amide",
        structure_code: Some(r#"destroy:linear:RC(=O)N"#),
        java_structure_code: None,
        translation_key: Some(r#"amide"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_amine",
        structure_code: Some(r#"destroy:linear:RC(R)(R)N"#),
        java_structure_code: None,
        translation_key: Some(r#"amine"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_borane",
        structure_code: Some(r#"destroy:linear:RC(R)(R)B(R)R"#),
        java_structure_code: None,
        translation_key: Some(r#"organic_borane"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_borate_ester",
        structure_code: Some(r#"destroy:linear:RC(R)(R)OB(R)R"#),
        java_structure_code: None,
        translation_key: Some(r#"borate_ester"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_carbonyl",
        structure_code: Some(r#"destroy:linear:RC=(R)O"#),
        java_structure_code: None,
        translation_key: Some(r#"carbonyl"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_carboxylic_acid",
        structure_code: Some(r#"destroy:linear:RC(=O)O"#),
        java_structure_code: None,
        translation_key: Some(r#"carboxylic_acid"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_chloride",
        structure_code: Some(r#"destroy:linear:RC(R)(R)Cl"#),
        java_structure_code: None,
        translation_key: Some(r#"organic_chloride"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_ester",
        structure_code: Some(r#"destroy:linear:RC(=O)OR"#),
        java_structure_code: None,
        translation_key: Some(r#"ester"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_isocyanate",
        structure_code: Some(r#"destroy:linear:RC(R)(R)N=C=O"#),
        java_structure_code: None,
        translation_key: Some(r#"isocyanate"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_nitrile",
        structure_code: Some(r#"destroy:linear:RC(R)(R)C#N"#),
        java_structure_code: None,
        translation_key: Some(r#"nitrile"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_nitro",
        structure_code: Some(r#"destroy:linear:RC(R)(R)N~(~O)O"#),
        java_structure_code: None,
        translation_key: Some(r#"nitro"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_organic_boric_acid",
        structure_code: Some(r#"destroy:linear:RB(R)OH"#),
        java_structure_code: None,
        translation_key: Some(r#"organic_boric_acid"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_primary_amine",
        structure_code: Some(r#"destroy:linear:RC(R)(R)N"#),
        java_structure_code: None,
        translation_key: Some(r#"non_tertiary_amine"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "generic_primary_borane",
        structure_code: Some(r#"destroy:linear:RC(R)(R)B"#),
        java_structure_code: None,
        translation_key: Some(r#"non_tertiary_borane"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "glycerol",
        structure_code: Some(r#"destroy:linear:OCC(O)CO"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(290.0),
        boiling_point_kelvin: None,
        density: Some(1261.0),
        molar_heat_capacity: Some(213.8),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "hexane_diisocyanate",
        structure_code: Some(r#"destroy:linear:O=C=NCCCCCCN=C=O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(255.0),
        boiling_point_kelvin: None,
        density: Some(1047.0),
        molar_heat_capacity: Some(222.7),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic"],
    },
    RawSubstance {
        id: "hexanediamine",
        structure_code: Some(r#"destroy:linear:NCCCCCCN"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(204.6),
        boiling_point_kelvin: None,
        density: Some(840.0),
        molar_heat_capacity: Some(250.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smelly"],
    },
    RawSubstance {
        id: "hydrazine",
        structure_code: Some(r#"destroy:linear:NN"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(114.0),
        boiling_point_kelvin: None,
        density: Some(1021.0),
        molar_heat_capacity: Some(70.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "carcinogen", "smelly", "smog"],
    },
    RawSubstance {
        id: "hydrochloric_acid",
        structure_code: Some(r#"destroy:linear:ClH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-85.05),
        boiling_point_kelvin: None,
        density: Some(1490.0),
        molar_heat_capacity: None,
        specific_heat_capacity: Some(798.1),
        color_argb: 0,
        tags: &["acutely_toxic", "ozone_depleter"],
    },
    RawSubstance {
        id: "hydrofluoric_acid",
        structure_code: Some(r#"destroy:linear:FH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(19.5),
        boiling_point_kelvin: None,
        density: Some(990.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "greenhouse"],
    },
    RawSubstance {
        id: "hydrogen",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.HYDROGEN).addAtom(LegacyElement.HYDROGEN)"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: Some(20.271),
        density: Some(70.85),
        molar_heat_capacity: Some(28.84),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "hydrogen_cyanide",
        structure_code: Some(r#"destroy:linear:N#C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(26.0),
        boiling_point_kelvin: None,
        density: Some(687.6),
        molar_heat_capacity: Some(35.9),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic"],
    },
    RawSubstance {
        id: "hydrogen_iodide",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.IODINE).addAtom(LegacyElement.HYDROGEN)"#,
        ),
        translation_key: None,
        boiling_point_celsius: Some(-35.36),
        boiling_point_kelvin: None,
        density: Some(2850.0),
        molar_heat_capacity: Some(29.2),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "hydrogen_peroxide",
        structure_code: Some(r#"destroy:linear:OO"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(150.2),
        boiling_point_kelvin: None,
        density: Some(1110.0),
        molar_heat_capacity: None,
        specific_heat_capacity: Some(2619.0),
        color_argb: 0x40C7F4FC,
        tags: &["acutely_toxic", "bleach"],
    },
    RawSubstance {
        id: "hydrogensulfate",
        structure_code: Some(r#"destroy:linear:O=S(=O)(OH)O^-1"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acid_rain"],
    },
    RawSubstance {
        id: "hydroxide",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.HYDROGEN).addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1))"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: Some(900.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "hypochlorous_acid",
        structure_code: Some(r#"destroy:linear:OCl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "bleach", "ozone_depleter", "bleach"],
    },
    RawSubstance {
        id: "hypochlorite",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.OXYGEN, -1).addAtom(LegacyElement.CHLORINE)"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "ozone_depleter"],
    },
    RawSubstance {
        id: "isoprene",
        structure_code: Some(r#"destroy:linear:C=C(C)C=C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(34.0),
        boiling_point_kelvin: None,
        density: Some(681.0),
        molar_heat_capacity: Some(102.69),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "iodide",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.IODINE, -1)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "iodine",
        structure_code: Some(r#"destroy:linear:II"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(184.3),
        boiling_point_kelvin: None,
        density: Some(4933.0),
        molar_heat_capacity: Some(54.44),
        specific_heat_capacity: None,
        color_argb: 0x80AA16A5,
        tags: &[],
    },
    RawSubstance {
        id: "iodomethane",
        structure_code: Some(r#"destroy:linear:CI"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(42.6),
        boiling_point_kelvin: None,
        density: Some(2280.0),
        molar_heat_capacity: Some(82.75),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smog"],
    },
    RawSubstance {
        id: "iron_ii",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.IRON, 2)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0x80A9C92A,
        tags: &[],
    },
    RawSubstance {
        id: "iron_iii",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.IRON, 3)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0x80F94939,
        tags: &[],
    },
    RawSubstance {
        id: "isopropanol",
        structure_code: Some(r#"destroy:linear:CC(O)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(82.6),
        boiling_point_kelvin: None,
        density: Some(786.0),
        molar_heat_capacity: None,
        specific_heat_capacity: Some(1.54),
        color_argb: 0,
        tags: &["solvent"],
    },
    RawSubstance {
        id: "lead_ii",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.LEAD, 2)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "mercury",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.MERCURY)"#),
        translation_key: None,
        boiling_point_celsius: Some(356.73),
        boiling_point_kelvin: None,
        density: Some(13534.0),
        molar_heat_capacity: Some(27.98),
        specific_heat_capacity: None,
        color_argb: 0xFFB3B3B3,
        tags: &["acutely_toxic"],
    },
    RawSubstance {
        id: "metaxylene",
        structure_code: Some(r#"destroy:benzene:C,,C,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(139.0),
        boiling_point_kelvin: None,
        density: Some(860.0),
        molar_heat_capacity: Some(181.5),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "methane",
        structure_code: Some(r#"destroy:linear:C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-161.5),
        boiling_point_kelvin: None,
        density: Some(424.0),
        molar_heat_capacity: Some(35.7),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["greenhouse", "smog"],
    },
    RawSubstance {
        id: "methanol",
        structure_code: Some(r#"destroy:linear:CO"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(65.0),
        boiling_point_kelvin: None,
        density: Some(792.0),
        molar_heat_capacity: Some(68.62),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smog", "solvent"],
    },
    RawSubstance {
        id: "methylamine",
        structure_code: Some(r#"destroy:linear:CN"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-6.3),
        boiling_point_kelvin: None,
        density: Some(656.2),
        molar_heat_capacity: Some(101.8),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "methyl_acetate",
        structure_code: Some(r#"destroy:linear:CC(=O)OC"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(56.9),
        boiling_point_kelvin: None,
        density: Some(932.0),
        molar_heat_capacity: Some(140.2),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smog", "solvent"],
    },
    RawSubstance {
        id: "methyl_methacrylate",
        structure_code: Some(r#"destroy:linear:CC(=C)C(=O)OC"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(101.0),
        boiling_point_kelvin: None,
        density: Some(940.0),
        molar_heat_capacity: Some(191.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "methyl_salicylate",
        structure_code: Some(r#"destroy:benzene:C(=O)OC,O,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(222.0),
        boiling_point_kelvin: None,
        density: Some(1174.0),
        molar_heat_capacity: Some(248.9),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["fragrant", "smog"],
    },
    RawSubstance {
        id: "nickel_ion",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.NICKEL, 2)"#),
        translation_key: Some(r#"nickel"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0xE09BEAAB,
        tags: &["carcinogen"],
    },
    RawSubstance {
        id: "nitrate",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.NITROGEN, 1) .addAtom(LegacyElement.OXYGEN, BondType.DOUBLE) .addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1)) .addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1))"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acid_rain"],
    },
    RawSubstance {
        id: "nitric_acid",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.NITROGEN) .addAtom(LegacyElement.OXYGEN, BondType.AROMATIC) .addAtom(LegacyElement.OXYGEN, BondType.AROMATIC) .addGroup(LegacyMolecularStructure.atom(LegacyElement.OXYGEN) .addAtom(LegacyElement.HYDROGEN) ) //TODO maybe add color (though this should come from a decomposition)"#,
        ),
        translation_key: None,
        boiling_point_celsius: Some(83.0),
        boiling_point_kelvin: None,
        density: Some(1510.0),
        molar_heat_capacity: Some(53.29),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "acid_rain"],
    },
    RawSubstance {
        id: "nitrogen",
        structure_code: Some(r#"destroy:linear:N#N"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: Some(77.355),
        density: Some(807.0),
        molar_heat_capacity: Some(29.12),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["abundant_in_air"],
    },
    RawSubstance {
        id: "nitrogen_dioxide",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.NITROGEN) .addAtom(LegacyElement.OXYGEN, BondType.AROMATIC) .addAtom(LegacyElement.OXYGEN, BondType.AROMATIC)"#,
        ),
        translation_key: None,
        boiling_point_celsius: Some(21.15),
        boiling_point_kelvin: None,
        density: Some(1880.0),
        molar_heat_capacity: Some(37.2),
        specific_heat_capacity: None,
        color_argb: 0xD089011A,
        tags: &["acid_rain", "acutely_toxic", "carcinogen"],
    },
    RawSubstance {
        id: "nitroglycerine",
        structure_code: Some(r#"destroy:linear:C(ON(~O)(~O))C(ON(~O)(~O))CON(~O)(~O)"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(50.0),
        boiling_point_kelvin: None,
        density: Some(1600.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "nitronium",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.NITROGEN, 1) .addAtom(LegacyElement.OXYGEN, BondType.DOUBLE) .addAtom(LegacyElement.OXYGEN, BondType.DOUBLE)"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "octasulfur",
        structure_code: Some(r#"destroy:octasulfur:hello"#),
        java_structure_code: None,
        translation_key: Some(r#"sulfur"#),
        boiling_point_celsius: Some(444.6),
        boiling_point_kelvin: None,
        density: Some(2070.0),
        molar_heat_capacity: Some(21.64),
        specific_heat_capacity: None,
        color_argb: 0xFFD00000,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "oleum",
        structure_code: Some(r#"destroy:linear:HOS(=O)(=O)OS(=O)(=O)OH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(10.0),
        boiling_point_kelvin: None,
        density: Some(1820.0),
        molar_heat_capacity: None,
        specific_heat_capacity: Some(2600.0),
        color_argb: 0,
        tags: &["acid_rain", "acutely_toxic"],
    },
    RawSubstance {
        id: "orthoxylene",
        structure_code: Some(r#"destroy:benzene:C,C,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(144.4),
        boiling_point_kelvin: None,
        density: Some(880.0),
        molar_heat_capacity: Some(187.1),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "oxide",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.OXYGEN, -2)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "oxygen",
        structure_code: Some(r#"destroy:linear:O=O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: Some(90.188),
        density: Some(1141.0),
        molar_heat_capacity: Some(29.378),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["abundant_in_air"],
    },
    RawSubstance {
        id: "paraxylene",
        structure_code: Some(r#"destroy:benzene:C,,,C,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(138.35),
        boiling_point_kelvin: None,
        density: Some(861.0),
        molar_heat_capacity: Some(182.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "phenol",
        structure_code: Some(r#"destroy:benzene:O,,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(181.7),
        boiling_point_kelvin: None,
        density: Some(1070.0),
        molar_heat_capacity: Some(220.9),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smelly", "smog"],
    },
    RawSubstance {
        id: "phenylacetic_acid",
        structure_code: Some(r#"destroy:benzene:CC(=O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(265.5),
        boiling_point_kelvin: None,
        density: Some(1080.9),
        molar_heat_capacity: Some(170.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["fragrant", "smog"],
    },
    RawSubstance {
        id: "phenylacetone",
        structure_code: Some(r#"destroy:benzene:CC(=O)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(215.0),
        boiling_point_kelvin: None,
        density: Some(1006.0),
        molar_heat_capacity: Some(250.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smog"],
    },
    RawSubstance {
        id: "phosgene",
        structure_code: Some(r#"destroy:linear:ClC(=O)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(8.3),
        boiling_point_kelvin: None,
        density: Some(1432.0),
        molar_heat_capacity: Some(101.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "ozone_depleter", "smelly"],
    },
    RawSubstance {
        id: "trimethylphosphine",
        structure_code: Some(r#"destroy:linear:CP(C)C"#),
        java_structure_code: None,
        translation_key: Some("destroy.chemical.trimethylphosphine"),
        boiling_point_celsius: Some(38.0),
        boiling_point_kelvin: None,
        density: Some(750.0),
        molar_heat_capacity: Some(180.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smelly"],
    },
    RawSubstance {
        id: "phthalic_anhydride",
        structure_code: Some(r#"destroy:isohydrobenzofuran:,,,O,O,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(295.0),
        boiling_point_kelvin: None,
        density: Some(1530.0),
        molar_heat_capacity: Some(161.8),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "picric_acid",
        structure_code: Some(r#"destroy:benzene:O,N~(~O)O,,N~(~O)O,,N~(~O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(300.0),
        boiling_point_kelvin: None,
        density: Some(1763.0),
        molar_heat_capacity: Some(200.0),
        specific_heat_capacity: None,
        color_argb: 0xC0ED7417,
        tags: &["smog"],
    },
    RawSubstance {
        id: "potassium_ion",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.POTASSIUM, 1)"#),
        translation_key: Some(r#"potassium"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "propene",
        structure_code: Some(r#"destroy:linear:C=CC"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-47.6),
        boiling_point_kelvin: None,
        density: Some(1810.0),
        molar_heat_capacity: Some(102.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "salicylic_acid",
        structure_code: Some(r#"destroy:benzene:C(=O)O,O,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(200.0),
        boiling_point_kelvin: None,
        density: Some(1443.0),
        molar_heat_capacity: Some(159.4),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "sodium_metal",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.SODIUM)"#),
        translation_key: Some(r#"sodium"#),
        boiling_point_celsius: Some(882.94),
        boiling_point_kelvin: None,
        density: Some(968.0),
        molar_heat_capacity: Some(28.23),
        specific_heat_capacity: None,
        color_argb: 0xFFB3B3B3,
        tags: &[],
    },
    RawSubstance {
        id: "sodium_ion",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.SODIUM, 1)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: Some(900.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "styrene",
        structure_code: Some(r#"destroy:benzene:C=C,,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(145.0),
        boiling_point_kelvin: None,
        density: Some(909.0),
        molar_heat_capacity: Some(183.2),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smog"],
    },
    RawSubstance {
        id: "sulfate",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.SULFUR, 2) .addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1)) .addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1)) .addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1)) .addAtom(new LegacyAtom(LegacyElement.OXYGEN, -1))"#,
        ),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acid_rain"],
    },
    RawSubstance {
        id: "sulfide",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.SULFUR, -2)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "sulfur_dioxide",
        structure_code: Some(r#"destroy:linear:O=S=O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-10.0),
        boiling_point_kelvin: None,
        density: Some(2628.8),
        molar_heat_capacity: Some(39.87),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "acid_rain"],
    },
    RawSubstance {
        id: "sulfuric_acid",
        structure_code: Some(r#"destroy:linear:OS(=O)(=O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(337.0),
        boiling_point_kelvin: None,
        density: Some(1830.2),
        molar_heat_capacity: Some(83.68),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "acid_rain"],
    },
    RawSubstance {
        id: "sulfur_trioxide",
        structure_code: Some(r#"destroy:linear:S=(=O)(=O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(45.0),
        boiling_point_kelvin: None,
        density: Some(1920.0),
        molar_heat_capacity: None,
        specific_heat_capacity: Some(50.6),
        color_argb: 0,
        tags: &["acutely_toxic", "acid_rain"],
    },
    RawSubstance {
        id: "tetraethyllead",
        structure_code: Some(r#"destroy:linear:CCPb(CC)(CC)CC"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(84.5),
        boiling_point_kelvin: None,
        density: Some(1653.0),
        molar_heat_capacity: Some(307.4),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["fuel_additive", "smog", "carcinogen", "acutely_toxic"],
    },
    RawSubstance {
        id: "tetrafluoroethene",
        structure_code: Some(r#"destroy:linear:FC=(F)C(F)F"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(-76.3),
        boiling_point_kelvin: None,
        density: Some(1519.0),
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "greenhouse", "smog"],
    },
    RawSubstance {
        id: "tetrahydroxyborate",
        structure_code: Some(r#"destroy:linear:OB^-1(O)(O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "tetrahydroxy_tetraborate",
        structure_code: Some(r#"destroy:tetraborate:O,O,O,O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "thionyl_chloride",
        structure_code: Some(r#"destroy:linear:S=(=O)(Cl)(Cl)"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(74.6),
        boiling_point_kelvin: None,
        density: Some(1638.0),
        molar_heat_capacity: Some(121.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "acutely_toxic"],
    },
    RawSubstance {
        id: "toluene",
        structure_code: Some(r#"destroy:benzene:C,,,,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(110.6),
        boiling_point_kelvin: None,
        density: Some(862.3),
        molar_heat_capacity: Some(157.3),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["carcinogen", "smog", "solvent"],
    },
    RawSubstance {
        id: "toluene_diisocyanate",
        structure_code: Some(r#"destroy:benzene:C,N=C=O,,N=C=O,,"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(251.0),
        boiling_point_kelvin: None,
        density: Some(1214.0),
        molar_heat_capacity: Some(222.7),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["carcinogen", "acutely_toxic"],
    },
    RawSubstance {
        id: "trichlorofluoromethane",
        structure_code: Some(r#"destroy:linear:FC(Cl)(Cl)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(23.77),
        boiling_point_kelvin: None,
        density: Some(1494.0),
        molar_heat_capacity: Some(122.5),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["greenhouse", "ozone_depleter", "refrigerant", "smog"],
    },
    RawSubstance {
        id: "trimethylamine",
        structure_code: Some(r#"destroy:linear:CN(C)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(5.0),
        boiling_point_kelvin: None,
        density: Some(670.0),
        molar_heat_capacity: Some(132.55),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "smog"],
    },
    RawSubstance {
        id: "trimethyl_borate",
        structure_code: Some(r#"destroy:linear:COB(OC)OC"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(68.5),
        boiling_point_kelvin: None,
        density: Some(932.0),
        molar_heat_capacity: Some(189.9),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "tnt",
        structure_code: Some(r#"destroy:benzene:C,N~(~O)O,,N~(~O)O,,N~(~O)O"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(240.0),
        boiling_point_kelvin: None,
        density: Some(1654.0),
        molar_heat_capacity: Some(243.3),
        specific_heat_capacity: None,
        color_argb: 0xD0FCF1E8,
        tags: &["acutely_toxic", "carcinogen", "smog"],
    },
    RawSubstance {
        id: "vinyl_acetate",
        structure_code: Some(r#"destroy:linear:C=COC(=O)C"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(72.7),
        boiling_point_kelvin: None,
        density: Some(934.0),
        molar_heat_capacity: Some(169.5),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["adhesive", "carcinogen", "fragrant", "smog"],
    },
    RawSubstance {
        id: "water",
        structure_code: Some(r#"destroy:linear:HOH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(100.0),
        boiling_point_kelvin: None,
        density: Some(1000.0),
        molar_heat_capacity: None,
        specific_heat_capacity: Some(4160.0),
        color_argb: 0,
        tags: &["greenhouse", "solvent"],
    },
    RawSubstance {
        id: "zinc_ion",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.ZINC, 2)"#),
        translation_key: Some(r#"zinc"#),
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "proton",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.HYDROGEN, 1)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "bromine",
        structure_code: Some(r#"destroy:linear:BrBr"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(58.8),
        boiling_point_kelvin: None,
        density: Some(3102.8),
        molar_heat_capacity: Some(75.69),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "hydrobromic_acid",
        structure_code: Some(r#"destroy:linear:BrH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(122.0),
        boiling_point_kelvin: None,
        density: Some(1490.0),
        molar_heat_capacity: Some(100.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "hydroiodic_acid",
        structure_code: Some(r#"destroy:linear:IH"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(127.0),
        boiling_point_kelvin: None,
        density: Some(1700.0),
        molar_heat_capacity: Some(100.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "ferric_chloride",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.IRON, 3)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "ferric_bromide",
        structure_code: None,
        java_structure_code: Some(r#"LegacyMolecularStructure.atom(LegacyElement.IRON, 3)"#),
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "aluminum_trichloride",
        structure_code: Some(r#"destroy:linear:Al(Cl)(Cl)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: None,
        boiling_point_kelvin: None,
        density: None,
        molar_heat_capacity: None,
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    RawSubstance {
        id: "acetyl_chloride",
        structure_code: Some(r#"destroy:linear:CC(=O)Cl"#),
        java_structure_code: None,
        translation_key: None,
        boiling_point_celsius: Some(52.0),
        boiling_point_kelvin: None,
        density: Some(1104.0),
        molar_heat_capacity: Some(100.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &[],
    },
    // Protecting group reagents
    // TMS-Cl: Me3SiCl - Si with 3 methyl groups and Cl
    RawSubstance {
        id: "trimethylsilyl_chloride",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.SILICON).addAtom(LegacyElement.CARBON).addAtom(LegacyElement.CARBON).addAtom(LegacyElement.CARBON).addAtom(LegacyElement.CHLORINE)"#,
        ),
        translation_key: Some(r#"trimethylsilyl_chloride"#),
        boiling_point_celsius: Some(57.0),
        boiling_point_kelvin: None,
        density: Some(856.0),
        molar_heat_capacity: Some(150.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly", "acutely_toxic"],
    },
    // TMS-F: Me3SiF - Si with 3 methyl groups and F
    RawSubstance {
        id: "trimethylsilyl_fluoride",
        structure_code: None,
        java_structure_code: Some(
            r#"LegacyMolecularStructure.atom(LegacyElement.SILICON).addAtom(LegacyElement.CARBON).addAtom(LegacyElement.CARBON).addAtom(LegacyElement.CARBON).addAtom(LegacyElement.FLUORINE)"#,
        ),
        translation_key: Some(r#"trimethylsilyl_fluoride"#),
        boiling_point_celsius: Some(16.0),
        boiling_point_kelvin: None,
        density: Some(793.0),
        molar_heat_capacity: Some(140.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic"],
    },
    RawSubstance {
        id: "di_tert_butyl_dicarbonate",
        structure_code: Some(r#"destroy:linear:CC(C)(C)OC(=O)OC(=O)OC(C)(C)C"#),
        java_structure_code: None,
        translation_key: Some(r#"di_tert_butyl_dicarbonate"#),
        boiling_point_celsius: Some(56.0),
        boiling_point_kelvin: None,
        density: Some(950.0),
        molar_heat_capacity: Some(280.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["smelly"],
    },
    RawSubstance {
        id: "benzyl_chloroformate",
        structure_code: Some(r#"destroy:benzene:COC(=O)Cl,,,,,"#),
        java_structure_code: None,
        translation_key: Some(r#"benzyl_chloroformate"#),
        boiling_point_celsius: Some(103.0),
        boiling_point_kelvin: None,
        density: Some(1212.0),
        molar_heat_capacity: Some(200.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["acutely_toxic", "smelly"],
    },
    RawSubstance {
        id: "tert_butanol",
        structure_code: Some(r#"destroy:linear:CC(C)(C)O"#),
        java_structure_code: None,
        translation_key: Some(r#"tert_butanol"#),
        boiling_point_celsius: Some(82.4),
        boiling_point_kelvin: None,
        density: Some(781.0),
        molar_heat_capacity: Some(220.0),
        specific_heat_capacity: None,
        color_argb: 0,
        tags: &["solvent"],
    },
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn destroy_catalog_builds_all_explicit_substances() {
        let registry = destroy_substances_registry_builder()
            .unwrap()
            .build()
            .unwrap();

        assert_eq!(DESTROY_SUBSTANCE_COUNT, 165);
        assert_eq!(registry.substance_count(), 211);
    }

    #[test]
    fn key_destroy_substances_match_legacy_properties() {
        let registry = destroy_substances_registry_builder()
            .unwrap()
            .build()
            .unwrap();

        let water = registry.substance(&"destroy:water".into()).unwrap();
        assert_eq!(water.charge, 0);
        assert!((water.molar_mass_grams - 18.02).abs() < 0.001);
        assert!((water.boiling_point_kelvin - 373.0).abs() < 0.001);
        assert!((water.molar_heat_capacity_j_per_mol_kelvin - (4160.0 / 18.02)).abs() < 0.001);

        let proton = registry.substance(&"destroy:proton".into()).unwrap();
        assert_eq!(proton.charge, 1);
        assert!((proton.molar_mass_grams - 1.01).abs() < 0.001);
        assert_eq!(proton.boiling_point_kelvin, f64::MAX);

        let hydroxide = registry.substance(&"destroy:hydroxide".into()).unwrap();
        assert_eq!(hydroxide.charge, -1);
        assert!((hydroxide.molar_mass_grams - 17.01).abs() < 0.001);
        assert_eq!(hydroxide.boiling_point_kelvin, f64::MAX);

        let ammonium = registry.substance(&"destroy:ammonium".into()).unwrap();
        assert_eq!(ammonium.charge, 1);
        assert!((ammonium.molar_mass_grams - 18.05).abs() < 0.001);

        let acetate = registry.substance(&"destroy:acetate".into()).unwrap();
        assert_eq!(acetate.charge, -1);
        assert!((acetate.molar_mass_grams - 59.04).abs() < 0.01);

        let chloride = registry.substance(&"destroy:chloride".into()).unwrap();
        assert_eq!(chloride.charge, -1);
        assert!((chloride.molar_mass_grams - 35.45).abs() < 0.001);

        let argon = registry.substance(&"destroy:argon".into()).unwrap();
        assert_eq!(argon.charge, 0);
        assert!((argon.molar_mass_grams - 39.95).abs() < 0.001);
        assert!((argon.boiling_point_kelvin - 87.302).abs() < 0.001);
    }

    #[test]
    fn metallurgy_substances_have_material_representations() {
        let registry = destroy_substances_registry_builder()
            .unwrap()
            .build()
            .unwrap();

        let iron = registry.substance(&"destroy:iron_metal".into()).unwrap();
        assert!(matches!(
            &iron.representation,
            SubstanceRepresentation::Metal { element_symbol } if element_symbol == "Fe"
        ));
        assert_eq!(iron.charge, 0);
        assert!(iron.phase_properties.can_precipitate);

        let hematite = registry
            .substance(&"destroy:iron_iii_oxide".into())
            .unwrap();
        assert!(matches!(
            &hematite.representation,
            SubstanceRepresentation::Oxide { formula_units }
                if formula_units
                    == &vec![
                        MaterialFormulaUnit::new("destroy:iron_iii", 2),
                        MaterialFormulaUnit::new("destroy:oxide", 3),
                    ]
        ));
        assert!(
            (hematite.molar_mass_grams
                - (element_mass("Fe").unwrap() * 2.0 + element_mass("O").unwrap() * 3.0))
                .abs()
                < 1.0e-9
        );

        let calcium_silicate = registry
            .substance(&"destroy:calcium_silicate".into())
            .unwrap();
        assert!(matches!(
            &calcium_silicate.representation,
            SubstanceRepresentation::IonicSolid { formula_units }
                if formula_units
                    == &vec![
                        MaterialFormulaUnit::new("destroy:calcium_ion", 2),
                        MaterialFormulaUnit::new("destroy:silicate", 1),
                    ]
        ));
        assert_eq!(calcium_silicate.charge, 0);
    }

    #[test]
    fn legacy_structure_rejects_unknown_element_and_bad_charge() {
        assert!(summarize_legacy_structure("destroy:linear:Xx").is_err());
        assert!(summarize_legacy_structure("destroy:linear:O^-nope").is_err());
    }

    #[test]
    fn destroy_catalog_populates_functional_groups() {
        let registry = destroy_substances_registry_builder()
            .unwrap()
            .build()
            .unwrap();

        let acetic_acid = registry.substance(&"destroy:acetic_acid".into()).unwrap();
        assert!(acetic_acid.functional_groups.iter().any(|group| {
            group.group_type == super::super::functional_group::FunctionalGroupType::CarboxylicAcid
        }));

        let ethanol = registry.substance(&"destroy:ethanol".into()).unwrap();
        assert!(ethanol.functional_groups.iter().any(|group| {
            group.group_type == super::super::functional_group::FunctionalGroupType::Alcohol
        }));

        let chloroethane = registry.substance(&"destroy:chloroethane".into()).unwrap();
        assert!(chloroethane.functional_groups.iter().any(|group| {
            group.group_type == super::super::functional_group::FunctionalGroupType::Halide
        }));
    }

    #[test]
    fn destroy_catalog_contains_inorganic_equilibrium_data() {
        let registry = destroy_substances_registry_builder()
            .unwrap()
            .build()
            .unwrap();

        let carbon_dioxide = registry
            .substance_index(&"destroy:carbon_dioxide".into())
            .unwrap();
        assert!(registry.gas_solubility(carbon_dioxide).is_some());

        assert!(registry
            .indexed_equilibria()
            .iter()
            .any(|equilibrium| equilibrium.spec.id == "destroy:boric_acid.hydrolysis"));

        assert!(registry
            .complex_specs()
            .any(|complex| complex.id == SubstanceId::from("destroy:zinc_tetraammine")));

        let dichromate = registry
            .redox_half_reaction("destroy:dichromate_to_chromium_iii")
            .unwrap();
        assert_eq!(dichromate.electron_count, 6);
        assert_eq!(dichromate.standard_potential_volts, Some(1.33));
    }
}
