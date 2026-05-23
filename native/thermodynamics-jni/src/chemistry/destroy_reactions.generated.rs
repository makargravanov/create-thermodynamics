use super::error::ChemistryResult;
use super::reaction::Reaction;
use super::registry::ChemistryRegistryBuilder;

pub const DESTROY_EXPLICIT_REACTION_COUNT: usize = 119;
pub const DESTROY_REGISTERED_REACTION_COUNT: usize = 149;

pub fn destroy_reactions_registry_builder(mut builder: ChemistryRegistryBuilder) -> ChemistryResult<ChemistryRegistryBuilder> {
    builder = builder.reaction(
        Reaction::builder("destroy:abs_copolymerization")
            .reactant("destroy:acrylonitrile", 1, 1)
            .reactant("destroy:butadiene", 1, 1)
            .reactant("destroy:styrene", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(1f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.ABS::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 1.0)
            .pre_exponential_factor(15.0)
            .activation_energy_kj_per_mol(20.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:acetylene_trimerization")
            .reactant("destroy:acetylene", 3, 3)
            .product("destroy:benzene", 1)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/nickel\"), 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:acrylonitrile_polymerization")
            .reactant("destroy:acrylonitrile", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(3f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.POLYACRYLONITRILE::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 3.0)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:aibn_synthesis")
            .reactant("destroy:acetone_cyanohydrin", 2, 2)
            .reactant("destroy:hydrazine", 1, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:aibn", 1)
            .product("destroy:water", 2)
            .product("destroy:hydrochloric_acid", 2)
            .reaction_result("withResult(0f, DestroyAdvancementTrigger.AIBN::asReactionResult)", 0.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:andrussow_process")
            .reactant("destroy:methane", 2, 1)
            .reactant("destroy:ammonia", 2, 1)
            .reactant("destroy:oxygen", 3, 1)
            .product("destroy:hydrogen_cyanide", 2)
            .product("destroy:water", 6)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/platinum\"), 1f)", 1.0)
            .reaction_result("withResult(0f, DestroyAdvancementTrigger.ANDRUSSOW_PROCESS::asReactionResult)", 0.0)
            .pre_exponential_factor(10000000000.0)
            .activation_energy_kj_per_mol(50.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:anthraquinone_process")
            .reactant("destroy:ethylanthrahydroquinone", 1, 1)
            .reactant("destroy:oxygen", 1, 1)
            .product("destroy:ethylanthraquinone", 1)
            .product("destroy:hydrogen_peroxide", 1)
            .reaction_result("withResult(0f, DestroyAdvancementTrigger.HYDROGEN_PEROXIDE::asReactionResult)", 0.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:anthraquinone_reduction")
            .reactant("destroy:ethylanthraquinone", 1, 1)
            .reactant("destroy:hydrogen", 1, 1)
            .product("destroy:ethylanthrahydroquinone", 1)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/palladium\"), 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:baby_blue_precipitation")
            .reactant("destroy:methyl_salicylate", 1, 1)
            .catalyst_order("destroy:sodium_ion", 0)
            .reaction_result("withResult(0.9f, PrecipitateReactionResult.of(DestroyItems.BABY_BLUE_CRYSTAL::asStack))", 0.9)
            .show_in_jei_condition("includeInJeiIf(DestroySubstancesConfigs::babyBlueEnabled)")
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:basic_diels_alder_reaction")
            .reactant("destroy:butadiene", 1, 1)
            .reactant("destroy:ethene", 1, 1)
            .product("destroy:cyclohexene", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:benzene_ethylation")
            .reactant("destroy:benzene", 1, 1)
            .reactant("destroy:ethene", 1, 1)
            .catalyst_order("destroy:proton", 1)
            .product("destroy:ethylbenzene", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:benzene_hydrogenation")
            .reactant("destroy:benzene", 1, 1)
            .reactant("destroy:hydrogen", 2, 1)
            .product("destroy:cyclohexene", 1)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/nickel\"), 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:borax_dissolution")
            .reactant("destroy:proton", 2, 1)
            .catalyst_order("destroy:chloride", 1)
            .product("destroy:sodium_ion", 2)
            .product("destroy:water", 5)
            .product("destroy:boric_acid", 4)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"raw_materials/borax\"), 15f)", 15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:borax_precipitation")
            .reactant("destroy:sodium_ion", 2, 1)
            .reactant("destroy:water", 4, 1)
            .reactant("destroy:tetrahydroxy_tetraborate", 1, 1)
            .reaction_result("withResult(15f, PrecipitateReactionResult.of(DestroyItems.BORAX::asStack))", 15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:boric_acid_neutralization")
            .reactant("destroy:boric_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:tetrahydroxyborate", 1)
            .display_as_reversible()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:borohydride_iodine_oxidation")
            .reactant("destroy:borohydride", 2, 2)
            .reactant("destroy:iodine", 1, 1)
            .product("destroy:diborane", 1)
            .product("destroy:iodide", 2)
            .product("destroy:hydrogen", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:brown_schlesinger_process")
            .reactant("destroy:trimethyl_borate", 1, 1)
            .product("destroy:borohydride", 1)
            .product("destroy:ethoxide", 3)
            .product("destroy:sodium_ion", 4)
            .external_reactant("addSimpleItemReactant(DestroyItems.SODIUM_HYDRIDE::get, 2.4f)", 2.4)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:butadiene_carbonylation")
            .reactant("destroy:butadiene", 1, 1)
            .reactant("destroy:carbon_monoxide", 2, 2)
            .reactant("destroy:water", 2, 1)
            .product("destroy:adipic_acid", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:carbide_hydrolysis")
            .reactant("destroy:water", 1, 1)
            .product("destroy:acetylene", 1)
            .product("destroy:calcium_ion", 1)
            .product("destroy:hydroxide", 2)
            .external_reactant("addSimpleItemReactant(DestroyItems.CALCIUM_CARBIDE, 2f)", 2.0)
            .activation_energy_kj_per_mol(1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:carbon_capture")
            .reactant("destroy:calcium_ion", 1, 1)
            .reactant("destroy:carbon_dioxide", 1, 1)
            .reactant("destroy:water", 1, 1)
            .product("destroy:proton", 2)
            .reaction_result("withResult(2f, PrecipitateReactionResult.of(DestroyItems.CHALK_DUST::asStack))", 2.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:carbon_monoxide_oxidation")
            .reactant("destroy:carbon_monoxide", 2, 1)
            .reactant("destroy:oxygen", 1, 1)
            .product("destroy:carbon_dioxide", 2)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:carbon_tetrachloride_fluorination")
            .reactant("destroy:carbon_tetrachloride", 2, 1)
            .reactant("destroy:hydrofluoric_acid", 3, 1)
            .product("destroy:dichlorodifluoromethane", 1)
            .product("destroy:trichlorofluoromethane", 1)
            .product("destroy:hydrochloric_acid", 3)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:cellulose_nitration")
            .reactant("destroy:nitronium", 1, 1)
            .product("destroy:proton", 1)
            .product("destroy:water", 1)
            .external_reactant("addSimpleItemReactant(AllItems.PULP, 2f)", 2.0)
            .reaction_result("withResult(2f, PrecipitateReactionResult.of(DestroyItems.NITROCELLULOSE::asStack))", 2.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chlorine_haloform_reaction")
            .reactant("destroy:hypochlorite", 3, 0)
            .reactant("destroy:acetone", 1, 1)
            .catalyst_order("destroy:hydroxide", 0)
            .product("destroy:acetate", 1)
            .product("destroy:chloroform", 1)
            .product("destroy:hydroxide", 2)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chlorine_solvation")
            .reactant("destroy:chlorine", 1, 1)
            .reactant("destroy:water", 1, 1)
            .product("destroy:hydrochloric_acid", 1)
            .product("destroy:hypochlorous_acid", 1)
            .requires_uv()
            .display_as_reversible()
            .pre_exponential_factor(2.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chlorodifluoromethane_pyrolysis")
            .reactant("destroy:chlorodifluoromethane", 2, 2)
            .product("destroy:hydrochloric_acid", 2)
            .product("destroy:tetrafluoroethene", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chloroethene_polymerization")
            .reactant("destroy:chloroethene", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(3f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.POLYVINYL_CHLORIDE::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 3.0)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chloroform_fluorination")
            .reactant("destroy:chloroform", 1, 1)
            .reactant("destroy:hydrofluoric_acid", 2, 2)
            .product("destroy:chlorodifluoromethane", 1)
            .product("destroy:hydrochloric_acid", 2)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chromium_dissolution")
            .reactant("destroy:proton", 6, 1)
            .product("destroy:hydrogen", 3)
            .product("destroy:chromium_iii", 2)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/chromium\"), 4.5f)", 4.5)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chromium_ore_dissolution")
            .reactant("destroy:proton", 6, 1)
            .product("destroy:hydrogen", 3)
            .product("destroy:chromium_iii", 2)
            .external_reactant("addSimpleItemReactant(DestroyItems.CRUSHED_RAW_CHROMIUM::get, 7.5f)", 7.5)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chromium_iii_oxidation")
            .reactant("destroy:chromium_iii", 2, 1)
            .reactant("destroy:hydrogen_peroxide", 3, 1)
            .reactant("destroy:hydroxide", 10, 1)
            .product("destroy:chromate", 2)
            .product("destroy:water", 8)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:cisplatin_synthesis")
            .reactant("destroy:chloride", 1, 1)
            .reactant("destroy:ammonia", 1, 1)
            .product("destroy:cisplatin", 1)
            .product("destroy:hydroxide", 1)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/platinum\"), 2f)", 2.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:contact_process")
            .reactant("destroy:sulfur_dioxide", 2, 1)
            .reactant("destroy:oxygen", 1, 1)
            .product("destroy:sulfur_trioxide", 2)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/platinum\"), 3f)", 3.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:cordite_precipitation")
            .reactant("destroy:acetone", 1, 1)
            .reactant("destroy:nitroglycerine", 1, 1)
            .external_reactant("addSimpleItemReactant(DestroyItems.NITROCELLULOSE::get, 1f)", 1.0)
            .reaction_result("withResult(2.99f, PrecipitateReactionResult.of(DestroyBlocks.CORDITE_BLOCK::asStack))", 2.99)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:crocoite_dissolution")
            .catalyst_order("destroy:nitrate", 1)
            .product("destroy:lead_ii", 1)
            .product("destroy:chromate", 1)
            .external_reactant("addSimpleItemReactant(DestroyItems.NETHER_CROCOITE, 15f)", 15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:cumene_process")
            .reactant("destroy:benzene", 1, 1)
            .reactant("destroy:propene", 1, 1)
            .reactant("destroy:oxygen", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .catalyst_order("destroy:proton", 1)
            .product("destroy:phenol", 1)
            .product("destroy:acetone", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:cyanamide_ion_hydrolysis")
            .reactant("destroy:cyanamide_ion", 1, 1)
            .reactant("destroy:water", 2, 1)
            .product("destroy:cyanamide", 1)
            .product("destroy:hydroxide", 2)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:cyclohexene_oxidative_cleavage")
            .reactant("destroy:cyclohexene", 1, 1)
            .reactant("destroy:hydrogen_peroxide", 3, 1)
            .product("destroy:adipic_acid", 1)
            .product("destroy:water", 2)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chromate_conversion")
            .reactant("destroy:chromate", 2, 1)
            .reactant("destroy:proton", 2, 1)
            .product("destroy:dichromate", 1)
            .product("destroy:water", 1)
            .display_as_reversible()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:copper_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:copper_ii", 1)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/copper\"), 9f)", 9.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:copper_ore_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:copper_ii", 1)
            .external_reactant("addSimpleItemReactant(AllItems.CRUSHED_COPPER::get, 15f)", 15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:creatine_precipitation")
            .reactant("destroy:creatine", 1, 1)
            .catalyst_order("destroy:water", 1)
            .reaction_result("withResult(10f, PrecipitateReactionResult.of(DestroyItems.CREATINE::asStack))", 10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:diborane_hydrolysis")
            .reactant("destroy:diborane", 1, 1)
            .reactant("destroy:water", 6, 1)
            .product("destroy:boric_acid", 2)
            .product("destroy:hydrogen", 6)
            .activation_energy_kj_per_mol(1.0)
            .enthalpy_change_kj_per_mol(-466.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:dinitrotoluene_nitration")
            .reactant("destroy:dinitrotoluene", 1, 1)
            .reactant("destroy:nitronium", 1, 1)
            .product("destroy:tnt", 1)
            .product("destroy:proton", 1)
            .pre_exponential_factor(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ethene_polymerization")
            .reactant("destroy:ethene", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(3f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.POLYETHENE::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 3.0)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ethylanthraquinone_synthesis")
            .reactant("destroy:phthalic_anhydride", 1, 1)
            .reactant("destroy:ethylbenzene", 1, 1)
            .product("destroy:water", 1)
            .product("destroy:ethylanthraquinone", 1)
            .external_catalyst("addSimpleItemCatalyst(DestroyItems.ZEOLITE::get, 1f)", 1.0)
            .reaction_result("withResult(0f, DestroyAdvancementTrigger.ETHYLANTHRAQUINONE::asReactionResult)", 0.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ethylbenzene_dehydrogenation")
            .reactant("destroy:ethylbenzene", 1, 1)
            .catalyst_order("destroy:water", 2)
            .catalyst_order("destroy:iron_iii", 1)
            .product("destroy:styrene", 1)
            .product("destroy:hydrogen", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ethylbenzene_transalkylation")
            .reactant("destroy:ethylbenzene", 3, 3)
            .product("destroy:metaxylene", 1)
            .product("destroy:orthoxylene", 1)
            .product("destroy:paraxylene", 1)
            .external_catalyst("addSimpleItemCatalyst(DestroyItems.ZEOLITE::get, 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:fluorite_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:calcium_ion", 1)
            .product("destroy:hydrofluoric_acid", 2)
            .external_reactant("addSimpleItemReactant(DestroyItems.FLUORITE::get, 5f)", 5.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:frank_caro_process")
            .reactant("destroy:nitrogen", 1, 1)
            .reactant("destroy:oxygen", 1, 1)
            .product("destroy:calcium_ion", 1)
            .product("destroy:cyanamide_ion", 1)
            .product("destroy:carbon_dioxide", 1)
            .external_reactant("addSimpleItemReactant(DestroyItems.CALCIUM_CARBIDE, 2f)", 2.0)
            .activation_energy_kj_per_mol(50.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:glycerol_nitration")
            .reactant("destroy:glycerol", 1, 1)
            .reactant("destroy:nitronium", 3, 3)
            .product("destroy:proton", 3)
            .product("destroy:nitroglycerine", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:gold_dissolution")
            .reactant("destroy:nitrate", 3, 1)
            .reactant("destroy:chloride", 4, 1)
            .reactant("destroy:proton", 6, 2)
            .product("destroy:chloroaurate", 1)
            .product("destroy:water", 3)
            .product("destroy:nitrogen_dioxide", 3)
            .external_reactant("addSimpleItemReactant(() -> Items.GOLDEN_CARROT, 10f)", 10.0)
            .pre_exponential_factor(1000.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:haber_process")
            .reactant("destroy:nitrogen", 1, 1)
            .reactant("destroy:hydrogen", 3, 0)
            .product("destroy:ammonia", 2)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/iron\"), 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_chloride_synthesis")
            .reactant("destroy:hydrogen", 1, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:hydrochloric_acid", 2)
            .requires_uv()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_combustion")
            .reactant("destroy:hydrogen", 2, 1)
            .reactant("destroy:oxygen", 1, 1)
            .product("destroy:water", 2)
            .pre_exponential_factor(10000000000.0)
            .activation_energy_kj_per_mol(100.0)
            .enthalpy_change_kj_per_mol(-500.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_cyanide_dissociation")
            .reactant("destroy:hydrogen_cyanide", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:water", 1)
            .product("destroy:cyanide", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_iodide_synthesis")
            .reactant("destroy:hydrazine", 1, 1)
            .reactant("destroy:iodine", 2, 2)
            .product("destroy:hydrogen_iodide", 4)
            .product("destroy:nitrogen", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydroxide_neutralization")
            .reactant("destroy:hydroxide", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:water", 1)
            .display_as_reversible()
            .pre_exponential_factor(130000000.0)
            .activation_energy_kj_per_mol(0.0)
            .enthalpy_change_kj_per_mol(-55.3745)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hypochlorite_formation")
            .reactant("destroy:chlorine", 1, 1)
            .reactant("destroy:hydroxide", 2, 1)
            .catalyst_order("destroy:sodium_ion", 1)
            .product("destroy:chloride", 1)
            .product("destroy:hypochlorite", 1)
            .product("destroy:water", 1)
            .requires_uv()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hypochlorous_acid_dissociation")
            .reactant("destroy:hypochlorous_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:water", 1)
            .product("destroy:hypochlorite", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:iodide_displacement")
            .reactant("destroy:iodide", 2, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:iodine", 1)
            .product("destroy:chloride", 2)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:iodine_dissolution")
            .catalyst_order("destroy:water", 0)
            .product("destroy:iodine", 1)
            .external_reactant("addSimpleItemReactant(DestroyItems.IODINE::get, 2f)", 2.0)
            .display_as_reversible()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:iron_dissolution")
            .reactant("destroy:proton", 6, 1)
            .product("destroy:hydrogen", 3)
            .product("destroy:iron_iii", 2)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/iron\"), 4.5f)", 4.5)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:iron_ore_dissolution")
            .reactant("destroy:proton", 6, 1)
            .product("destroy:hydrogen", 3)
            .product("destroy:iron_iii", 2)
            .external_reactant("addSimpleItemReactant(AllItems.CRUSHED_IRON::get, 7.5f)", 7.5)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:isoprene_polymerization")
            .reactant("destroy:isoprene", 1, 1)
            .reaction_result("withResult(3f, PrecipitateReactionResult.of(DestroyItems.POLYISOPRENE::asStack))", 3.0)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:iron_iii_reduction")
            .reactant("destroy:iron_iii", 1, 1)
            .product("destroy:iron_ii", 1)
            .display_as_reversible()
            .allow_mass_imbalance()
            .allow_charge_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:kelp_dissolution")
            .catalyst_order("destroy:ethanol", 0)
            .product("destroy:potassium_ion", 1)
            .product("destroy:iodide", 1)
            .external_reactant("addSimpleItemReactant(() -> Items.DRIED_KELP, 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:kolbe_schmitt_reaction")
            .reactant("destroy:carbon_dioxide", 1, 1)
            .reactant("destroy:phenol", 1, 1)
            .catalyst_order("destroy:sodium_ion", 1)
            .catalyst_order("destroy:proton", 1)
            .product("destroy:salicylic_acid", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:lead_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:lead_ii", 1)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/lead\"), 9f)", 9.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:lead_ore_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:lead_ii", 1)
            .external_reactant("addSimpleItemReactant(AllItems.CRUSHED_LEAD::get, 15f)", 15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:lime_slaking")
            .reactant("destroy:water", 1, 1)
            .product("destroy:calcium_ion", 1)
            .product("destroy:hydroxide", 2)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/lime\"), 2f)", 2.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:mercury_fulmination")
            .reactant("destroy:mercury", 3, 1)
            .reactant("destroy:nitrate", 12, 2)
            .reactant("destroy:ethanol", 4, 1)
            .product("destroy:carbon_dioxide", 2)
            .product("destroy:hydroxide", 12)
            .product("destroy:water", 6)
            .product("destroy:nitrogen_dioxide", 6)
            .reaction_result("withResult(5f, PrecipitateReactionResult.of(DestroyItems.FULMINATED_MERCURY::asStack))", 5.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:metaxylene_transalkylation")
            .reactant("destroy:metaxylene", 3, 3)
            .product("destroy:orthoxylene", 1)
            .product("destroy:paraxylene", 1)
            .product("destroy:ethylbenzene", 1)
            .external_catalyst("addSimpleItemCatalyst(DestroyItems.ZEOLITE::get, 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:methanol_synthesis")
            .reactant("destroy:carbon_monoxide", 1, 1)
            .reactant("destroy:hydrogen", 2, 1)
            .product("destroy:methanol", 1)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/copper\"), 1f)", 1.0)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/zinc\"), 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:methyl_acetate_carbonylation")
            .reactant("destroy:methanol", 1, 1)
            .reactant("destroy:carbon_monoxide", 1, 1)
            .product("destroy:acetic_acid", 1)
            .external_catalyst("addSimpleItemCatalyst(DestroyItems.SILICA::get, 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:methyl_methacrylate_polymerization")
            .reactant("destroy:methyl_methacrylate", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(3f, PrecipitateReactionResult.of(DestroyItems.POLYMETHYL_METHACRYLATE::asStack))", 3.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:naughty_reaction")
            .reactant("destroy:phenylacetone", 1, 1)
            .reactant("destroy:methylamine", 1, 1)
            .reaction_result("withResult(0f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(ExplosionReactionResult::small)\n            .with(DestroyAdvancementTrigger.TRY_TO_MAKE_METH::asReactionResult))", 0.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nhn_synthesis")
            .reactant("destroy:nickel_ion", 1, 1)
            .reactant("destroy:nitrate", 2, 0)
            .reactant("destroy:hydrazine", 3, 3)
            .reaction_result("withResult(3f, PrecipitateReactionResult.of(DestroyItems.NICKEL_HYDRAZINE_NITRATE::asStack))", 3.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nickel_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:nickel_ion", 1)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/nickel\"), 9f)", 9.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nickel_ore_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:nickel_ion", 1)
            .external_reactant("addSimpleItemReactant(AllItems.CRUSHED_NICKEL::get, 15f)", 15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nitronium_formation")
            .reactant("destroy:nitric_acid", 1, 1)
            .reactant("destroy:sulfuric_acid", 1, 1)
            .product("destroy:nitronium", 1)
            .product("destroy:water", 1)
            .product("destroy:hydrogensulfate", 1)
            .display_as_reversible()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nylon_polymerisation")
            .reactant("destroy:adipic_acid", 1, 1)
            .reactant("destroy:hexanediamine", 1, 1)
            .reaction_result("withResult(3f, PrecipitateReactionResult.of(DestroyItems.NYLON::asStack))", 3.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:oleum_formation")
            .reactant("destroy:sulfuric_acid", 1, 1)
            .reactant("destroy:sulfur_trioxide", 1, 1)
            .product("destroy:oleum", 1)
            .display_as_reversible()
            .pre_exponential_factor(20000.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:oleum_hydration")
            .reactant("destroy:oleum", 1, 1)
            .reactant("destroy:water", 1, 1)
            .product("destroy:sulfuric_acid", 2)
            .pre_exponential_factor(10000000000.0)
            .activation_energy_kj_per_mol(2.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:orthoxylene_oxidation")
            .reactant("destroy:orthoxylene", 1, 1)
            .reactant("destroy:oxygen", 3, 1)
            .catalyst_order("destroy:mercury", 1)
            .product("destroy:phthalic_anhydride", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:orthoxylene_transalkylation")
            .reactant("destroy:orthoxylene", 3, 3)
            .product("destroy:metaxylene", 1)
            .product("destroy:paraxylene", 1)
            .product("destroy:ethylbenzene", 1)
            .external_catalyst("addSimpleItemCatalyst(DestroyItems.ZEOLITE::get, 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ostwald_process")
            .reactant("destroy:ammonia", 1, 1)
            .reactant("destroy:oxygen", 2, 2)
            .product("destroy:water", 1)
            .product("destroy:nitric_acid", 1)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/rhodium\"), 1f)", 1.0)
            .reaction_result("withResult(0f, DestroyAdvancementTrigger.OSTWALD_PROCESS::asReactionResult)", 0.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:paraxylene_transalkylation")
            .reactant("destroy:paraxylene", 3, 3)
            .product("destroy:metaxylene", 1)
            .product("destroy:orthoxylene", 1)
            .product("destroy:ethylbenzene", 1)
            .external_catalyst("addSimpleItemCatalyst(DestroyItems.ZEOLITE::get, 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:peroxide_process")
            .reactant("destroy:hydrogen_peroxide", 1, 1)
            .reactant("destroy:ammonia", 2, 1)
            .catalyst_order("destroy:acetone", 1)
            .catalyst_order("destroy:proton", 0)
            .product("destroy:hydrazine", 1)
            .product("destroy:water", 2)
            .reaction_result("withResult(0f, DestroyAdvancementTrigger.HYDRAZINE::asReactionResult)", 0.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:phenol_nitration")
            .reactant("destroy:phenol", 1, 1)
            .reactant("destroy:nitronium", 3, 1)
            .product("destroy:picric_acid", 1)
            .product("destroy:proton", 3)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:phosgene_formation")
            .reactant("destroy:carbon_monoxide", 1, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:phosgene", 1)
            .enthalpy_change_kj_per_mol(-107.6)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:propene_polymerization")
            .reactant("destroy:propene", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(3f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.POLYPROPENE::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 3.0)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sodium_amalgamization")
            .catalyst_order("destroy:mercury", 1)
            .product("destroy:sodium_metal", 1)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"ingots/sodium\"), 9.9f)", 9.9)
            .display_as_reversible()
            .activation_energy_kj_per_mol(1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sodium_dissolution")
            .reactant("destroy:sodium_metal", 2, 1)
            .reactant("destroy:water", 2, 1)
            .product("destroy:sodium_ion", 2)
            .product("destroy:hydroxide", 2)
            .product("destroy:hydrogen", 1)
            .activation_energy_kj_per_mol(1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sodium_hydride_formation")
            .reactant("destroy:sodium_metal", 2, 2)
            .reactant("destroy:hydrogen", 1, 1)
            .reaction_result("withResult(10f, PrecipitateReactionResult.of(DestroyItems.SODIUM_HYDRIDE::asStack))", 10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sodium_hydride_hydrolysis")
            .reactant("destroy:water", 1, 1)
            .product("destroy:sodium_ion", 1)
            .product("destroy:hydroxide", 1)
            .product("destroy:hydrogen", 1)
            .external_reactant("addSimpleItemReactant(DestroyItems.SODIUM_HYDRIDE::get, 10.1f)", 10.1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sodium_ingot_dissolution")
            .reactant("destroy:water", 2, 1)
            .product("destroy:sodium_ion", 2)
            .product("destroy:hydroxide", 2)
            .product("destroy:hydrogen", 1)
            .external_reactant("addSimpleItemReactant(DestroyItems.SODIUM_INGOT::get, 4.9f)", 4.9)
            .activation_energy_kj_per_mol(1.0)
            .enthalpy_change_kj_per_mol(-370.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sodium_oxide_dissolution")
            .reactant("destroy:water", 1, 1)
            .product("destroy:sodium_ion", 2)
            .product("destroy:hydroxide", 2)
            .external_reactant("addSimpleItemReactant(DestroyItems.OXIDIZED_SODIUM_INGOT, 4.9f)", 4.9)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:steam_reformation")
            .reactant("destroy:water", 1, 1)
            .reactant("destroy:methane", 1, 1)
            .product("destroy:carbon_monoxide", 1)
            .product("destroy:hydrogen", 3)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/nickel\"), 1f)", 1.0)
            .display_as_reversible()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:styrene_butadiene_copolymerization")
            .reactant("destroy:styrene", 1, 1)
            .reactant("destroy:butadiene", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(1.5f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.POLYSTYRENE_BUTADIENE::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 1.5)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:styrene_polymerization")
            .reactant("destroy:styrene", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(3f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.POLYSTYRENE::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 3.0)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sulfur_oxidation")
            .reactant("destroy:octasulfur", 1, 1)
            .reactant("destroy:oxygen", 8, 1)
            .product("destroy:sulfur_dioxide", 8)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sulfur_trioxide_hydration")
            .reactant("destroy:sulfur_trioxide", 1, 1)
            .reactant("destroy:water", 1, 1)
            .product("destroy:sulfuric_acid", 1)
            .display_as_reversible()
            .activation_energy_kj_per_mol(10.0)
            .enthalpy_change_kj_per_mol(-200.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:tatp")
            .reactant("destroy:acetone", 1, 1)
            .reactant("destroy:hydrogen_peroxide", 1, 1)
            .catalyst_order("destroy:proton", 1)
            .reaction_result("withResult(3f, PrecipitateReactionResult.of(DestroyItems.ACETONE_PEROXIDE::asStack))", 3.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:tetraborate_equilibrium")
            .reactant("destroy:boric_acid", 2, 1)
            .reactant("destroy:tetrahydroxyborate", 2, 1)
            .product("destroy:tetrahydroxy_tetraborate", 1)
            .product("destroy:water", 5)
            .display_as_reversible()
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:tetraethyllead_synthesis")
            .reactant("destroy:sodium_metal", 4, 4)
            .reactant("destroy:chloroethane", 4, 4)
            .product("destroy:tetraethyllead", 1)
            .product("destroy:sodium_ion", 4)
            .product("destroy:chloride", 4)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/lead\"), 2.5f)", 2.5)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:tetrafluoroethene_polymerization")
            .reactant("destroy:tetrafluoroethene", 1, 1)
            .catalyst_order("destroy:aibn", 0)
            .reaction_result("withResult(3f, (m, r) -> new CombinedReactionResult(m, r)\n            .with(PrecipitateReactionResult.of(DestroyItems.POLYTETRAFLUOROETHENE::asStack))\n            .with(DestroyAdvancementTrigger.ADDITION_POLYMER::asReactionResult))", 3.0)
            .pre_exponential_factor(10.0)
            .activation_energy_kj_per_mol(10.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:thionyl_chloride_synthesis")
            .reactant("destroy:sulfur_dioxide", 1, 1)
            .reactant("destroy:phosgene", 1, 1)
            .product("destroy:thionyl_chloride", 1)
            .product("destroy:carbon_dioxide", 1)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:toluene_nitration")
            .reactant("destroy:toluene", 1, 1)
            .reactant("destroy:nitronium", 2, 1)
            .product("destroy:dinitrotoluene", 1)
            .product("destroy:proton", 2)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:toluene_transalkylation")
            .reactant("destroy:toluene", 8, 8)
            .product("destroy:benzene", 4)
            .product("destroy:metaxylene", 1)
            .product("destroy:orthoxylene", 1)
            .product("destroy:paraxylene", 1)
            .product("destroy:ethylbenzene", 1)
            .external_catalyst("addSimpleItemCatalyst(DestroyItems.ZEOLITE::get, 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:touch_powder_synthesis")
            .reactant("destroy:ammonia", 1, 1)
            .external_reactant("addSimpleItemReactant(DestroyItems.IODINE::get, 3f)", 3.0)
            .reaction_result("withResult(3f, PrecipitateReactionResult.of(DestroyItems.TOUCH_POWDER::asStack))", 3.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:urethane_hdi_polymerization")
            .reactant("destroy:glycerol", 1, 1)
            .reactant("destroy:hexane_diisocyanate", 1, 1)
            .reaction_result("withResult(1f, PrecipitateReactionResult.of(DestroyItems.POLYURETHANE::asStack))", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:urethane_tdi_polymerization")
            .reactant("destroy:glycerol", 1, 1)
            .reactant("destroy:toluene_diisocyanate", 1, 1)
            .reaction_result("withResult(1f, PrecipitateReactionResult.of(DestroyItems.POLYURETHANE::asStack))", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:vinyl_acetate_synthesis")
            .reactant("destroy:ethene", 2, 1)
            .reactant("destroy:acetic_acid", 2, 1)
            .reactant("destroy:oxygen", 1, 1)
            .product("destroy:vinyl_acetate", 2)
            .product("destroy:water", 2)
            .external_catalyst("addSimpleItemTagCatalyst(AllTags.forgeItemTag(\"dusts/palladium\"), 1f)", 1.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:zinc_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:zinc_ion", 1)
            .external_reactant("addSimpleItemTagReactant(AllTags.forgeItemTag(\"dusts/zinc\"), 9f)", 9.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:zinc_ore_dissolution")
            .reactant("destroy:proton", 2, 1)
            .product("destroy:hydrogen", 1)
            .product("destroy:zinc_ion", 1)
            .external_reactant("addSimpleItemReactant(AllItems.CRUSHED_ZINC::get, 15f)", 15.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:methane_uv_chlorination")
            .reactant("destroy:methane", 1, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:chloromethane", 1)
            .product("destroy:hydrochloric_acid", 1)
            .requires_uv()
            .activation_energy_kj_per_mol(22.5)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chloromethane_uv_chlorination")
            .reactant("destroy:chloromethane", 1, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:dichloromethane", 1)
            .product("destroy:hydrochloric_acid", 1)
            .requires_uv()
            .activation_energy_kj_per_mol(25.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:dichloromethane_uv_chlorination")
            .reactant("destroy:dichloromethane", 1, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:chloroform", 1)
            .product("destroy:hydrochloric_acid", 1)
            .requires_uv()
            .activation_energy_kj_per_mol(27.5)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:chloroform_uv_chlorination")
            .reactant("destroy:chloroform", 1, 1)
            .reactant("destroy:chlorine", 1, 1)
            .product("destroy:carbon_tetrachloride", 1)
            .product("destroy:hydrochloric_acid", 1)
            .requires_uv()
            .activation_energy_kj_per_mol(30.0)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:acetic_acid.dissociation")
            .reactant("destroy:acetic_acid", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:acetate", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(8.689004143746882e-6)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:acetic_acid.neutralization")
            .reactant("destroy:acetic_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:acetate", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(8.689004143746882e-6)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:acetic_acid.association")
            .reactant("destroy:acetate", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:acetic_acid", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ammonium.dissociation")
            .reactant("destroy:ammonium", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:ammonia", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(2.8117066259517455e-10)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ammonium.neutralization")
            .reactant("destroy:ammonium", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:ammonia", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(2.8117066259517455e-10)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:ammonium.association")
            .reactant("destroy:ammonia", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:ammonium", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrochloric_acid.dissociation")
            .reactant("destroy:hydrochloric_acid", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:chloride", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(997631.1574844394)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrochloric_acid.neutralization")
            .reactant("destroy:hydrochloric_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:chloride", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(997631.1574844394)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrochloric_acid.association")
            .reactant("destroy:chloride", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:hydrochloric_acid", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrofluoric_acid.dissociation")
            .reactant("destroy:hydrofluoric_acid", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:fluoride", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(0.00033804148769599093)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrofluoric_acid.neutralization")
            .reactant("destroy:hydrofluoric_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:fluoride", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(0.00033804148769599093)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrofluoric_acid.association")
            .reactant("destroy:fluoride", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:hydrofluoric_acid", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_cyanide.dissociation")
            .reactant("destroy:hydrogen_cyanide", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:cyanide", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(3.154786722400971e-10)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_cyanide.neutralization")
            .reactant("destroy:hydrogen_cyanide", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:cyanide", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(3.154786722400971e-10)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_cyanide.association")
            .reactant("destroy:cyanide", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:hydrogen_cyanide", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_iodide.dissociation")
            .reactant("destroy:hydrogen_iodide", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:iodide", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(997631157.4844414)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_iodide.neutralization")
            .reactant("destroy:hydrogen_iodide", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:iodide", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(997631157.4844414)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogen_iodide.association")
            .reactant("destroy:iodide", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:hydrogen_iodide", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogensulfate.dissociation")
            .reactant("destroy:hydrogensulfate", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:sulfate", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(0.005116464961403771)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogensulfate.neutralization")
            .reactant("destroy:hydrogensulfate", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:sulfate", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(0.005116464961403771)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hydrogensulfate.association")
            .reactant("destroy:sulfate", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:hydrogensulfate", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hypochlorous_acid.dissociation")
            .reactant("destroy:hypochlorous_acid", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:hypochlorite", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.475604613333192e-8)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hypochlorous_acid.neutralization")
            .reactant("destroy:hypochlorous_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:hypochlorite", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.475604613333192e-8)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:hypochlorous_acid.association")
            .reactant("destroy:hypochlorite", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:hypochlorous_acid", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nitric_acid.dissociation")
            .reactant("destroy:nitric_acid", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:nitrate", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(9.976311574844399)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nitric_acid.neutralization")
            .reactant("destroy:nitric_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:nitrate", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(9.976311574844399)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:nitric_acid.association")
            .reactant("destroy:nitrate", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:nitric_acid", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sulfuric_acid.dissociation")
            .reactant("destroy:sulfuric_acid", 1, 1)
            .catalyst_order("destroy:water", 1)
            .product("destroy:proton", 1)
            .product("destroy:hydrogensulfate", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(75.67806242181044)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sulfuric_acid.neutralization")
            .reactant("destroy:sulfuric_acid", 1, 1)
            .reactant("destroy:hydroxide", 1, 1)
            .product("destroy:hydrogensulfate", 1)
            .product("destroy:water", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(75.67806242181044)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    builder = builder.reaction(
        Reaction::builder("destroy:sulfuric_acid.association")
            .reactant("destroy:hydrogensulfate", 1, 1)
            .reactant("destroy:proton", 1, 1)
            .product("destroy:sulfuric_acid", 1)
            .activation_energy_kj_per_mol(2.477709860209665)
            .pre_exponential_factor(1.0)
            .show_in_jei(false)
            .allow_mass_imbalance()
            .build(),
    );
    Ok(builder)
}
