//! Core types for kinetic selectivity system

use std::collections::HashMap;

use crate::chemistry::error::ChemistryResult;
use crate::chemistry::mixture::{
    LiquidPhaseSnapshot, Mixture, MixturePhase, STANDARD_PRESSURE_PASCAL,
    TRACE_CONCENTRATION_MOL_PER_BUCKET,
};
use crate::chemistry::redox::RedoxRole;
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::simulation::ReactionContext;
use crate::chemistry::substance::{LiquidPhasePreference, SolventRole, SubstanceId};

/// Substitution degree at a reactive center
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SubstitutionDegree {
    /// Primary: attached to ≤1 carbon (e.g., CH3OH, R-CH2-OH, R-CHO)
    Primary,
    /// Secondary: attached to 2 carbons (e.g., R2CH-OH, R2C=O ketone)
    Secondary,
    /// Tertiary: attached to 3 carbons (e.g., R3C-OH)
    Tertiary,
    /// Benzylic: attached to benzene ring (Ph-CH2-X)
    Benzylic,
    /// Allylic: attached to alkene carbon (C=C-CH2-X)
    Allylic,
}

impl SubstitutionDegree {
    /// Base steric score from degree alone (0.0 = unhindered, 1.0 = blocked)
    pub fn base_steric_score(&self) -> f64 {
        match self {
            SubstitutionDegree::Primary => 0.0,
            SubstitutionDegree::Secondary => 0.3,
            SubstitutionDegree::Tertiary => 0.6,
            SubstitutionDegree::Benzylic => 0.2, // phenyl is flat, less hindrance
            SubstitutionDegree::Allylic => 0.25,
        }
    }
}

/// Electronic environment around a reactive site
#[derive(Debug, Clone, Default)]
pub struct ElectronicEnvironment {
    /// Count of electron-donating groups (+I, +M effects) in conjugation/position
    pub electron_donating_groups: u32,
    /// Count of electron-withdrawing groups (-I, -M effects)
    pub electron_withdrawing_groups: u32,
    /// Whether site is resonance-stabilized (conjugated)
    pub resonance_stabilization: bool,
    /// Whether site is part of aromatic system
    pub aromatic: bool,
}

impl ElectronicEnvironment {
    /// Net electronic effect: positive = activating, negative = deactivating
    pub fn net_effect(&self) -> f64 {
        let edg_factor = self.electron_donating_groups as f64 * 0.1;
        let ewg_factor = self.electron_withdrawing_groups as f64 * 0.1;
        let resonance_bonus = if self.resonance_stabilization {
            0.15
        } else {
            0.0
        };
        edg_factor - ewg_factor + resonance_bonus
    }
}

/// Complete descriptor for selectivity evaluation
#[derive(Debug, Clone)]
pub struct SiteDescriptor {
    /// Type of reactive site
    pub site_kind: crate::chemistry::reactive_site::ReactiveSiteKind,
    /// Substitution degree at the reactive center
    pub degree: SubstitutionDegree,
    /// Electronic environment
    pub electronics: ElectronicEnvironment,
    /// Steric hindrance score (0.0 = open, 1.0 = fully hindered)
    pub steric_score: f64,
    /// Whether β-hydrogen is available (for elimination reactions)
    pub has_beta_hydrogen: bool,
}

impl SiteDescriptor {
    /// Create descriptor with steric score calculated from degree
    pub fn new(
        site_kind: crate::chemistry::reactive_site::ReactiveSiteKind,
        degree: SubstitutionDegree,
        electronics: ElectronicEnvironment,
        bulky_substituents: u32,
        has_beta_hydrogen: bool,
    ) -> Self {
        let base_steric = degree.base_steric_score();
        let additional_steric = (bulky_substituents as f64 * 0.1).min(0.4);
        let steric_score = (base_steric + additional_steric).min(1.0);

        Self {
            site_kind,
            degree,
            electronics,
            steric_score,
            has_beta_hydrogen,
        }
    }

    /// Steric accessibility factor (1.0 = fully accessible, 0.0 = blocked)
    pub fn steric_accessibility(&self) -> f64 {
        (1.0 - self.steric_score).max(0.1)
    }
}

/// Selectivity evaluation conditions
///
/// This is separate from simulation::ReactionContext which tracks
/// actual reaction state (external reactants, catalysts, etc.)
#[derive(Debug, Clone)]
pub struct SelectivityContext {
    /// Temperature in Kelvin
    pub temperature: f64,
    /// pH if applicable (None for non-aqueous/neutral)
    pub ph: Option<f64>,
    /// Physical reaction medium inferred from actual mixture phases.
    pub medium: ReactionMedium,
    /// Solvent type
    pub solvent_type: SolventType,
    pub water_activity: f64,
    pub fluoride_mol_per_bucket: f64,
    pub hydrogen_mol_per_bucket: f64,
    pub hydrogen_partial_pressure_pascal: f64,
    pub oxygen_activity: f64,
    pub oxygen_partial_pressure_pascal: f64,
    pub total_gas_pressure_pascal: f64,
    pub uv_power: f64,
    pub total_light_power: f64,
    pub oxidizing_strength: f64,
    pub reducing_strength: f64,
    pub free_complexable_metal_activity: f64,
    pub complexed_metal_mol_per_bucket: f64,
    pub available_surface_sites_mol_per_bucket: f64,
    pub palladium_available: bool,
}

impl Default for SelectivityContext {
    fn default() -> Self {
        Self {
            temperature: 298.15, // 25°C
            ph: None,
            medium: ReactionMedium::Vacuum,
            solvent_type: SolventType::Neutral,
            water_activity: 0.0,
            fluoride_mol_per_bucket: 0.0,
            hydrogen_mol_per_bucket: 0.0,
            hydrogen_partial_pressure_pascal: 0.0,
            oxygen_activity: 0.0,
            oxygen_partial_pressure_pascal: 0.0,
            total_gas_pressure_pascal: 0.0,
            uv_power: 0.0,
            total_light_power: 0.0,
            oxidizing_strength: 0.0,
            reducing_strength: 0.0,
            free_complexable_metal_activity: 0.0,
            complexed_metal_mol_per_bucket: 0.0,
            available_surface_sites_mol_per_bucket: 0.0,
            palladium_available: false,
        }
    }
}

impl SelectivityContext {
    /// Build selectivity conditions from the actual simulated mixture.
    ///
    /// This is intentionally conservative: if the phase composition cannot be
    /// classified confidently, neutral selectivity rules are used.
    pub fn from_mixture(
        registry: &ChemistryRegistry,
        mixture: &Mixture,
        reaction_context: &ReactionContext,
    ) -> ChemistryResult<Self> {
        let ph = mixture.ph(registry)?;
        let water = SubstanceId::from("destroy:water");
        let water_activity = max_known_activity(
            registry,
            mixture,
            &water,
            &[MixturePhase::Aqueous, MixturePhase::Organic],
        )?;
        let fluoride = SubstanceId::from("destroy:fluoride");
        let fluoride_mol_per_bucket = registry
            .substance(&fluoride)
            .ok()
            .map(|_| mixture.concentration_of(&fluoride))
            .unwrap_or(0.0);
        let hydrogen = SubstanceId::from("destroy:hydrogen");
        let hydrogen_mol_per_bucket = registry
            .substance(&hydrogen)
            .ok()
            .map(|_| mixture.concentration_of(&hydrogen))
            .unwrap_or(0.0);
        let hydrogen_partial_pressure_pascal = known_partial_pressure(registry, mixture, &hydrogen);
        let oxygen = SubstanceId::from("destroy:oxygen");
        let oxygen_activity = max_known_activity(
            registry,
            mixture,
            &oxygen,
            &[
                MixturePhase::Aqueous,
                MixturePhase::Organic,
                MixturePhase::Gas,
            ],
        )?;
        let oxygen_partial_pressure_pascal = known_partial_pressure(registry, mixture, &oxygen);
        let (oxidizing_strength, reducing_strength) = redox_role_strengths(registry, mixture)?;
        let (free_complexable_metal_activity, complexed_metal_mol_per_bucket) =
            metal_availability(registry, mixture)?;
        let total_light_power = reaction_context
            .light_power_by_band
            .values()
            .copied()
            .sum::<f64>()
            .max(reaction_context.uv_power);
        let available_surface_sites_mol_per_bucket = reaction_context
            .surfaces
            .values()
            .map(|surface| surface.free_sites())
            .sum::<f64>();
        let palladium_available =
            reaction_context
                .external_catalysts
                .iter()
                .any(|(description, amount)| {
                    *amount > TRACE_CONCENTRATION_MOL_PER_BUCKET
                        && description.to_ascii_lowercase().contains("palladium")
                })
                || reaction_context
                    .surfaces
                    .iter()
                    .any(|(surface_id, surface)| {
                        surface.free_sites() > TRACE_CONCENTRATION_MOL_PER_BUCKET
                            && surface_id
                                .as_str()
                                .to_ascii_lowercase()
                                .contains("palladium")
                    });
        let medium = ReactionMedium::from_mixture(registry, mixture, reaction_context)?;
        let mut context = Self {
            temperature: mixture.temperature_kelvin(),
            ph,
            medium,
            solvent_type: SolventType::Neutral,
            water_activity,
            fluoride_mol_per_bucket,
            hydrogen_mol_per_bucket,
            hydrogen_partial_pressure_pascal,
            oxygen_activity,
            oxygen_partial_pressure_pascal,
            total_gas_pressure_pascal: mixture.gas_pressure_pascal(),
            uv_power: reaction_context.uv_power,
            total_light_power,
            oxidizing_strength,
            reducing_strength,
            free_complexable_metal_activity,
            complexed_metal_mol_per_bucket,
            available_surface_sites_mol_per_bucket,
            palladium_available,
        };
        context.solvent_type = if context.is_basic() {
            SolventType::Basic
        } else if context.is_acidic() {
            SolventType::Acidic
        } else {
            context.medium.solvent_type()
        };
        Ok(context)
    }

    /// Create context for specific temperature
    pub fn at_temperature(kelvin: f64) -> Self {
        Self {
            temperature: kelvin,
            ..Default::default()
        }
    }

    /// Check if conditions are acidic (pH < 6)
    pub fn is_acidic(&self) -> bool {
        self.ph.map(|p| p < 6.0).unwrap_or(false)
    }

    /// Check if conditions are basic (pH > 8)
    pub fn is_basic(&self) -> bool {
        self.ph.map(|p| p > 8.0).unwrap_or(false)
    }

    /// Check if temperature is high (> 80°C)
    pub fn is_high_temperature(&self) -> bool {
        self.temperature > 353.15 // 80°C
    }

    /// Check if temperature is very high (> 150°C)
    pub fn is_very_high_temperature(&self) -> bool {
        self.temperature > 423.15 // 150°C
    }
    pub fn is_water_rich(&self) -> bool {
        self.water_activity >= 0.35
    }

    pub fn is_water_poor(&self) -> bool {
        self.water_activity <= 0.2
    }

    pub fn has_fluoride(&self) -> bool {
        self.fluoride_mol_per_bucket > TRACE_CONCENTRATION_MOL_PER_BUCKET
    }

    pub fn has_hydrogen(&self) -> bool {
        self.hydrogen_mol_per_bucket > TRACE_CONCENTRATION_MOL_PER_BUCKET
            || self.hydrogen_partial_pressure_pascal > STANDARD_PRESSURE_PASCAL * 1.0e-6
    }

    pub fn has_uv(&self) -> bool {
        self.uv_power > 0.0
    }

    pub fn is_oxygen_rich(&self) -> bool {
        self.oxygen_activity > TRACE_CONCENTRATION_MOL_PER_BUCKET
            || self.oxygen_partial_pressure_pascal > STANDARD_PRESSURE_PASCAL * 1.0e-4
    }

    pub fn is_oxidizing(&self) -> bool {
        self.oxidizing_strength > self.reducing_strength + TRACE_CONCENTRATION_MOL_PER_BUCKET
    }

    pub fn is_reducing(&self) -> bool {
        self.reducing_strength > self.oxidizing_strength + TRACE_CONCENTRATION_MOL_PER_BUCKET
    }

    pub fn has_available_surface(&self) -> bool {
        self.available_surface_sites_mol_per_bucket > TRACE_CONCENTRATION_MOL_PER_BUCKET
    }

    pub fn has_free_complexable_metal(&self) -> bool {
        self.free_complexable_metal_activity > TRACE_CONCENTRATION_MOL_PER_BUCKET
    }

    pub fn has_complexed_metal(&self) -> bool {
        self.complexed_metal_mol_per_bucket > TRACE_CONCENTRATION_MOL_PER_BUCKET
    }
}

/// Physical medium used by organic selectivity rules.
///
/// This is a derived view of `Mixture`; it must not duplicate or own phase state.
#[derive(Debug, Clone, PartialEq)]
pub enum ReactionMedium {
    Vacuum,
    Gas {
        total_mol_per_bucket: f64,
        pressure_pascal: f64,
    },
    Liquid {
        coarse_phase: MixturePhase,
        representative_solvent_id: SubstanceId,
        solvent_type: SolventType,
        total_solvent_mol_per_bucket: f64,
        water_mole_fraction: f64,
    },
    GasLiquidInterface {
        liquid: Box<ReactionMedium>,
        gas_mol_per_bucket: f64,
        pressure_pascal: f64,
    },
    SupercriticalFluid {
        carrier_id: Option<SubstanceId>,
        total_mol_per_bucket: f64,
        pressure_pascal: f64,
        solvent_type: SolventType,
    },
    MoltenMetal {
        total_mol_per_bucket: f64,
    },
    MoltenSlag {
        total_mol_per_bucket: f64,
    },
    SolidSurface {
        solid_mol_per_bucket: f64,
        surface_sites_mol_per_bucket: f64,
    },
    Heterogeneous {
        dominant_phase: MixturePhase,
        total_mol_per_bucket: f64,
    },
}

impl ReactionMedium {
    pub fn from_mixture(
        registry: &ChemistryRegistry,
        mixture: &Mixture,
        reaction_context: &ReactionContext,
    ) -> ChemistryResult<Self> {
        let gas_amount = mixture.total_concentration_in_phase(MixturePhase::Gas);
        let supercritical_amount =
            mixture.total_concentration_in_phase(MixturePhase::SupercriticalFluid);
        let molten_metal = mixture.total_concentration_in_phase(MixturePhase::MoltenMetal);
        let molten_slag = mixture.total_concentration_in_phase(MixturePhase::MoltenSlag);
        let solid = mixture.total_concentration_in_phase(MixturePhase::Solid);
        let surface_sites = reaction_context
            .surfaces
            .values()
            .map(|surface| surface.free_sites())
            .sum::<f64>();

        if supercritical_amount > TRACE_CONCENTRATION_MOL_PER_BUCKET
            && supercritical_amount >= gas_amount
            && supercritical_amount >= molten_metal
            && supercritical_amount >= molten_slag
        {
            return Ok(Self::SupercriticalFluid {
                carrier_id: dominant_substance_in_phase(
                    registry,
                    mixture,
                    MixturePhase::SupercriticalFluid,
                )?,
                total_mol_per_bucket: supercritical_amount,
                pressure_pascal: mixture.gas_pressure_pascal(),
                solvent_type: supercritical_solvent_type(registry, mixture)?,
            });
        }

        if molten_metal > TRACE_CONCENTRATION_MOL_PER_BUCKET && molten_metal >= molten_slag {
            return Ok(Self::MoltenMetal {
                total_mol_per_bucket: molten_metal,
            });
        }
        if molten_slag > TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Ok(Self::MoltenSlag {
                total_mol_per_bucket: molten_slag,
            });
        }

        if let Some(liquid) = dominant_liquid_medium(registry, mixture)? {
            if gas_amount > TRACE_CONCENTRATION_MOL_PER_BUCKET {
                return Ok(Self::GasLiquidInterface {
                    liquid: Box::new(liquid),
                    gas_mol_per_bucket: gas_amount,
                    pressure_pascal: mixture.gas_pressure_pascal(),
                });
            }
            return Ok(liquid);
        }

        if surface_sites > TRACE_CONCENTRATION_MOL_PER_BUCKET
            || solid > TRACE_CONCENTRATION_MOL_PER_BUCKET
        {
            return Ok(Self::SolidSurface {
                solid_mol_per_bucket: solid,
                surface_sites_mol_per_bucket: surface_sites,
            });
        }

        if gas_amount > TRACE_CONCENTRATION_MOL_PER_BUCKET {
            return Ok(Self::Gas {
                total_mol_per_bucket: gas_amount,
                pressure_pascal: mixture.gas_pressure_pascal(),
            });
        }

        let (dominant_phase, total_mol_per_bucket) = dominant_phase_by_amount(mixture);
        if total_mol_per_bucket > TRACE_CONCENTRATION_MOL_PER_BUCKET {
            Ok(Self::Heterogeneous {
                dominant_phase,
                total_mol_per_bucket,
            })
        } else {
            Ok(Self::Vacuum)
        }
    }

    pub fn solvent_type(&self) -> SolventType {
        match self {
            ReactionMedium::Liquid { solvent_type, .. }
            | ReactionMedium::SupercriticalFluid { solvent_type, .. } => *solvent_type,
            ReactionMedium::GasLiquidInterface { liquid, .. } => liquid.solvent_type(),
            ReactionMedium::Vacuum
            | ReactionMedium::Gas { .. }
            | ReactionMedium::MoltenMetal { .. }
            | ReactionMedium::MoltenSlag { .. }
            | ReactionMedium::SolidSurface { .. }
            | ReactionMedium::Heterogeneous { .. } => SolventType::Neutral,
        }
    }

    pub fn is_supercritical(&self) -> bool {
        matches!(self, ReactionMedium::SupercriticalFluid { .. })
    }

    pub fn has_gas_liquid_interface(&self) -> bool {
        matches!(self, ReactionMedium::GasLiquidInterface { .. })
    }
}

/// Solvent classification for selectivity
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolventType {
    /// Protic polar (water, alcohols) - favors SN1, E1
    Protic,
    /// Aprotic polar (DMSO, DMF) - favors SN2
    AproticPolar,
    /// Non-polar (hexane, toluene)
    NonPolar,
    /// Basic conditions
    Basic,
    /// Acidic conditions
    Acidic,
    /// Neutral/default
    Neutral,
}

/// Type of organic reaction mechanism
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReactionType {
    /// SN2 bimolecular substitution
    SN2,
    /// SN1 unimolecular substitution
    SN1,
    /// E2 bimolecular elimination
    E2,
    /// E1 unimolecular elimination
    E1,
    /// Acid-catalyzed esterification (Fischer)
    FischerEsterification,
    /// Nucleophilic addition to carbonyl
    CarbonylAddition,
    /// Hydride reduction of aldehydes and ketones to alcohols
    CarbonylReduction,
    /// Oxygen-transfer or dehydrogenating oxidation of an organic reactive site
    OrganicOxidation,
    /// Electrophilic addition to alkene/alkyne
    ElectrophilicAddition,
    /// Halogenation at the alpha carbon of an enolizable carbonyl
    AlphaHalogenation,
    /// Carbon-carbon bond formation between an enol/enolate and carbonyl
    AldolAddition,
    /// Dehydration of beta-hydroxy carbonyls into enones/enals
    AldolDehydration,
    /// Enamine formation from an enolizable carbonyl and secondary amine
    EnamineFormation,
    /// Carbon-carbon bond formation between an enolate and alkyl halide
    EnolateAlkylation,
    /// Conjugate addition of an enolate to an activated alkene
    MichaelAddition,
    /// Condensation of an ester enolate with an ester
    ClaisenCondensation,
    /// Formation of a phosphonium salt from phosphine and an alkyl halide
    PhosphoniumSaltFormation,
    /// Base-induced phosphonium ylide formation
    PhosphoniumYlideFormation,
    /// Carbonyl olefination through a phosphonium ylide
    WittigOlefination,
    /// Carbonyl olefination through a phosphonate carbanion
    HornerWadsworthEmmonsOlefination,
    /// Carbonyl olefination through a sulfone carbanion
    JuliaOlefination,
    /// Protection of an alcohol as a silyl ether
    SilylEtherFormation,
    /// Deprotection of a silyl ether
    SilylEtherCleavage,
    /// Formation of an acetal or ketal from a carbonyl and alcohol
    AcetalFormation,
    /// Hydrolysis of an acetal or ketal
    AcetalHydrolysis,
    /// Formation of a Boc or Cbz carbamate
    CarbamateFormation,
    /// Acidic or hydrogenolytic carbamate cleavage
    CarbamateCleavage,
    /// Formation of an ester protecting group
    EsterProtection,
    /// Hydrolysis or cleavage of an ester protecting group
    EsterHydrolysis,
    /// Nucleophilic acyl substitution of acid chlorides, anhydrides, esters and amides
    AcylSubstitution,
    /// Intramolecular esterification closing a lactone ring
    Lactonization,
    /// Intramolecular amide formation closing a lactam ring
    Lactamization,
    /// Dehydrative condensation closing an aromatic heterocycle
    /// (Paal–Knorr family: 1,4-dicarbonyl + heteroatom donor → pyrrole/furan/thiophene)
    HeterocycleCondensation,
    /// [4+2] cycloaddition of a conjugated diene and a dienophile
    DielsAlder,
    /// Thermal cycloreversion of a cyclohexene-like Diels-Alder adduct
    RetroDielsAlder,
    /// Light-driven double-bond geometric isomerization
    PhotochemicalIsomerization,
    /// Alkylation of an amine or amide nitrogen by an alkyl halide
    NAlkylation,
    /// Radical-chain substitution of an abstractable C-H bond
    RadicalHalogenation,
    /// Thermal or catalytic C-C bond scission in hydrocarbons.
    HydrocarbonCracking,
    /// Thermal dehydrogenation of hydrocarbons into more unsaturated products.
    HydrocarbonPyrolysis,
    /// Intramolecular migration with sigma-bond reorganization.
    SkeletalRearrangement,
    /// Chain-growth (addition) polymerization of an alkene into a repeat unit.
    ChainGrowthPolymerization,
}

/// Nucleophile strength classification for selectivity models.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleophileStrength {
    /// Grignard, organolithium and comparable very strong nucleophiles.
    VeryStrong,
    /// Borohydride-like or strong anionic nucleophiles.
    Strong,
    /// Hydroxide, cyanide, amines.
    Moderate,
    /// Water, alcohols and other weak neutral nucleophiles.
    Weak,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectivitySuppressionPolicy {
    SuppressWhenDisfavored,
    NeverSuppress,
}

#[derive(Debug, Clone)]
pub struct SelectivityProfile {
    pub mechanism: ReactionType,
    pub primary_site: SiteDescriptor,
    pub secondary_site: Option<SiteDescriptor>,
    pub nucleophile_strength: Option<NucleophileStrength>,
    pub suppression_policy: SelectivitySuppressionPolicy,
}

#[derive(Debug, Clone)]
pub struct SelectivityRuntimeEffect {
    pub rate_multiplier: f64,
    pub activation_delta_kj_per_mol: f64,
    pub pre_exp_multiplier: f64,
    pub suppressed: bool,
    pub reason: String,
}

impl SelectivityProfile {
    pub fn new(mechanism: ReactionType, primary_site: SiteDescriptor) -> Self {
        Self {
            mechanism,
            primary_site,
            secondary_site: None,
            nucleophile_strength: None,
            suppression_policy: SelectivitySuppressionPolicy::SuppressWhenDisfavored,
        }
    }

    pub fn with_secondary_site(mut self, site: SiteDescriptor) -> Self {
        self.secondary_site = Some(site);
        self
    }

    pub fn with_nucleophile_strength(mut self, strength: NucleophileStrength) -> Self {
        self.nucleophile_strength = Some(strength);
        self
    }

    pub fn never_suppress(mut self) -> Self {
        self.suppression_policy = SelectivitySuppressionPolicy::NeverSuppress;
        self
    }
}

/// Reactivity score with metadata
#[derive(Debug, Clone)]
pub struct ReactivityScore {
    /// Primary score (0.0 to 1.0+, >1.0 means faster than reference)
    pub value: f64,
    /// Activation energy modification in kJ/mol (negative = faster)
    pub activation_delta: f64,
    /// Pre-exponential factor multiplier
    pub pre_exp_multiplier: f64,
    /// Competing mechanisms with their scores
    pub competing: Vec<(CompetingMechanism, f64)>,
    /// Reasoning for this score
    pub reason: String,
}

impl ReactivityScore {
    /// Create new score with default competing mechanisms
    pub fn new(value: f64, reason: impl Into<String>) -> Self {
        Self {
            value: value.max(0.0),
            activation_delta: Self::value_to_activation_delta(value),
            pre_exp_multiplier: 1.0,
            competing: Vec::new(),
            reason: reason.into(),
        }
    }

    /// Create score from explicit activation energy delta
    pub fn with_activation_delta(delta: f64, reason: impl Into<String>) -> Self {
        Self {
            value: Self::activation_delta_to_value(delta),
            activation_delta: delta,
            pre_exp_multiplier: 1.0,
            competing: Vec::new(),
            reason: reason.into(),
        }
    }

    /// Add competing mechanism
    pub fn with_competing(mut self, mechanism: CompetingMechanism, score: f64) -> Self {
        self.competing.push((mechanism, score));
        self
    }

    /// Set pre-exponential multiplier
    pub fn with_pre_exp_multiplier(mut self, multiplier: f64) -> Self {
        self.pre_exp_multiplier = multiplier;
        self
    }

    /// Convert relative value to activation energy delta
    /// Uses approximation: ΔEa ≈ -RT ln(k/k0) at 298K
    fn value_to_activation_delta(value: f64) -> f64 {
        const R: f64 = 8.314; // J/(mol·K)
        const T: f64 = 298.15; // K
        if value <= 0.0 {
            return 50.0; // Very slow
        }
        -R * T * value.ln() / 1000.0 // Convert to kJ/mol
    }

    /// Convert activation energy delta to relative value
    fn activation_delta_to_value(delta: f64) -> f64 {
        const R: f64 = 8.314;
        const T: f64 = 298.15;
        (-delta * 1000.0 / (R * T)).exp()
    }
}

/// Competing reaction mechanisms
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompetingMechanism {
    None,
    SN1,
    SN2,
    E1,
    E2,
    /// General elimination (E1/E2 unspecified)
    Elimination,
    /// General substitution
    Substitution,
    Rearrangement,
}

/// Selectivity evaluation recommendation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectivityRecommendation {
    /// Reaction proceeds exclusively (ratio > 10:1)
    Exclusive,
    /// Reaction is preferred but mixture expected (ratio 3:1 to 10:1)
    Preferred,
    /// Significant mixture of products (ratio 1:3 to 3:1)
    Mixed,
    /// Reaction is suppressed (ratio < 1:3)
    Suppressed,
    /// Reaction does not occur (score effectively zero)
    None,
}

/// Complete result of selectivity evaluation
#[derive(Debug, Clone)]
pub struct SelectivityResult {
    /// Primary mechanism score
    pub primary: ReactivityScore,
    /// All competing mechanisms with scores
    pub all_scores: HashMap<ReactionType, ReactivityScore>,
    /// Final recommendation
    pub recommendation: SelectivityRecommendation,
    /// Dominant competing mechanism (if any)
    pub dominant_competitor: Option<CompetingMechanism>,
}

impl SelectivityResult {
    /// Create result from primary score alone (no competition)
    pub fn exclusive(score: ReactivityScore) -> Self {
        let mut all_scores = HashMap::new();
        all_scores.insert(ReactionType::SN2, score.clone());

        Self {
            primary: score,
            all_scores,
            recommendation: SelectivityRecommendation::Exclusive,
            dominant_competitor: None,
        }
    }

    /// Determine recommendation from competing scores
    pub fn from_scores(
        primary_type: ReactionType,
        primary: ReactivityScore,
        competitors: Vec<(ReactionType, ReactivityScore)>,
    ) -> Self {
        let mut all_scores = HashMap::new();
        all_scores.insert(primary_type, primary.clone());

        let mut max_competitor_score = 0.0;
        let mut dominant_competitor = None;

        for (mech, score) in &competitors {
            all_scores.insert(*mech, score.clone());
            if score.value > max_competitor_score {
                max_competitor_score = score.value;
                dominant_competitor = CompetingMechanism::from_reaction_type(*mech);
            }
        }

        let ratio = if max_competitor_score > 0.0 {
            primary.value / max_competitor_score
        } else {
            f64::INFINITY
        };

        let recommendation = if primary.value < 0.01 {
            SelectivityRecommendation::None
        } else if ratio >= 10.0 {
            SelectivityRecommendation::Exclusive
        } else if ratio >= 3.0 {
            SelectivityRecommendation::Preferred
        } else if ratio >= 0.5 {
            SelectivityRecommendation::Mixed
        } else if ratio >= 0.1 {
            SelectivityRecommendation::Suppressed
        } else {
            SelectivityRecommendation::None
        };

        Self {
            primary,
            all_scores,
            recommendation,
            dominant_competitor,
        }
    }
}

impl CompetingMechanism {
    fn from_reaction_type(rt: ReactionType) -> Option<Self> {
        match rt {
            ReactionType::SN1 => Some(Self::SN1),
            ReactionType::SN2 => Some(Self::SN2),
            ReactionType::E1 => Some(Self::E1),
            ReactionType::E2 => Some(Self::E2),
            _ => None,
        }
    }
}

fn dominant_liquid_medium(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
) -> ChemistryResult<Option<ReactionMedium>> {
    let mut selected = None;
    for phase in mixture.liquid_phase_snapshots(registry)? {
        if phase.total_solvent_mol_per_bucket <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            continue;
        }
        if matches!(
            phase.coarse_phase,
            MixturePhase::MoltenMetal | MixturePhase::MoltenSlag
        ) {
            continue;
        }
        if selected
            .as_ref()
            .is_some_and(|(_, amount)| *amount >= phase.total_solvent_mol_per_bucket)
        {
            continue;
        }
        let water_mole_fraction = phase
            .solvents
            .iter()
            .find(|amount| amount.substance_id.as_str() == "destroy:water")
            .map(|amount| amount.mole_fraction)
            .unwrap_or(0.0);
        let total_solvent_mol_per_bucket = phase.total_solvent_mol_per_bucket;
        selected = Some((
            ReactionMedium::Liquid {
                coarse_phase: phase.coarse_phase,
                representative_solvent_id: phase.representative_solvent_id.clone(),
                solvent_type: solvent_type_from_liquid_phase(registry, &phase)?,
                total_solvent_mol_per_bucket,
                water_mole_fraction,
            },
            total_solvent_mol_per_bucket,
        ));
    }
    Ok(selected.map(|(medium, _)| medium))
}

fn solvent_type_from_liquid_phase(
    registry: &ChemistryRegistry,
    phase: &LiquidPhaseSnapshot,
) -> ChemistryResult<SolventType> {
    if phase.coarse_phase == MixturePhase::Aqueous {
        return Ok(SolventType::Protic);
    }
    let mut protic = 0.0;
    let mut polar_aprotic = 0.0;
    let mut nonpolar = 0.0;
    for solvent in &phase.solvents {
        let substance = registry.substance(&solvent.substance_id)?;
        if substance.phase_properties.solvent_role == SolventRole::NotSolvent {
            continue;
        }
        match substance.phase_properties.preferred_liquid_phase {
            LiquidPhasePreference::Aqueous => protic += solvent.concentration_mol_per_bucket,
            LiquidPhasePreference::Organic => {
                if organic_solvent_looks_protic(solvent.substance_id.as_str()) {
                    protic += solvent.concentration_mol_per_bucket;
                } else if organic_solvent_looks_polar_aprotic(solvent.substance_id.as_str()) {
                    polar_aprotic += solvent.concentration_mol_per_bucket;
                } else {
                    nonpolar += solvent.concentration_mol_per_bucket;
                }
            }
            LiquidPhasePreference::MoltenMetal | LiquidPhasePreference::MoltenSlag => {}
        }
    }
    if protic >= polar_aprotic && protic >= nonpolar && protic > TRACE_CONCENTRATION_MOL_PER_BUCKET
    {
        Ok(SolventType::Protic)
    } else if polar_aprotic >= nonpolar && polar_aprotic > TRACE_CONCENTRATION_MOL_PER_BUCKET {
        Ok(SolventType::AproticPolar)
    } else if nonpolar > TRACE_CONCENTRATION_MOL_PER_BUCKET {
        Ok(SolventType::NonPolar)
    } else {
        Ok(SolventType::Neutral)
    }
}

fn dominant_substance_in_phase(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    phase: MixturePhase,
) -> ChemistryResult<Option<SubstanceId>> {
    let mut selected = None;
    for substance_id in mixture.substances() {
        registry.substance(substance_id)?;
        let amount = mixture.concentration_in_phase(substance_id, phase);
        if amount <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            continue;
        }
        if selected
            .as_ref()
            .is_none_or(|(_, current)| amount > *current)
        {
            selected = Some((substance_id.clone(), amount));
        }
    }
    Ok(selected.map(|(id, _)| id))
}

fn supercritical_solvent_type(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
) -> ChemistryResult<SolventType> {
    let Some(carrier) =
        dominant_substance_in_phase(registry, mixture, MixturePhase::SupercriticalFluid)?
    else {
        return Ok(SolventType::Neutral);
    };
    if carrier.as_str() == "destroy:water" {
        return Ok(SolventType::Protic);
    }
    if matches!(
        carrier.as_str(),
        "destroy:carbon_dioxide" | "destroy:nitrogen" | "destroy:methane"
    ) {
        return Ok(SolventType::NonPolar);
    }
    let substance = registry.substance(&carrier)?;
    match substance.phase_properties.preferred_liquid_phase {
        LiquidPhasePreference::Aqueous => Ok(SolventType::Protic),
        LiquidPhasePreference::Organic => {
            if organic_solvent_looks_polar_aprotic(carrier.as_str()) {
                Ok(SolventType::AproticPolar)
            } else if organic_solvent_looks_protic(carrier.as_str()) {
                Ok(SolventType::Protic)
            } else {
                Ok(SolventType::NonPolar)
            }
        }
        LiquidPhasePreference::MoltenMetal | LiquidPhasePreference::MoltenSlag => {
            Ok(SolventType::Neutral)
        }
    }
}

fn dominant_phase_by_amount(mixture: &Mixture) -> (MixturePhase, f64) {
    MixturePhase::ALL
        .iter()
        .copied()
        .map(|phase| (phase, mixture.total_concentration_in_phase(phase)))
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .unwrap_or((MixturePhase::Aqueous, 0.0))
}

fn max_known_activity(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    substance_id: &SubstanceId,
    phases: &[MixturePhase],
) -> ChemistryResult<f64> {
    if registry.substance(substance_id).is_err() {
        return Ok(0.0);
    }
    phases.iter().try_fold(0.0_f64, |current, phase| {
        Ok(current.max(mixture.activity_of(registry, substance_id, *phase)?))
    })
}

fn known_partial_pressure(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
    substance_id: &SubstanceId,
) -> f64 {
    if registry.substance(substance_id).is_err() {
        0.0
    } else {
        mixture.gas_partial_pressure_pascal(substance_id)
    }
}

fn redox_role_strengths(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
) -> ChemistryResult<(f64, f64)> {
    let mut oxidizing_strength = 0.0;
    let mut reducing_strength = 0.0;
    for substance_id in mixture.substances() {
        let activity = max_known_activity(
            registry,
            mixture,
            substance_id,
            &[
                MixturePhase::Aqueous,
                MixturePhase::Organic,
                MixturePhase::Gas,
            ],
        )?;
        if activity <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
            continue;
        }
        let substance = registry.substance(substance_id)?;
        for role in &substance.redox_roles {
            match role {
                RedoxRole::Oxidant => oxidizing_strength += activity,
                RedoxRole::Reductant => reducing_strength += activity,
                RedoxRole::OxidantAndReductant => {
                    oxidizing_strength += activity;
                    reducing_strength += activity;
                }
            }
        }
    }
    Ok((oxidizing_strength, reducing_strength))
}

fn metal_availability(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
) -> ChemistryResult<(f64, f64)> {
    let mut free_metal_activity = 0.0_f64;
    let mut complexed_metal = 0.0_f64;
    let mut seen_central_ions = std::collections::BTreeSet::new();
    for spec in registry.complex_specs() {
        if seen_central_ions.insert(spec.central_ion.clone()) {
            free_metal_activity = free_metal_activity.max(max_known_activity(
                registry,
                mixture,
                &spec.central_ion,
                &[MixturePhase::Aqueous, MixturePhase::Organic],
            )?);
        }
        complexed_metal += mixture.concentration_of(&spec.id);
    }
    Ok((free_metal_activity, complexed_metal))
}

fn organic_solvent_looks_polar_aprotic(id: &str) -> bool {
    matches!(
        id,
        "destroy:acetone"
            | "destroy:acetonitrile"
            | "destroy:dimethyl_sulfoxide"
            | "destroy:dmf"
            | "destroy:ethyl_acetate"
    ) || id.contains("acetone")
        || id.contains("acetonitrile")
        || id.contains("sulfoxide")
        || id.contains("dmf")
        || id.contains("ether")
        || id.contains("ester")
}

fn organic_solvent_looks_protic(id: &str) -> bool {
    id.contains("alcohol")
        || id.contains("ethanol")
        || id.contains("methanol")
        || id.contains("propanol")
        || id.contains("butanol")
        || id.contains("phenol")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::complex::{ComplexLigand, ComplexSpec};
    use crate::chemistry::mixture::Mixture;
    use crate::chemistry::mixture::STANDARD_PRESSURE_PASCAL;
    use crate::chemistry::redox::RedoxRole;
    use crate::chemistry::simulation::ReactionContext;
    use crate::chemistry::substance::{
        LiquidPhasePreference, SolventRole, Substance, SubstancePhaseProperties,
    };
    use crate::chemistry::ChemistryRegistryBuilder;

    #[test]
    fn substitution_degree_base_scores() {
        assert!(
            SubstitutionDegree::Primary.base_steric_score()
                < SubstitutionDegree::Secondary.base_steric_score()
        );
        assert!(
            SubstitutionDegree::Secondary.base_steric_score()
                < SubstitutionDegree::Tertiary.base_steric_score()
        );
        assert!(
            SubstitutionDegree::Benzylic.base_steric_score()
                < SubstitutionDegree::Allylic.base_steric_score()
        );
    }

    #[test]
    fn site_descriptor_steric_calculation() {
        let desc = SiteDescriptor::new(
            crate::chemistry::reactive_site::ReactiveSiteKind::Alcohol,
            SubstitutionDegree::Secondary,
            ElectronicEnvironment::default(),
            2, // two bulky groups
            true,
        );
        // Base 0.3 + 0.2 = 0.5
        assert!((desc.steric_score - 0.5).abs() < 0.01);
        assert!((desc.steric_accessibility() - 0.5).abs() < 0.01);
    }

    #[test]
    fn reactivity_score_activation_conversion() {
        let score = ReactivityScore::new(1.0, "reference");
        assert!(score.activation_delta.abs() < 0.1); // ~0 for k=k0

        let fast = ReactivityScore::new(2.0, "2x faster");
        assert!(fast.activation_delta < 0.0); // negative = faster

        let slow = ReactivityScore::new(0.5, "half speed");
        assert!(slow.activation_delta > 0.0); // positive = slower
    }

    #[test]
    fn selectivity_result_ratios() {
        let primary = ReactivityScore::new(1.0, "primary");
        let competitor = ReactivityScore::new(0.1, "competitor");

        let result = SelectivityResult::from_scores(
            ReactionType::SN2,
            primary,
            vec![(ReactionType::E2, competitor)],
        );

        assert_eq!(result.recommendation, SelectivityRecommendation::Exclusive);
    }

    #[test]
    fn selectivity_context_uses_mixture_temperature_ph_and_solvent() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1000.0, 373.15, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(Substance::new(
                "destroy:proton",
                1,
                1.0,
                1000.0,
                f64::MAX,
                0.0,
                0.0,
            ))
            .build()
            .unwrap();
        let mut mixture = Mixture::new(320.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:proton", 0.1)
            .unwrap();

        let context =
            SelectivityContext::from_mixture(&registry, &mixture, &ReactionContext::default())
                .unwrap();
        assert!(
            (context.temperature - 320.0).abs() < 5.0,
            "temperature should be near 320K, got {}",
            context.temperature
        );
        assert!(context.is_acidic());
        assert_eq!(context.solvent_type, SolventType::Acidic);
    }

    #[test]
    fn selectivity_context_exposes_gas_liquid_reaction_medium() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1000.0, 373.15, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new("destroy:oxygen", 0, 32.0, 1_140.0, 90.0, 29.4, 6_820.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: Some(0.0),
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: false,
                        can_form_liquid_phase: false,
                        solvent_role: SolventRole::NotSolvent,
                    }),
            )
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.15).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:oxygen", 0.5)
            .unwrap();

        let context =
            SelectivityContext::from_mixture(&registry, &mixture, &ReactionContext::default())
                .unwrap();

        assert!(context.medium.has_gas_liquid_interface());
        assert_eq!(context.solvent_type, SolventType::Protic);
    }

    #[test]
    fn selectivity_context_exposes_supercritical_medium() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new(
                    "destroy:carbon_dioxide",
                    0,
                    44.0,
                    1_100.0,
                    194.7,
                    37.0,
                    16_000.0,
                )
                .with_critical_point(304.1, 7_377_000.0)
                .with_phase_properties(SubstancePhaseProperties {
                    preferred_liquid_phase: LiquidPhasePreference::Organic,
                    aqueous_solubility_mol_per_bucket: Some(0.0),
                    organic_solubility_mol_per_bucket: Some(0.0),
                    can_precipitate: false,
                    can_form_liquid_phase: false,
                    solvent_role: SolventRole::NotSolvent,
                }),
            )
            .build()
            .unwrap();
        let mut mixture = Mixture::new(310.0).unwrap();
        mixture
            .add_substance(&registry, "destroy:carbon_dioxide", 3.0)
            .unwrap();

        let context =
            SelectivityContext::from_mixture(&registry, &mixture, &ReactionContext::default())
                .unwrap();

        assert!(context.medium.is_supercritical());
        assert_eq!(context.solvent_type, SolventType::NonPolar);
    }

    #[test]
    fn selectivity_context_uses_gases_light_surfaces_and_redox_roles() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1000.0, 373.15, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(
                Substance::new("destroy:hydrogen", 0, 2.0, 70.0, 20.0, 28.0, 900.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: Some(0.0),
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: false,
                        can_form_liquid_phase: false,
                        solvent_role: SolventRole::NotSolvent,
                    })
                    .with_redox_roles(vec![RedoxRole::Reductant]),
            )
            .substance(
                Substance::new("destroy:oxygen", 0, 32.0, 1_140.0, 90.0, 29.4, 6_820.0)
                    .with_phase_properties(SubstancePhaseProperties {
                        preferred_liquid_phase: LiquidPhasePreference::Aqueous,
                        aqueous_solubility_mol_per_bucket: Some(0.0),
                        organic_solubility_mol_per_bucket: Some(0.0),
                        can_precipitate: false,
                        can_form_liquid_phase: false,
                        solvent_role: SolventRole::NotSolvent,
                    })
                    .with_redox_roles(vec![RedoxRole::Oxidant]),
            )
            .build()
            .unwrap();
        let mut mixture = Mixture::new(298.15).unwrap();
        mixture
            .add_substance(&registry, "destroy:water", 1.0)
            .unwrap();
        mixture
            .exchange_gases_with_atmosphere(
                &registry,
                &[
                    (SubstanceId::from("destroy:hydrogen"), 0.8),
                    (SubstanceId::from("destroy:oxygen"), 0.2),
                ],
                STANDARD_PRESSURE_PASCAL,
                10.0,
                1.0,
            )
            .unwrap();
        let context = ReactionContext::default()
            .with_uv_power(0.5)
            .unwrap()
            .with_open_atmosphere(
                [(SubstanceId::from("destroy:hydrogen"), 1.0)],
                STANDARD_PRESSURE_PASCAL,
                1.0,
            )
            .unwrap();
        let mut context = context;
        context.add_external_catalyst("palladium", 0.1).unwrap();

        let selectivity = SelectivityContext::from_mixture(&registry, &mixture, &context).unwrap();

        assert!(selectivity.has_hydrogen());
        assert!(selectivity.is_oxygen_rich());
        assert!(selectivity.has_uv());
        assert!(selectivity.palladium_available);
        assert!(selectivity.has_available_surface());
        assert!(selectivity.total_gas_pressure_pascal > 0.0);
        assert!(selectivity.total_light_power >= 0.5);
        assert!(selectivity.oxidizing_strength > 0.0);
        assert!(selectivity.reducing_strength > 0.0);
    }

    #[test]
    fn selectivity_context_distinguishes_free_and_complexed_metals() {
        let registry = ChemistryRegistryBuilder::new()
            .substance(
                Substance::new("destroy:water", 0, 18.0, 1000.0, 373.15, 75.0, 40_000.0)
                    .with_phase_properties(SubstancePhaseProperties::aqueous_solvent()),
            )
            .substance(Substance::new(
                "metal:copper_ion",
                2,
                63.55,
                1000.0,
                f64::MAX,
                50.0,
                0.0,
            ))
            .substance(Substance::new(
                "ligand:ammonia",
                0,
                17.0,
                680.0,
                240.0,
                35.0,
                20_000.0,
            ))
            .complex_spec(ComplexSpec::new(
                "complex:copper_tetraammine",
                "metal:copper_ion",
                [ComplexLigand::new("ligand:ammonia", 4)],
                2,
                1.0e8,
            ))
            .build()
            .unwrap();

        let mut free = Mixture::new(298.15).unwrap();
        free.add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        free.add_substance(&registry, "metal:copper_ion", 0.05)
            .unwrap();
        let free_context =
            SelectivityContext::from_mixture(&registry, &free, &ReactionContext::default())
                .unwrap();

        let mut bound = Mixture::new(298.15).unwrap();
        bound
            .add_substance(&registry, "destroy:water", 55.5)
            .unwrap();
        bound
            .add_substance(&registry, "complex:copper_tetraammine", 0.05)
            .unwrap();
        let bound_context =
            SelectivityContext::from_mixture(&registry, &bound, &ReactionContext::default())
                .unwrap();

        assert!(free_context.has_free_complexable_metal());
        assert!(!free_context.has_complexed_metal());
        assert!(!bound_context.has_free_complexable_metal());
        assert!(bound_context.has_complexed_metal());
    }
}
