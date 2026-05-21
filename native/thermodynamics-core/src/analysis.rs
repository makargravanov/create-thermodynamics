use crate::activity::{davies_log10_gamma, DAVIES_MAX_IONIC_STRENGTH_MOLAL};
use crate::chemistry::{PhaseKind, SpeciesAmount, SpeciesId};
use crate::equilibrium::EquilibriumResult;
use crate::registry::SpeciesRegistry;

const WATER_MOLAR_MASS_KILOGRAM_PER_MOL: f64 = 0.018_015_28;
const MIN_ACTIVITY: f64 = 1.0e-300;

#[derive(Debug, Clone, PartialEq)]
pub struct AqueousEquilibriumSummary {
    pub solvent_water_mass_kg: f64,
    pub ionic_strength_molal: f64,
    pub hydrogen_activity: f64,
    pub ph: f64,
    pub aqueous_species: Vec<AqueousSpeciesSummary>,
    pub solid_species_amounts_mol: Vec<SpeciesAmount>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AqueousSpeciesSummary {
    pub species_id: SpeciesId,
    pub amount_mol: f64,
    pub molality_mol_per_kg_water: f64,
    pub activity_coefficient: f64,
    pub activity: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum EquilibriumAnalysisError {
    MissingSpeciesData(SpeciesId),
    MissingWaterSolvent(SpeciesId),
    MissingHydrogenIon(SpeciesId),
    NonPositiveWaterAmount {
        species_id: SpeciesId,
        amount_mol: f64,
    },
    NonAqueousWaterSpecies {
        species_id: SpeciesId,
    },
    NonAqueousHydrogenIonSpecies {
        species_id: SpeciesId,
    },
    DaviesModelOutOfRange {
        ionic_strength_molal: f64,
        max_valid_ionic_strength_molal: f64,
    },
}

pub fn analyze_aqueous_equilibrium(
    registry: &SpeciesRegistry,
    result: &EquilibriumResult,
    water_species_id: SpeciesId,
    hydrogen_ion_species_id: SpeciesId,
) -> Result<AqueousEquilibriumSummary, EquilibriumAnalysisError> {
    let water_species =
        registry
            .species(water_species_id)
            .ok_or(EquilibriumAnalysisError::MissingWaterSolvent(
                water_species_id,
            ))?;
    if water_species.phase != PhaseKind::Aqueous {
        return Err(EquilibriumAnalysisError::NonAqueousWaterSpecies {
            species_id: water_species_id,
        });
    }

    let hydrogen_species = registry.species(hydrogen_ion_species_id).ok_or(
        EquilibriumAnalysisError::MissingHydrogenIon(hydrogen_ion_species_id),
    )?;
    if hydrogen_species.phase != PhaseKind::Aqueous {
        return Err(EquilibriumAnalysisError::NonAqueousHydrogenIonSpecies {
            species_id: hydrogen_ion_species_id,
        });
    }

    let water_amount_mol = result_amount(result, water_species_id);
    if !water_amount_mol.is_finite() || water_amount_mol <= 0.0 {
        return Err(EquilibriumAnalysisError::NonPositiveWaterAmount {
            species_id: water_species_id,
            amount_mol: water_amount_mol,
        });
    }

    let solvent_water_mass_kg = water_amount_mol * WATER_MOLAR_MASS_KILOGRAM_PER_MOL;
    let ionic_strength_molal = ionic_strength_molal(registry, result, solvent_water_mass_kg)?;
    if ionic_strength_molal > DAVIES_MAX_IONIC_STRENGTH_MOLAL {
        return Err(EquilibriumAnalysisError::DaviesModelOutOfRange {
            ionic_strength_molal,
            max_valid_ionic_strength_molal: DAVIES_MAX_IONIC_STRENGTH_MOLAL,
        });
    }

    let mut aqueous_species = Vec::new();
    let mut solid_species_amounts_mol = Vec::new();
    for amount in &result.species_amounts_mol {
        let species = registry.species(amount.species_id).ok_or(
            EquilibriumAnalysisError::MissingSpeciesData(amount.species_id),
        )?;
        match species.phase {
            PhaseKind::Aqueous => {
                if amount.species_id == water_species_id {
                    continue;
                }
                let molality_mol_per_kg_water = amount.amount_mol.max(0.0) / solvent_water_mass_kg;
                let log10_gamma = davies_log10_gamma(species.charge_number, ionic_strength_molal)
                    .expect("ionic strength range checked before activity calculation");
                let activity_coefficient = 10.0_f64.powf(log10_gamma);
                aqueous_species.push(AqueousSpeciesSummary {
                    species_id: amount.species_id,
                    amount_mol: amount.amount_mol,
                    molality_mol_per_kg_water,
                    activity_coefficient,
                    activity: (molality_mol_per_kg_water * activity_coefficient).max(MIN_ACTIVITY),
                });
            }
            PhaseKind::Solid => solid_species_amounts_mol.push(*amount),
            PhaseKind::Gas => {}
        }
    }

    let hydrogen_activity = aqueous_species
        .iter()
        .find(|summary| summary.species_id == hydrogen_ion_species_id)
        .map(|summary| summary.activity)
        .unwrap_or(MIN_ACTIVITY);

    Ok(AqueousEquilibriumSummary {
        solvent_water_mass_kg,
        ionic_strength_molal,
        hydrogen_activity,
        ph: -hydrogen_activity.log10(),
        aqueous_species,
        solid_species_amounts_mol,
    })
}

fn result_amount(result: &EquilibriumResult, species_id: SpeciesId) -> f64 {
    result
        .species_amounts_mol
        .iter()
        .find(|amount| amount.species_id == species_id)
        .map(|amount| amount.amount_mol)
        .unwrap_or_default()
}

fn ionic_strength_molal(
    registry: &SpeciesRegistry,
    result: &EquilibriumResult,
    solvent_water_mass_kg: f64,
) -> Result<f64, EquilibriumAnalysisError> {
    let mut ionic_strength = 0.0;
    for amount in &result.species_amounts_mol {
        let species = registry.species(amount.species_id).ok_or(
            EquilibriumAnalysisError::MissingSpeciesData(amount.species_id),
        )?;
        if species.phase != PhaseKind::Aqueous {
            continue;
        }
        let charge = f64::from(species.charge_number);
        let molality = amount.amount_mol.max(0.0) / solvent_water_mass_kg;
        ionic_strength += 0.5 * molality * charge * charge;
    }
    Ok(ionic_strength)
}
