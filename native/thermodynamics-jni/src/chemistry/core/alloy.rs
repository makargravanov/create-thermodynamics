use super::error::{ChemistryError, ChemistryResult};
use super::mixture::{LiquidPhaseId, Mixture, MixturePhase, TRACE_CONCENTRATION_MOL_PER_BUCKET};
use super::registry::ChemistryRegistry;
use super::substance::{Substance, SubstanceId, SubstanceRepresentation};

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum AlloyConstituentRole {
    Metal,
    DissolvedMaterial,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlloyConstituent {
    pub substance_id: SubstanceId,
    pub metallurgical_component_id: String,
    pub role: AlloyConstituentRole,
    pub concentration_mol_per_bucket: f64,
    pub mole_fraction: f64,
    pub mass_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct AlloyPhaseSnapshot {
    pub phase_id: LiquidPhaseId,
    pub representative_substance_id: SubstanceId,
    pub temperature_kelvin: f64,
    pub total_mol_per_bucket: f64,
    pub total_mass_grams_per_bucket: f64,
    pub density_grams_per_bucket: f64,
    pub average_molar_mass_grams: f64,
    pub molar_heat_capacity_j_per_mol_kelvin: f64,
    pub estimated_solidus_kelvin: f64,
    pub estimated_liquidus_kelvin: f64,
    pub constituents: Vec<AlloyConstituent>,
}

pub fn alloy_phase_snapshots(
    registry: &ChemistryRegistry,
    mixture: &Mixture,
) -> ChemistryResult<Vec<AlloyPhaseSnapshot>> {
    let phases = mixture
        .liquid_phase_snapshots(registry)?
        .into_iter()
        .filter(|phase| phase.coarse_phase == MixturePhase::MoltenMetal)
        .collect::<Vec<_>>();
    let mut snapshots = Vec::new();
    for phase in phases {
        let mut raw_constituents = Vec::new();
        for substance in registry.substances() {
            let amount = mixture
                .liquid_phase_amounts_of(registry, &substance.id)?
                .into_iter()
                .find(|amount| amount.phase_id == phase.id)
                .map(|amount| amount.concentration_mol_per_bucket)
                .unwrap_or(0.0);
            if amount <= TRACE_CONCENTRATION_MOL_PER_BUCKET {
                continue;
            }
            raw_constituents.push((substance, amount));
        }
        if raw_constituents.is_empty() {
            continue;
        }
        snapshots.push(build_alloy_snapshot(
            phase.id,
            phase.representative_solvent_id,
            mixture.temperature_kelvin(),
            raw_constituents,
        )?);
    }
    Ok(snapshots)
}

fn build_alloy_snapshot(
    phase_id: LiquidPhaseId,
    representative_substance_id: SubstanceId,
    temperature_kelvin: f64,
    raw_constituents: Vec<(&Substance, f64)>,
) -> ChemistryResult<AlloyPhaseSnapshot> {
    let total_mol_per_bucket = raw_constituents
        .iter()
        .map(|(_, amount)| *amount)
        .sum::<f64>();
    validate_positive_finite(total_mol_per_bucket, "alloy total amount")?;

    let mut total_mass_grams_per_bucket = 0.0;
    let mut volume_buckets = 0.0;
    let mut molar_heat_capacity = 0.0;
    let mut estimated_solidus_kelvin = f64::INFINITY;
    let mut estimated_liquidus_kelvin = 0.0_f64;
    for (substance, amount) in &raw_constituents {
        validate_positive_finite(*amount, "alloy constituent amount")?;
        validate_positive_finite(substance.molar_mass_grams, "alloy constituent molar mass")?;
        validate_positive_finite(
            substance.liquid_density_grams_per_bucket,
            "alloy constituent liquid density",
        )?;
        validate_positive_finite(
            substance.molar_heat_capacity_j_per_mol_kelvin,
            "alloy constituent heat capacity",
        )?;
        validate_non_negative_finite(
            substance.melting_point_kelvin,
            "alloy constituent melting point",
        )?;
        let mass = amount * substance.molar_mass_grams;
        total_mass_grams_per_bucket += mass;
        volume_buckets += mass / substance.liquid_density_grams_per_bucket;
        molar_heat_capacity +=
            amount / total_mol_per_bucket * substance.molar_heat_capacity_j_per_mol_kelvin;
        estimated_solidus_kelvin = estimated_solidus_kelvin.min(substance.melting_point_kelvin);
        estimated_liquidus_kelvin = estimated_liquidus_kelvin.max(substance.melting_point_kelvin);
    }
    validate_positive_finite(total_mass_grams_per_bucket, "alloy total mass")?;
    validate_positive_finite(volume_buckets, "alloy volume")?;
    validate_positive_finite(molar_heat_capacity, "alloy heat capacity")?;
    validate_non_negative_finite(estimated_solidus_kelvin, "alloy solidus estimate")?;
    validate_non_negative_finite(estimated_liquidus_kelvin, "alloy liquidus estimate")?;

    let constituents = raw_constituents
        .into_iter()
        .map(|(substance, amount)| {
            let mass = amount * substance.molar_mass_grams;
            AlloyConstituent {
                substance_id: substance.id.clone(),
                metallurgical_component_id: alloy_component_id(substance),
                role: alloy_constituent_role(substance),
                concentration_mol_per_bucket: amount,
                mole_fraction: amount / total_mol_per_bucket,
                mass_fraction: mass / total_mass_grams_per_bucket,
            }
        })
        .collect::<Vec<_>>();

    Ok(AlloyPhaseSnapshot {
        phase_id,
        representative_substance_id,
        temperature_kelvin,
        total_mol_per_bucket,
        total_mass_grams_per_bucket,
        density_grams_per_bucket: total_mass_grams_per_bucket / volume_buckets,
        average_molar_mass_grams: total_mass_grams_per_bucket / total_mol_per_bucket,
        molar_heat_capacity_j_per_mol_kelvin: molar_heat_capacity,
        estimated_solidus_kelvin,
        estimated_liquidus_kelvin,
        constituents,
    })
}

fn alloy_component_id(substance: &Substance) -> String {
    match &substance.representation {
        SubstanceRepresentation::Metal { element_symbol } => element_symbol.clone(),
        SubstanceRepresentation::MetallurgicalSolute { component_id } => component_id.clone(),
        _ => substance.id.to_string(),
    }
}

fn alloy_constituent_role(substance: &Substance) -> AlloyConstituentRole {
    match &substance.representation {
        SubstanceRepresentation::Metal { .. } => AlloyConstituentRole::Metal,
        SubstanceRepresentation::MetallurgicalSolute { .. } => {
            AlloyConstituentRole::DissolvedMaterial
        }
        _ => AlloyConstituentRole::DissolvedMaterial,
    }
}

fn validate_positive_finite(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be positive and finite"
        )));
    }
    Ok(())
}

fn validate_non_negative_finite(value: f64, label: &str) -> ChemistryResult<()> {
    if !value.is_finite() || value < 0.0 {
        return Err(ChemistryError::InvalidMixtureState(format!(
            "{label} must be non-negative and finite"
        )));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chemistry::mixture::Mixture;
    use crate::chemistry::registry::ChemistryRegistryBuilder;
    use crate::chemistry::substance::{
        LiquidPhasePreference, SolventRole, SubstancePhaseProperties,
    };

    #[test]
    fn molten_metals_form_one_alloy_phase_without_named_alloy_data() {
        let registry = alloy_registry().build().unwrap();
        let mut mixture = Mixture::new(2000.0).unwrap();

        mixture
            .add_substance(&registry, "destroy:test_iron", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:test_copper", 3.0)
            .unwrap();

        assert_eq!(mixture.liquid_phase_count(&registry).unwrap(), 1);
        let alloys = alloy_phase_snapshots(&registry, &mixture).unwrap();
        assert_eq!(alloys.len(), 1);
        assert_eq!(alloys[0].constituents.len(), 2);
        assert!((alloys[0].total_mol_per_bucket - 4.0).abs() < 1.0e-12);
        assert!(alloys[0]
            .constituents
            .iter()
            .any(
                |component| component.substance_id == SubstanceId::from("destroy:test_iron")
                    && (component.mole_fraction - 0.25).abs() < 1.0e-12
            ));
        assert!(alloys[0]
            .constituents
            .iter()
            .any(
                |component| component.substance_id == SubstanceId::from("destroy:test_copper")
                    && (component.mole_fraction - 0.75).abs() < 1.0e-12
            ));
        assert!(alloys[0].density_grams_per_bucket > 7_800.0);
        assert!(alloys[0].density_grams_per_bucket < 8_960.0);
    }

    #[test]
    fn molten_metal_and_slag_remain_separate_phases() {
        let registry = alloy_registry().build().unwrap();
        let mut mixture = Mixture::new(2000.0).unwrap();

        mixture
            .add_substance(&registry, "destroy:test_iron", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:test_slag", 1.0)
            .unwrap();

        assert_eq!(mixture.liquid_phase_count(&registry).unwrap(), 2);
        let alloys = alloy_phase_snapshots(&registry, &mixture).unwrap();
        assert_eq!(alloys.len(), 1);
        assert_eq!(alloys[0].constituents.len(), 1);
        assert_eq!(
            alloys[0].constituents[0].substance_id,
            SubstanceId::from("destroy:test_iron")
        );
    }

    #[test]
    fn solid_metals_do_not_create_alloy_phase_below_melting_point() {
        let registry = alloy_registry().build().unwrap();
        let mut mixture = Mixture::new(298.0).unwrap();

        mixture
            .add_substance(&registry, "destroy:test_iron", 1.0)
            .unwrap();
        mixture
            .add_substance(&registry, "destroy:test_copper", 1.0)
            .unwrap();

        assert_eq!(mixture.liquid_phase_count(&registry).unwrap(), 0);
        assert_eq!(
            mixture.concentration_in_phase(
                &SubstanceId::from("destroy:test_iron"),
                MixturePhase::Solid
            ),
            1.0
        );
        assert!(alloy_phase_snapshots(&registry, &mixture)
            .unwrap()
            .is_empty());
    }

    fn alloy_registry() -> ChemistryRegistryBuilder {
        ChemistryRegistryBuilder::new()
            .substance(test_metal(
                "destroy:test_iron",
                "Fe",
                55.845,
                7_874.0,
                1811.0,
            ))
            .substance(test_metal(
                "destroy:test_copper",
                "Cu",
                63.546,
                8_960.0,
                1357.77,
            ))
            .substance(
                Substance::new("destroy:test_slag", 0, 60.084, 2_650.0, 3200.0, 45.0, 0.0)
                    .with_solid_density_grams_per_bucket(2_650.0)
                    .with_melting_point_kelvin(1200.0)
                    .with_phase_properties(molten_phase_properties(
                        LiquidPhasePreference::MoltenSlag,
                    )),
            )
    }

    fn test_metal(
        id: &'static str,
        element: &'static str,
        molar_mass: f64,
        density: f64,
        melting_point: f64,
    ) -> Substance {
        Substance::new(id, 0, molar_mass, density, 3200.0, 25.0, 0.0)
            .with_solid_density_grams_per_bucket(density)
            .with_melting_point_kelvin(melting_point)
            .with_phase_properties(molten_phase_properties(LiquidPhasePreference::MoltenMetal))
            .with_representation(SubstanceRepresentation::Metal {
                element_symbol: element.to_string(),
            })
    }

    fn molten_phase_properties(
        preferred_liquid_phase: LiquidPhasePreference,
    ) -> SubstancePhaseProperties {
        SubstancePhaseProperties {
            preferred_liquid_phase,
            aqueous_solubility_mol_per_bucket: Some(0.0),
            organic_solubility_mol_per_bucket: Some(0.0),
            can_precipitate: true,
            can_form_liquid_phase: true,
            solvent_role: SolventRole::NotSolvent,
        }
    }
}
