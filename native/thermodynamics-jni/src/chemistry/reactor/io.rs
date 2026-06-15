use crate::chemistry::error::ChemistryResult;
use crate::chemistry::mixture::MixturePhase;
use crate::chemistry::registry::ChemistryRegistry;
use crate::chemistry::substance::SubstanceId;

use super::zone::ReactorZone;

#[derive(Debug, Clone)]
pub struct MixtureSnapshot {
    pub substances: Vec<SubstanceComponent>,
    pub temperature_kelvin: f64,
    pub pressure_pascal: f64,
    pub volume_cubic_meters: f64,
}

#[derive(Debug, Clone)]
pub struct SubstanceComponent {
    pub id: SubstanceId,
    pub total_mol_per_bucket: f64,
    pub phase_amounts: Vec<PhaseAmount>,
}

#[derive(Debug, Clone)]
pub struct PhaseAmount {
    pub phase: MixturePhase,
    pub mol_per_bucket: f64,
}

#[derive(Debug, Clone)]
pub struct ExtractedStream {
    pub id: SubstanceId,
    pub mol_per_bucket: f64,
    pub thermal_energy_j_per_bucket: f64,
}

pub fn mixture_snapshot(zone: &ReactorZone) -> MixtureSnapshot {
    let mixture = zone.mixture();
    let substances: Vec<SubstanceComponent> = mixture
        .substances()
        .map(|id| {
            let total = mixture.concentration_of(id);
            let phase_amounts = MixturePhase::ALL
                .iter()
                .filter_map(|&phase| {
                    let mol = mixture.concentration_in_phase(id, phase);
                    if mol > 0.0 {
                        Some(PhaseAmount {
                            phase,
                            mol_per_bucket: mol,
                        })
                    } else {
                        None
                    }
                })
                .collect();
            SubstanceComponent {
                id: id.clone(),
                total_mol_per_bucket: total,
                phase_amounts,
            }
        })
        .collect();

    MixtureSnapshot {
        substances,
        temperature_kelvin: zone.temperature_kelvin(),
        pressure_pascal: zone.pressure_pascal(),
        volume_cubic_meters: zone.volume_cubic_meters(),
    }
}

pub fn insert_substance(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    substance_id: &SubstanceId,
    mol_per_bucket: f64,
) -> ChemistryResult<()> {
    zone.mixture_mut()
        .add_substance(registry, substance_id.clone(), mol_per_bucket)
}

pub fn extract_substance(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    substance_id: &SubstanceId,
    max_mol_per_bucket: f64,
) -> ChemistryResult<f64> {
    let available = zone.concentration_of(substance_id);
    let amount = available.min(max_mol_per_bucket);
    if amount <= 0.0 {
        return Ok(0.0);
    }
    zone.mixture_mut()
        .change_concentration(registry, substance_id, -amount)?;
    Ok(amount)
}

pub fn extract_substance_stream(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    substance_id: &SubstanceId,
    max_mol_per_bucket: f64,
) -> ChemistryResult<Option<ExtractedStream>> {
    let available = zone.concentration_of(substance_id);
    let amount = available.min(max_mol_per_bucket);
    if amount <= 0.0 {
        return Ok(None);
    }
    let thermal_energy_j_per_bucket =
        stream_energy_for_substance(zone, registry, substance_id, amount)?;
    let actual = extract_substance(zone, registry, substance_id, amount)?;
    if actual <= 0.0 {
        return Ok(None);
    }
    Ok(Some(ExtractedStream {
        id: substance_id.clone(),
        mol_per_bucket: actual,
        thermal_energy_j_per_bucket: thermal_energy_j_per_bucket * actual / amount,
    }))
}

pub fn extract_from_phase(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    phase: MixturePhase,
    max_mol_per_bucket: f64,
) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
    let snapshot = mixture_snapshot(zone);
    let total_in_phase: f64 = snapshot
        .substances
        .iter()
        .filter_map(|s| {
            s.phase_amounts
                .iter()
                .find(|p| p.phase == phase)
                .map(|p| p.mol_per_bucket)
        })
        .sum();

    if total_in_phase <= 0.0 || max_mol_per_bucket <= 0.0 {
        return Ok(Vec::new());
    }

    let scale = (max_mol_per_bucket / total_in_phase).min(1.0);
    let mut extracted = Vec::new();

    for component in &snapshot.substances {
        let mol_in_phase = component
            .phase_amounts
            .iter()
            .find(|p| p.phase == phase)
            .map(|p| p.mol_per_bucket)
            .unwrap_or(0.0);

        if mol_in_phase <= 0.0 {
            continue;
        }

        let take = mol_in_phase * scale;
        let Some(substance_index) = registry.substance_index(&component.id) else {
            continue;
        };
        let actual = zone.mixture_mut().extract_from_phase_by_index(
            registry,
            substance_index,
            phase,
            take,
        )?;
        if actual > 0.0 {
            extracted.push((component.id.clone(), actual));
        }
    }

    Ok(extracted)
}

pub fn extract_from_phase_streams(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    phase: MixturePhase,
    max_mol_per_bucket: f64,
) -> ChemistryResult<Vec<ExtractedStream>> {
    let snapshot = mixture_snapshot(zone);
    let total_in_phase: f64 = snapshot
        .substances
        .iter()
        .filter_map(|s| {
            s.phase_amounts
                .iter()
                .find(|p| p.phase == phase)
                .map(|p| p.mol_per_bucket)
        })
        .sum();

    if total_in_phase <= 0.0 || max_mol_per_bucket <= 0.0 {
        return Ok(Vec::new());
    }

    let scale = (max_mol_per_bucket / total_in_phase).min(1.0);
    let mut extracted = Vec::new();

    for component in &snapshot.substances {
        let mol_in_phase = component
            .phase_amounts
            .iter()
            .find(|p| p.phase == phase)
            .map(|p| p.mol_per_bucket)
            .unwrap_or(0.0);

        if mol_in_phase <= 0.0 {
            continue;
        }

        let take = mol_in_phase * scale;
        let Some(substance_index) = registry.substance_index(&component.id) else {
            continue;
        };
        let thermal_energy_j_per_bucket =
            stream_energy_for_phase(zone, registry, &component.id, phase, take)?;
        let actual = zone.mixture_mut().extract_from_phase_by_index(
            registry,
            substance_index,
            phase,
            take,
        )?;
        if actual > 0.0 {
            extracted.push(ExtractedStream {
                id: component.id.clone(),
                mol_per_bucket: actual,
                thermal_energy_j_per_bucket: thermal_energy_j_per_bucket * actual / take,
            });
        }
    }

    Ok(extracted)
}

pub fn extract_all(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
    let snapshot = mixture_snapshot(zone);
    let mut extracted = Vec::new();

    for component in &snapshot.substances {
        if component.total_mol_per_bucket <= 0.0 {
            continue;
        }
        zone.mixture_mut().change_concentration(
            registry,
            &component.id,
            -component.total_mol_per_bucket,
        )?;
        extracted.push((component.id.clone(), component.total_mol_per_bucket));
    }

    Ok(extracted)
}

pub fn insert_stream(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    stream: &ExtractedStream,
) -> ChemistryResult<()> {
    let target_temperature = zone.temperature_kelvin();
    let target_energy = baseline_energy_for_substance(
        registry,
        &stream.id,
        stream.mol_per_bucket,
        target_temperature,
    )?;
    insert_substance(zone, registry, &stream.id, stream.mol_per_bucket)?;
    let energy_delta = stream.thermal_energy_j_per_bucket - target_energy;
    if energy_delta.abs() > 1.0e-12 {
        zone.mixture_mut().heat(registry, energy_delta)?;
    }
    Ok(())
}

pub fn transfer_all(
    from: &mut ReactorZone,
    to: &mut ReactorZone,
    registry: &ChemistryRegistry,
) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
    let snapshot = mixture_snapshot(from);
    let mut extracted = Vec::new();
    for component in &snapshot.substances {
        if component.total_mol_per_bucket <= 0.0 {
            continue;
        }
        if let Some(stream) = extract_substance_stream(
            from,
            registry,
            &component.id,
            component.total_mol_per_bucket,
        )? {
            insert_stream(to, registry, &stream)?;
            extracted.push((stream.id, stream.mol_per_bucket));
        }
    }
    Ok(extracted)
}

fn stream_energy_for_substance(
    zone: &ReactorZone,
    registry: &ChemistryRegistry,
    substance_id: &SubstanceId,
    amount: f64,
) -> ChemistryResult<f64> {
    let total = zone.concentration_of(substance_id);
    if total <= 0.0 || amount <= 0.0 {
        return Ok(0.0);
    }
    let mut energy = 0.0;
    for phase in MixturePhase::ALL {
        let phase_amount = zone.mixture().concentration_in_phase(substance_id, phase);
        if phase_amount <= 0.0 {
            continue;
        }
        energy += stream_energy_for_phase(
            zone,
            registry,
            substance_id,
            phase,
            amount * phase_amount / total,
        )?;
    }
    Ok(energy)
}

fn stream_energy_for_phase(
    zone: &ReactorZone,
    registry: &ChemistryRegistry,
    substance_id: &SubstanceId,
    phase: MixturePhase,
    amount: f64,
) -> ChemistryResult<f64> {
    let mut energy =
        baseline_energy_for_substance(registry, substance_id, amount, zone.temperature_kelvin())?;
    if matches!(phase, MixturePhase::Gas | MixturePhase::SupercriticalFluid) {
        let substance = registry.substance(substance_id)?;
        energy += substance.latent_heat_j_per_mol * amount;
    }
    Ok(energy)
}

fn baseline_energy_for_substance(
    registry: &ChemistryRegistry,
    substance_id: &SubstanceId,
    amount: f64,
    temperature_kelvin: f64,
) -> ChemistryResult<f64> {
    let substance = registry.substance(substance_id)?;
    Ok(substance.molar_heat_capacity_j_per_mol_kelvin * amount * temperature_kelvin)
}
