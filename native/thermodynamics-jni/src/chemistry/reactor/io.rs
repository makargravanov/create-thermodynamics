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
        zone.mixture_mut()
            .change_concentration(registry, &component.id, -take)?;
        extracted.push((component.id.clone(), take));
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

pub fn transfer_all(
    from: &mut ReactorZone,
    to: &mut ReactorZone,
    registry: &ChemistryRegistry,
) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
    let extracted = extract_all(from, registry)?;
    for (id, amount) in &extracted {
        to.mixture_mut()
            .add_substance(registry, id.clone(), *amount)?;
    }
    Ok(extracted)
}

pub fn insert_from_phase(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    phase: MixturePhase,
    max_mol_per_bucket: f64,
) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
    let mut inserted = Vec::new();
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

        let amount = mol_in_phase * scale;
        zone.mixture_mut()
            .add_substance(registry, component.id.clone(), amount)?;
        inserted.push((component.id.clone(), amount));
    }

    Ok(inserted)
}

pub fn insert_all(
    zone: &mut ReactorZone,
    registry: &ChemistryRegistry,
    max_mol_per_bucket: f64,
) -> ChemistryResult<Vec<(SubstanceId, f64)>> {
    let snapshot = mixture_snapshot(zone);
    let mut inserted = Vec::new();

    let total: f64 = snapshot
        .substances
        .iter()
        .map(|s| s.total_mol_per_bucket)
        .sum();
    if total <= 0.0 || max_mol_per_bucket <= 0.0 {
        return Ok(Vec::new());
    }

    let scale = (max_mol_per_bucket / total).min(1.0);

    for component in &snapshot.substances {
        if component.total_mol_per_bucket <= 0.0 {
            continue;
        }
        let amount = component.total_mol_per_bucket * scale;
        zone.mixture_mut()
            .add_substance(registry, component.id.clone(), amount)?;
        inserted.push((component.id.clone(), amount));
    }

    Ok(inserted)
}
