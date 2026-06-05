use std::collections::BTreeSet;

use crate::chemistry::error::{ChemistryError, ChemistryResult};
use crate::chemistry::metallurgy::{
    MetallurgicalComponentId, MetallurgicalCompoundPhaseData, MetallurgicalElementData,
    MetallurgicalPairInteractionData,
};

pub(crate) fn validate_element_data(
    element_data: Vec<MetallurgicalElementData>,
) -> ChemistryResult<Vec<MetallurgicalElementData>> {
    let mut seen = BTreeSet::new();
    for data in &element_data {
        data.validate()?;
        if !seen.insert(data.component.clone()) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "duplicate metallurgical element data '{}'",
                data.component.as_str()
            )));
        }
    }
    Ok(element_data)
}

pub(crate) fn validate_pair_interactions(
    pair_interactions: Vec<MetallurgicalPairInteractionData>,
) -> ChemistryResult<Vec<MetallurgicalPairInteractionData>> {
    let mut seen = BTreeSet::new();
    for interaction in &pair_interactions {
        interaction.validate()?;
        let key = ordered_pair_key(&interaction.first, &interaction.second);
        if !seen.insert(key.clone()) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "duplicate metallurgical pair interaction '{}:{}'",
                key.0.as_str(),
                key.1.as_str()
            )));
        }
    }
    Ok(pair_interactions)
}

pub(crate) fn validate_compound_phases(
    compound_phases: Vec<MetallurgicalCompoundPhaseData>,
) -> ChemistryResult<Vec<MetallurgicalCompoundPhaseData>> {
    let mut seen = BTreeSet::new();
    for phase in &compound_phases {
        phase.validate()?;
        if !seen.insert(phase.id.clone()) {
            return Err(ChemistryError::InvalidMixtureState(format!(
                "duplicate metallurgical compound phase '{}'",
                phase.id
            )));
        }
    }
    Ok(compound_phases)
}

fn ordered_pair_key(
    left: &MetallurgicalComponentId,
    right: &MetallurgicalComponentId,
) -> (MetallurgicalComponentId, MetallurgicalComponentId) {
    if left <= right {
        (left.clone(), right.clone())
    } else {
        (right.clone(), left.clone())
    }
}
