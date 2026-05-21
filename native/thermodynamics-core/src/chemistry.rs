use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ElementId(pub u16);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpeciesId(pub u16);

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PhaseId(pub u16);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseKind {
    Aqueous,
    Solid,
    Gas,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Element {
    pub id: ElementId,
    pub atomic_number: u8,
    pub symbol: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StandardThermo {
    pub standard_gibbs_energy_joule_per_mol_298_15: f64,
    pub valid_min_temperature_kelvin: f64,
    pub valid_max_temperature_kelvin: f64,
    pub provenance: &'static str,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Species {
    pub id: SpeciesId,
    pub symbol: &'static str,
    pub composition: BTreeMap<ElementId, u16>,
    pub charge_number: i8,
    pub phase: PhaseKind,
    pub thermo: StandardThermo,
}

impl Species {
    pub fn element_count(&self, element_id: ElementId) -> f64 {
        self.composition
            .get(&element_id)
            .copied()
            .unwrap_or_default() as f64
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpeciesAmount {
    pub species_id: SpeciesId,
    pub amount_mol: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_ids_do_not_compare_across_domains() {
        let element = ElementId(1);
        let species = SpeciesId(1);
        let phase = PhaseId(1);

        assert_eq!(element, ElementId(1));
        assert_eq!(species, SpeciesId(1));
        assert_eq!(phase, PhaseId(1));
    }
}
