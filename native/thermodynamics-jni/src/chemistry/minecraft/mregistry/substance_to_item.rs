use std::collections::HashMap;

use crate::chemistry::substance::SubstanceId;

use super::item_to_substance::MinecraftId;

#[derive(Debug, Clone)]
pub struct SubstanceItemPair {
    pub item_id: MinecraftId,
    pub mol_per_item: f64,
}

impl SubstanceItemPair {
    pub fn new(item_id: MinecraftId, mol_per_item: f64) -> Self {
        Self {
            item_id,
            mol_per_item,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct SubstanceToItemMappingRegistry {
    entries: HashMap<SubstanceId, Vec<SubstanceItemPair>>,
}

impl SubstanceToItemMappingRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register(&mut self, substance_id: SubstanceId, pair: SubstanceItemPair) {
        self.entries.entry(substance_id).or_default().push(pair);
    }

    pub fn lookup(&self, substance_id: &SubstanceId) -> Option<&[SubstanceItemPair]> {
        self.entries.get(substance_id).map(|v| v.as_slice())
    }

    pub fn contains(&self, substance_id: &SubstanceId) -> bool {
        self.entries.contains_key(substance_id)
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn items(&self) -> impl Iterator<Item = (&SubstanceId, &[SubstanceItemPair])> {
        self.entries.iter().map(|(k, v)| (k, v.as_slice()))
    }

    pub fn merge(&mut self, other: SubstanceToItemMappingRegistry) {
        for (substance_id, pairs) in other.entries {
            self.entries.entry(substance_id).or_default().extend(pairs);
        }
    }
}
