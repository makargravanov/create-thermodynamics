use std::borrow::Borrow;
use std::collections::HashMap;

use crate::chemistry::substance::SubstanceId;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MinecraftId(String);

impl MinecraftId {
    pub fn new(id: String) -> Self {
        Self(id)
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn into_string(self) -> String {
        self.0
    }
}

impl From<&str> for MinecraftId {
    fn from(s: &str) -> Self {
        Self(s.to_string())
    }
}

impl Borrow<str> for MinecraftId {
    fn borrow(&self) -> &str {
        &self.0
    }
}

#[derive(Debug, Clone)]
pub struct ItemSubstancePair {
    pub substance_id: SubstanceId,
    pub mol_per_item: f64,
}

impl ItemSubstancePair {
    pub fn new(substance_id: SubstanceId, mol_per_item: f64) -> Self {
        Self {
            substance_id,
            mol_per_item,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct ItemToSubstanceMappingRegistry {
    entries: HashMap<MinecraftId, ItemSubstancePair>,
}

impl ItemToSubstanceMappingRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register(&mut self, item_id: MinecraftId, mapping: ItemSubstancePair) {
        self.entries.insert(item_id, mapping);
    }

    pub fn lookup(&self, item_id: &str) -> Option<&ItemSubstancePair> {
        self.entries.get(item_id)
    }

    pub fn contains(&self, item_id: &str) -> bool {
        self.entries.contains_key(item_id)
    }

    pub fn len(&self) -> usize {
        self.entries.len()
    }

    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    pub fn items(&self) -> impl Iterator<Item = (&MinecraftId, &ItemSubstancePair)> {
        self.entries.iter()
    }

    pub fn merge(&mut self, other: ItemToSubstanceMappingRegistry) {
        self.entries.extend(other.entries);
    }
}
