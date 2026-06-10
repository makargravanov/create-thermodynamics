use std::collections::HashMap;
use std::borrow::Borrow;

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
pub struct SubstanceMapping {
    pub substance_id: SubstanceId,
    pub mol_per_item: f64,
}

impl SubstanceMapping {
    pub fn new(substance_id: SubstanceId, mol_per_item: f64) -> Self {
        Self {
            substance_id,
            mol_per_item,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct MappingRegistry {
    entries: HashMap<MinecraftId, Vec<SubstanceMapping>>,
}

impl MappingRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register(&mut self, item_id: MinecraftId, mapping: SubstanceMapping) {
        self.entries.entry(item_id).or_default().push(mapping);
    }

    pub fn register_many(&mut self, item_id: MinecraftId, mappings: Vec<SubstanceMapping>) {
        if !mappings.is_empty() {
            self.entries.entry(item_id).or_default().extend(mappings);
        }
    }

    pub fn lookup(&self, item_id: &str) -> Option<&[SubstanceMapping]> {
        self.entries.get(item_id).map(|v| v.as_slice())
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

    pub fn items(&self) -> impl Iterator<Item = (&MinecraftId, &[SubstanceMapping])> {
        self.entries.iter().map(|(k, v)| (k, v.as_slice()))
    }

    pub fn substance_ids(&self) -> impl Iterator<Item = &SubstanceId> {
        self.entries
            .values()
            .flat_map(|mappings| mappings.iter())
            .map(|m| &m.substance_id)
    }

    pub fn merge(&mut self, other: MappingRegistry) {
        for (item_id, mappings) in other.entries {
            self.entries.entry(item_id).or_default().extend(mappings);
        }
    }
}
