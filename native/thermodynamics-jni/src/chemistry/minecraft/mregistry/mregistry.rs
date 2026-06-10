use super::item_to_substance::{ItemSubstancePair, ItemToSubstanceMappingRegistry, MinecraftId};
use super::substance_to_item::{SubstanceItemPair, SubstanceToItemMappingRegistry};
use crate::chemistry::substance::SubstanceId;

#[derive(Debug, Clone, Default)]
pub struct MinecraftChemicalRegistry {
    item_to_substance: ItemToSubstanceMappingRegistry,
    substance_to_item: SubstanceToItemMappingRegistry,
}

impl MinecraftChemicalRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn register(
        &mut self,
        item_id: MinecraftId,
        substance_id: SubstanceId,
        mol_per_item: f64,
    ) {
        self.item_to_substance.register(
            item_id.clone(),
            ItemSubstancePair::new(substance_id.clone(), mol_per_item),
        );
        self.substance_to_item.register(
            substance_id,
            SubstanceItemPair::new(item_id, mol_per_item),
        );
    }

    pub fn lookup_by_item(&self, item_id: &str) -> Option<&ItemSubstancePair> {
        self.item_to_substance.lookup(item_id)
    }

    pub fn lookup_by_substance(
        &self,
        substance_id: &SubstanceId,
    ) -> Option<&[SubstanceItemPair]> {
        self.substance_to_item.lookup(substance_id)
    }

    pub fn contains_item(&self, item_id: &str) -> bool {
        self.item_to_substance.contains(item_id)
    }

    pub fn contains_substance(&self, substance_id: &SubstanceId) -> bool {
        self.substance_to_item.contains(substance_id)
    }

    pub fn item_count(&self) -> usize {
        self.item_to_substance.len()
    }

    pub fn substance_count(&self) -> usize {
        self.substance_to_item.len()
    }

    pub fn is_empty(&self) -> bool {
        self.item_to_substance.is_empty()
    }

    pub fn items(&self) -> impl Iterator<Item = (&MinecraftId, &ItemSubstancePair)> {
        self.item_to_substance.items()
    }

    pub fn substances(&self) -> impl Iterator<Item = (&SubstanceId, &[SubstanceItemPair])> {
        self.substance_to_item.items()
    }
}
