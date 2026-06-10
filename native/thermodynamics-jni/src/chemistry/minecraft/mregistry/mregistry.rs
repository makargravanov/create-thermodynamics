use super::item_to_substance::{ItemSubstancePair, ItemToSubstanceMappingRegistry, MinecraftId};
use super::substance_to_item::{SubstanceItemPair, SubstanceToItemMappingRegistry};
use crate::chemistry::substance::SubstanceId;
#[derive(Debug, Clone)]
pub struct DuplicateItemError {
    pub item_id: MinecraftId,
    pub existing_substance_id: SubstanceId,
    pub new_substance_id: SubstanceId,
}

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
    ) -> Result<(), DuplicateItemError> {
        if let Some(existing) = self.item_to_substance.lookup(item_id.as_str()) {
            return Err(DuplicateItemError {
                item_id,
                existing_substance_id: existing.substance_id.clone(),
                new_substance_id: substance_id,
            });
        }
        self.item_to_substance.register(
            item_id.clone(),
            ItemSubstancePair::new(substance_id.clone(), mol_per_item),
        );
        self.substance_to_item
            .register(substance_id, SubstanceItemPair::new(item_id, mol_per_item));
        Ok(())
    }

    pub fn lookup_by_item(&self, item_id: &str) -> Option<&ItemSubstancePair> {
        self.item_to_substance.lookup(item_id)
    }

    pub fn lookup_by_substance(&self, substance_id: &SubstanceId) -> Option<&[SubstanceItemPair]> {
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

#[cfg(test)]
mod tests {
    use crate::chemistry::minecraft::mregistry::{
        item_to_substance::MinecraftId, mregistry::MinecraftChemicalRegistry,
    };
    use crate::chemistry::substance::SubstanceId;

    fn fe_id() -> SubstanceId {
        SubstanceId::from("destroy:iron_iii")
    }

    fn fe_block_id() -> SubstanceId {
        SubstanceId::from("destroy:iron_block")
    }

    fn cu_id() -> SubstanceId {
        SubstanceId::from("destroy:copper_ii")
    }

    #[test]
    fn same_substance_can_map_to_multiple_items() {
        let mut registry = MinecraftChemicalRegistry::new();
        registry
            .register(MinecraftId::from("minecraft:iron_ore"), fe_id(), 10.0)
            .unwrap();
        registry
            .register(MinecraftId::from("minecraft:iron_block"), fe_id(), 90.0)
            .unwrap();

        assert_eq!(registry.item_count(), 2);
        assert_eq!(registry.substance_count(), 1);
        let items = registry.lookup_by_substance(&fe_id()).unwrap();
        assert_eq!(items.len(), 2);
    }

    #[test]
    fn duplicate_item_id_is_rejected() {
        let mut registry = MinecraftChemicalRegistry::new();
        registry
            .register(MinecraftId::from("minecraft:iron_ore"), fe_id(), 10.0)
            .unwrap();
        let err = registry
            .register(MinecraftId::from("minecraft:iron_ore"), cu_id(), 5.0)
            .unwrap_err();
        assert_eq!(err.existing_substance_id, fe_id());
        assert_eq!(err.new_substance_id, cu_id());
    }

    #[test]
    fn lookup_by_item_returns_correct_substance() {
        let mut registry = MinecraftChemicalRegistry::new();
        registry
            .register(MinecraftId::from("minecraft:iron_ore"), fe_id(), 10.0)
            .unwrap();
        registry
            .register(MinecraftId::from("minecraft:copper_ore"), cu_id(), 8.0)
            .unwrap();

        let iron = registry.lookup_by_item("minecraft:iron_ore").unwrap();
        assert_eq!(iron.substance_id, fe_id());
        assert_eq!(iron.mol_per_item, 10.0);

        let copper = registry.lookup_by_item("minecraft:copper_ore").unwrap();
        assert_eq!(copper.substance_id, cu_id());
        assert_eq!(copper.mol_per_item, 8.0);
    }

    #[test]
    fn lookup_by_substance_returns_all_items() {
        let mut registry = MinecraftChemicalRegistry::new();
        registry
            .register(MinecraftId::from("minecraft:iron_ore"), fe_id(), 10.0)
            .unwrap();
        registry
            .register(MinecraftId::from("minecraft:iron_ingot"), fe_id(), 1.0)
            .unwrap();
        registry
            .register(MinecraftId::from("minecraft:iron_block"), fe_id(), 90.0)
            .unwrap();

        let items = registry.lookup_by_substance(&fe_id()).unwrap();
        assert_eq!(items.len(), 3);
        assert_eq!(items[0].item_id.as_str(), "minecraft:iron_ore");
        assert_eq!(items[0].mol_per_item, 10.0);
        assert_eq!(items[1].item_id.as_str(), "minecraft:iron_ingot");
        assert_eq!(items[1].mol_per_item, 1.0);
        assert_eq!(items[2].item_id.as_str(), "minecraft:iron_block");
        assert_eq!(items[2].mol_per_item, 90.0);
    }

    #[test]
    fn contains_checks_both_directions() {
        let mut registry = MinecraftChemicalRegistry::new();
        registry
            .register(
                MinecraftId::from("minecraft:iron_ore"),
                SubstanceId::from("destroy:iron_iii"),
                10.0,
            )
            .unwrap();

        assert!(registry.contains_item("minecraft:iron_ore"));
        assert!(registry.contains_substance(&SubstanceId::from("destroy:iron_iii")));
        assert!(!registry.contains_item("minecraft:copper_ore"));
        assert!(!registry.contains_substance(&SubstanceId::from("destroy:copper_ii")));
    }

    #[test]
    fn items_and_substances_iterators() {
        let mut registry = MinecraftChemicalRegistry::new();
        registry
            .register(MinecraftId::from("minecraft:iron_ore"), fe_id(), 10.0)
            .unwrap();
        registry
            .register(MinecraftId::from("minecraft:copper_ore"), cu_id(), 8.0)
            .unwrap();

        let item_ids: Vec<&str> = registry.items().map(|(id, _)| id.as_str()).collect();
        assert!(item_ids.contains(&"minecraft:iron_ore"));
        assert!(item_ids.contains(&"minecraft:copper_ore"));

        let substance_ids: Vec<&SubstanceId> = registry
            .substances()
            .map(|(id, _)| id)
            .collect();
        assert!(substance_ids.contains(&&SubstanceId::from("destroy:iron_iii")));
        assert!(substance_ids.contains(&&SubstanceId::from("destroy:copper_ii")));
    }

    #[test]
    fn empty_registry() {
        let registry = MinecraftChemicalRegistry::new();
        assert!(registry.is_empty());
        assert_eq!(registry.item_count(), 0);
        assert_eq!(registry.substance_count(), 0);
        assert!(registry.lookup_by_item("anything").is_none());
        assert!(registry
            .lookup_by_substance(&SubstanceId::from("destroy:nonexistent"))
            .is_none());
    }

    #[test]
    fn duplicate_error_does_not_corrupt_registry() {
        let mut registry = MinecraftChemicalRegistry::new();
        registry
            .register(MinecraftId::from("minecraft:iron_ore"), fe_id(), 10.0)
            .unwrap();
        let _err = registry
            .register(MinecraftId::from("minecraft:iron_ore"), cu_id(), 5.0)
            .unwrap_err();

        assert!(registry.contains_item("minecraft:iron_ore"));
        let pair = registry.lookup_by_item("minecraft:iron_ore").unwrap();
        assert_eq!(pair.substance_id, fe_id());
        assert_eq!(pair.mol_per_item, 10.0);
        assert_eq!(registry.item_count(), 1);
        assert_eq!(registry.substance_count(), 1);
    }
}
