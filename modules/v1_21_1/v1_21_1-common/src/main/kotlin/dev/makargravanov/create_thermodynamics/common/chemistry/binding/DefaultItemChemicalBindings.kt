package dev.makargravanov.create_thermodynamics.common.chemistry.binding

object DefaultItemChemicalBindings {
    val bindings: List<ItemChemicalBinding> = listOf(
        binding("minecraft:water_bucket", "destroy:water", 1.0),
        binding("minecraft:ice", "destroy:water", 1.0),
        binding("minecraft:snowball", "destroy:water", 0.25),
        binding("minecraft:iron_ingot", "destroy:iron_metal", 1.0),
        binding("minecraft:copper_ingot", "destroy:copper_metal", 1.0),
    )

    private fun binding(
        itemId: String,
        substanceId: String,
        molPerItem: Double,
    ): ItemChemicalBinding =
        ItemChemicalBinding(
            itemId = itemId,
            substanceId = substanceId,
            molPerItem = molPerItem,
        )
}
