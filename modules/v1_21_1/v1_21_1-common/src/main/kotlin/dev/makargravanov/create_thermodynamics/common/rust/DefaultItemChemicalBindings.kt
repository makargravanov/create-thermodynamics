package dev.makargravanov.create_thermodynamics.common.rust

object DefaultItemChemicalBindings {
    val bindings: List<ThermodynamicsNative.ItemChemicalBinding> = listOf(
        binding("minecraft:water_bucket", "destroy:water", 1.0),
        binding("minecraft:ice", "destroy:water", 1.0),
        binding("minecraft:snowball", "destroy:water", 0.25),
        binding("minecraft:coal", "destroy:carbon", 1.0),
        binding("minecraft:charcoal", "destroy:carbon", 1.0),
        binding("minecraft:iron_ingot", "destroy:iron", 1.0),
        binding("minecraft:copper_ingot", "destroy:copper", 1.0),
    )

    private fun binding(
        itemId: String,
        substanceId: String,
        molPerItem: Double,
    ): ThermodynamicsNative.ItemChemicalBinding =
        ThermodynamicsNative.ItemChemicalBinding(
            itemId = itemId,
            substanceId = substanceId,
            molPerItem = molPerItem,
        )
}
