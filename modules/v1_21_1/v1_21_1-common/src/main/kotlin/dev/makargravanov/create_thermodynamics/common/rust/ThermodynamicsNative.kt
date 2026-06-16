package dev.makargravanov.create_thermodynamics.common.rust

object ThermodynamicsNative {
    init {
        NativeLibraryLoader.load()
    }

    data class ItemChemicalBinding(
        val itemId: String,
        val substanceId: String,
        val molPerItem: Double,
    )

    fun idealGasPressure(moles: Double, temperatureKelvin: Double, volumeCubicMeters: Double): Double =
        nativeIdealGasPressure(moles, temperatureKelvin, volumeCubicMeters)

    fun abiVersion(): Int = nativeAbiVersion()

    fun configureItemChemicalBindings(bindings: Collection<ItemChemicalBinding>) {
        val itemIds = Array(bindings.size) { "" }
        val substanceIds = Array(bindings.size) { "" }
        val molPerItems = DoubleArray(bindings.size)
        bindings.forEachIndexed { index, binding ->
            itemIds[index] = binding.itemId
            substanceIds[index] = binding.substanceId
            molPerItems[index] = binding.molPerItem
        }
        nativeReplaceMinecraftItemChemicalBindings(itemIds, substanceIds, molPerItems)
    }

    fun clearItemChemicalBindings() {
        nativeClearMinecraftItemChemicalBindings()
    }

    fun itemChemicalBindingCount(): Int = nativeMinecraftItemChemicalBindingCount()

    fun hasItemChemicalBinding(itemId: String): Boolean =
        nativeHasMinecraftItemChemicalBinding(itemId)

    @JvmStatic
    private external fun nativeIdealGasPressure(
        moles: Double,
        temperatureKelvin: Double,
        volumeCubicMeters: Double,
    ): Double

    @JvmStatic
    private external fun nativeAbiVersion(): Int

    @JvmStatic
    private external fun nativeReplaceMinecraftItemChemicalBindings(
        itemIds: Array<String>,
        substanceIds: Array<String>,
        molPerItems: DoubleArray,
    )

    @JvmStatic
    private external fun nativeClearMinecraftItemChemicalBindings()

    @JvmStatic
    private external fun nativeMinecraftItemChemicalBindingCount(): Int

    @JvmStatic
    private external fun nativeHasMinecraftItemChemicalBinding(itemId: String): Boolean
}
