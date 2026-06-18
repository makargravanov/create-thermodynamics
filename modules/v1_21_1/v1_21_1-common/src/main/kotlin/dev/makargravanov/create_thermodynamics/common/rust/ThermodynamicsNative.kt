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

    @JvmInline
    value class NativeReactorId(val value: Long)

    fun idealGasPressure(moles: Double, temperatureKelvin: Double, volumeCubicMeters: Double): Double =
        nativeIdealGasPressure(moles, temperatureKelvin, volumeCubicMeters)

    fun abiVersion(): Int = nativeAbiVersion()

    fun configureDefaultItemChemicalBindings() {
        configureItemChemicalBindings(DefaultItemChemicalBindings.bindings)
    }

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

    fun staticSubstanceIds(): List<String> =
        nativeStaticSubstanceIds().asList()

    fun createSingleZoneReactor(
        volumeCubicMeters: Double,
        itemInputCount: Int,
        itemOutputCount: Int,
        fluidInputCount: Int,
        fluidOutputCount: Int,
    ): NativeReactorId =
        NativeReactorId(
            nativeCreateSingleZoneReactor(
                volumeCubicMeters,
                itemInputCount,
                itemOutputCount,
                fluidInputCount,
                fluidOutputCount,
            ),
        )

    fun removeReactor(reactorId: NativeReactorId) {
        nativeRemoveReactor(reactorId.value)
    }

    fun reactorCount(): Int = nativeReactorCount()

    fun tickReactor(reactorId: NativeReactorId, dtSeconds: Double) {
        nativeTickReactor(reactorId.value, dtSeconds)
    }

    fun insertItemStackToReactorInput(
        reactorId: NativeReactorId,
        inputIndex: Int,
        itemId: String,
        itemCount: Int,
    ): Double =
        nativeInsertItemStackToReactorInput(reactorId.value, inputIndex, itemId, itemCount)

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

    @JvmStatic
    private external fun nativeStaticSubstanceIds(): Array<String>

    @JvmStatic
    private external fun nativeCreateSingleZoneReactor(
        volumeCubicMeters: Double,
        itemInputCount: Int,
        itemOutputCount: Int,
        fluidInputCount: Int,
        fluidOutputCount: Int,
    ): Long

    @JvmStatic
    private external fun nativeRemoveReactor(reactorId: Long)

    @JvmStatic
    private external fun nativeReactorCount(): Int

    @JvmStatic
    private external fun nativeTickReactor(reactorId: Long, dtSeconds: Double)

    @JvmStatic
    private external fun nativeInsertItemStackToReactorInput(
        reactorId: Long,
        inputIndex: Int,
        itemId: String,
        itemCount: Int,
    ): Double
}
