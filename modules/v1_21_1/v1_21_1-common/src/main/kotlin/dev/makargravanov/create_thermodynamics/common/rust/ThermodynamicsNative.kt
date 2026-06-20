package dev.makargravanov.create_thermodynamics.common.rust

import dev.makargravanov.create_thermodynamics.common.chemistry.binding.ItemChemicalBinding

object ThermodynamicsNative {
    init {
        NativeLibraryLoader.load()
    }

    @JvmInline
    value class NativeReactorId(val value: Long)

    data class NativeReactorZoneSubstanceSnapshot(
        val substanceId: String,
        val concentrationMolPerBucket: Double,
    ) {
        init {
            require(substanceId.isNotBlank()) { "substanceId must not be blank" }
            require(concentrationMolPerBucket.isFinite() && concentrationMolPerBucket >= 0.0) {
                "concentrationMolPerBucket must be non-negative and finite"
            }
        }
    }

    data class NativeReactorZoneSnapshot(
        val temperatureKelvin: Double,
        val pressurePascal: Double,
        val substances: List<NativeReactorZoneSubstanceSnapshot>,
    ) {
        init {
            require(temperatureKelvin.isFinite() && temperatureKelvin >= 0.0) {
                "temperatureKelvin must be non-negative and finite"
            }
            require(pressurePascal.isFinite() && pressurePascal >= 0.0) {
                "pressurePascal must be non-negative and finite"
            }
        }
    }

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

    fun reactorZoneSnapshot(
        reactorId: NativeReactorId,
        zoneIndex: Int,
    ): NativeReactorZoneSnapshot {
        require(zoneIndex >= 0) { "zoneIndex must be non-negative" }
        val substanceIds = nativeReactorZoneSubstanceIds(reactorId.value, zoneIndex)
        val concentrations = nativeReactorZoneConcentrationsMolPerBucket(reactorId.value, zoneIndex)
        require(substanceIds.size == concentrations.size) {
            "native reactor zone snapshot returned ${substanceIds.size} substance ids and ${concentrations.size} concentrations"
        }
        return NativeReactorZoneSnapshot(
            temperatureKelvin = nativeReactorZoneTemperatureKelvin(reactorId.value, zoneIndex),
            pressurePascal = nativeReactorZonePressurePascal(reactorId.value, zoneIndex),
            substances = substanceIds.indices.map { index ->
                NativeReactorZoneSubstanceSnapshot(
                    substanceId = substanceIds[index],
                    concentrationMolPerBucket = concentrations[index],
                )
            },
        )
    }

    fun insertItemStackToReactorInput(
        reactorId: NativeReactorId,
        inputIndex: Int,
        itemId: String,
        itemCount: Int,
    ): Int =
        nativeInsertItemStackToReactorInput(reactorId.value, inputIndex, itemId, itemCount)

    fun extractItemStackFromReactorOutput(
        reactorId: NativeReactorId,
        outputIndex: Int,
        itemId: String,
        maxItemCount: Int,
        dtSeconds: Double,
    ): Int =
        nativeExtractItemStackFromReactorOutput(
            reactorId.value,
            outputIndex,
            itemId,
            maxItemCount,
            dtSeconds,
        )

    fun exportCatalogCheckpoint(contentVersion: Long): ByteArray =
        nativeExportCatalogCheckpoint(contentVersion)

    fun importCatalogCheckpoint(encoded: ByteArray) {
        nativeImportCatalogCheckpoint(encoded)
    }

    fun exportReactorCheckpoint(
        reactorId: NativeReactorId,
        contentVersion: Long,
    ): ByteArray =
        nativeExportReactorCheckpoint(reactorId.value, contentVersion)

    fun createReactorFromCheckpoint(encoded: ByteArray): NativeReactorId =
        NativeReactorId(nativeCreateReactorFromCheckpoint(encoded))

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
    private external fun nativeReactorZoneTemperatureKelvin(reactorId: Long, zoneIndex: Int): Double

    @JvmStatic
    private external fun nativeReactorZonePressurePascal(reactorId: Long, zoneIndex: Int): Double

    @JvmStatic
    private external fun nativeReactorZoneSubstanceIds(reactorId: Long, zoneIndex: Int): Array<String>

    @JvmStatic
    private external fun nativeReactorZoneConcentrationsMolPerBucket(reactorId: Long, zoneIndex: Int): DoubleArray

    @JvmStatic
    private external fun nativeInsertItemStackToReactorInput(
        reactorId: Long,
        inputIndex: Int,
        itemId: String,
        itemCount: Int,
    ): Int

    @JvmStatic
    private external fun nativeExtractItemStackFromReactorOutput(
        reactorId: Long,
        outputIndex: Int,
        itemId: String,
        maxItemCount: Int,
        dtSeconds: Double,
    ): Int

    @JvmStatic
    private external fun nativeExportCatalogCheckpoint(contentVersion: Long): ByteArray

    @JvmStatic
    private external fun nativeImportCatalogCheckpoint(encoded: ByteArray)

    @JvmStatic
    private external fun nativeExportReactorCheckpoint(
        reactorId: Long,
        contentVersion: Long,
    ): ByteArray

    @JvmStatic
    private external fun nativeCreateReactorFromCheckpoint(encoded: ByteArray): Long
}
