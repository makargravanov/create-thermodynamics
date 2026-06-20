package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative

data class NativeReactorPortBinding(
    val port: ReactorPortDescriptor,
    val nativePortIndex: Int,
)

data class NativeReactorMultiblockBinding(
    val structureId: ReactorStructureId,
    val nativeReactorId: ThermodynamicsNative.NativeReactorId,
    val itemInputs: List<NativeReactorPortBinding>,
    val itemOutputs: List<NativeReactorPortBinding>,
    val fluidInputs: List<NativeReactorPortBinding>,
    val fluidOutputs: List<NativeReactorPortBinding>,
) {
    fun allPortBindings(): List<NativeReactorPortBinding> =
        itemInputs + itemOutputs + fluidInputs + fluidOutputs
}

object ReactorMultiblockNativeBridge : NativeReactorBridge {
    override fun createNativeReactor(definition: ReactorMultiblockDefinition): NativeReactorMultiblockBinding {
        val reactorId = ThermodynamicsNative.createSingleZoneReactor(
            volumeCubicMeters = definition.totalVolumeCubicMeters,
            itemInputCount = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).size,
            itemOutputCount = definition.portsOfKind(ReactorPortKind.ITEM_OUTPUT).size,
            fluidInputCount = definition.portsOfKind(ReactorPortKind.FLUID_INPUT).size,
            fluidOutputCount = definition.portsOfKind(ReactorPortKind.FLUID_OUTPUT).size,
        )
        return definition.toNativeBinding(reactorId)
    }

    override fun createNativeReactorFromCheckpoint(
        definition: ReactorMultiblockDefinition,
        encodedCheckpoint: ByteArray,
    ): NativeReactorMultiblockBinding =
        definition.toNativeBinding(ThermodynamicsNative.createReactorFromCheckpoint(encodedCheckpoint))

    override fun exportReactorCheckpoint(
        binding: NativeReactorMultiblockBinding,
        contentVersion: Long,
    ): ByteArray =
        ThermodynamicsNative.exportReactorCheckpoint(binding.nativeReactorId, contentVersion)

    private fun ReactorMultiblockDefinition.toNativeBinding(
        reactorId: ThermodynamicsNative.NativeReactorId,
    ): NativeReactorMultiblockBinding {
        val itemInputs = portsOfKind(ReactorPortKind.ITEM_INPUT)
        val itemOutputs = portsOfKind(ReactorPortKind.ITEM_OUTPUT)
        val fluidInputs = portsOfKind(ReactorPortKind.FLUID_INPUT)
        val fluidOutputs = portsOfKind(ReactorPortKind.FLUID_OUTPUT)

        return NativeReactorMultiblockBinding(
            structureId = structureId,
            nativeReactorId = reactorId,
            itemInputs = itemInputs.toBindings(startIndex = 0),
            itemOutputs = itemOutputs.toBindings(startIndex = 0),
            fluidInputs = fluidInputs.toBindings(startIndex = itemInputs.size),
            fluidOutputs = fluidOutputs.toBindings(startIndex = itemOutputs.size),
        )
    }

    override fun removeNativeReactor(binding: NativeReactorMultiblockBinding) {
        ThermodynamicsNative.removeReactor(binding.nativeReactorId)
    }

    override fun tickNativeReactor(binding: NativeReactorMultiblockBinding, dtSeconds: Double) {
        ThermodynamicsNative.tickReactor(binding.nativeReactorId, dtSeconds)
    }

    override fun readZoneMetrics(
        binding: NativeReactorMultiblockBinding,
        zoneIndex: Int,
        simulatedSeconds: Double,
    ): ReactorTickMetrics {
        val snapshot = ThermodynamicsNative.reactorZoneSnapshot(binding.nativeReactorId, zoneIndex)
        return ReactorTickMetrics(
            simulatedSeconds = simulatedSeconds,
            temperatureKelvin = snapshot.temperatureKelvin,
            pressurePascal = snapshot.pressurePascal,
            substances = snapshot.substances.map { substance ->
                ReactorMixtureSubstanceMetric(
                    substanceId = substance.substanceId,
                    concentrationMolPerBucket = substance.concentrationMolPerBucket,
                )
            },
        )
    }

    override fun insertItemStack(
        binding: NativeReactorMultiblockBinding,
        itemInputPort: ReactorPortDescriptor,
        itemId: String,
        itemCount: Int,
    ): Int {
        val nativeInputIndex = binding.itemInputs
            .singleOrNull { it.port == itemInputPort }
            ?.nativePortIndex
            ?: throw IllegalArgumentException("port ${itemInputPort.position} is not an item input of this reactor")
        return ThermodynamicsNative.insertItemStackToReactorInput(
            reactorId = binding.nativeReactorId,
            inputIndex = nativeInputIndex,
            itemId = itemId,
            itemCount = itemCount,
        )
    }

    override fun extractItemStack(
        binding: NativeReactorMultiblockBinding,
        itemOutputPort: ReactorPortDescriptor,
        itemId: String,
        maxItemCount: Int,
        dtSeconds: Double,
    ): Int {
        val nativeOutputIndex = binding.itemOutputs
            .singleOrNull { it.port == itemOutputPort }
            ?.nativePortIndex
            ?: throw IllegalArgumentException("port ${itemOutputPort.position} is not an item output of this reactor")
        return ThermodynamicsNative.extractItemStackFromReactorOutput(
            reactorId = binding.nativeReactorId,
            outputIndex = nativeOutputIndex,
            itemId = itemId,
            maxItemCount = maxItemCount,
            dtSeconds = dtSeconds,
        )
    }

    private fun List<ReactorPortDescriptor>.toBindings(startIndex: Int): List<NativeReactorPortBinding> =
        mapIndexed { offset, port ->
            NativeReactorPortBinding(
                port = port,
                nativePortIndex = startIndex + offset,
            )
        }
}
