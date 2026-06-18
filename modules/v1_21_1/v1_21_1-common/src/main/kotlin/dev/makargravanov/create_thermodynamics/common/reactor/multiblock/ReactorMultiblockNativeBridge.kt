package dev.makargravanov.create_thermodynamics.common.reactor.multiblock

import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative
import java.util.UUID

data class NativeReactorPortBinding(
    val port: ReactorPortDescriptor,
    val nativePortIndex: Int,
)

data class NativeReactorMultiblockBinding(
    val structureId: UUID,
    val nativeReactorId: ThermodynamicsNative.NativeReactorId,
    val itemInputs: List<NativeReactorPortBinding>,
    val itemOutputs: List<NativeReactorPortBinding>,
    val fluidInputs: List<NativeReactorPortBinding>,
    val fluidOutputs: List<NativeReactorPortBinding>,
) {
    fun allPortBindings(): List<NativeReactorPortBinding> =
        itemInputs + itemOutputs + fluidInputs + fluidOutputs
}

object ReactorMultiblockNativeBridge {
    fun createNativeReactor(definition: ReactorMultiblockDefinition): NativeReactorMultiblockBinding {
        val itemInputs = definition.portsOfKind(ReactorPortKind.ITEM_INPUT)
        val itemOutputs = definition.portsOfKind(ReactorPortKind.ITEM_OUTPUT)
        val fluidInputs = definition.portsOfKind(ReactorPortKind.FLUID_INPUT)
        val fluidOutputs = definition.portsOfKind(ReactorPortKind.FLUID_OUTPUT)

        val reactorId = ThermodynamicsNative.createSingleZoneReactor(
            volumeCubicMeters = definition.totalVolumeCubicMeters,
            itemInputCount = itemInputs.size,
            itemOutputCount = itemOutputs.size,
            fluidInputCount = fluidInputs.size,
            fluidOutputCount = fluidOutputs.size,
        )

        return NativeReactorMultiblockBinding(
            structureId = definition.structureId,
            nativeReactorId = reactorId,
            itemInputs = itemInputs.toBindings(startIndex = 0),
            itemOutputs = itemOutputs.toBindings(startIndex = 0),
            fluidInputs = fluidInputs.toBindings(startIndex = itemInputs.size),
            fluidOutputs = fluidOutputs.toBindings(startIndex = itemOutputs.size),
        )
    }

    fun removeNativeReactor(binding: NativeReactorMultiblockBinding) {
        ThermodynamicsNative.removeReactor(binding.nativeReactorId)
    }

    fun tickNativeReactor(binding: NativeReactorMultiblockBinding, dtSeconds: Double) {
        ThermodynamicsNative.tickReactor(binding.nativeReactorId, dtSeconds)
    }

    fun insertItemStack(
        binding: NativeReactorMultiblockBinding,
        itemInputPort: ReactorPortDescriptor,
        itemId: String,
        itemCount: Int,
    ): Double {
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

    private fun List<ReactorPortDescriptor>.toBindings(startIndex: Int): List<NativeReactorPortBinding> =
        mapIndexed { offset, port ->
            NativeReactorPortBinding(
                port = port,
                nativePortIndex = startIndex + offset,
            )
        }
}
