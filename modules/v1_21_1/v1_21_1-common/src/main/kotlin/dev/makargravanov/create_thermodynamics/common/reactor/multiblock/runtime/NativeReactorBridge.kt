package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor

interface NativeReactorBridge {
    fun createNativeReactor(definition: ReactorMultiblockDefinition): NativeReactorMultiblockBinding

    fun removeNativeReactor(binding: NativeReactorMultiblockBinding)

    fun tickNativeReactor(binding: NativeReactorMultiblockBinding, dtSeconds: Double)

    fun insertItemStack(
        binding: NativeReactorMultiblockBinding,
        itemInputPort: ReactorPortDescriptor,
        itemId: String,
        itemCount: Int,
    ): Double
}
