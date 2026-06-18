package dev.makargravanov.create_thermodynamics.common.reactor.multiblock

import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative
import java.util.UUID
import kotlin.test.Test
import kotlin.test.assertEquals

class ReactorMultiblockNativeBridgeTest {
    @Test
    fun `creates native reactor and maps typed ports to native indexes`() {
        val definition = ReactorMultiblockAssembler(
            ReactorMultiblockRules(chamberVolumeCubicMeters = 0.001),
        ).assemble(
            structureId = UUID.fromString("96913175-87c7-4bf7-ab1d-f2f9e5d12924"),
            blocks = listOf(
                block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                block(0, -1, 0, ReactorMultiblockBlockKind.FLUID_INPUT_PORT),
                block(1, 0, 0, ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT),
                block(0, 0, 1, ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT),
            ),
        )
        val countBefore = ThermodynamicsNative.reactorCount()

        val binding = ReactorMultiblockNativeBridge.createNativeReactor(definition)
        try {
            assertEquals(countBefore + 1, ThermodynamicsNative.reactorCount())
            assertEquals(definition.structureId, binding.structureId)
            assertEquals(0, binding.itemInputs.single().nativePortIndex)
            assertEquals(1, binding.fluidInputs.single().nativePortIndex)
            assertEquals(0, binding.itemOutputs.single().nativePortIndex)
            assertEquals(1, binding.fluidOutputs.single().nativePortIndex)
        } finally {
            ReactorMultiblockNativeBridge.removeNativeReactor(binding)
        }

        assertEquals(countBefore, ThermodynamicsNative.reactorCount())
    }

    @Test
    fun `inserts item stack through bound item input port`() {
        ThermodynamicsNative.configureItemChemicalBindings(
            listOf(
                ThermodynamicsNative.ItemChemicalBinding(
                    itemId = "minecraft:water_bucket",
                    substanceId = "destroy:water",
                    molPerItem = 1.0,
                ),
            ),
        )
        val definition = ReactorMultiblockAssembler(
            ReactorMultiblockRules(chamberVolumeCubicMeters = 0.001),
        ).assemble(
            structureId = UUID.fromString("3ad9bb6d-e98b-45f7-8e54-8d04ba264e7a"),
            blocks = listOf(
                block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
            ),
        )

        val binding = ReactorMultiblockNativeBridge.createNativeReactor(definition)
        try {
            val molInserted = ReactorMultiblockNativeBridge.insertItemStack(
                binding = binding,
                itemInputPort = binding.itemInputs.single().port,
                itemId = "minecraft:water_bucket",
                itemCount = 2,
            )

            assertEquals(2.0, molInserted)
        } finally {
            ReactorMultiblockNativeBridge.removeNativeReactor(binding)
        }
    }

    private fun block(
        x: Int,
        y: Int,
        z: Int,
        kind: ReactorMultiblockBlockKind,
    ): ReactorMultiblockBlock =
        ReactorMultiblockBlock(ReactorBlockPosition(x, y, z), kind)
}
