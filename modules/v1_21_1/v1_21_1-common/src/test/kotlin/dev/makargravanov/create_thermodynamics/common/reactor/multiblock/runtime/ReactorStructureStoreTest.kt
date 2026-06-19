package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.chemistry.binding.ItemChemicalBinding
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationRejection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockAssembler
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockRules
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockDirection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlock
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobKind
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage
import java.util.UUID
import java.nio.ByteBuffer
import java.nio.ByteOrder
import java.nio.file.Files
import java.security.MessageDigest
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertIs

class ReactorStructureStoreTest {
    @Test
    fun `registers active structure and creates native binding`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()

        val record = store.register(definition)

        assertEquals(ReactorStructureState.ACTIVE, record.state)
        assertEquals(definition.structureId, record.structureId)
        assertEquals(1, nativeBridge.created.size)
        assertEquals(listOf(record), store.activeRecords())
    }

    @Test
    fun `rejects repeated registration of active structure`() {
        val store = ReactorStructureStore(FakeNativeReactorBridge())
        val definition = testDefinition()
        store.register(definition)

        assertFailsWith<IllegalStateException> {
            store.register(definition)
        }
    }

    @Test
    fun `removes structure and skips it during tick all`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)

        assertEquals(ReactorOperationResult.Completed, store.remove(definition.structureId))
        assertEquals(ReactorStructureState.REMOVED, store.record(definition.structureId)?.state)
        assertEquals(1, nativeBridge.removed.size)
        assertEquals(emptyList(), store.tickAll(0.05))
    }

    @Test
    fun `suspends active structure to native checkpoint`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val storage = NativeBlobStorage(Files.createTempDirectory("ct-reactor-checkpoint-test"))

        val result = store.suspendToCheckpoint(
            structureId = definition.structureId,
            blobStorage = storage,
            contentVersion = 7,
        )

        assertIs<ReactorOperationResult.ReactorSuspended>(result)
        val record = store.record(definition.structureId)!!
        assertEquals(ReactorStructureState.SUSPENDED_UNLOADED, record.state)
        assertEquals(NativeBlobKind.ReactorSnapshot, record.reactorCheckpoint?.kind)
        assertEquals(1, nativeBridge.exportedCheckpoints.size)
        assertEquals(1, nativeBridge.removed.size)
        assertEquals(emptyList(), store.tickAll(0.05))
    }

    @Test
    fun `resumes suspended structure from native checkpoint`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val storage = NativeBlobStorage(Files.createTempDirectory("ct-reactor-checkpoint-test"))
        store.suspendToCheckpoint(definition.structureId, storage, contentVersion = 8)

        val result = store.resumeFromCheckpoint(definition.structureId, storage)

        assertIs<ReactorOperationResult.ReactorResumed>(result)
        val record = store.record(definition.structureId)!!
        assertEquals(ReactorStructureState.ACTIVE, record.state)
        assertEquals(1, nativeBridge.restoredCheckpoints.size)
        assertEquals(listOf(ReactorOperationResult.Completed), store.tickAll(0.05))
    }

    @Test
    fun `rejects resume of active structure`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val storage = NativeBlobStorage(Files.createTempDirectory("ct-reactor-checkpoint-test"))

        val result = store.resumeFromCheckpoint(definition.structureId, storage)

        val rejected = assertIs<ReactorOperationResult.Rejected>(result)
        assertEquals(ReactorOperationRejection.STRUCTURE_NOT_SUSPENDED, rejected.reason)
    }

    @Test
    fun `ticks only active structures`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val first = testDefinition("96913175-87c7-4bf7-ab1d-f2f9e5d12924")
        val second = testDefinition("3ad9bb6d-e98b-45f7-8e54-8d04ba264e7a")
        store.register(first)
        store.register(second)
        store.remove(first.structureId)

        assertEquals(listOf(ReactorOperationResult.Completed), store.tickAll(0.05))
        assertEquals(listOf(second.structureId), nativeBridge.ticked)
    }

    @Test
    fun `inserts item through item input port`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val itemInput = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single()

        val result = store.insertItem(
            structureId = definition.structureId,
            portPosition = itemInput.position,
            itemId = "minecraft:water_bucket",
            itemCount = 2,
        )

        assertEquals(ReactorOperationResult.ItemInserted(2), result)
        assertEquals(listOf(itemInput.position), nativeBridge.itemInsertPorts)
    }

    @Test
    fun `native insertion failure is reported without hidden item buffering`() {
        val nativeBridge = FakeNativeReactorBridge()
        nativeBridge.failItemInsertion = true
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val itemInput = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single()

        val result = store.insertItem(
            structureId = definition.structureId,
            portPosition = itemInput.position,
            itemId = "minecraft:water_bucket",
            itemCount = 2,
        )

        assertIs<ReactorOperationResult.Failed>(result)
    }

    @Test
    fun `rejects invalid item count before touching input buffer`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val itemInput = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single()

        val result = store.insertItem(
            structureId = definition.structureId,
            portPosition = itemInput.position,
            itemId = "minecraft:water_bucket",
            itemCount = 0,
        )

        val rejected = assertIs<ReactorOperationResult.Rejected>(result)
        assertEquals(ReactorOperationRejection.INVALID_ITEM_COUNT, rejected.reason)
        assertEquals(0, nativeBridge.itemInsertPorts.size)
    }

    @Test
    fun `rejects blank item id before touching input buffer`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val itemInput = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single()

        val result = store.insertItem(
            structureId = definition.structureId,
            portPosition = itemInput.position,
            itemId = " ",
            itemCount = 1,
        )

        val rejected = assertIs<ReactorOperationResult.Rejected>(result)
        assertEquals(ReactorOperationRejection.INVALID_ITEM_ID, rejected.reason)
        assertEquals(0, nativeBridge.itemInsertPorts.size)
    }

    @Test
    fun `rejects item insertion through non item input port`() {
        val store = ReactorStructureStore(FakeNativeReactorBridge())
        val definition = testDefinition()
        store.register(definition)
        val fluidInput = definition.portsOfKind(ReactorPortKind.FLUID_INPUT).single()

        val result = store.insertItem(
            structureId = definition.structureId,
            portPosition = fluidInput.position,
            itemId = "minecraft:water_bucket",
            itemCount = 1,
        )

        val rejected = assertIs<ReactorOperationResult.Rejected>(result)
        assertEquals(ReactorOperationRejection.WRONG_PORT_KIND, rejected.reason)
    }

    @Test
    fun `fluid insertion validates fluid input port before reporting unsupported native operation`() {
        val store = ReactorStructureStore(FakeNativeReactorBridge())
        val definition = testDefinition()
        store.register(definition)
        val fluidInput = definition.portsOfKind(ReactorPortKind.FLUID_INPUT).single()

        val result = store.insertFluid(
            structureId = definition.structureId,
            portPosition = fluidInput.position,
            fluidId = "minecraft:water",
            millibuckets = 1000,
        )

        val rejected = assertIs<ReactorOperationResult.Rejected>(result)
        assertEquals(ReactorOperationRejection.OPERATION_NOT_SUPPORTED, rejected.reason)
    }

    @Test
    fun `item extraction rejects wrong port kind before unsupported native operation`() {
        val store = ReactorStructureStore(FakeNativeReactorBridge())
        val definition = testDefinition()
        store.register(definition)
        val itemInput = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single()

        val result = store.extractItem(
            structureId = definition.structureId,
            portPosition = itemInput.position,
            itemId = "minecraft:water_bucket",
            maxItemCount = 1,
            dtSeconds = 1.0,
        )

        val rejected = assertIs<ReactorOperationResult.Rejected>(result)
        assertEquals(ReactorOperationRejection.WRONG_PORT_KIND, rejected.reason)
    }

    @Test
    fun `extracts item through item output port`() {
        val nativeBridge = FakeNativeReactorBridge()
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val itemOutput = definition.portsOfKind(ReactorPortKind.ITEM_OUTPUT).single()

        val result = store.extractItem(
            structureId = definition.structureId,
            portPosition = itemOutput.position,
            itemId = "minecraft:water_bucket",
            maxItemCount = 3,
            dtSeconds = 1.0,
        )

        assertEquals(ReactorOperationResult.ItemExtracted("minecraft:water_bucket", 3), result)
        assertEquals(listOf(itemOutput.position), nativeBridge.itemExtractPorts)
    }

    @Test
    fun `native extraction failure is reported without changing minecraft inventory`() {
        val nativeBridge = FakeNativeReactorBridge()
        nativeBridge.failItemExtraction = true
        val store = ReactorStructureStore(nativeBridge)
        val definition = testDefinition()
        store.register(definition)
        val itemOutput = definition.portsOfKind(ReactorPortKind.ITEM_OUTPUT).single()

        val result = store.extractItem(
            structureId = definition.structureId,
            portPosition = itemOutput.position,
            itemId = "minecraft:water_bucket",
            maxItemCount = 3,
            dtSeconds = 1.0,
        )

        assertIs<ReactorOperationResult.Failed>(result)
    }

    @Test
    fun `rejects port from another structure`() {
        val store = ReactorStructureStore(FakeNativeReactorBridge())
        val definition = testDefinition()
        store.register(definition)

        val result = store.insertItem(
            structureId = definition.structureId,
            portPosition = ReactorBlockPosition(100, 0, 0),
            itemId = "minecraft:water_bucket",
            itemCount = 1,
        )

        val rejected = assertIs<ReactorOperationResult.Rejected>(result)
        assertEquals(ReactorOperationRejection.PORT_NOT_FOUND, rejected.reason)
    }

    @Test
    fun `native bridge maps typed ports to native indexes`() {
        val definition = testDefinition()
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
    fun `native item insertion is available through structure store`() {
        ThermodynamicsNative.configureItemChemicalBindings(
            listOf(
                ItemChemicalBinding(
                    itemId = "minecraft:water_bucket",
                    substanceId = "destroy:water",
                    molPerItem = 1.0,
                ),
            ),
        )
        val store = ReactorStructureStore()
        val definition = testDefinition()
        store.register(definition)
        val itemInput = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single()

        try {
            val result = store.insertItem(
                structureId = definition.structureId,
                portPosition = itemInput.position,
                itemId = "minecraft:water_bucket",
                itemCount = 2,
            )

            assertEquals(ReactorOperationResult.ItemInserted(2), result)
        } finally {
            store.remove(definition.structureId)
        }
    }

    private fun testDefinition(
        structureId: String = "96913175-87c7-4bf7-ab1d-f2f9e5d12924",
    ): ReactorMultiblockDefinition =
        ReactorMultiblockAssembler(
            ReactorMultiblockRules(chamberVolumeCubicMeters = 0.001),
        ).assemble(
            structureId = UUID.fromString(structureId),
            blocks = listOf(
                block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT, ReactorBlockDirection.UP),
                block(0, -1, 0, ReactorMultiblockBlockKind.FLUID_INPUT_PORT, ReactorBlockDirection.DOWN),
                block(1, 0, 0, ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT),
                block(0, 0, 1, ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT),
            ),
        )

    private fun block(
        x: Int,
        y: Int,
        z: Int,
        kind: ReactorMultiblockBlockKind,
        facing: ReactorBlockDirection? = null,
    ): ReactorMultiblockBlock =
        ReactorMultiblockBlock(ReactorBlockPosition(x, y, z), kind, facing)

    private class FakeNativeReactorBridge : NativeReactorBridge {
        val created = mutableListOf<ReactorStructureId>()
        val removed = mutableListOf<ReactorStructureId>()
        val ticked = mutableListOf<ReactorStructureId>()
        val itemInsertPorts = mutableListOf<ReactorBlockPosition>()
        val itemExtractPorts = mutableListOf<ReactorBlockPosition>()
        val exportedCheckpoints = mutableListOf<ReactorStructureId>()
        val restoredCheckpoints = mutableListOf<ReactorStructureId>()
        var failItemInsertion = false
        var failItemExtraction = false

        override fun createNativeReactor(definition: ReactorMultiblockDefinition): NativeReactorMultiblockBinding {
            created += definition.structureId
            return bindingFor(definition)
        }

        override fun createNativeReactorFromCheckpoint(
            definition: ReactorMultiblockDefinition,
            encodedCheckpoint: ByteArray,
        ): NativeReactorMultiblockBinding {
            restoredCheckpoints += definition.structureId
            return bindingFor(definition)
        }

        private fun bindingFor(definition: ReactorMultiblockDefinition): NativeReactorMultiblockBinding {
            val itemInputs = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).toBindings(0)
            val itemOutputs = definition.portsOfKind(ReactorPortKind.ITEM_OUTPUT).toBindings(0)
            val fluidInputs = definition.portsOfKind(ReactorPortKind.FLUID_INPUT).toBindings(itemInputs.size)
            val fluidOutputs = definition.portsOfKind(ReactorPortKind.FLUID_OUTPUT).toBindings(itemOutputs.size)
            return NativeReactorMultiblockBinding(
                structureId = definition.structureId,
                nativeReactorId = ThermodynamicsNative.NativeReactorId((created.size + restoredCheckpoints.size).toLong()),
                itemInputs = itemInputs,
                itemOutputs = itemOutputs,
                fluidInputs = fluidInputs,
                fluidOutputs = fluidOutputs,
            )
        }

        override fun removeNativeReactor(binding: NativeReactorMultiblockBinding) {
            removed += binding.structureId
        }

        override fun tickNativeReactor(binding: NativeReactorMultiblockBinding, dtSeconds: Double) {
            ticked += binding.structureId
        }

        override fun exportReactorCheckpoint(
            binding: NativeReactorMultiblockBinding,
            contentVersion: Long,
        ): ByteArray {
            exportedCheckpoints += binding.structureId
            return nativeBlobBytes(NativeBlobKind.ReactorSnapshot, contentVersion)
        }

        override fun insertItemStack(
            binding: NativeReactorMultiblockBinding,
            itemInputPort: ReactorPortDescriptor,
            itemId: String,
            itemCount: Int,
        ): Int {
            if (failItemInsertion) {
                throw IllegalStateException("configured item insertion failure")
            }
            itemInsertPorts += itemInputPort.position
            return itemCount
        }

        override fun extractItemStack(
            binding: NativeReactorMultiblockBinding,
            itemOutputPort: ReactorPortDescriptor,
            itemId: String,
            maxItemCount: Int,
            dtSeconds: Double,
        ): Int {
            if (failItemExtraction) {
                throw IllegalStateException("configured item extraction failure")
            }
            itemExtractPorts += itemOutputPort.position
            return maxItemCount
        }

        private fun List<ReactorPortDescriptor>.toBindings(startIndex: Int): List<NativeReactorPortBinding> =
            mapIndexed { offset, port ->
                NativeReactorPortBinding(
                    port = port,
                    nativePortIndex = startIndex + offset,
                )
            }

        private fun nativeBlobBytes(
            kind: NativeBlobKind,
            contentVersion: Long,
        ): ByteArray {
            val modelVersion = "test-model".encodeToByteArray()
            val payload = byteArrayOf(1, 2, 3, 4)
            val payloadHash = MessageDigest.getInstance("SHA-256").digest(payload)
            val header = ByteBuffer.allocate(8 + 2 + 2 + 2 + 8 + 8 + 8 + 32 + modelVersion.size)
                .order(ByteOrder.LITTLE_ENDIAN)
                .put(byteArrayOf(0x43, 0x54, 0x4e, 0x42, 0x4c, 0x42, 0x31, 0x00))
                .putShort(1.toShort())
                .putShort(kind.wireId.toShort())
                .putShort(modelVersion.size.toShort())
                .putLong(contentVersion)
                .putLong(payload.size.toLong())
                .putLong(payload.size.toLong())
                .put(payloadHash)
                .put(modelVersion)
                .array()
            return header + payload
        }
    }
}
