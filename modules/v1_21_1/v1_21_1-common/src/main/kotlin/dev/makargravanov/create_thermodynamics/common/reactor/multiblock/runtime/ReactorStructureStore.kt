package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationRejection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorPortAccess
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import java.util.LinkedHashMap

class ReactorStructureStore(
    private val nativeBridge: NativeReactorBridge = ReactorMultiblockNativeBridge,
) : ReactorPortAccess {
    private val records = LinkedHashMap<ReactorStructureId, ReactorStructureRecord>()

    fun register(definition: ReactorMultiblockDefinition): ReactorStructureRecord {
        val existing = records[definition.structureId]
        if (existing != null && existing.state != ReactorStructureState.REMOVED) {
            throw IllegalStateException("reactor structure ${definition.structureId.value} is already registered")
        }

        val binding = nativeBridge.createNativeReactor(definition)
        val record = ReactorStructureRecord(
            structureId = definition.structureId,
            definition = definition,
            nativeBinding = binding,
            state = ReactorStructureState.ACTIVE,
        )
        records[definition.structureId] = record
        return record
    }

    fun record(structureId: ReactorStructureId): ReactorStructureRecord? =
        records[structureId]

    fun activeRecords(): List<ReactorStructureRecord> =
        records.values.filter { it.state == ReactorStructureState.ACTIVE }

    fun remove(structureId: ReactorStructureId): ReactorOperationResult {
        val record = records[structureId]
            ?: return rejected(ReactorOperationRejection.STRUCTURE_NOT_FOUND, "reactor structure ${structureId.value} is not registered")
        if (record.state != ReactorStructureState.ACTIVE) {
            return rejected(ReactorOperationRejection.STRUCTURE_NOT_ACTIVE, "reactor structure ${structureId.value} is ${record.state}")
        }

        return try {
            nativeBridge.removeNativeReactor(record.nativeBinding)
            records[structureId] = record.copy(state = ReactorStructureState.REMOVED)
            ReactorOperationResult.Completed
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed("failed to remove native reactor for structure ${structureId.value}: ${error.message}")
        }
    }

    fun tick(structureId: ReactorStructureId, dtSeconds: Double): ReactorOperationResult {
        if (!dtSeconds.isFinite() || dtSeconds < 0.0) {
            return ReactorOperationResult.Failed("dtSeconds must be non-negative and finite")
        }
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)

        return try {
            nativeBridge.tickNativeReactor(record.nativeBinding, dtSeconds)
            ReactorOperationResult.Completed
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed("failed to tick native reactor for structure ${structureId.value}: ${error.message}")
        }
    }

    fun tickAll(dtSeconds: Double): List<ReactorOperationResult> =
        activeRecords().map { tick(it.structureId, dtSeconds) }

    override fun insertItem(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        itemId: String,
        itemCount: Int,
    ): ReactorOperationResult {
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)
        val port = record.inputPort(portPosition, ReactorPortKind.ITEM_INPUT) ?: return portRejected(record, portPosition, ReactorPortKind.ITEM_INPUT)

        return try {
            ReactorOperationResult.ItemInserted(
                nativeBridge.insertItemStack(
                    binding = record.nativeBinding,
                    itemInputPort = port,
                    itemId = itemId,
                    itemCount = itemCount,
                ),
            )
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed("failed to insert item stack into reactor ${structureId.value}: ${error.message}")
        }
    }

    override fun extractItem(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        itemId: String,
        maxItemCount: Int,
    ): ReactorOperationResult {
        validatePort(structureId, portPosition, ReactorPortKind.ITEM_OUTPUT)?.let { return it }
        return rejected(ReactorOperationRejection.OPERATION_NOT_SUPPORTED, "item extraction is not wired to the native reactor yet")
    }

    override fun insertFluid(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        fluidId: String,
        millibuckets: Int,
    ): ReactorOperationResult {
        validatePort(structureId, portPosition, ReactorPortKind.FLUID_INPUT)?.let { return it }
        return rejected(ReactorOperationRejection.OPERATION_NOT_SUPPORTED, "fluid insertion is not wired to the native reactor yet")
    }

    override fun extractFluid(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        fluidId: String,
        maxMillibuckets: Int,
    ): ReactorOperationResult {
        validatePort(structureId, portPosition, ReactorPortKind.FLUID_OUTPUT)?.let { return it }
        return rejected(ReactorOperationRejection.OPERATION_NOT_SUPPORTED, "fluid extraction is not wired to the native reactor yet")
    }

    private fun activeRecord(structureId: ReactorStructureId): ReactorStructureRecord? =
        records[structureId]?.takeIf { it.state == ReactorStructureState.ACTIVE }

    private fun inactiveOrMissing(structureId: ReactorStructureId): ReactorOperationResult {
        val record = records[structureId]
        return if (record == null) {
            rejected(ReactorOperationRejection.STRUCTURE_NOT_FOUND, "reactor structure ${structureId.value} is not registered")
        } else {
            rejected(ReactorOperationRejection.STRUCTURE_NOT_ACTIVE, "reactor structure ${structureId.value} is ${record.state}")
        }
    }

    private fun ReactorStructureRecord.inputPort(
        portPosition: ReactorBlockPosition,
        kind: ReactorPortKind,
    ): ReactorPortDescriptor? =
        portOfKind(portPosition, kind)

    private fun validatePort(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        expectedKind: ReactorPortKind,
    ): ReactorOperationResult? {
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)
        if (record.portOfKind(portPosition, expectedKind) == null) {
            return portRejected(record, portPosition, expectedKind)
        }
        return null
    }

    private fun portRejected(
        record: ReactorStructureRecord,
        portPosition: ReactorBlockPosition,
        expectedKind: ReactorPortKind,
    ): ReactorOperationResult {
        val actual = record.portAt(portPosition)
        return if (actual == null) {
            rejected(ReactorOperationRejection.PORT_NOT_FOUND, "reactor port $portPosition is not part of structure ${record.structureId.value}")
        } else {
            rejected(
                ReactorOperationRejection.WRONG_PORT_KIND,
                "reactor port $portPosition is ${actual.kind}, expected $expectedKind",
            )
        }
    }

    private fun rejected(
        reason: ReactorOperationRejection,
        message: String,
    ): ReactorOperationResult.Rejected =
        ReactorOperationResult.Rejected(reason, message)
}
