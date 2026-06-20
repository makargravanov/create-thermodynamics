package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationRejection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorPortAccess
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobKind
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorageResult
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

    fun reconcileFreshStructures(
        definitions: Collection<ReactorMultiblockDefinition>,
        removeMissingStructureIds: Set<ReactorStructureId> = emptySet(),
    ): List<ReactorOperationResult> {
        val results = mutableListOf<ReactorOperationResult>()
        val nextDefinitions = definitions.associateBy { it.structureId }
        val removedIds = removeMissingStructureIds
            .filter { it !in nextDefinitions }

        for (structureId in removedIds) {
            results += remove(structureId)
        }

        for (definition in definitions) {
            val existing = records[definition.structureId]
            if (existing != null && existing.state != ReactorStructureState.REMOVED) {
                if (existing.definition == definition && existing.state == ReactorStructureState.ACTIVE) {
                    continue
                }
                results += remove(definition.structureId)
            }
            register(definition)
            results += ReactorOperationResult.Completed
        }

        return results
    }

    fun removeStructures(structureIds: Collection<ReactorStructureId>): List<ReactorOperationResult> =
        structureIds.map { remove(it) }

    fun record(structureId: ReactorStructureId): ReactorStructureRecord? =
        records[structureId]

    fun activeRecords(): List<ReactorStructureRecord> =
        records.values.filter { it.state == ReactorStructureState.ACTIVE }

    fun suspendToCheckpoint(
        structureId: ReactorStructureId,
        blobStorage: NativeBlobStorage,
        contentVersion: Long,
    ): ReactorOperationResult {
        if (contentVersion < 0) {
            return rejected(ReactorOperationRejection.INVALID_CONTENT_VERSION, "contentVersion must be non-negative")
        }
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)
        val binding = record.activeBinding() ?: return inactiveOrMissing(structureId)

        val encoded = try {
            nativeBridge.exportReactorCheckpoint(binding, contentVersion)
        } catch (error: RuntimeException) {
            return ReactorOperationResult.Failed("failed to export reactor checkpoint for structure ${structureId.value}: ${error.message}")
        }

        val stored = when (val result = blobStorage.store(reactorCheckpointStorageKey(structureId, contentVersion), encoded)) {
            is NativeBlobStorageResult.Stored -> result
            is NativeBlobStorageResult.Rejected -> {
                return rejected(
                    ReactorOperationRejection.SNAPSHOT_STORAGE_REJECTED,
                    "failed to store reactor checkpoint for structure ${structureId.value}: ${result.message}",
                )
            }
            is NativeBlobStorageResult.Loaded -> {
                return ReactorOperationResult.Failed("native blob storage returned Loaded while storing reactor checkpoint")
            }
        }
        if (stored.ref.kind != NativeBlobKind.ReactorSnapshot) {
            return ReactorOperationResult.Failed("native checkpoint for structure ${structureId.value} is ${stored.ref.kind}, expected ${NativeBlobKind.ReactorSnapshot}")
        }

        return try {
            nativeBridge.removeNativeReactor(binding)
            records[structureId] = record.copy(
                nativeBinding = null,
                reactorCheckpoint = stored.ref,
                snapshotVersion = ReactorSnapshotVersion(contentVersion),
                state = ReactorStructureState.SUSPENDED_UNLOADED,
            )
            ReactorOperationResult.ReactorSuspended("reactor structure ${structureId.value} was suspended to ${stored.ref.storageKey}")
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed("failed to remove native reactor after checkpoint for structure ${structureId.value}: ${error.message}")
        }
    }

    fun exportCheckpoint(
        structureId: ReactorStructureId,
        blobStorage: NativeBlobStorage,
        contentVersion: Long,
    ): ReactorOperationResult {
        if (contentVersion < 0) {
            return rejected(ReactorOperationRejection.INVALID_CONTENT_VERSION, "contentVersion must be non-negative")
        }
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)
        val binding = record.activeBinding() ?: return inactiveOrMissing(structureId)

        val encoded = try {
            nativeBridge.exportReactorCheckpoint(binding, contentVersion)
        } catch (error: RuntimeException) {
            return ReactorOperationResult.Failed("failed to export reactor checkpoint for structure ${structureId.value}: ${error.message}")
        }

        val stored = when (val result = blobStorage.store(reactorCheckpointStorageKey(structureId, contentVersion), encoded)) {
            is NativeBlobStorageResult.Stored -> result
            is NativeBlobStorageResult.Rejected -> {
                return rejected(
                    ReactorOperationRejection.SNAPSHOT_STORAGE_REJECTED,
                    "failed to store reactor checkpoint for structure ${structureId.value}: ${result.message}",
                )
            }
            is NativeBlobStorageResult.Loaded -> {
                return ReactorOperationResult.Failed("native blob storage returned Loaded while storing reactor checkpoint")
            }
        }
        if (stored.ref.kind != NativeBlobKind.ReactorSnapshot) {
            return ReactorOperationResult.Failed("native checkpoint for structure ${structureId.value} is ${stored.ref.kind}, expected ${NativeBlobKind.ReactorSnapshot}")
        }

        return ReactorOperationResult.ReactorCheckpointExported(stored.ref)
    }

    fun resumeFromCheckpoint(
        structureId: ReactorStructureId,
        blobStorage: NativeBlobStorage,
    ): ReactorOperationResult {
        val record = records[structureId]
            ?: return rejected(ReactorOperationRejection.STRUCTURE_NOT_FOUND, "reactor structure ${structureId.value} is not registered")
        if (record.state != ReactorStructureState.SUSPENDED_UNLOADED) {
            return rejected(
                ReactorOperationRejection.STRUCTURE_NOT_SUSPENDED,
                "reactor structure ${structureId.value} is ${record.state}, expected ${ReactorStructureState.SUSPENDED_UNLOADED}",
            )
        }
        val checkpoint = record.reactorCheckpoint
            ?: return ReactorOperationResult.Failed("reactor structure ${structureId.value} is suspended without a checkpoint reference")

        val bytes = when (val result = blobStorage.load(checkpoint)) {
            is NativeBlobStorageResult.Loaded -> result.bytes
            is NativeBlobStorageResult.Rejected -> {
                return rejected(
                    ReactorOperationRejection.SNAPSHOT_STORAGE_REJECTED,
                    "failed to load reactor checkpoint for structure ${structureId.value}: ${result.message}",
                )
            }
            is NativeBlobStorageResult.Stored -> {
                return ReactorOperationResult.Failed("native blob storage returned Stored while loading reactor checkpoint")
            }
        }

        return try {
            val binding = nativeBridge.createNativeReactorFromCheckpoint(record.definition, bytes)
            records[structureId] = record.copy(
                nativeBinding = binding,
                snapshotVersion = ReactorSnapshotVersion(checkpoint.contentVersion),
                state = ReactorStructureState.ACTIVE,
            )
            ReactorOperationResult.ReactorResumed("reactor structure ${structureId.value} was resumed from ${checkpoint.storageKey}")
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed("failed to create native reactor from checkpoint for structure ${structureId.value}: ${error.message}")
        }
    }

    fun remove(structureId: ReactorStructureId): ReactorOperationResult {
        val record = records[structureId]
            ?: return rejected(ReactorOperationRejection.STRUCTURE_NOT_FOUND, "reactor structure ${structureId.value} is not registered")
        if (record.state == ReactorStructureState.REMOVED) {
            return rejected(ReactorOperationRejection.STRUCTURE_NOT_ACTIVE, "reactor structure ${structureId.value} is ${record.state}")
        }

        return try {
            record.nativeBinding?.let { nativeBridge.removeNativeReactor(it) }
            records[structureId] = record.copy(
                nativeBinding = null,
                reactorCheckpoint = null,
                state = ReactorStructureState.REMOVED,
            )
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
        val binding = record.activeBinding() ?: return inactiveOrMissing(structureId)

        return try {
            nativeBridge.tickNativeReactor(binding, dtSeconds)
            ReactorOperationResult.Completed
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed("failed to tick native reactor for structure ${structureId.value}: ${error.message}")
        }
    }

    fun readZoneMetrics(
        structureId: ReactorStructureId,
        zoneIndex: Int,
        simulatedSeconds: Double,
    ): ReactorOperationResult {
        if (zoneIndex < 0) {
            return ReactorOperationResult.Failed("zoneIndex must be non-negative")
        }
        if (!simulatedSeconds.isFinite() || simulatedSeconds < 0.0) {
            return ReactorOperationResult.Failed("simulatedSeconds must be non-negative and finite")
        }
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)
        val binding = record.activeBinding() ?: return inactiveOrMissing(structureId)
        return try {
            ReactorOperationResult.ReactorMetricsRead(
                nativeBridge.readZoneMetrics(binding, zoneIndex, simulatedSeconds),
            )
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed("failed to read native reactor metrics for structure ${structureId.value}: ${error.message}")
        }
    }

    fun tickAll(dtSeconds: Double): List<ReactorOperationResult> =
        activeRecords().map { tick(it.structureId, dtSeconds) }

    fun validateItemInputPort(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
    ): ReactorOperationResult.Rejected? =
        validatePort(structureId, portPosition, ReactorPortKind.ITEM_INPUT) as? ReactorOperationResult.Rejected

    fun validateItemOutputPort(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
    ): ReactorOperationResult.Rejected? =
        validatePort(structureId, portPosition, ReactorPortKind.ITEM_OUTPUT) as? ReactorOperationResult.Rejected

    fun applyReport(report: ReactorReport): ReactorOperationResult {
        val record = records[report.structureId]
            ?: return rejected(
                ReactorOperationRejection.STRUCTURE_NOT_FOUND,
                "reactor structure ${report.structureId.value} is not registered",
            )
        if (report.snapshotVersion.value < record.snapshotVersion.value) {
            return rejected(
                ReactorOperationRejection.STALE_REPORT,
                "reactor report ${report.reportId.value} has snapshot ${report.snapshotVersion.value}, current snapshot is ${record.snapshotVersion.value}",
            )
        }

        records[report.structureId] = when (report) {
            is ReactorReport.SnapshotReady -> record.copy(
                reactorCheckpoint = report.checkpoint,
                snapshotVersion = report.snapshotVersion,
            )
            is ReactorReport.ReactorFailed -> record.copy(
                snapshotVersion = report.snapshotVersion,
                state = ReactorStructureState.INVALID,
            )
            else -> record.copy(snapshotVersion = report.snapshotVersion)
        }
        return ReactorOperationResult.Completed
    }

    override fun insertItem(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        itemId: String,
        itemCount: Int,
    ): ReactorOperationResult {
        if (itemId.isBlank()) {
            return rejected(ReactorOperationRejection.INVALID_ITEM_ID, "itemId must not be blank")
        }
        if (itemCount <= 0) {
            return rejected(ReactorOperationRejection.INVALID_ITEM_COUNT, "itemCount must be positive")
        }
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)
        val port = record.inputPort(portPosition, ReactorPortKind.ITEM_INPUT) ?: return portRejected(record, portPosition, ReactorPortKind.ITEM_INPUT)
        return try {
            val inserted = nativeBridge.insertItemStack(
                binding = record.activeBinding() ?: return inactiveOrMissing(structureId),
                itemInputPort = port,
                itemId = itemId,
                itemCount = itemCount,
            )
            ReactorOperationResult.ItemInserted(inserted)
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed(
                "native insertion failed for reactor input $portPosition: ${error.message}",
            )
        }
    }

    override fun extractItem(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        itemId: String,
        maxItemCount: Int,
        dtSeconds: Double,
    ): ReactorOperationResult {
        if (itemId.isBlank()) {
            return rejected(ReactorOperationRejection.INVALID_ITEM_ID, "itemId must not be blank")
        }
        if (maxItemCount <= 0) {
            return rejected(ReactorOperationRejection.INVALID_ITEM_COUNT, "maxItemCount must be positive")
        }
        if (!dtSeconds.isFinite() || dtSeconds <= 0.0) {
            return ReactorOperationResult.Failed("dtSeconds must be positive and finite")
        }
        val record = activeRecord(structureId) ?: return inactiveOrMissing(structureId)
        val port = record.inputPort(portPosition, ReactorPortKind.ITEM_OUTPUT) ?: return portRejected(record, portPosition, ReactorPortKind.ITEM_OUTPUT)
        return try {
            val extracted = nativeBridge.extractItemStack(
                binding = record.activeBinding() ?: return inactiveOrMissing(structureId),
                itemOutputPort = port,
                itemId = itemId,
                maxItemCount = maxItemCount,
                dtSeconds = dtSeconds,
            )
            ReactorOperationResult.ItemExtracted(
                itemId = itemId,
                itemCount = extracted,
            )
        } catch (error: RuntimeException) {
            ReactorOperationResult.Failed(
                "native extraction failed for reactor output $portPosition: ${error.message}",
            )
        }
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

    private fun ReactorStructureRecord.activeBinding(): NativeReactorMultiblockBinding? =
        nativeBinding?.takeIf { state == ReactorStructureState.ACTIVE }

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

    private fun reactorCheckpointStorageKey(
        structureId: ReactorStructureId,
        contentVersion: Long,
    ): String =
        "reactors/${structureId.value}/checkpoint_${contentVersion.toString().padStart(12, '0')}.bin.zst"
}
