package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationRejection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobKind
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobRef
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorageResult

data class ReactorCatalogRuntimeState(
    val checkpoint: NativeBlobRef? = null,
) {
    val catalogVersion: Long
        get() = checkpoint?.contentVersion ?: 0L
}

enum class ReactorWorldRuntimeRejection {
    STRUCTURE_NOT_FOUND,
    STRUCTURE_NOT_ACTIVE,
    QUEUE_FULL,
    WRONG_BLOB_KIND,
    INVALID_PORT_OPERATION,
    CATALOG_EXPORT_FAILED,
    BLOB_STORAGE_REJECTED,
}

sealed interface ReactorWorldRuntimeResult {
    data class CommandQueued(
        val commandId: ReactorCommandId,
        val queueSize: Int,
    ) : ReactorWorldRuntimeResult

    data class ReportQueued(
        val queueSize: Int,
    ) : ReactorWorldRuntimeResult

    data class ReportsApplied(
        val results: List<ReactorOperationResult>,
    ) : ReactorWorldRuntimeResult

    data class CommandsSubmitted(
        val commandCount: Int,
        val reportCount: Int,
    ) : ReactorWorldRuntimeResult

    data class CatalogUpdated(
        val catalog: ReactorCatalogRuntimeState,
    ) : ReactorWorldRuntimeResult

    data class BatchQueued(
        val commandCount: Int,
        val queueSize: Int,
    ) : ReactorWorldRuntimeResult

    data class Rejected(
        val reason: ReactorWorldRuntimeRejection,
        val message: String,
    ) : ReactorWorldRuntimeResult
}

class ReactorWorldRuntime(
    val structures: ReactorStructureStore = ReactorStructureStore(),
    val blobStorage: NativeBlobStorage,
    private val catalogBridge: NativeCatalogBridge = ThermodynamicsNativeCatalogBridge,
    val commandOutbox: ReactorCommandOutbox = ReactorCommandOutbox(),
    val reportInbox: ReactorReportInbox = ReactorReportInbox(),
    initialCatalog: ReactorCatalogRuntimeState = ReactorCatalogRuntimeState(),
) {
    private var nextCommandId = 0L

    var catalog: ReactorCatalogRuntimeState = initialCatalog
        private set

    fun registerStructure(definition: ReactorMultiblockDefinition): ReactorStructureRecord =
        structures.register(definition)

    fun suspendStructureToCheckpoint(
        structureId: ReactorStructureId,
        contentVersion: Long,
    ): ReactorOperationResult =
        structures.suspendToCheckpoint(structureId, blobStorage, contentVersion)

    fun resumeStructureFromCheckpoint(structureId: ReactorStructureId): ReactorOperationResult =
        structures.resumeFromCheckpoint(structureId, blobStorage)

    fun installCatalogCheckpoint(checkpoint: NativeBlobRef): ReactorWorldRuntimeResult =
        if (checkpoint.kind == NativeBlobKind.DynamicCatalogCheckpoint) {
            val updated = catalog.copy(checkpoint = checkpoint)
            catalog = updated
            ReactorWorldRuntimeResult.CatalogUpdated(updated)
        } else {
            rejected(
                ReactorWorldRuntimeRejection.WRONG_BLOB_KIND,
                "catalog checkpoint blob is ${checkpoint.kind}, expected ${NativeBlobKind.DynamicCatalogCheckpoint}",
            )
        }

    fun exportCatalogCheckpoint(contentVersion: Long): ReactorWorldRuntimeResult {
        if (contentVersion < 0) {
            return rejected(
                ReactorWorldRuntimeRejection.CATALOG_EXPORT_FAILED,
                "catalog checkpoint contentVersion must be non-negative",
            )
        }
        val encoded = try {
            catalogBridge.exportCatalogCheckpoint(contentVersion)
        } catch (error: RuntimeException) {
            return rejected(
                ReactorWorldRuntimeRejection.CATALOG_EXPORT_FAILED,
                "failed to export catalog checkpoint $contentVersion: ${error.message}",
            )
        }
        val stored = when (val result = blobStorage.store(catalogCheckpointStorageKey(contentVersion), encoded)) {
            is NativeBlobStorageResult.Stored -> result
            is NativeBlobStorageResult.Rejected -> {
                return rejected(
                    ReactorWorldRuntimeRejection.BLOB_STORAGE_REJECTED,
                    "failed to store catalog checkpoint $contentVersion: ${result.message}",
                )
            }
            is NativeBlobStorageResult.Loaded -> {
                return rejected(
                    ReactorWorldRuntimeRejection.BLOB_STORAGE_REJECTED,
                    "native blob storage returned Loaded while storing catalog checkpoint",
                )
            }
        }
        return installCatalogCheckpoint(stored.ref)
    }

    fun queueTick(
        structureId: ReactorStructureId,
        dtSeconds: Double,
    ): ReactorWorldRuntimeResult {
        val record = structures.record(structureId)
            ?: return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_FOUND,
                "reactor structure ${structureId.value} is not registered",
            )
        if (record.state != ReactorStructureState.ACTIVE) {
            return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_ACTIVE,
                "reactor structure ${structureId.value} is ${record.state}",
            )
        }
        return enqueueCommand(
            ReactorCommand.Tick(
                commandId = nextCommandId(),
                structureId = structureId,
                expectedSnapshotVersion = record.snapshotVersion,
                dtSeconds = dtSeconds,
            ),
        )
    }

    fun queueTickForActiveStructures(dtSeconds: Double): ReactorWorldRuntimeResult {
        if (!dtSeconds.isFinite() || dtSeconds < 0.0) {
            return rejected(
                ReactorWorldRuntimeRejection.INVALID_PORT_OPERATION,
                "reactor tick duration must be non-negative and finite",
            )
        }
        val records = structures.activeRecords()
        if (records.size > commandOutbox.remainingCapacity) {
            return rejected(
                ReactorWorldRuntimeRejection.QUEUE_FULL,
                "reactor command queue cannot accept ${records.size} tick commands; remaining capacity is ${commandOutbox.remainingCapacity}",
            )
        }
        val commands = records.map { record ->
            ReactorCommand.Tick(
                commandId = nextCommandId(),
                structureId = record.structureId,
                expectedSnapshotVersion = record.snapshotVersion,
                dtSeconds = dtSeconds,
            )
        }
        return enqueueCommands(commands)
    }

    fun queueEnsureLoaded(structureId: ReactorStructureId): ReactorWorldRuntimeResult {
        val record = structures.record(structureId)
            ?: return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_FOUND,
                "reactor structure ${structureId.value} is not registered",
            )
        if (record.state != ReactorStructureState.ACTIVE && record.state != ReactorStructureState.SUSPENDED_UNLOADED) {
            return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_ACTIVE,
                "reactor structure ${structureId.value} is ${record.state}, expected ${ReactorStructureState.ACTIVE} or ${ReactorStructureState.SUSPENDED_UNLOADED}",
            )
        }
        return enqueueCommand(
            ReactorCommand.EnsureLoaded(
                commandId = nextCommandId(),
                structureId = structureId,
                expectedSnapshotVersion = record.snapshotVersion,
                checkpoint = record.reactorCheckpoint,
            ),
        )
    }

    fun queueUnloadAndExportSnapshot(
        structureId: ReactorStructureId,
        reason: String,
    ): ReactorWorldRuntimeResult {
        val record = structures.record(structureId)
            ?: return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_FOUND,
                "reactor structure ${structureId.value} is not registered",
            )
        if (record.state != ReactorStructureState.ACTIVE) {
            return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_ACTIVE,
                "reactor structure ${structureId.value} is ${record.state}, expected ${ReactorStructureState.ACTIVE}",
            )
        }
        return enqueueCommand(
            ReactorCommand.UnloadAndExportSnapshot(
                commandId = nextCommandId(),
                structureId = structureId,
                expectedSnapshotVersion = record.snapshotVersion,
                reason = reason,
            ),
        )
    }

    fun queueRemoveReactor(
        structureId: ReactorStructureId,
        reason: String,
    ): ReactorWorldRuntimeResult {
        val record = structures.record(structureId)
            ?: return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_FOUND,
                "reactor structure ${structureId.value} is not registered",
            )
        if (record.state == ReactorStructureState.REMOVED) {
            return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_ACTIVE,
                "reactor structure ${structureId.value} is ${record.state}",
            )
        }
        return enqueueCommand(
            ReactorCommand.RemoveReactor(
                commandId = nextCommandId(),
                structureId = structureId,
                expectedSnapshotVersion = record.snapshotVersion,
                reason = reason,
            ),
        )
    }

    fun queueInsertItem(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        itemId: String,
        itemCount: Int,
    ): ReactorWorldRuntimeResult {
        val record = structures.record(structureId)
            ?: return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_FOUND,
                "reactor structure ${structureId.value} is not registered",
            )
        if (record.state != ReactorStructureState.ACTIVE) {
            return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_ACTIVE,
                "reactor structure ${structureId.value} is ${record.state}, expected ${ReactorStructureState.ACTIVE}",
            )
        }
        return enqueueCommand(
            ReactorCommand.InsertItem(
                commandId = nextCommandId(),
                structureId = structureId,
                expectedSnapshotVersion = record.snapshotVersion,
                portPosition = portPosition,
                itemId = itemId,
                itemCount = itemCount,
            ),
        )
    }

    fun queueSnapshotRequest(
        structureId: ReactorStructureId,
        reason: String,
    ): ReactorWorldRuntimeResult {
        val record = structures.record(structureId)
            ?: return rejected(
                ReactorWorldRuntimeRejection.STRUCTURE_NOT_FOUND,
                "reactor structure ${structureId.value} is not registered",
            )
        return enqueueCommand(
            ReactorCommand.RequestSnapshot(
                commandId = nextCommandId(),
                structureId = structureId,
                expectedSnapshotVersion = record.snapshotVersion,
                reason = reason,
            ),
        )
    }

    fun receiveReport(report: ReactorReport): ReactorWorldRuntimeResult =
        when (val result = reportInbox.enqueue(report)) {
            is ReactorReportInboxResult.Enqueued -> ReactorWorldRuntimeResult.ReportQueued(result.size)
            is ReactorReportInboxResult.Rejected -> rejected(
                ReactorWorldRuntimeRejection.QUEUE_FULL,
                result.message,
            )
        }

    fun applyReadyReports(maxReports: Int): ReactorWorldRuntimeResult {
        val reports = reportInbox.drain(maxReports)
        return ReactorWorldRuntimeResult.ReportsApplied(
            reports.map { structures.applyReport(it) },
        )
    }

    fun drainCommands(maxCommands: Int): List<ReactorCommand> =
        commandOutbox.drain(maxCommands)

    fun submitQueuedCommands(
        nativeSession: ReactorNativeSession,
        maxCommands: Int,
    ): ReactorWorldRuntimeResult {
        require(maxCommands > 0) { "maxCommands must be positive" }
        val commandCount = commandOutbox.drainableCount(maxCommands)
        if (commandCount == 0) {
            return ReactorWorldRuntimeResult.CommandsSubmitted(commandCount = 0, reportCount = 0)
        }
        if (reportInbox.remainingCapacity < commandCount) {
            return rejected(
                ReactorWorldRuntimeRejection.QUEUE_FULL,
                "reactor report queue cannot accept $commandCount command reports; remaining capacity is ${reportInbox.remainingCapacity}",
            )
        }

        val commands = commandOutbox.drain(commandCount)
        val reports = nativeSession.submit(commands)
        for (report in reports) {
            when (val result = reportInbox.enqueue(report)) {
                is ReactorReportInboxResult.Enqueued -> Unit
                is ReactorReportInboxResult.Rejected -> {
                    return rejected(
                        ReactorWorldRuntimeRejection.QUEUE_FULL,
                        "reactor report queue rejected report ${report.reportId.value}: ${result.message}",
                    )
                }
            }
        }
        return ReactorWorldRuntimeResult.CommandsSubmitted(
            commandCount = commands.size,
            reportCount = reports.size,
        )
    }

    private fun enqueueCommand(command: ReactorCommand): ReactorWorldRuntimeResult =
        when (val result = commandOutbox.enqueue(command)) {
            is ReactorCommandOutboxResult.Enqueued -> ReactorWorldRuntimeResult.CommandQueued(command.commandId, result.size)
            is ReactorCommandOutboxResult.Rejected -> rejected(
                ReactorWorldRuntimeRejection.QUEUE_FULL,
                result.message,
            )
        }

    private fun enqueueCommands(commands: List<ReactorCommand>): ReactorWorldRuntimeResult {
        if (commands.size > commandOutbox.remainingCapacity) {
            return rejected(
                ReactorWorldRuntimeRejection.QUEUE_FULL,
                "reactor command queue cannot accept ${commands.size} commands; remaining capacity is ${commandOutbox.remainingCapacity}",
            )
        }
        var queueSize = commandOutbox.size
        for (command in commands) {
            when (val result = commandOutbox.enqueue(command)) {
                is ReactorCommandOutboxResult.Enqueued -> queueSize = result.size
                is ReactorCommandOutboxResult.Rejected -> {
                    return rejected(
                        ReactorWorldRuntimeRejection.QUEUE_FULL,
                        result.message,
                    )
                }
            }
        }
        return ReactorWorldRuntimeResult.BatchQueued(
            commandCount = commands.size,
            queueSize = queueSize,
        )
    }

    private fun nextCommandId(): ReactorCommandId =
        ReactorCommandId(nextCommandId++)

    private fun rejected(
        reason: ReactorWorldRuntimeRejection,
        message: String,
    ): ReactorWorldRuntimeResult.Rejected =
        ReactorWorldRuntimeResult.Rejected(reason, message)

    private fun catalogCheckpointStorageKey(contentVersion: Long): String =
        "catalog/checkpoint_${contentVersion.toString().padStart(12, '0')}.bin.zst"
}
