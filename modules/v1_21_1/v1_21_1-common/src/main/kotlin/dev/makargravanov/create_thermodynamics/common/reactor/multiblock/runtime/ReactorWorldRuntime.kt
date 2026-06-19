package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationRejection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobKind
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobRef
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage

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

    data class CatalogUpdated(
        val catalog: ReactorCatalogRuntimeState,
    ) : ReactorWorldRuntimeResult

    data class Rejected(
        val reason: ReactorWorldRuntimeRejection,
        val message: String,
    ) : ReactorWorldRuntimeResult
}

class ReactorWorldRuntime(
    val structures: ReactorStructureStore = ReactorStructureStore(),
    val blobStorage: NativeBlobStorage,
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

    private fun enqueueCommand(command: ReactorCommand): ReactorWorldRuntimeResult =
        when (val result = commandOutbox.enqueue(command)) {
            is ReactorCommandOutboxResult.Enqueued -> ReactorWorldRuntimeResult.CommandQueued(command.commandId, result.size)
            is ReactorCommandOutboxResult.Rejected -> rejected(
                ReactorWorldRuntimeRejection.QUEUE_FULL,
                result.message,
            )
        }

    private fun nextCommandId(): ReactorCommandId =
        ReactorCommandId(nextCommandId++)

    private fun rejected(
        reason: ReactorWorldRuntimeRejection,
        message: String,
    ): ReactorWorldRuntimeResult.Rejected =
        ReactorWorldRuntimeResult.Rejected(reason, message)
}
