package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage

class ReactorNativeSession(
    private val structures: ReactorStructureStore,
    private val blobStorage: NativeBlobStorage,
) {
    private var nextReportId = 0L

    fun submit(commands: List<ReactorCommand>): List<ReactorReport> {
        val batchVersions = LinkedHashMap<ReactorStructureId, ReactorSnapshotVersion>()
        return buildList(commands.size) {
            for (command in commands) {
                add(execute(command, batchVersions))
            }
        }
    }

    private fun execute(
        command: ReactorCommand,
        batchVersions: MutableMap<ReactorStructureId, ReactorSnapshotVersion>,
    ): ReactorReport {
        val stored = structures.record(command.structureId)?.snapshotVersion
        if (stored == null) {
            return rejected(command, ReactorSnapshotVersion(0), "reactor structure ${command.structureId.value} is not registered")
        }
        val current = batchVersions[command.structureId] ?: stored
        val expected = command.expectedSnapshotVersion
        if (expected != null && expected != stored) {
            return rejected(
                command,
                current,
                "reactor command ${command.commandId.value} expected snapshot ${expected.value}, stored snapshot is ${stored.value}",
            )
        }

        return when (command) {
            is ReactorCommand.EnsureLoaded -> ensureLoaded(command, current)
            is ReactorCommand.UnloadAndExportSnapshot -> unloadAndExport(command, current, batchVersions)
            is ReactorCommand.RemoveReactor -> remove(command, current)
            is ReactorCommand.Tick -> tick(command, current, batchVersions)
            is ReactorCommand.InsertItem -> insertItem(command, current, batchVersions)
            is ReactorCommand.ExtractItem -> extractItem(command, current, batchVersions)
            is ReactorCommand.RequestSnapshot -> exportSnapshot(command, current, batchVersions)
            is ReactorCommand.RequestMetrics -> accepted(command, current)
        }
    }

    private fun ensureLoaded(
        command: ReactorCommand.EnsureLoaded,
        current: ReactorSnapshotVersion,
    ): ReactorReport {
        val record = structures.record(command.structureId)
            ?: return rejected(command, current, "reactor structure ${command.structureId.value} is not registered")
        if (record.state == ReactorStructureState.ACTIVE) {
            return accepted(command, current)
        }
        if (record.state != ReactorStructureState.SUSPENDED_UNLOADED) {
            return rejected(command, current, "reactor structure ${command.structureId.value} is ${record.state}")
        }
        if (command.checkpoint != null && record.reactorCheckpoint != command.checkpoint) {
            return rejected(command, current, "reactor structure ${command.structureId.value} checkpoint does not match the command")
        }
        return when (val result = structures.resumeFromCheckpoint(command.structureId, blobStorage)) {
            is ReactorOperationResult.ReactorResumed -> accepted(command, current)
            else -> rejected(command, current, result.message())
        }
    }

    private fun unloadAndExport(
        command: ReactorCommand.UnloadAndExportSnapshot,
        current: ReactorSnapshotVersion,
        batchVersions: MutableMap<ReactorStructureId, ReactorSnapshotVersion>,
    ): ReactorReport {
        val next = current.next()
        return when (val result = structures.suspendToCheckpoint(command.structureId, blobStorage, next.value)) {
            is ReactorOperationResult.ReactorSuspended -> {
                batchVersions[command.structureId] = next
                val record = structures.record(command.structureId)
                val checkpoint = record?.reactorCheckpoint
                    ?: return rejected(command, current, "reactor structure ${command.structureId.value} was suspended without checkpoint")
                ReactorReport.SnapshotReady(
                    reportId = nextReportId(),
                    commandId = command.commandId,
                    structureId = command.structureId,
                    snapshotVersion = next,
                    checkpoint = checkpoint,
                )
            }
            else -> rejected(command, current, result.message())
        }
    }

    private fun remove(
        command: ReactorCommand.RemoveReactor,
        current: ReactorSnapshotVersion,
    ): ReactorReport =
        when (val result = structures.remove(command.structureId)) {
            is ReactorOperationResult.Completed -> accepted(command, current)
            else -> rejected(command, current, result.message())
        }

    private fun tick(
        command: ReactorCommand.Tick,
        current: ReactorSnapshotVersion,
        batchVersions: MutableMap<ReactorStructureId, ReactorSnapshotVersion>,
    ): ReactorReport {
        val next = current.next()
        return when (val result = structures.tick(command.structureId, command.dtSeconds)) {
            is ReactorOperationResult.Completed -> {
                batchVersions[command.structureId] = next
                ReactorReport.TickCompleted(
                    reportId = nextReportId(),
                    commandId = command.commandId,
                    structureId = command.structureId,
                    snapshotVersion = next,
                    metrics = ReactorTickMetrics(
                        simulatedSeconds = command.dtSeconds,
                        temperatureKelvin = null,
                        pressurePascal = null,
                    ),
                )
            }
            else -> rejected(command, current, result.message())
        }
    }

    private fun insertItem(
        command: ReactorCommand.InsertItem,
        current: ReactorSnapshotVersion,
        batchVersions: MutableMap<ReactorStructureId, ReactorSnapshotVersion>,
    ): ReactorReport {
        val next = current.next()
        return when (val result = structures.insertItem(
            structureId = command.structureId,
            portPosition = command.portPosition,
            itemId = command.itemId,
            itemCount = command.itemCount,
        )) {
            is ReactorOperationResult.ItemInserted -> {
                batchVersions[command.structureId] = next
                ReactorReport.PortInputAccepted(
                    reportId = nextReportId(),
                    commandId = command.commandId,
                    structureId = command.structureId,
                    snapshotVersion = next,
                    portPosition = command.portPosition,
                    itemId = command.itemId,
                    acceptedCount = result.itemCount,
                )
            }
            else -> ReactorReport.PortInputRejected(
                reportId = nextReportId(),
                commandId = command.commandId,
                structureId = command.structureId,
                snapshotVersion = current,
                portPosition = command.portPosition,
                reason = result.message(),
            )
        }
    }

    private fun extractItem(
        command: ReactorCommand.ExtractItem,
        current: ReactorSnapshotVersion,
        batchVersions: MutableMap<ReactorStructureId, ReactorSnapshotVersion>,
    ): ReactorReport {
        val next = current.next()
        return when (val result = structures.extractItem(
            structureId = command.structureId,
            portPosition = command.portPosition,
            itemId = command.itemId,
            maxItemCount = command.maxItemCount,
            dtSeconds = command.dtSeconds,
        )) {
            is ReactorOperationResult.ItemExtracted -> {
                batchVersions[command.structureId] = next
                ReactorReport.PortOutputAccepted(
                    reportId = nextReportId(),
                    commandId = command.commandId,
                    structureId = command.structureId,
                    snapshotVersion = next,
                    portPosition = command.portPosition,
                    itemId = result.itemId,
                    extractedCount = result.itemCount,
                )
            }
            else -> ReactorReport.PortOutputRejected(
                reportId = nextReportId(),
                commandId = command.commandId,
                structureId = command.structureId,
                snapshotVersion = current,
                portPosition = command.portPosition,
                reason = result.message(),
            )
        }
    }

    private fun exportSnapshot(
        command: ReactorCommand.RequestSnapshot,
        current: ReactorSnapshotVersion,
        batchVersions: MutableMap<ReactorStructureId, ReactorSnapshotVersion>,
    ): ReactorReport {
        val next = current.next()
        return when (val result = structures.exportCheckpoint(command.structureId, blobStorage, next.value)) {
            is ReactorOperationResult.ReactorCheckpointExported -> {
                batchVersions[command.structureId] = next
                ReactorReport.SnapshotReady(
                    reportId = nextReportId(),
                    commandId = command.commandId,
                    structureId = command.structureId,
                    snapshotVersion = next,
                    checkpoint = result.checkpoint,
                )
            }
            else -> rejected(command, current, result.message())
        }
    }

    private fun accepted(
        command: ReactorCommand,
        snapshotVersion: ReactorSnapshotVersion,
    ): ReactorReport.CommandAccepted =
        ReactorReport.CommandAccepted(
            reportId = nextReportId(),
            commandId = command.commandId,
            structureId = command.structureId,
            snapshotVersion = snapshotVersion,
        )

    private fun rejected(
        command: ReactorCommand,
        snapshotVersion: ReactorSnapshotVersion,
        reason: String,
    ): ReactorReport.CommandRejected =
        ReactorReport.CommandRejected(
            reportId = nextReportId(),
            commandId = command.commandId,
            structureId = command.structureId,
            snapshotVersion = snapshotVersion,
            reason = reason,
        )

    private fun nextReportId(): ReactorReportId =
        ReactorReportId(nextReportId++)

    private fun ReactorSnapshotVersion.next(): ReactorSnapshotVersion =
        ReactorSnapshotVersion(value + 1)

    private fun ReactorOperationResult.message(): String =
        when (this) {
            ReactorOperationResult.Completed -> "operation completed"
            is ReactorOperationResult.ItemInserted -> "inserted $itemCount items"
            is ReactorOperationResult.ItemExtracted -> "extracted $itemCount of $itemId"
            is ReactorOperationResult.ReactorSuspended -> message
            is ReactorOperationResult.ReactorResumed -> message
            is ReactorOperationResult.ReactorCheckpointExported -> "reactor checkpoint exported to ${checkpoint.storageKey}"
            is ReactorOperationResult.Rejected -> message
            is ReactorOperationResult.Failed -> message
        }

}
