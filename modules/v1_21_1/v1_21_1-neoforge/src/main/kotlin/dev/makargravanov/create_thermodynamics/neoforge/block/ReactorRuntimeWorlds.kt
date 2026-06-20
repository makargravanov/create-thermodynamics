package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.common.chemistry.binding.DefaultItemChemicalBindings
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorAppliedReport
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorNativeSession
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorPortInputTransfer
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorPortOutputTransfer
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorPortTransferCycle
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorReport
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorWorldRuntimeResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorWorldRuntime
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage
import net.minecraft.core.BlockPos
import net.minecraft.server.level.ServerLevel
import net.minecraft.world.level.storage.LevelResource
import java.util.WeakHashMap

object ReactorRuntimeWorlds {
    private const val TRANSFER_INTERVAL_TICKS = 20L
    private const val TRANSFER_DT_SECONDS = 1.0
    private const val MAX_SUBMITTED_COMMANDS = 256
    private const val MAX_APPLIED_REPORTS = 256
    private const val MAX_OUTPUT_ITEMS_PER_PORT = 64 * 26

    private val knownBoundItemIds = DefaultItemChemicalBindings.bindings
        .mapTo(linkedSetOf()) { it.itemId }
    private val contextsByLevel = WeakHashMap<ServerLevel, ReactorRuntimeWorldContext>()

    fun runtime(level: ServerLevel): ReactorWorldRuntime =
        context(level).runtime

    fun tickLevel(level: ServerLevel) {
        val context = context(level)
        applyReactorReports(level, context.runtime, context.runtime.applyReadyReports(MAX_APPLIED_REPORTS))
        context.tickIndex += 1
        if (context.tickIndex % TRANSFER_INTERVAL_TICKS != 0L) {
            return
        }

        queuePortTransfers(level, context.runtime)
        context.runtime.queueTickForActiveStructures(TRANSFER_DT_SECONDS)
        context.runtime.submitQueuedCommands(context.nativeSession, MAX_SUBMITTED_COMMANDS)
        applyReactorReports(level, context.runtime, context.runtime.applyReadyReports(MAX_APPLIED_REPORTS))
    }

    private fun context(level: ServerLevel): ReactorRuntimeWorldContext =
        contextsByLevel.getOrPut(level) {
            val blobStorage = NativeBlobStorage(
                level.server.getWorldPath(LevelResource.ROOT)
                    .resolve("create_thermodynamics")
                    .resolve("native"),
            )
            val runtime = ReactorWorldRuntime(blobStorage = blobStorage)
            ReactorRuntimeWorldContext(
                runtime = runtime,
                nativeSession = ReactorNativeSession(runtime.structures, blobStorage),
            )
        }

    private fun queuePortTransfers(level: ServerLevel, runtime: ReactorWorldRuntime) {
        for (record in runtime.structures.activeRecords()) {
            val outputs = record.definition.ports
                .filter { it.kind == ReactorPortKind.ITEM_OUTPUT }
                .flatMap { port ->
                    val blockEntity = reactorBlockEntity(level, port.position) ?: return@flatMap emptyList()
                    buildList {
                        for (itemId in knownBoundItemIds) {
                            if (!blockEntity.portFilterAllowsItemId(itemId)) {
                                continue
                            }
                            val maxItemCount = blockEntity.insertablePortOutputCount(itemId, MAX_OUTPUT_ITEMS_PER_PORT)
                            if (maxItemCount > 0) {
                                add(
                                    ReactorPortOutputTransfer(
                                        portPosition = port.position,
                                        itemId = itemId,
                                        maxItemCount = maxItemCount,
                                        dtSeconds = TRANSFER_DT_SECONDS,
                                    ),
                                )
                            }
                        }
                    }
                }
            val inputs = record.definition.ports
                .filter { it.kind == ReactorPortKind.ITEM_INPUT }
                .mapNotNull { port ->
                    val blockEntity = reactorBlockEntity(level, port.position) ?: return@mapNotNull null
                    val stack = blockEntity.firstPortInputStack() ?: return@mapNotNull null
                    if (stack.itemId !in knownBoundItemIds) {
                        return@mapNotNull null
                    }
                    ReactorPortInputTransfer(
                        portPosition = port.position,
                        itemId = stack.itemId,
                        itemCount = stack.count,
                    )
                }
            if (outputs.isNotEmpty() || inputs.isNotEmpty()) {
                runtime.queuePortTransferCycle(
                    ReactorPortTransferCycle(
                        structureId = record.structureId,
                        outputs = outputs,
                        inputs = inputs,
                    ),
                )
            }
        }
    }

    private fun applyReactorReports(
        level: ServerLevel,
        runtime: ReactorWorldRuntime,
        result: ReactorWorldRuntimeResult,
    ) {
        if (result !is ReactorWorldRuntimeResult.ReportsApplied) {
            return
        }
        for (applied in result.reports) {
            applyReactorReport(level, runtime, applied)
        }
    }

    private fun applyReactorReport(
        level: ServerLevel,
        runtime: ReactorWorldRuntime,
        applied: ReactorAppliedReport,
    ) {
        when (val report = applied.report) {
            is ReactorReport.TickCompleted -> {
                val controllerPosition = runtime.structures.record(report.structureId)
                    ?.definition
                    ?.controllerPosition
                    ?: return
                reactorBlockEntity(level, controllerPosition)
                    ?.applyNativeMetrics(report.metrics)
            }
            is ReactorReport.PortInputAccepted -> {
                reactorBlockEntity(level, report.portPosition)
                    ?.removeConfirmedPortInput(report.itemId, report.acceptedCount)
            }
            is ReactorReport.PortOutputAccepted -> {
                if (report.extractedCount > 0) {
                    reactorBlockEntity(level, report.portPosition)
                        ?.insertConfirmedPortOutput(report.itemId, report.extractedCount)
                }
            }
            else -> Unit
        }
    }

    private fun reactorBlockEntity(
        level: ServerLevel,
        position: ReactorBlockPosition,
    ): ReactorMultiblockBlockEntity? {
        val pos = BlockPos(position.x, position.y, position.z)
        if (!level.isLoaded(pos)) {
            return null
        }
        return level.getBlockEntity(pos) as? ReactorMultiblockBlockEntity
    }
}

private data class ReactorRuntimeWorldContext(
    val runtime: ReactorWorldRuntime,
    val nativeSession: ReactorNativeSession,
    var tickIndex: Long = 0L,
)
