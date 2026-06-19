package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

data class ReactorRuntimeTickSchedulerConfig(
    val simulationIntervalTicks: Int = 20,
    val dtSeconds: Double = 1.0,
    val maxReportsAppliedPerTick: Int = 128,
    val maxCommandsSubmittedPerTick: Int = 128,
) {
    init {
        require(simulationIntervalTicks > 0) { "simulationIntervalTicks must be positive" }
        require(dtSeconds.isFinite() && dtSeconds >= 0.0) {
            "dtSeconds must be non-negative and finite"
        }
        require(maxReportsAppliedPerTick > 0) { "maxReportsAppliedPerTick must be positive" }
        require(maxCommandsSubmittedPerTick > 0) { "maxCommandsSubmittedPerTick must be positive" }
    }
}

data class ReactorRuntimeTickCycle(
    val tickIndex: Long,
    val reportsBeforeSubmit: ReactorWorldRuntimeResult,
    val tickQueueResult: ReactorWorldRuntimeResult?,
    val commandSubmitResult: ReactorWorldRuntimeResult,
    val reportsAfterSubmit: ReactorWorldRuntimeResult,
)

class ReactorRuntimeTickScheduler(
    private val runtime: ReactorWorldRuntime,
    private val nativeSession: ReactorNativeSession,
    private val config: ReactorRuntimeTickSchedulerConfig = ReactorRuntimeTickSchedulerConfig(),
) {
    private var tickIndex: Long = 0

    fun tick(): ReactorRuntimeTickCycle {
        tickIndex += 1

        val reportsBeforeSubmit = runtime.applyReadyReports(config.maxReportsAppliedPerTick)
        val tickQueueResult = if (tickIndex % config.simulationIntervalTicks == 0L) {
            runtime.queueTickForActiveStructures(config.dtSeconds)
        } else {
            null
        }
        val commandSubmitResult = runtime.submitQueuedCommands(
            nativeSession = nativeSession,
            maxCommands = config.maxCommandsSubmittedPerTick,
        )
        val reportsAfterSubmit = runtime.applyReadyReports(config.maxReportsAppliedPerTick)

        return ReactorRuntimeTickCycle(
            tickIndex = tickIndex,
            reportsBeforeSubmit = reportsBeforeSubmit,
            tickQueueResult = tickQueueResult,
            commandSubmitResult = commandSubmitResult,
            reportsAfterSubmit = reportsAfterSubmit,
        )
    }
}
