package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

data class ReactorReportInboxLimits(
    val maxReports: Int = 1024,
    val maxDrainBatch: Int = 256,
) {
    init {
        require(maxReports > 0) { "maxReports must be positive" }
        require(maxDrainBatch > 0) { "maxDrainBatch must be positive" }
    }
}

sealed interface ReactorReportInboxResult {
    data class Enqueued(val size: Int) : ReactorReportInboxResult
    data class Rejected(val reason: ReactorQueueRejection, val message: String) : ReactorReportInboxResult
}

class ReactorReportInbox(
    private val limits: ReactorReportInboxLimits = ReactorReportInboxLimits(),
) {
    private val reports = ArrayDeque<ReactorReport>()

    val size: Int
        get() = reports.size

    fun isEmpty(): Boolean =
        reports.isEmpty()

    fun enqueue(report: ReactorReport): ReactorReportInboxResult {
        coalesceReport(report)
        if (reports.size >= limits.maxReports) {
            return ReactorReportInboxResult.Rejected(
                ReactorQueueRejection.QUEUE_FULL,
                "reactor report queue is full: ${reports.size}/${limits.maxReports}",
            )
        }
        reports.addLast(report)
        return ReactorReportInboxResult.Enqueued(reports.size)
    }

    fun drain(maxReports: Int = limits.maxDrainBatch): List<ReactorReport> {
        require(maxReports > 0) { "drain maxReports must be positive" }
        val count = minOf(maxReports, limits.maxDrainBatch, reports.size)
        return buildList(count) {
            repeat(count) {
                add(reports.removeFirst())
            }
        }
    }

    private fun coalesceReport(report: ReactorReport) {
        if (report !is ReactorReport.TickCompleted) {
            return
        }

        val kept = ArrayDeque<ReactorReport>(reports.size)
        while (reports.isNotEmpty()) {
            val existing = reports.removeFirst()
            if (existing !is ReactorReport.TickCompleted || existing.structureId != report.structureId) {
                kept.addLast(existing)
            }
        }
        reports.addAll(kept)
    }
}
