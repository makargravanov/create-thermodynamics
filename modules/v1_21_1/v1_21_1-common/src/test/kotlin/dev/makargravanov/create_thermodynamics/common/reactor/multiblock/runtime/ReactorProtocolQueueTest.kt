package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import java.util.UUID
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertIs
import kotlin.test.assertTrue

class ReactorProtocolQueueTest {
    @Test
    fun `command outbox preserves order and rejects overflow`() {
        val outbox = ReactorCommandOutbox(
            ReactorCommandOutboxLimits(maxCommands = 2, maxDrainBatch = 8),
        )
        val first = tickCommand(1, snapshotVersion = 10, dtSeconds = 0.05)
        val second = tickCommand(2, snapshotVersion = 11, dtSeconds = 0.05)
        val third = tickCommand(3, snapshotVersion = 12, dtSeconds = 0.05)

        assertIs<ReactorCommandOutboxResult.Enqueued>(outbox.enqueue(first))
        assertIs<ReactorCommandOutboxResult.Enqueued>(outbox.enqueue(second))
        val rejected = assertIs<ReactorCommandOutboxResult.Rejected>(outbox.enqueue(third))

        assertEquals(ReactorQueueRejection.QUEUE_FULL, rejected.reason)
        assertEquals(listOf(first, second), outbox.drain(maxCommands = 10))
        assertTrue(outbox.isEmpty())
    }

    @Test
    fun `command outbox drains at configured batch limit`() {
        val outbox = ReactorCommandOutbox(
            ReactorCommandOutboxLimits(maxCommands = 4, maxDrainBatch = 2),
        )
        val first = tickCommand(1, snapshotVersion = 1, dtSeconds = 0.05)
        val second = tickCommand(2, snapshotVersion = 2, dtSeconds = 0.05)
        val third = tickCommand(3, snapshotVersion = 3, dtSeconds = 0.05)
        outbox.enqueue(first)
        outbox.enqueue(second)
        outbox.enqueue(third)

        assertEquals(listOf(first, second), outbox.drain(maxCommands = 4))
        assertEquals(listOf(third), outbox.drain(maxCommands = 4))
    }

    @Test
    fun `report inbox coalesces tick metrics for the same structure`() {
        val inbox = ReactorReportInbox(
            ReactorReportInboxLimits(maxReports = 1, maxDrainBatch = 8),
        )
        val first = tickReport(1, snapshotVersion = 1, simulatedSeconds = 0.05)
        val second = tickReport(2, snapshotVersion = 2, simulatedSeconds = 0.10)

        assertIs<ReactorReportInboxResult.Enqueued>(inbox.enqueue(first))
        assertIs<ReactorReportInboxResult.Enqueued>(inbox.enqueue(second))

        assertEquals(listOf(second), inbox.drain(maxReports = 8))
    }

    @Test
    fun `report inbox rejects required report when full without discarding existing report`() {
        val inbox = ReactorReportInbox(
            ReactorReportInboxLimits(maxReports = 1, maxDrainBatch = 8),
        )
        val accepted = ReactorReport.CommandAccepted(
            reportId = ReactorReportId(1),
            commandId = ReactorCommandId(1),
            structureId = STRUCTURE_ID,
            snapshotVersion = ReactorSnapshotVersion(1),
        )
        val rejectedReport = ReactorReport.CommandRejected(
            reportId = ReactorReportId(2),
            commandId = ReactorCommandId(2),
            structureId = STRUCTURE_ID,
            snapshotVersion = ReactorSnapshotVersion(1),
            reason = "queue pressure test",
        )

        assertIs<ReactorReportInboxResult.Enqueued>(inbox.enqueue(accepted))
        val rejected = assertIs<ReactorReportInboxResult.Rejected>(inbox.enqueue(rejectedReport))

        assertEquals(ReactorQueueRejection.QUEUE_FULL, rejected.reason)
        assertEquals(listOf(accepted), inbox.drain(maxReports = 8))
    }

    @Test
    fun `report inbox does not coalesce metrics across different structures`() {
        val inbox = ReactorReportInbox(
            ReactorReportInboxLimits(maxReports = 2, maxDrainBatch = 8),
        )
        val first = tickReport(1, snapshotVersion = 1, simulatedSeconds = 0.05)
        val second = tickReport(
            reportId = 2,
            structureId = ReactorStructureId(UUID.fromString("85e18c43-8dec-4485-a8f3-f4657959bcc8")),
            snapshotVersion = 1,
            simulatedSeconds = 0.05,
        )

        assertIs<ReactorReportInboxResult.Enqueued>(inbox.enqueue(first))
        assertIs<ReactorReportInboxResult.Enqueued>(inbox.enqueue(second))

        assertEquals(listOf(first, second), inbox.drain(maxReports = 8))
    }

    @Test
    fun `protocol ids and payloads reject invalid values`() {
        assertFailsWith<IllegalArgumentException> { ReactorCommandId(-1) }
        assertFailsWith<IllegalArgumentException> { ReactorReportId(-1) }
        assertFailsWith<IllegalArgumentException> { ReactorSnapshotVersion(-1) }
        assertFailsWith<IllegalArgumentException> {
            tickCommand(1, snapshotVersion = 1, dtSeconds = Double.NaN)
        }
        assertFailsWith<IllegalArgumentException> {
            ReactorCommand.InsertItem(
                commandId = ReactorCommandId(1),
                structureId = STRUCTURE_ID,
                expectedSnapshotVersion = ReactorSnapshotVersion(1),
                portPosition = ReactorBlockPosition(0, 0, 0),
                itemId = " ",
                itemCount = 1,
            )
        }
    }

    private fun tickCommand(
        commandId: Long,
        snapshotVersion: Long,
        dtSeconds: Double,
    ): ReactorCommand.Tick =
        ReactorCommand.Tick(
            commandId = ReactorCommandId(commandId),
            structureId = STRUCTURE_ID,
            expectedSnapshotVersion = ReactorSnapshotVersion(snapshotVersion),
            dtSeconds = dtSeconds,
        )

    private fun tickReport(
        reportId: Long,
        structureId: ReactorStructureId = STRUCTURE_ID,
        snapshotVersion: Long,
        simulatedSeconds: Double,
    ): ReactorReport.TickCompleted =
        ReactorReport.TickCompleted(
            reportId = ReactorReportId(reportId),
            commandId = null,
            structureId = structureId,
            snapshotVersion = ReactorSnapshotVersion(snapshotVersion),
            metrics = ReactorTickMetrics(
                simulatedSeconds = simulatedSeconds,
                temperatureKelvin = 300.0,
                pressurePascal = 101_325.0,
            ),
        )

    private companion object {
        val STRUCTURE_ID: ReactorStructureId =
            ReactorStructureId(UUID.fromString("b8aefdf9-9ddc-46d0-a212-615c4b8e7ddb"))
    }
}
