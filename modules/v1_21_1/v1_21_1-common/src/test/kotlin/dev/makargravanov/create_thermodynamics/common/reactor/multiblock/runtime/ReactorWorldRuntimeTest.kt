package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationRejection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access.ReactorOperationResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockDirection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorZoneDescriptor
import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobHash
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobKind
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobRef
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage
import java.nio.file.Files
import java.util.UUID
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs

class ReactorWorldRuntimeTest {
    @Test
    fun `queued tick uses current snapshot version from structure record`() {
        val runtime = testRuntime()
        val definition = testDefinition()
        runtime.registerStructure(definition)
        runtime.receiveReport(
            tickReport(
                reportId = 1,
                structureId = definition.structureId,
                snapshotVersion = 4,
            ),
        )
        assertIs<ReactorWorldRuntimeResult.ReportsApplied>(runtime.applyReadyReports(8))

        val queued = assertIs<ReactorWorldRuntimeResult.CommandQueued>(
            runtime.queueTick(definition.structureId, 0.05),
        )
        val command = assertIs<ReactorCommand.Tick>(runtime.drainCommands(8).single())

        assertEquals(ReactorCommandId(0), queued.commandId)
        assertEquals(ReactorSnapshotVersion(4), command.expectedSnapshotVersion)
        assertEquals(0.05, command.dtSeconds)
    }

    @Test
    fun `stale report is rejected and does not roll snapshot version back`() {
        val runtime = testRuntime()
        val definition = testDefinition()
        runtime.registerStructure(definition)
        runtime.receiveReport(tickReport(1, definition.structureId, snapshotVersion = 4))
        runtime.applyReadyReports(8)
        runtime.receiveReport(tickReport(2, definition.structureId, snapshotVersion = 3))

        val applied = assertIs<ReactorWorldRuntimeResult.ReportsApplied>(runtime.applyReadyReports(8))
        val rejected = assertIs<ReactorOperationResult.Rejected>(applied.results.single())

        assertEquals(ReactorOperationRejection.STALE_REPORT, rejected.reason)
        assertEquals(ReactorSnapshotVersion(4), runtime.structures.record(definition.structureId)?.snapshotVersion)
    }

    @Test
    fun `snapshot ready report stores reactor checkpoint reference`() {
        val runtime = testRuntime()
        val definition = testDefinition()
        runtime.registerStructure(definition)
        val checkpoint = blobRef(NativeBlobKind.ReactorSnapshot, contentVersion = 9)
        runtime.receiveReport(
            ReactorReport.SnapshotReady(
                reportId = ReactorReportId(1),
                commandId = null,
                structureId = definition.structureId,
                snapshotVersion = ReactorSnapshotVersion(9),
                checkpoint = checkpoint,
            ),
        )

        val applied = assertIs<ReactorWorldRuntimeResult.ReportsApplied>(runtime.applyReadyReports(8))

        assertEquals(listOf(ReactorOperationResult.Completed), applied.results)
        val record = runtime.structures.record(definition.structureId)
        assertEquals(ReactorSnapshotVersion(9), record?.snapshotVersion)
        assertEquals(checkpoint, record?.reactorCheckpoint)
    }

    @Test
    fun `catalog checkpoint accepts only dynamic catalog checkpoint blobs`() {
        val runtime = testRuntime()
        val accepted = assertIs<ReactorWorldRuntimeResult.CatalogUpdated>(
            runtime.installCatalogCheckpoint(blobRef(NativeBlobKind.DynamicCatalogCheckpoint, contentVersion = 12)),
        )
        val rejected = assertIs<ReactorWorldRuntimeResult.Rejected>(
            runtime.installCatalogCheckpoint(blobRef(NativeBlobKind.ReactorSnapshot, contentVersion = 13)),
        )

        assertEquals(12, accepted.catalog.catalogVersion)
        assertEquals(12, runtime.catalog.catalogVersion)
        assertEquals(ReactorWorldRuntimeRejection.WRONG_BLOB_KIND, rejected.reason)
        assertEquals(12, runtime.catalog.catalogVersion)
    }

    @Test
    fun `command queue overflow is reported by world runtime`() {
        val runtime = testRuntime(
            commandOutbox = ReactorCommandOutbox(ReactorCommandOutboxLimits(maxCommands = 1, maxDrainBatch = 8)),
        )
        val definition = testDefinition()
        runtime.registerStructure(definition)

        assertIs<ReactorWorldRuntimeResult.CommandQueued>(runtime.queueTick(definition.structureId, 0.05))
        val rejected = assertIs<ReactorWorldRuntimeResult.Rejected>(runtime.queueSnapshotRequest(definition.structureId, "test"))

        assertEquals(ReactorWorldRuntimeRejection.QUEUE_FULL, rejected.reason)
    }

    @Test
    fun `queued tick can be submitted through native session and applied as report`() {
        val bridge = FakeNativeReactorBridge()
        val runtime = testRuntime(nativeBridge = bridge)
        val definition = testDefinition()
        runtime.registerStructure(definition)
        assertIs<ReactorWorldRuntimeResult.CommandQueued>(runtime.queueTick(definition.structureId, 0.05))

        val submitted = assertIs<ReactorWorldRuntimeResult.CommandsSubmitted>(
            runtime.submitQueuedCommands(ReactorNativeSession(runtime.structures, runtime.blobStorage), maxCommands = 8),
        )
        val applied = assertIs<ReactorWorldRuntimeResult.ReportsApplied>(runtime.applyReadyReports(8))

        assertEquals(1, submitted.commandCount)
        assertEquals(1, submitted.reportCount)
        assertEquals(listOf(ReactorOperationResult.Completed), applied.results)
        assertEquals(1, bridge.tickCount)
        assertEquals(ReactorSnapshotVersion(1), runtime.structures.record(definition.structureId)?.snapshotVersion)
    }

    @Test
    fun `snapshot command exports checkpoint without suspending active reactor`() {
        val runtime = testRuntime()
        val definition = testDefinition()
        runtime.registerStructure(definition)
        assertIs<ReactorWorldRuntimeResult.CommandQueued>(
            runtime.queueSnapshotRequest(definition.structureId, "test snapshot"),
        )

        assertIs<ReactorWorldRuntimeResult.CommandsSubmitted>(
            runtime.submitQueuedCommands(ReactorNativeSession(runtime.structures, runtime.blobStorage), maxCommands = 8),
        )
        assertIs<ReactorWorldRuntimeResult.ReportsApplied>(runtime.applyReadyReports(8))

        val record = runtime.structures.record(definition.structureId)
        assertEquals(ReactorStructureState.ACTIVE, record?.state)
        assertEquals(ReactorSnapshotVersion(1), record?.snapshotVersion)
        assertEquals(NativeBlobKind.ReactorSnapshot, record?.reactorCheckpoint?.kind)
    }

    @Test
    fun `world runtime does not execute commands when report inbox has no capacity`() {
        val bridge = FakeNativeReactorBridge()
        val runtime = testRuntime(
            nativeBridge = bridge,
            reportInbox = ReactorReportInbox(ReactorReportInboxLimits(maxReports = 1, maxDrainBatch = 8)),
        )
        val definition = testDefinition()
        runtime.registerStructure(definition)
        assertIs<ReactorWorldRuntimeResult.ReportQueued>(runtime.receiveReport(tickReport(1, definition.structureId, snapshotVersion = 1)))
        assertIs<ReactorWorldRuntimeResult.CommandQueued>(runtime.queueTick(definition.structureId, 0.05))

        val rejected = assertIs<ReactorWorldRuntimeResult.Rejected>(
            runtime.submitQueuedCommands(ReactorNativeSession(runtime.structures, runtime.blobStorage), maxCommands = 8),
        )

        assertEquals(ReactorWorldRuntimeRejection.QUEUE_FULL, rejected.reason)
        assertEquals(0, bridge.tickCount)
        assertIs<ReactorCommand.Tick>(runtime.drainCommands(8).single())
    }

    @Test
    fun `native session rejects second command when batch snapshot version has advanced`() {
        val bridge = FakeNativeReactorBridge()
        val runtime = testRuntime(nativeBridge = bridge)
        val definition = testDefinition()
        runtime.registerStructure(definition)
        val session = ReactorNativeSession(runtime.structures, runtime.blobStorage)
        val first = ReactorCommand.Tick(
            commandId = ReactorCommandId(1),
            structureId = definition.structureId,
            expectedSnapshotVersion = ReactorSnapshotVersion(0),
            dtSeconds = 0.05,
        )
        val second = first.copy(commandId = ReactorCommandId(2))

        val reports = session.submit(listOf(first, second))

        assertIs<ReactorReport.TickCompleted>(reports[0])
        val rejected = assertIs<ReactorReport.CommandRejected>(reports[1])
        assertEquals(ReactorSnapshotVersion(1), rejected.snapshotVersion)
        assertEquals(1, bridge.tickCount)
    }

    private fun testRuntime(
        commandOutbox: ReactorCommandOutbox = ReactorCommandOutbox(),
        reportInbox: ReactorReportInbox = ReactorReportInbox(),
        nativeBridge: NativeReactorBridge = FakeNativeReactorBridge(),
    ): ReactorWorldRuntime =
        ReactorWorldRuntime(
            structures = ReactorStructureStore(nativeBridge),
            blobStorage = NativeBlobStorage(Files.createTempDirectory("ct-reactor-world-runtime-test")),
            commandOutbox = commandOutbox,
            reportInbox = reportInbox,
        )

    private fun tickReport(
        reportId: Long,
        structureId: ReactorStructureId,
        snapshotVersion: Long,
    ): ReactorReport.TickCompleted =
        ReactorReport.TickCompleted(
            reportId = ReactorReportId(reportId),
            commandId = null,
            structureId = structureId,
            snapshotVersion = ReactorSnapshotVersion(snapshotVersion),
            metrics = ReactorTickMetrics(
                simulatedSeconds = 0.05,
                temperatureKelvin = 300.0,
                pressurePascal = 101_325.0,
            ),
        )

    private fun testDefinition(): ReactorMultiblockDefinition {
        val structureId = ReactorStructureId(UUID.fromString("50352d37-1e3d-45c1-beb4-885f5bd83ba1"))
        val chamber = ReactorBlockPosition(0, 0, 0)
        val itemInput = ReactorBlockPosition(1, 0, 0)
        return ReactorMultiblockDefinition(
            structureId = structureId,
            controllerPosition = ReactorBlockPosition(-1, 0, 0),
            controllerContactDirection = ReactorBlockDirection.EAST,
            zone = ReactorZoneDescriptor(
                zoneIndex = 0,
                volumePositions = setOf(chamber),
                plainChamberPositions = setOf(chamber),
                volumeCubicMeters = 1.0,
            ),
            ports = listOf(
                ReactorPortDescriptor(
                    portIndex = 0,
                    kind = ReactorPortKind.ITEM_INPUT,
                    position = itemInput,
                    zoneIndex = 0,
                    attachedChamberPosition = chamber,
                    contactDirection = ReactorBlockDirection.WEST,
                ),
            ),
        )
    }

    private fun blobRef(
        kind: NativeBlobKind,
        contentVersion: Long,
    ): NativeBlobRef =
        NativeBlobRef(
            kind = kind,
            schemaVersion = 1,
            modelVersion = "test-model",
            contentVersion = contentVersion,
            payloadHash = ZERO_HASH,
            encodedHash = ZERO_HASH,
            encodedByteSize = 64,
            compressedPayloadByteSize = 16,
            uncompressedByteSize = 16,
            storageKey = "test/${kind.name.lowercase()}_$contentVersion.bin.zst",
        )

    private class FakeNativeReactorBridge : NativeReactorBridge {
        var tickCount: Int = 0
            private set

        override fun createNativeReactor(definition: ReactorMultiblockDefinition): NativeReactorMultiblockBinding =
            NativeReactorMultiblockBinding(
                structureId = definition.structureId,
                nativeReactorId = ThermodynamicsNative.NativeReactorId(1),
                itemInputs = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).mapIndexed { index, port ->
                    NativeReactorPortBinding(port, index)
                },
                itemOutputs = emptyList(),
                fluidInputs = emptyList(),
                fluidOutputs = emptyList(),
            )

        override fun createNativeReactorFromCheckpoint(
            definition: ReactorMultiblockDefinition,
            encodedCheckpoint: ByteArray,
        ): NativeReactorMultiblockBinding =
            createNativeReactor(definition)

        override fun removeNativeReactor(binding: NativeReactorMultiblockBinding) = Unit

        override fun tickNativeReactor(binding: NativeReactorMultiblockBinding, dtSeconds: Double) {
            tickCount += 1
        }

        override fun exportReactorCheckpoint(
            binding: NativeReactorMultiblockBinding,
            contentVersion: Long,
        ): ByteArray =
            nativeBlobBytes(
                kind = NativeBlobKind.ReactorSnapshot,
                contentVersion = contentVersion,
            )

        override fun insertItemStack(
            binding: NativeReactorMultiblockBinding,
            itemInputPort: ReactorPortDescriptor,
            itemId: String,
            itemCount: Int,
        ): Double =
            itemCount.toDouble()
    }

    private companion object {
        val ZERO_HASH = NativeBlobHash("0".repeat(64))

        fun nativeBlobBytes(
            kind: NativeBlobKind,
            contentVersion: Long = 1,
            modelVersion: String = "test-model",
            uncompressedByteSize: Long = 16,
            compressedPayload: ByteArray = byteArrayOf(1),
        ): ByteArray =
            buildList<Int> {
                addAll(byteArrayOf(0x43, 0x54, 0x4e, 0x42, 0x4c, 0x42, 0x31, 0x00).map { it.toInt() and 0xff })
                addU16(1)
                addU16(kind.wireId)
                addU16(modelVersion.encodeToByteArray().size)
                addU64(contentVersion)
                addU64(uncompressedByteSize)
                addU64(compressedPayload.size.toLong())
                repeat(32) { add(0) }
                addAll(modelVersion.encodeToByteArray().map { it.toInt() and 0xff })
                addAll(compressedPayload.map { it.toInt() and 0xff })
            }.map { it.toByte() }.toByteArray()

        private fun MutableList<Int>.addU16(value: Int) {
            add(value and 0xff)
            add((value ushr 8) and 0xff)
        }

        private fun MutableList<Int>.addU64(value: Long) {
            for (index in 0 until 8) {
                add(((value ushr (index * 8)) and 0xff).toInt())
            }
        }
    }
}
