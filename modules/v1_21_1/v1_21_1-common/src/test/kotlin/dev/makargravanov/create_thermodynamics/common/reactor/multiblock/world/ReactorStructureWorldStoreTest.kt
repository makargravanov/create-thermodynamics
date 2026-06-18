package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorAssemblyResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockAssembler
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockRules
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockDirection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlock
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import java.util.UUID
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertNotNull
import kotlin.test.assertNull

class ReactorStructureWorldStoreTest {
    private val rules = ReactorMultiblockRules(chamberVolumeCubicMeters = 1.0)
    private val planner = ReactorAssemblyPlanner(rules)

    @Test
    fun `assembler returns rejected result instead of using exception as normal path`() {
        val assembler = ReactorMultiblockAssembler(rules)

        val result = assembler.tryAssemble(
            structureId = ReactorStructureId(UUID.fromString("c962f353-6fe0-44d2-ad22-1fe4674d66d7")),
            blocks = listOf(
                ReactorMultiblockBlock(pos(0, 0, 0), ReactorMultiblockBlockKind.CONTROLLER),
            ),
        )

        val rejected = assertIs<ReactorAssemblyResult.Rejected>(result)
        assertEquals(true, rejected.errors.any { it.contains("chamber blocks") })
    }

    @Test
    fun `formed controller view does not depend on controller being active volume block`() {
        val store = ReactorStructureWorldStore()
        val scan = scanOf(
            snapshot(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
            snapshot(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )
        val plan = planner.buildPlan(scan)

        store.applyPlan(
            scannedPositions = scan.loadedReactorPositions,
            plan = plan,
            removeMissingStructures = true,
        )

        val controller = pos(-1, 0, 0)
        val viewState = store.controllerViewState(controller)
        assertEquals(ReactorControllerFormationState.FORMED, viewState.formationState)
        assertEquals(1, viewState.zoneCount)
        assertEquals(1, viewState.chamberBlockCount)
        assertEquals(false, store.membershipAt(controller)?.activeVolumeBlock)
    }

    @Test
    fun `rejected controller keeps diagnostics in view state`() {
        val store = ReactorStructureWorldStore()
        val scan = scanOf(
            snapshot(0, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )
        val plan = planner.buildPlan(scan)

        store.applyPlan(
            scannedPositions = scan.loadedReactorPositions,
            plan = plan,
            removeMissingStructures = true,
        )

        val viewState = store.controllerViewState(pos(0, 0, 0))
        assertEquals(ReactorControllerFormationState.NOT_FORMED, viewState.formationState)
        assertNotNull(viewState.diagnostic)
        assertEquals(true, viewState.diagnostic!!.contains("chamber blocks"))
    }

    @Test
    fun `unknown boundary does not destroy previously known structure`() {
        val store = ReactorStructureWorldStore()
        val formedScan = scanOf(
            snapshot(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
            snapshot(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )
        val formedPlan = planner.buildPlan(formedScan)
        store.applyPlan(
            scannedPositions = formedScan.loadedReactorPositions,
            plan = formedPlan,
            removeMissingStructures = true,
        )

        val unknownScan = ReactorAssemblyScan(
            blocks = setOf(
                snapshot(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                ReactorWorldBlockSnapshot(pos(-2, 0, 0), kind = null, facing = null, loaded = false),
            ),
            hasUnknownBoundary = true,
        )
        val unknownPlan = planner.buildPlan(unknownScan)
        store.applyPlan(
            scannedPositions = unknownScan.loadedReactorPositions,
            plan = unknownPlan,
            removeMissingStructures = false,
        )

        assertNotNull(store.membershipAt(pos(0, 0, 0)))
        assertEquals(ReactorControllerFormationState.UNKNOWN, store.controllerViewState(pos(-1, 0, 0)).formationState)
    }

    @Test
    fun `removing a position removes only its owning structure`() {
        val store = ReactorStructureWorldStore()
        val scan = scanOf(
            snapshot(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
            snapshot(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
            snapshot(3, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
            snapshot(2, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )
        val plan = planner.buildPlan(scan)
        store.applyPlan(
            scannedPositions = scan.loadedReactorPositions,
            plan = plan,
            removeMissingStructures = true,
        )

        store.removePosition(pos(0, 0, 0))

        assertNull(store.membershipAt(pos(-1, 0, 0)))
        assertNotNull(store.membershipAt(pos(2, 0, 0)))
        assertNotNull(store.membershipAt(pos(3, 0, 0)))
    }

    @Test
    fun `adjacent maximum columns remain separate structures`() {
        val store = ReactorStructureWorldStore()
        val scan = scanOf(
            *verticalColumn(x = 0, z = 0, height = 4).toTypedArray(),
            snapshot(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
            *verticalColumn(x = 1, z = 0, height = 4).toTypedArray(),
            snapshot(2, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )
        val plan = planner.buildPlan(scan)

        store.applyPlan(
            scannedPositions = scan.loadedReactorPositions,
            plan = plan,
            removeMissingStructures = true,
        )

        val first = store.controllerViewState(pos(-1, 0, 0))
        val second = store.controllerViewState(pos(2, 0, 0))
        assertEquals(ReactorControllerFormationState.FORMED, first.formationState)
        assertEquals(ReactorControllerFormationState.FORMED, second.formationState)
        assertEquals(4, first.chamberBlockCount)
        assertEquals(4, second.chamberBlockCount)
        assertEquals(first.structureId, store.membershipAt(pos(0, 0, 0))?.structureId)
        assertEquals(second.structureId, store.membershipAt(pos(1, 0, 0))?.structureId)
    }

    @Test
    fun `extra chamber next to maximum structure stays outside membership`() {
        val store = ReactorStructureWorldStore()
        val scan = scanOf(
            *verticalColumn(x = 0, z = 0, height = 4).toTypedArray(),
            snapshot(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
            snapshot(1, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
        )
        val plan = planner.buildPlan(scan)

        store.applyPlan(
            scannedPositions = scan.loadedReactorPositions,
            plan = plan,
            removeMissingStructures = true,
        )

        assertEquals(ReactorControllerFormationState.FORMED, store.controllerViewState(pos(-1, 0, 0)).formationState)
        assertNull(store.membershipAt(pos(1, 0, 0)))
    }

    private fun scanOf(vararg blocks: ReactorWorldBlockSnapshot): ReactorAssemblyScan =
        ReactorAssemblyScan(blocks = blocks.toSet(), hasUnknownBoundary = false)

    private fun verticalColumn(x: Int, z: Int, height: Int): List<ReactorWorldBlockSnapshot> =
        (0 until height).map { y -> snapshot(x, y, z, ReactorMultiblockBlockKind.CHAMBER) }

    private fun snapshot(
        x: Int,
        y: Int,
        z: Int,
        kind: ReactorMultiblockBlockKind,
        facing: ReactorBlockDirection? = null,
    ): ReactorWorldBlockSnapshot =
        ReactorWorldBlockSnapshot(
            position = pos(x, y, z),
            kind = kind,
            facing = facing,
            loaded = true,
        )

    private fun pos(x: Int, y: Int, z: Int): ReactorBlockPosition =
        ReactorBlockPosition(x, y, z)
}
