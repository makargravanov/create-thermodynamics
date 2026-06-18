package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlock
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockValidationException
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import java.util.UUID
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith

class ReactorMultiblockAssemblerTest {
    private val assembler = ReactorMultiblockAssembler(
        ReactorMultiblockRules(chamberVolumeCubicMeters = 0.002),
    )

    @Test
    fun `assembles one connected chamber zone with typed ports`() {
        val structureId = UUID.fromString("2f10b8a4-fd3a-4d6e-b2b1-61c7ef0cdb6b")

        val definition = assembler.assemble(
            structureId = structureId,
            blocks = listOf(
                block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(1, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                block(1, 1, 0, ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT),
            ),
        )

        assertEquals(ReactorStructureId(structureId), definition.structureId)
        assertEquals(pos(-1, 0, 0), definition.controllerPosition)
        assertEquals(0.004, definition.totalVolumeCubicMeters)
        assertEquals(2, definition.zone.chamberPositions.size)
        assertEquals(
            listOf(ReactorPortKind.ITEM_INPUT, ReactorPortKind.FLUID_OUTPUT),
            definition.ports.map { it.kind },
        )
        assertEquals(0, definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single().portIndex)
        assertEquals(0, definition.portsOfKind(ReactorPortKind.FLUID_OUTPUT).single().portIndex)
    }

    @Test
    fun `rejects duplicate positions`() {
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            assembler.assemble(
                structureId = UUID.randomUUID(),
                blocks = listOf(
                    block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(0, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                ),
            )
        }

        assertEquals(true, error.validationErrors.any { it.contains("duplicate reactor multiblock block") })
    }

    @Test
    fun `rejects disconnected chamber blocks`() {
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            assembler.assemble(
                structureId = UUID.randomUUID(),
                blocks = listOf(
                    block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(10, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                ),
            )
        }

        assertEquals(true, error.validationErrors.any { it.contains("one face-connected zone") })
    }

    @Test
    fun `rejects controller that does not touch chamber`() {
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            assembler.assemble(
                structureId = UUID.randomUUID(),
                blocks = listOf(
                    block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(3, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                ),
            )
        }

        assertEquals(true, error.validationErrors.any { it.contains("controller") && it.contains("must touch") })
    }

    @Test
    fun `rejects port that does not touch chamber`() {
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            assembler.assemble(
                structureId = UUID.randomUUID(),
                blocks = listOf(
                    block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                    block(3, 0, 0, ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT),
                ),
            )
        }

        assertEquals(true, error.validationErrors.any { it.contains("ITEM_OUTPUT") && it.contains("must touch") })
    }

    @Test
    fun `assigns port indexes deterministically within each port kind`() {
        val definition = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = listOf(
                block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(1, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(2, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(0, -1, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(2, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                block(1, 1, 0, ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT),
            ),
        )

        val itemInputs = definition.portsOfKind(ReactorPortKind.ITEM_INPUT)
        assertEquals(listOf(pos(0, 1, 0), pos(2, 1, 0)), itemInputs.map { it.position })
        assertEquals(listOf(0, 1), itemInputs.map { it.portIndex })
        assertEquals(0, definition.portsOfKind(ReactorPortKind.ITEM_OUTPUT).single().portIndex)
    }

    private fun block(
        x: Int,
        y: Int,
        z: Int,
        kind: ReactorMultiblockBlockKind,
    ): ReactorMultiblockBlock =
        ReactorMultiblockBlock(pos(x, y, z), kind)

    private fun pos(x: Int, y: Int, z: Int): ReactorBlockPosition =
        ReactorBlockPosition(x, y, z)
}
