package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockDirection
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
            blocks = squareTank(baseSize = 2, height = 1) + listOf(
                block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                block(1, 1, 0, ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT),
            ),
        )

        assertEquals(ReactorStructureId(structureId), definition.structureId)
        assertEquals(pos(-1, 0, 0), definition.controllerPosition)
        assertEquals(0.008, definition.totalVolumeCubicMeters)
        assertEquals(4, definition.zone.plainChamberPositions.size)
        assertEquals(4, definition.zone.volumePositions.size)
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
    fun `freeform strategy rejects disconnected chamber blocks`() {
        val freeformAssembler = ReactorMultiblockAssembler(
            ReactorMultiblockRules(
                chamberVolumeCubicMeters = 0.002,
                chamberShapeStrategy = FreeformChamberShapeStrategy,
            ),
        )
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            freeformAssembler.assemble(
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
    fun `vertical tank strategy does not absorb non rectangular adjacent chamber`() {
        val definition = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = listOf(
                block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(1, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
            ),
        )

        assertEquals(setOf(pos(0, 0, 0)), definition.zone.plainChamberPositions)
    }

    @Test
    fun `rejects controller that does not touch chamber`() {
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            assembler.assemble(
                structureId = UUID.randomUUID(),
                blocks = listOf(
                    block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(0, 5, 0, ReactorMultiblockBlockKind.CONTROLLER),
                ),
            )
        }

        assertEquals(true, error.validationErrors.any { it.contains("complete vertical tank") })
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
    fun `rejects embedded port without explicit facing`() {
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            assembler.assemble(
                structureId = UUID.randomUUID(),
                blocks = listOf(
                    block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                    block(1, 0, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                    block(0, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(0, 0, 1, ReactorMultiblockBlockKind.CHAMBER),
                    block(1, 0, 1, ReactorMultiblockBlockKind.CHAMBER),
                ),
            )
        }

        assertEquals(true, error.validationErrors.any { it.contains("embedded reactor port") && it.contains("must have explicit facing") })
    }

    @Test
    fun `assigns port indexes deterministically within each port kind`() {
        val definition = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = squareTank(baseSize = 2, height = 1) + listOf(
                block(0, -1, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(2, 0, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                block(1, 1, 0, ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT),
            ),
        )

        val itemInputs = definition.portsOfKind(ReactorPortKind.ITEM_INPUT)
        assertEquals(listOf(pos(0, 1, 0), pos(2, 0, 0)), itemInputs.map { it.position })
        assertEquals(listOf(0, 1), itemInputs.map { it.portIndex })
        assertEquals(0, definition.portsOfKind(ReactorPortKind.ITEM_OUTPUT).single().portIndex)
    }

    @Test
    fun `vertical tank allows configured maximum heights by base size`() {
        val oneByOne = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = squareTank(baseSize = 1, height = 4) + block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )
        val twoByTwo = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = squareTank(baseSize = 2, height = 8) + block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )
        val threeByThree = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = squareTank(baseSize = 3, height = 16) + block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )

        assertEquals(4, oneByOne.zone.volumePositions.size)
        assertEquals(32, twoByTwo.zone.volumePositions.size)
        assertEquals(144, threeByThree.zone.volumePositions.size)
    }

    @Test
    fun `vertical tank does not absorb extra chamber next to maximum structure`() {
        val blocks = squareTank(baseSize = 2, height = 8) + listOf(
            block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
            block(2, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
        )

        val definition = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = blocks,
        )

        assertEquals(32, definition.zone.volumePositions.size)
        assertEquals(setOf(pos(2, 0, 0)), definition.inactiveChamberPositions)
        assertEquals(false, pos(2, 0, 0) in definition.zone.volumePositions)
    }

    @Test
    fun `maximum chamber blocks chooses smaller valid tank instead of rejecting larger candidates`() {
        val limitedAssembler = ReactorMultiblockAssembler(
            ReactorMultiblockRules(
                chamberVolumeCubicMeters = 0.002,
                maximumVolumeBlocks = 4,
            ),
        )

        val definition = limitedAssembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = squareTank(baseSize = 2, height = 2) + block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
        )

        assertEquals(4, definition.zone.volumePositions.size)
        assertEquals(4, definition.inactiveChamberPositions.size)
    }

    @Test
    fun `embedded controller and ports count as reactor volume without being chamber contacts`() {
        val definition = assembler.assemble(
            structureId = UUID.randomUUID(),
            blocks = listOf(
                block(0, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                block(0, 1, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT, ReactorBlockDirection.NORTH),
                block(0, 2, 0, ReactorMultiblockBlockKind.CHAMBER),
            ),
        )

        assertEquals(pos(0, 0, 0), definition.controllerPosition)
        assertEquals(setOf(pos(0, 2, 0)), definition.zone.plainChamberPositions)
        assertEquals(setOf(pos(0, 0, 0), pos(0, 1, 0), pos(0, 2, 0)), definition.zone.volumePositions)
        assertEquals(0.006, definition.totalVolumeCubicMeters)
        val inputPort = definition.portsOfKind(ReactorPortKind.ITEM_INPUT).single()
        assertEquals(pos(0, 1, 0), inputPort.attachedChamberPosition)
        assertEquals(ReactorBlockDirection.NORTH, inputPort.contactDirection)
    }

    @Test
    fun `port attached only to extra chamber is rejected`() {
        val error = assertFailsWith<ReactorMultiblockValidationException> {
            assembler.assemble(
                structureId = UUID.randomUUID(),
                blocks = squareTank(baseSize = 2, height = 8) + listOf(
                    block(-1, 0, 0, ReactorMultiblockBlockKind.CONTROLLER),
                    block(2, 0, 0, ReactorMultiblockBlockKind.CHAMBER),
                    block(3, 0, 0, ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
                ),
            )
        }

        assertEquals(true, error.validationErrors.any { it.contains("ITEM_INPUT") && it.contains("must touch") })
    }

    private fun block(
        x: Int,
        y: Int,
        z: Int,
        kind: ReactorMultiblockBlockKind,
        facing: ReactorBlockDirection? = null,
    ): ReactorMultiblockBlock =
        ReactorMultiblockBlock(pos(x, y, z), kind, facing)

    private fun squareTank(baseSize: Int, height: Int): List<ReactorMultiblockBlock> =
        buildList {
            for (x in 0 until baseSize) {
                for (y in 0 until height) {
                    for (z in 0 until baseSize) {
                        add(block(x, y, z, ReactorMultiblockBlockKind.CHAMBER))
                    }
                }
            }
        }

    private fun pos(x: Int, y: Int, z: Int): ReactorBlockPosition =
        ReactorBlockPosition(x, y, z)
}
