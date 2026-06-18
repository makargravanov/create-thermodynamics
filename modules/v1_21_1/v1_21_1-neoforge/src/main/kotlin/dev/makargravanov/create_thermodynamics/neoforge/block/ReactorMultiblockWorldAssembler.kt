package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockAssembler
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockRules
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlock as ModelReactorMultiblockBlock
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockValidationException
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import net.minecraft.core.BlockPos
import net.minecraft.core.Direction
import net.minecraft.server.level.ServerLevel
import net.minecraft.world.level.block.Block
import java.util.UUID

object ReactorMultiblockWorldAssembler {
    private const val SCAN_LIMIT = 512

    private val assembler = ReactorMultiblockAssembler(
        ReactorMultiblockRules(chamberVolumeCubicMeters = 1.0),
    )

    fun rebuildAround(level: ServerLevel, origin: BlockPos) {
        val positions = collectConnectedReactorBlocks(level, origin)
        if (positions.isEmpty()) {
            clearNeighbouringDanglingMembership(level, origin)
            return
        }

        val assignments = linkedMapOf<BlockPos, StructureAssignment>()
        val used = linkedSetOf<BlockPos>()

        val controllers = positions
            .filter { (level.getBlockState(it).block as? ReactorMultiblockBlock)?.kind == ReactorMultiblockKind.CONTROLLER }
            .sortedWith(compareBy<BlockPos> { it.x }.thenBy { it.y }.thenBy { it.z })

        for (controller in controllers) {
            if (controller in used) {
                continue
            }
            val candidatePositions = collectCandidateForController(level, positions, controller, used)
            val structureId = membership(level, controller)?.structureId ?: controllerStructureId(controller)
            val definition = assembleCandidate(level, structureId, candidatePositions)
            if (definition != null) {
                assignDefinition(definition, assignments, used)
            }
        }

        for (pos in positions) {
            val assignment = assignments[pos]
            setMembership(level, pos, assignment?.structureId, assignment?.activeVolumeBlock ?: false)
        }
        notifyVisibleShape(level, positions)
    }

    fun clearMembership(level: ServerLevel, origin: BlockPos) {
        clearNeighbouringDanglingMembership(level, origin)
        notifyVisibleShape(level, Direction.entries.mapTo(linkedSetOf()) { origin.relative(it) })
    }

    private fun collectConnectedReactorBlocks(level: ServerLevel, origin: BlockPos): Set<BlockPos> {
        val visited = linkedSetOf<BlockPos>()
        val queue = ArrayDeque<BlockPos>()

        fun enqueue(pos: BlockPos) {
            if (pos !in visited && level.getBlockState(pos).block is ReactorMultiblockBlock) {
                queue += pos.immutable()
            }
        }

        enqueue(origin)
        for (direction in Direction.entries) {
            enqueue(origin.relative(direction))
        }

        while (queue.isNotEmpty()) {
            val pos = queue.removeFirst()
            if (!visited.add(pos)) {
                continue
            }
            check(visited.size <= SCAN_LIMIT) {
                "reactor multiblock scan exceeded $SCAN_LIMIT blocks near $origin"
            }
            for (direction in Direction.entries) {
                enqueue(pos.relative(direction))
            }
        }
        return visited
    }

    private fun clearNeighbouringDanglingMembership(level: ServerLevel, origin: BlockPos) {
        for (direction in Direction.entries) {
            val neighbour = origin.relative(direction)
            val blockEntity = level.getBlockEntity(neighbour) as? ReactorMultiblockBlockEntity ?: continue
            blockEntity.setStructureMembership(null, false)
        }
    }

    private fun assembleCandidate(
        level: ServerLevel,
        structureId: UUID,
        positions: Set<BlockPos>,
    ): ReactorMultiblockDefinition? {
        val blocks = positions.mapNotNull { pos ->
            blockKind(level, pos)?.let { ModelReactorMultiblockBlock(toReactorPosition(pos), it.modelKind) }
        }
        return try {
            assembler.assemble(ReactorStructureId(structureId), blocks)
        } catch (_: ReactorMultiblockValidationException) {
            null
        }
    }

    private fun assignDefinition(
        definition: ReactorMultiblockDefinition,
        assignments: MutableMap<BlockPos, StructureAssignment>,
        used: MutableSet<BlockPos>,
    ) {
        val volumePositions = definition.zone.volumePositions.mapTo(linkedSetOf(), ::toBlockPos)
        val structurePositions = linkedSetOf<BlockPos>()
        structurePositions += volumePositions
        structurePositions += toBlockPos(definition.controllerPosition)
        structurePositions += definition.ports.map { toBlockPos(it.position) }

        for (pos in structurePositions) {
            if (pos in used) {
                continue
            }
            val assignment = StructureAssignment(
                structureId = definition.structureId.value,
                activeVolumeBlock = pos in volumePositions,
            )
            assignments[pos] = assignment
            used += pos
        }
    }

    private fun collectCandidateForController(
        level: ServerLevel,
        componentPositions: Set<BlockPos>,
        controller: BlockPos,
        alreadyUsed: Set<BlockPos>,
    ): Set<BlockPos> {
        val visited = linkedSetOf<BlockPos>()
        val queue = ArrayDeque<BlockPos>()
        queue += controller

        while (queue.isNotEmpty()) {
            val pos = queue.removeFirst()
            if (!visited.add(pos)) {
                continue
            }
            for (direction in Direction.entries) {
                val next = pos.relative(direction)
                if (
                    next in componentPositions &&
                    next !in alreadyUsed &&
                    (next == controller || blockKind(level, next) != ReactorMultiblockKind.CONTROLLER)
                ) {
                    queue += next
                }
            }
        }
        return visited
    }

    private fun setMembership(level: ServerLevel, pos: BlockPos, structureId: UUID?, active: Boolean) {
        val blockEntity = level.getBlockEntity(pos) as? ReactorMultiblockBlockEntity ?: return
        blockEntity.setStructureMembership(structureId, active)
    }

    private fun notifyVisibleShape(level: ServerLevel, positions: Set<BlockPos>) {
        for (pos in positions) {
            val state = level.getBlockState(pos)
            if (state.block is ReactorMultiblockBlock) {
                level.sendBlockUpdated(pos, state, state, Block.UPDATE_CLIENTS)
            }
        }
    }

    private fun controllerStructureId(controller: BlockPos): UUID =
        UUID.nameUUIDFromBytes("reactor-controller:${controller.x},${controller.y},${controller.z}".toByteArray(Charsets.UTF_8))

    private fun toReactorPosition(pos: BlockPos): ReactorBlockPosition =
        ReactorBlockPosition(pos.x, pos.y, pos.z)

    private fun toBlockPos(pos: ReactorBlockPosition): BlockPos =
        BlockPos(pos.x, pos.y, pos.z)

    private fun membership(level: ServerLevel, pos: BlockPos): ReactorMultiblockBlockEntity? =
        level.getBlockEntity(pos) as? ReactorMultiblockBlockEntity

    private fun blockKind(level: ServerLevel, pos: BlockPos): ReactorMultiblockKind? =
        (level.getBlockState(pos).block as? ReactorMultiblockBlock)?.kind

    private data class StructureAssignment(
        val structureId: UUID,
        val activeVolumeBlock: Boolean,
    )
}
