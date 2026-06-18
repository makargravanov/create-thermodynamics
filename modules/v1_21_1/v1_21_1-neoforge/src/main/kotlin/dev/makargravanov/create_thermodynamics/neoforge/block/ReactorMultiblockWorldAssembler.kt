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
import java.util.WeakHashMap

object ReactorMultiblockWorldAssembler {
    private const val SCAN_LIMIT = 512

    private val assembler = ReactorMultiblockAssembler(
        ReactorMultiblockRules(chamberVolumeCubicMeters = 1.0),
    )
    private val statesByLevel = WeakHashMap<ServerLevel, ReactorMultiblockWorldState>()

    fun rebuildAround(level: ServerLevel, origin: BlockPos) {
        val scan = scanConnectedReactorBlocks(level, origin)
        if (scan.positions.isEmpty()) {
            if (!scan.hasUnknownBoundary) {
                clearNeighbouringDanglingMembership(level, origin)
            }
            return
        }

        val snapshot = buildAssemblySnapshot(level, scan.positions)
        val affectedPositions = worldState(level).applySnapshot(
            scannedPositions = scan.positions,
            definitions = snapshot.definitions,
            removeMissingStructures = !scan.hasUnknownBoundary,
        )

        applyMembership(level, affectedPositions)
        notifyVisibleShape(level, affectedPositions)
    }

    fun clearMembership(level: ServerLevel, origin: BlockPos) {
        val affectedPositions = worldState(level).removePosition(origin) +
            Direction.entries.mapTo(linkedSetOf()) { origin.relative(it) }
        clearNeighbouringDanglingMembership(level, origin)
        applyMembership(level, affectedPositions)
        notifyVisibleShape(level, affectedPositions)
    }

    private fun buildAssemblySnapshot(
        level: ServerLevel,
        positions: Set<BlockPos>,
    ): ReactorWorldAssemblySnapshot {
        val definitions = mutableListOf<ReactorMultiblockDefinition>()
        val used = linkedSetOf<BlockPos>()

        val controllers = positions
            .filter { (level.getBlockState(it).block as? ReactorMultiblockBlock)?.kind == ReactorMultiblockKind.CONTROLLER }
            .sortedWith(compareBy<BlockPos> { it.x }.thenBy { it.y }.thenBy { it.z })

        for (controller in controllers) {
            if (controller in used) {
                continue
            }
            val candidatePositions = collectCandidateForController(level, positions, controller, used)
            val structureId = controllerStructureId(controller)
            val definition = assembleCandidate(level, structureId, candidatePositions)
            if (definition != null) {
                definitions += definition
                used += definition.structurePositions()
            }
        }

        return ReactorWorldAssemblySnapshot(definitions)
    }

    private fun scanConnectedReactorBlocks(level: ServerLevel, origin: BlockPos): ReactorWorldScan {
        val visited = linkedSetOf<BlockPos>()
        val queue = ArrayDeque<BlockPos>()
        var hasUnknownBoundary = false

        fun enqueue(pos: BlockPos) {
            if (!level.isLoaded(pos)) {
                hasUnknownBoundary = true
                return
            }
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
        return ReactorWorldScan(visited, hasUnknownBoundary)
    }

    private fun clearNeighbouringDanglingMembership(level: ServerLevel, origin: BlockPos) {
        for (direction in Direction.entries) {
            val neighbour = origin.relative(direction)
            val blockEntity = reactorBlockEntity(level, neighbour) ?: continue
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

    private fun applyMembership(level: ServerLevel, positions: Set<BlockPos>) {
        val state = worldState(level)
        for (pos in positions) {
            if (!level.isLoaded(pos) || level.getBlockState(pos).block !is ReactorMultiblockBlock) {
                continue
            }
            val assignment = state.assignmentAt(pos)
            requireNotNull(reactorBlockEntity(level, pos)) {
                "cannot assign reactor membership to non-reactor block at $pos"
            }.setStructureMembership(assignment?.structureId, assignment?.activeVolumeBlock ?: false)
        }
    }

    private fun notifyVisibleShape(level: ServerLevel, positions: Set<BlockPos>) {
        for (pos in positions) {
            val state = level.getBlockState(pos)
            if (level.isLoaded(pos) && state.block is ReactorMultiblockBlock) {
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

    private fun blockKind(level: ServerLevel, pos: BlockPos): ReactorMultiblockKind? =
        (level.getBlockState(pos).block as? ReactorMultiblockBlock)?.kind

    private fun reactorBlockEntity(level: ServerLevel, pos: BlockPos): ReactorMultiblockBlockEntity? {
        if (level.getBlockState(pos).block !is ReactorMultiblockBlock) {
            return null
        }
        return requireNotNull(level.getBlockEntity(pos) as? ReactorMultiblockBlockEntity) {
            "reactor block at $pos must have ReactorMultiblockBlockEntity"
        }
    }

    private fun worldState(level: ServerLevel): ReactorMultiblockWorldState =
        statesByLevel.getOrPut(level) { ReactorMultiblockWorldState() }

    private fun ReactorMultiblockDefinition.structurePositions(): Set<BlockPos> =
        linkedSetOf<BlockPos>().also { positions ->
            positions += zone.volumePositions.map(::toBlockPos)
            positions += toBlockPos(controllerPosition)
            positions += ports.map { toBlockPos(it.position) }
        }

    private fun ReactorMultiblockDefinition.assignmentAt(pos: BlockPos): StructureAssignment? {
        if (pos !in structurePositions()) {
            return null
        }
        return StructureAssignment(
            structureId = structureId.value,
            activeVolumeBlock = toReactorPosition(pos) in zone.volumePositions,
        )
    }

    private data class ReactorWorldScan(
        val positions: Set<BlockPos>,
        val hasUnknownBoundary: Boolean,
    )

    private data class ReactorWorldAssemblySnapshot(
        val definitions: List<ReactorMultiblockDefinition>,
    )

    private class ReactorMultiblockWorldState {
        private val definitions = linkedMapOf<UUID, ReactorMultiblockDefinition>()
        private val positionsByStructure = linkedMapOf<UUID, Set<BlockPos>>()
        private val structureByPosition = linkedMapOf<BlockPos, UUID>()

        fun applySnapshot(
            scannedPositions: Set<BlockPos>,
            definitions: List<ReactorMultiblockDefinition>,
            removeMissingStructures: Boolean,
        ): Set<BlockPos> {
            val affectedPositions = linkedSetOf<BlockPos>()
            val newIds = definitions.mapTo(linkedSetOf()) { it.structureId.value }

            val intersectingIds = scannedPositions
                .mapNotNullTo(linkedSetOf()) { structureByPosition[it] }

            if (removeMissingStructures) {
                for (structureId in intersectingIds - newIds) {
                    affectedPositions += removeStructure(structureId)
                }
            }

            for (definition in definitions) {
                val structureId = definition.structureId.value
                affectedPositions += removeStructure(structureId)
                val positions = definition.structurePositions()
                this.definitions[structureId] = definition
                positionsByStructure[structureId] = positions
                for (pos in positions) {
                    structureByPosition[pos] = structureId
                }
                affectedPositions += positions
            }

            if (removeMissingStructures) {
                affectedPositions += scannedPositions
            }

            return affectedPositions
        }

        fun removePosition(pos: BlockPos): Set<BlockPos> {
            val structureId = structureByPosition[pos] ?: return setOf(pos)
            return removeStructure(structureId) + pos
        }

        fun assignmentAt(pos: BlockPos): StructureAssignment? {
            val structureId = structureByPosition[pos] ?: return null
            return definitions[structureId]?.assignmentAt(pos)
        }

        private fun removeStructure(structureId: UUID): Set<BlockPos> {
            definitions.remove(structureId)
            val positions = positionsByStructure.remove(structureId).orEmpty()
            for (pos in positions) {
                structureByPosition.remove(pos)
            }
            return positions
        }
    }

    private data class StructureAssignment(
        val structureId: UUID,
        val activeVolumeBlock: Boolean,
    )
}
