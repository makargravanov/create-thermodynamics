package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockRules
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockDirection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorAssemblyPlanner
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorAssemblyScan
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorStructureWorldStore
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorWorldBlockSnapshot
import net.minecraft.core.BlockPos
import net.minecraft.core.Direction
import net.minecraft.server.level.ServerLevel
import net.minecraft.world.level.block.Block
import java.util.WeakHashMap

object ReactorMultiblockWorldAssembler {
    private const val SCAN_LIMIT = 512

    private val planner = ReactorAssemblyPlanner(
        ReactorMultiblockRules(chamberVolumeCubicMeters = 1.0),
    )
    private val storesByLevel = WeakHashMap<ServerLevel, ReactorStructureWorldStore>()

    fun rebuildAround(level: ServerLevel, origin: BlockPos) {
        val scan = scanWorld(level, origin)
        val store = worldStore(level)
        if (scan.loadedReactorPositions.isEmpty()) {
            if (!scan.hasUnknownBoundary) {
                val affected = (setOf(origin) + Direction.entries.map { origin.relative(it) })
                    .mapTo(linkedSetOf()) { it.toReactorPosition() }
                    .flatMapTo(linkedSetOf()) { store.removePosition(it) }
                    .mapTo(linkedSetOf(), ::toBlockPos)
                applyPlan(level, affected)
                syncChangedBlocks(level, affected)
            }
            return
        }

        val plan = planner.buildPlan(scan)
        val affected = store.applyPlan(
            scannedPositions = scan.loadedReactorPositions,
            plan = plan,
            removeMissingStructures = !scan.hasUnknownBoundary,
        ).mapTo(linkedSetOf(), ::toBlockPos)

        applyPlan(level, affected)
        syncChangedBlocks(level, affected)
    }

    fun clearMembership(level: ServerLevel, origin: BlockPos) {
        val store = worldStore(level)
        val affected = (store.removePosition(origin.toReactorPosition()) +
            Direction.entries.mapTo(linkedSetOf()) { origin.relative(it).toReactorPosition() })
            .mapTo(linkedSetOf(), ::toBlockPos)
        applyPlan(level, affected)
        syncChangedBlocks(level, affected)
    }

    private fun scanWorld(level: ServerLevel, origin: BlockPos): ReactorAssemblyScan {
        val snapshots = linkedMapOf<ReactorBlockPosition, ReactorWorldBlockSnapshot>()
        val visited = linkedSetOf<BlockPos>()
        val queue = ArrayDeque<BlockPos>()
        var hasUnknownBoundary = false

        fun enqueue(pos: BlockPos) {
            if (!level.isLoaded(pos)) {
                hasUnknownBoundary = true
                snapshots[pos.toReactorPosition()] = ReactorWorldBlockSnapshot(
                    position = pos.toReactorPosition(),
                    kind = null,
                    facing = null,
                    loaded = false,
                )
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
            val state = level.getBlockState(pos)
            val block = state.block as? ReactorMultiblockBlock
            snapshots[pos.toReactorPosition()] = ReactorWorldBlockSnapshot(
                position = pos.toReactorPosition(),
                kind = block?.kind?.modelKind,
                facing = block?.modelFacing(state)?.toReactorDirection(),
                loaded = true,
            )
            for (direction in Direction.entries) {
                enqueue(pos.relative(direction))
            }
        }

        return ReactorAssemblyScan(
            blocks = snapshots.values.toSet(),
            hasUnknownBoundary = hasUnknownBoundary,
        )
    }

    private fun applyPlan(level: ServerLevel, positions: Set<BlockPos>) {
        val store = worldStore(level)
        for (pos in positions) {
            if (!level.isLoaded(pos) || level.getBlockState(pos).block !is ReactorMultiblockBlock) {
                continue
            }
            val blockEntity = reactorBlockEntity(level, pos)
            val membership = store.membershipAt(pos.toReactorPosition())
            val controllerViewState = if ((level.getBlockState(pos).block as ReactorMultiblockBlock).kind == ReactorMultiblockKind.CONTROLLER) {
                store.controllerViewState(pos.toReactorPosition())
            } else {
                null
            }
            blockEntity.applyWorldProjection(membership, controllerViewState)
        }
    }

    private fun syncChangedBlocks(level: ServerLevel, positions: Set<BlockPos>) {
        for (pos in positions) {
            val state = level.getBlockState(pos)
            if (level.isLoaded(pos) && state.block is ReactorMultiblockBlock) {
                level.sendBlockUpdated(pos, state, state, Block.UPDATE_CLIENTS)
            }
        }
    }

    private fun reactorBlockEntity(level: ServerLevel, pos: BlockPos): ReactorMultiblockBlockEntity =
        requireNotNull(level.getBlockEntity(pos) as? ReactorMultiblockBlockEntity) {
            "reactor block at $pos must have ReactorMultiblockBlockEntity"
        }

    private fun worldStore(level: ServerLevel): ReactorStructureWorldStore =
        storesByLevel.getOrPut(level) { ReactorStructureWorldStore() }

    private fun BlockPos.toReactorPosition(): ReactorBlockPosition =
        ReactorBlockPosition(x, y, z)

    private fun toBlockPos(pos: ReactorBlockPosition): BlockPos =
        BlockPos(pos.x, pos.y, pos.z)

    private fun Direction.toReactorDirection(): ReactorBlockDirection =
        when (this) {
            Direction.EAST -> ReactorBlockDirection.EAST
            Direction.WEST -> ReactorBlockDirection.WEST
            Direction.UP -> ReactorBlockDirection.UP
            Direction.DOWN -> ReactorBlockDirection.DOWN
            Direction.SOUTH -> ReactorBlockDirection.SOUTH
            Direction.NORTH -> ReactorBlockDirection.NORTH
        }
}
