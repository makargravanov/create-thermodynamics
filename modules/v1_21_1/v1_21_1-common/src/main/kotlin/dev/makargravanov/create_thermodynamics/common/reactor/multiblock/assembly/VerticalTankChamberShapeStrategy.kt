package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorZoneDescriptor

object VerticalTankChamberShapeStrategy : ReactorChamberShapeStrategy {
    private val maxHeightByBaseSize = linkedMapOf(
        1 to 4,
        2 to 8,
        3 to 16,
    )

    override fun buildZone(
        chambers: Set<ReactorBlockPosition>,
        volumeCapableBlocks: Map<ReactorBlockPosition, ReactorMultiblockBlockKind>,
        controller: ReactorBlockPosition?,
        chamberVolumeCubicMeters: Double,
        maximumVolumeBlocks: Int?,
    ): ReactorChamberShapeResult {
        if (chambers.isEmpty()) {
            return ReactorChamberShapeResult(null, errors = listOf("reactor multiblock must contain chamber blocks"))
        }
        if (controller == null) {
            return ReactorChamberShapeResult(null, errors = listOf("reactor vertical tank chamber requires a controller"))
        }

        val selected = findBestTank(chambers, volumeCapableBlocks, controller, maximumVolumeBlocks)
        if (selected == null) {
            return ReactorChamberShapeResult(
                zone = null,
                inactiveChamberPositions = emptySet(),
                errors = listOf(
                    "reactor chamber must be a complete vertical tank: square base 1x1, 2x2 or 3x3; max heights are 4, 8 and 16",
                ),
            )
        }

        val plainChamberPositions = selected.positions.toSortedSet()
        val inactiveChambers = (chambers - selected.positions).toSortedSet()
        return ReactorChamberShapeResult(
            zone = ReactorZoneDescriptor(
                zoneIndex = 0,
                volumePositions = selected.volumePositions.toSortedSet(),
                plainChamberPositions = plainChamberPositions,
                volumeCubicMeters = selected.volumePositions.size * chamberVolumeCubicMeters,
            ),
            inactiveChamberPositions = inactiveChambers,
            errors = emptyList(),
        )
    }

    private fun findBestTank(
        chambers: Set<ReactorBlockPosition>,
        volumeCapableBlocks: Map<ReactorBlockPosition, ReactorMultiblockBlockKind>,
        controller: ReactorBlockPosition,
        maximumVolumeBlocks: Int?,
    ): TankCandidate? {
        val candidates = mutableListOf<TankCandidate>()
        val ys = volumeCapableBlocks.keys.map { it.y }.toSet()

        for ((baseSize, maxHeight) in maxHeightByBaseSize) {
            val possibleMinX = volumeCapableBlocks.keys.flatMap { block -> (block.x - baseSize + 1)..block.x }.toSet()
            val possibleMinZ = volumeCapableBlocks.keys.flatMap { block -> (block.z - baseSize + 1)..block.z }.toSet()
            for (minX in possibleMinX) {
                for (minZ in possibleMinZ) {
                    for (minY in ys) {
                        for (height in 1..maxHeight) {
                            val volumePositions = tankPositions(minX, minY, minZ, baseSize, height)
                            val plainChamberPositions = volumePositions.filterTo(linkedSetOf()) { volumeCapableBlocks[it] == ReactorMultiblockBlockKind.CHAMBER }
                            if (
                                volumePositions.all { it in volumeCapableBlocks } &&
                                plainChamberPositions.isNotEmpty() &&
                                (maximumVolumeBlocks == null || volumePositions.size <= maximumVolumeBlocks) &&
                                (controller in volumePositions || plainChamberPositions.any { it in controller.faceNeighbours() })
                            ) {
                                candidates += TankCandidate(
                                    baseSize = baseSize,
                                    height = height,
                                    minPosition = ReactorBlockPosition(minX, minY, minZ),
                                    positions = plainChamberPositions,
                                    volumePositions = volumePositions,
                                )
                            }
                        }
                    }
                }
            }
        }

        return candidates.maxWithOrNull(
            compareBy<TankCandidate> { it.positions.size }
                .thenBy { it.baseSize }
                .thenBy { it.height }
                .thenByDescending { it.minPosition },
        )
    }

    private fun tankPositions(
        minX: Int,
        minY: Int,
        minZ: Int,
        baseSize: Int,
        height: Int,
    ): Set<ReactorBlockPosition> {
        val positions = linkedSetOf<ReactorBlockPosition>()
        for (x in minX until minX + baseSize) {
            for (y in minY until minY + height) {
                for (z in minZ until minZ + baseSize) {
                    positions += ReactorBlockPosition(x, y, z)
                }
            }
        }
        return positions
    }

    private data class TankCandidate(
        val baseSize: Int,
        val height: Int,
        val minPosition: ReactorBlockPosition,
        val positions: Set<ReactorBlockPosition>,
        val volumePositions: Set<ReactorBlockPosition>,
    )
}
