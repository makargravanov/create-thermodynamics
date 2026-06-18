package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorZoneDescriptor

object FreeformChamberShapeStrategy : ReactorChamberShapeStrategy {
    override fun buildZone(
        chambers: Set<ReactorBlockPosition>,
        volumeCapableBlocks: Map<ReactorBlockPosition, ReactorMultiblockBlockKind>,
        controller: ReactorBlockPosition?,
        chamberVolumeCubicMeters: Double,
        maximumChamberBlocks: Int?,
    ): ReactorChamberShapeResult {
        if (chambers.isEmpty()) {
            return ReactorChamberShapeResult(null, errors = listOf("reactor multiblock must contain chamber blocks"))
        }
        if (maximumChamberBlocks != null && chambers.size > maximumChamberBlocks) {
            return ReactorChamberShapeResult(
                zone = null,
                inactiveChamberPositions = emptySet(),
                errors = listOf("reactor multiblock may contain at most $maximumChamberBlocks chamber blocks, got ${chambers.size}"),
            )
        }

        val connectedChambers = connectedChamberComponent(chambers)
        if (connectedChambers.size != chambers.size) {
            val disconnected = (chambers - connectedChambers).sorted()
            return ReactorChamberShapeResult(
                zone = null,
                inactiveChamberPositions = emptySet(),
                errors = listOf("reactor chamber blocks must form one face-connected zone; disconnected chambers: $disconnected"),
            )
        }

        val chamberPositions = chambers.toSortedSet()
        return ReactorChamberShapeResult(
            zone = ReactorZoneDescriptor(
                zoneIndex = 0,
                chamberPositions = chamberPositions,
                volumePositions = chamberPositions,
                volumeCubicMeters = chamberPositions.size * chamberVolumeCubicMeters,
            ),
            inactiveChamberPositions = emptySet(),
            errors = emptyList(),
        )
    }

    private fun connectedChamberComponent(chambers: Set<ReactorBlockPosition>): Set<ReactorBlockPosition> {
        val start = chambers.minOrNull() ?: return emptySet()
        val visited = mutableSetOf<ReactorBlockPosition>()
        val queue = ArrayDeque<ReactorBlockPosition>()
        queue += start
        while (queue.isNotEmpty()) {
            val current = queue.removeFirst()
            if (!visited.add(current)) {
                continue
            }
            for (neighbour in current.faceNeighbours()) {
                if (neighbour in chambers && neighbour !in visited) {
                    queue += neighbour
                }
            }
        }
        return visited
    }
}
