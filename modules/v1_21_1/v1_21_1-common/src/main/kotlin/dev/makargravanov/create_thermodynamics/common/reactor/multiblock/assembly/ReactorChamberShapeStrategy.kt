package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorZoneDescriptor

interface ReactorChamberShapeStrategy {
    fun buildZone(
        chambers: Set<ReactorBlockPosition>,
        volumeCapableBlocks: Map<ReactorBlockPosition, ReactorMultiblockBlockKind>,
        controller: ReactorBlockPosition?,
        chamberVolumeCubicMeters: Double,
        maximumChamberBlocks: Int?,
    ): ReactorChamberShapeResult
}

data class ReactorChamberShapeResult(
    val zone: ReactorZoneDescriptor?,
    val inactiveChamberPositions: Set<ReactorBlockPosition> = emptySet(),
    val errors: List<String>,
)
