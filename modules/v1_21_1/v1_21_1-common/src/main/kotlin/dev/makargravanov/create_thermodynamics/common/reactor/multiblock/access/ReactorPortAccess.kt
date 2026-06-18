package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId

interface ReactorPortAccess {
    fun insertItem(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        itemId: String,
        itemCount: Int,
    ): ReactorOperationResult

    fun extractItem(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        itemId: String,
        maxItemCount: Int,
    ): ReactorOperationResult

    fun insertFluid(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        fluidId: String,
        millibuckets: Int,
    ): ReactorOperationResult

    fun extractFluid(
        structureId: ReactorStructureId,
        portPosition: ReactorBlockPosition,
        fluidId: String,
        maxMillibuckets: Int,
    ): ReactorOperationResult
}
