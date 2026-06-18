package dev.makargravanov.create_thermodynamics.neoforge.block

import net.minecraft.core.BlockPos
import net.minecraft.core.Direction
import net.minecraft.world.level.BlockGetter
import net.minecraft.world.level.block.state.BlockBehaviour
import net.minecraft.world.level.block.state.BlockState

class ReactorChamberBlock(properties: BlockBehaviour.Properties) : ReactorMultiblockBlock(properties, ReactorMultiblockKind.CHAMBER) {
    override fun skipRendering(state: BlockState, adjacentBlockState: BlockState, direction: Direction): Boolean =
        super.skipRendering(state, adjacentBlockState, direction)

    override fun useShapeForLightOcclusion(state: BlockState): Boolean = false

    override fun propagatesSkylightDown(state: BlockState, level: BlockGetter, pos: BlockPos): Boolean = false

}
