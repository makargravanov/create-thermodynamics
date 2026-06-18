package dev.makargravanov.create_thermodynamics.neoforge.block

import net.minecraft.core.BlockPos
import net.minecraft.core.Direction
import net.minecraft.world.item.context.BlockPlaceContext
import net.minecraft.world.level.BlockGetter
import net.minecraft.world.level.LevelAccessor
import net.minecraft.world.level.block.state.BlockBehaviour
import net.minecraft.world.level.block.state.BlockState
import net.minecraft.world.level.block.state.StateDefinition
import net.minecraft.world.level.block.state.properties.BlockStateProperties
import net.minecraft.world.level.block.state.properties.BooleanProperty

class ReactorChamberBlock(properties: BlockBehaviour.Properties) : ReactorMultiblockBlock(properties) {
    init {
        registerDefaultState(
            stateDefinition.any()
                .setValue(NORTH, false)
                .setValue(EAST, false)
                .setValue(SOUTH, false)
                .setValue(WEST, false)
                .setValue(UP, false)
                .setValue(DOWN, false),
        )
    }

    override fun getStateForPlacement(context: BlockPlaceContext): BlockState =
        connectedState(defaultBlockState(), context.level, context.clickedPos)

    override fun updateShape(
        state: BlockState,
        direction: Direction,
        neighborState: BlockState,
        level: LevelAccessor,
        currentPos: BlockPos,
        neighborPos: BlockPos,
    ): BlockState =
        state.setValue(propertyFor(direction), canConnectTo(neighborState))

    override fun createBlockStateDefinition(builder: StateDefinition.Builder<net.minecraft.world.level.block.Block, BlockState>) {
        builder.add(NORTH, EAST, SOUTH, WEST, UP, DOWN)
    }

    private fun connectedState(state: BlockState, level: BlockGetter, pos: BlockPos): BlockState {
        var connectedState = state
        for (direction in Direction.entries) {
            connectedState = connectedState.setValue(
                propertyFor(direction),
                canConnectTo(level.getBlockState(pos.relative(direction))),
            )
        }
        return connectedState
    }

    private fun canConnectTo(state: BlockState): Boolean =
        state.block is ReactorMultiblockBlock

    private fun propertyFor(direction: Direction): BooleanProperty =
        when (direction) {
            Direction.NORTH -> NORTH
            Direction.EAST -> EAST
            Direction.SOUTH -> SOUTH
            Direction.WEST -> WEST
            Direction.UP -> UP
            Direction.DOWN -> DOWN
        }

    companion object {
        val NORTH: BooleanProperty = BlockStateProperties.NORTH
        val EAST: BooleanProperty = BlockStateProperties.EAST
        val SOUTH: BooleanProperty = BlockStateProperties.SOUTH
        val WEST: BooleanProperty = BlockStateProperties.WEST
        val UP: BooleanProperty = BlockStateProperties.UP
        val DOWN: BooleanProperty = BlockStateProperties.DOWN
    }
}
