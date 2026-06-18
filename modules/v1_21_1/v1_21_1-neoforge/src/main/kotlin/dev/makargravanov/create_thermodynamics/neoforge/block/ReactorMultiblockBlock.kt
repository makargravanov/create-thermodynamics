package dev.makargravanov.create_thermodynamics.neoforge.block

import net.minecraft.core.BlockPos
import net.minecraft.core.Direction
import net.minecraft.server.level.ServerLevel
import net.minecraft.world.InteractionResult
import net.minecraft.world.entity.LivingEntity
import net.minecraft.world.entity.player.Player
import net.minecraft.world.item.ItemStack
import net.minecraft.world.item.context.BlockPlaceContext
import net.minecraft.world.level.Level
import net.minecraft.world.level.block.Block
import net.minecraft.world.level.block.EntityBlock
import net.minecraft.world.level.block.entity.BlockEntity
import net.minecraft.world.level.block.state.BlockBehaviour
import net.minecraft.world.level.block.state.BlockState
import net.minecraft.world.level.block.state.StateDefinition
import net.minecraft.world.level.block.state.properties.BlockStateProperties
import net.minecraft.world.level.block.state.properties.DirectionProperty
import net.minecraft.world.phys.BlockHitResult

class ReactorMultiblockBlock(
    properties: BlockBehaviour.Properties,
    val kind: ReactorMultiblockKind,
) : Block(properties), EntityBlock {
    init {
        registerDefaultState(stateDefinition.any().setValue(FACING, Direction.NORTH))
    }

    override fun newBlockEntity(pos: BlockPos, state: BlockState): BlockEntity =
        ReactorMultiblockBlockEntity(pos, state)

    override fun getStateForPlacement(context: BlockPlaceContext): BlockState =
        defaultBlockState().setValue(FACING, if (kind.hasFacing) context.clickedFace else Direction.NORTH)

    fun modelFacing(state: BlockState): Direction? =
        state.getValue(FACING).takeIf { kind.hasFacing }

    override fun useWithoutItem(
        state: BlockState,
        level: Level,
        pos: BlockPos,
        player: Player,
        hitResult: BlockHitResult,
    ): InteractionResult {
        if (!kind.isPort) {
            return super.useWithoutItem(state, level, pos, player, hitResult)
        }
        if (!level.isClientSide) {
            val blockEntity = level.getBlockEntity(pos) as? ReactorMultiblockBlockEntity
                ?: error("reactor port at $pos must have ReactorMultiblockBlockEntity")
            player.openMenu(blockEntity)
        }
        return InteractionResult.sidedSuccess(level.isClientSide)
    }

    override fun onPlace(state: BlockState, level: Level, pos: BlockPos, oldState: BlockState, movedByPiston: Boolean) {
        super.onPlace(state, level, pos, oldState, movedByPiston)
        if (level is ServerLevel && state.block != oldState.block) {
            rebuildIfBlockEntityIsReady(level, pos)
        }
    }

    override fun setPlacedBy(level: Level, pos: BlockPos, state: BlockState, placer: LivingEntity?, stack: ItemStack) {
        super.setPlacedBy(level, pos, state, placer, stack)
        if (level is ServerLevel) {
            rebuildIfBlockEntityIsReady(level, pos)
        }
    }

    override fun onRemove(state: BlockState, level: Level, pos: BlockPos, newState: BlockState, movedByPiston: Boolean) {
        super.onRemove(state, level, pos, newState, movedByPiston)
        if (level is ServerLevel && state.block != newState.block) {
            ReactorMultiblockWorldAssembler.clearMembership(level, pos)
            for (direction in Direction.entries) {
                ReactorMultiblockWorldAssembler.rebuildAround(level, pos.relative(direction))
            }
        }
    }

    override fun neighborChanged(
        state: BlockState,
        level: Level,
        pos: BlockPos,
        neighborBlock: Block,
        neighborPos: BlockPos,
        movedByPiston: Boolean,
    ) {
        super.neighborChanged(state, level, pos, neighborBlock, neighborPos, movedByPiston)
        if (level is ServerLevel) {
            ReactorMultiblockWorldAssembler.rebuildAround(level, pos)
        }
    }

    private fun rebuildIfBlockEntityIsReady(level: ServerLevel, pos: BlockPos) {
        if (level.getBlockEntity(pos) is ReactorMultiblockBlockEntity) {
            ReactorMultiblockWorldAssembler.rebuildAround(level, pos)
        }
    }

    override fun createBlockStateDefinition(builder: StateDefinition.Builder<Block, BlockState>) {
        builder.add(FACING)
    }

    companion object {
        val FACING: DirectionProperty = BlockStateProperties.FACING
    }
}
