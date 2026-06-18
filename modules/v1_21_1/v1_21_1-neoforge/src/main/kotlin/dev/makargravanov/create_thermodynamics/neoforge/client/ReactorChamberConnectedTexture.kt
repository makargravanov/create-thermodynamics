package dev.makargravanov.create_thermodynamics.neoforge.client

import com.simibubi.create.foundation.block.connected.AllCTTypes
import com.simibubi.create.foundation.block.connected.CTModel
import com.simibubi.create.foundation.block.connected.CTSpriteShiftEntry
import com.simibubi.create.foundation.block.connected.CTSpriteShifter
import com.simibubi.create.foundation.block.connected.CTType
import com.simibubi.create.foundation.block.connected.ConnectedTextureBehaviour
import dev.makargravanov.create_thermodynamics.neoforge.CreateThermodynamicsMod
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlock
import net.minecraft.client.renderer.texture.TextureAtlasSprite
import net.minecraft.client.resources.model.ModelResourceLocation
import net.minecraft.core.BlockPos
import net.minecraft.core.Direction
import net.minecraft.resources.ResourceLocation
import net.minecraft.world.level.BlockAndTintGetter
import net.minecraft.world.level.block.Block
import net.minecraft.world.level.block.state.BlockState
import net.neoforged.neoforge.client.event.ModelEvent

object ReactorChamberConnectedTexture {
    private val blockId = id("reactor_chamber")
    private val modelLocation = ModelResourceLocation(blockId, "")
    private val behaviour = ReactorChamberConnectedTextureBehaviour()

    fun onModifyBakingResult(event: ModelEvent.ModifyBakingResult) {
        val models = event.models
        val originalModel = models[modelLocation]
            ?: error("Missing baked model for reactor chamber: $modelLocation")
        models[modelLocation] = CTModel(originalModel, behaviour)
    }

    private fun id(path: String): ResourceLocation =
        ResourceLocation.fromNamespaceAndPath(CreateThermodynamicsMod.MOD_ID, path)

    private class ReactorChamberConnectedTextureBehaviour : ConnectedTextureBehaviour() {
        private val shifts: Map<Direction, CTSpriteShiftEntry> =
            Direction.entries.associateWith { direction ->
                val face = direction.getSerializedName()
                CTSpriteShifter.getCT(
                    AllCTTypes.RECTANGLE,
                    id("block/reactor_chamber_$face"),
                    id("block/reactor_chamber_connected"),
                )
            }

        override fun getShift(
            state: BlockState,
            direction: Direction,
            sprite: TextureAtlasSprite,
        ): CTSpriteShiftEntry? =
            shifts.getValue(direction)

        override fun getDataType(
            world: BlockAndTintGetter,
            pos: BlockPos,
            state: BlockState,
            direction: Direction,
        ): CTType =
            AllCTTypes.RECTANGLE

        override fun connectsTo(
            state: BlockState,
            other: BlockState,
            reader: BlockAndTintGetter,
            pos: BlockPos,
            otherPos: BlockPos,
            face: Direction,
        ): Boolean {
            if (other.block !is ReactorMultiblockBlock) {
                return false
            }

            val nextPosInFrontOfOther = otherPos.relative(face)
            return Block.shouldRenderFace(other, reader, otherPos, face, nextPosInFrontOfOther)
        }
    }
}
