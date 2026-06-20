package dev.makargravanov.create_thermodynamics.neoforge.client

import com.mojang.blaze3d.vertex.PoseStack
import com.simibubi.create.foundation.blockEntity.behaviour.filtering.FilteringRenderer
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlockEntity
import net.minecraft.client.renderer.LightTexture
import net.minecraft.client.renderer.MultiBufferSource
import net.minecraft.client.renderer.blockentity.BlockEntityRenderer
import net.minecraft.client.renderer.blockentity.BlockEntityRendererProvider

class ReactorMultiblockBlockEntityRenderer(
    @Suppress("UNUSED_PARAMETER")
    context: BlockEntityRendererProvider.Context,
) : BlockEntityRenderer<ReactorMultiblockBlockEntity> {
    override fun render(
        blockEntity: ReactorMultiblockBlockEntity,
        partialTick: Float,
        poseStack: PoseStack,
        bufferSource: MultiBufferSource,
        packedLight: Int,
        packedOverlay: Int,
    ) {
        FilteringRenderer.renderOnBlockEntity(
            blockEntity,
            partialTick,
            poseStack,
            bufferSource,
            LightTexture.FULL_BRIGHT,
            packedOverlay,
        )
    }
}
