package dev.makargravanov.create_thermodynamics.neoforge.client

import com.simibubi.create.foundation.gui.menu.AbstractSimiContainerScreen
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorPortMenu
import net.minecraft.client.gui.GuiGraphics
import net.minecraft.network.chat.Component
import net.minecraft.resources.ResourceLocation
import net.minecraft.world.entity.player.Inventory

class ReactorPortScreen(
    menu: ReactorPortMenu,
    playerInventory: Inventory,
    title: Component,
) : AbstractSimiContainerScreen<ReactorPortMenu>(menu, playerInventory, title) {
    init {
        setWindowSize(176, 168)
    }

    override fun renderBg(graphics: GuiGraphics, partialTick: Float, mouseX: Int, mouseY: Int) {
        val x = leftPos
        val y = topPos
        val topSectionHeight = PortRows * 18 + 17
        graphics.blit(ContainerTexture, x, y, 0, 0, imageWidth, topSectionHeight)
        graphics.blit(ContainerTexture, x, y + topSectionHeight, 0, 126, imageWidth, 96)
        graphics.drawString(font, title, x + 8, y + 6, LabelColor, false)
        graphics.drawString(
            font,
            playerInventoryTitle,
            x + 8,
            y + imageHeight - 94,
            LabelColor,
            false,
        )
    }

    companion object {
        private const val PortRows = 3
        private const val LabelColor = 0x404040
        private val ContainerTexture: ResourceLocation =
            ResourceLocation.fromNamespaceAndPath("minecraft", "textures/gui/container/generic_54.png")
    }
}
