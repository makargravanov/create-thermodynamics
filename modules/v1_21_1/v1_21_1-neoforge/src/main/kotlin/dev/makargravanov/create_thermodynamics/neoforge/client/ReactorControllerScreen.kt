package dev.makargravanov.create_thermodynamics.neoforge.client

import com.simibubi.create.foundation.gui.menu.AbstractSimiContainerScreen
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerFormationState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerViewState
import dev.makargravanov.create_thermodynamics.neoforge.CreateThermodynamicsMod
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorControllerMenu
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlockEntity
import dev.makargravanov.create_thermodynamics.neoforge.registry.CreateThermodynamicsRegistries
import net.createmod.catnip.gui.element.GuiGameElement
import net.minecraft.client.gui.GuiGraphics
import net.minecraft.network.chat.Component
import net.minecraft.resources.ResourceLocation
import net.minecraft.world.entity.player.Inventory
import net.minecraft.world.item.ItemStack
import java.util.Locale

class ReactorControllerScreen(
    menu: ReactorControllerMenu,
    playerInventory: Inventory,
    title: Component,
) : AbstractSimiContainerScreen<ReactorControllerMenu>(menu, playerInventory, title) {
    override fun init() {
        setWindowSize(BackgroundWidth, BackgroundHeight)
        setWindowOffset(-16, 0)
        super.init()
        clearWidgets()
    }

    override fun renderBg(graphics: GuiGraphics, partialTick: Float, mouseX: Int, mouseY: Int) {
        val x = leftPos
        val y = topPos
        val state = controllerState()

        graphics.blit(BackgroundTexture, x, y, 0.0f, 0.0f, BackgroundWidth, BackgroundHeight, BackgroundWidth, BackgroundHeight)
        graphics.drawString(font, title, x + BackgroundWidth / 2 - font.width(title) / 2, y + 10, CreateTitleColor, false)

        renderStatus(graphics, x, y, state)
        renderStructure(graphics, x, y, state)
        renderZone(graphics, x, y, state)
        renderMixture(graphics, x, y, state)

        GuiGameElement.of(ItemStack(CreateThermodynamicsRegistries.reactorControllerItem.get()))
            .scale(3.0)
            .render(graphics, x + BackgroundWidth + 16, y + BackgroundHeight - 34)
    }

    private fun renderStatus(graphics: GuiGraphics, x: Int, y: Int, state: ReactorControllerViewState) {
        graphics.drawString(font, "State", x + 18, y + 33, CreateMutedColor, false)
        graphics.drawString(
            font,
            state.formationState.label(),
            x + 18,
            y + 43,
            CreateTextColor,
            false,
        )

        graphics.drawString(font, "Native", x + 120, y + 33, CreateMutedColor, false)
        graphics.drawString(font, state.nativeBinding, x + 120, y + 43, CreateTextColor, false)
    }

    private fun renderStructure(graphics: GuiGraphics, x: Int, y: Int, state: ReactorControllerViewState) {
        graphics.drawString(font, "Structure", x + 18, y + 72, CreateMutedColor, false)
        graphics.drawString(font, "zones ${state.zoneCount}", x + 18, y + 83, CreateTextColor, false)
        graphics.drawString(font, "blocks ${state.chamberBlockCount}", x + 18, y + 94, CreateTextColor, false)
        graphics.drawString(font, "ports ${state.portCount}", x + 18, y + 105, CreateTextColor, false)
    }

    private fun renderZone(graphics: GuiGraphics, x: Int, y: Int, state: ReactorControllerViewState) {
        graphics.drawString(font, "Zone 0", x + 92, y + 72, CreateMutedColor, false)
        graphics.drawString(font, formatTemperature(state.temperatureKelvin), x + 92, y + 83, CreateTextColor, false)
        graphics.drawString(font, formatPressure(state.pressurePascal), x + 92, y + 94, CreateTextColor, false)
        graphics.drawString(font, state.structureId?.value?.toString()?.take(8) ?: "no structure", x + 92, y + 105, CreateMutedColor, false)
    }

    private fun renderMixture(graphics: GuiGraphics, x: Int, y: Int, state: ReactorControllerViewState) {
        graphics.drawString(font, "Mixture", x + 154, y + 72, CreateMutedColor, false)
        if (state.mixture.isEmpty()) {
            graphics.drawString(font, "empty", x + 154, y + 83, CreateTextColor, false)
            return
        }

        for ((index, entry) in state.mixture.take(MaxVisibleMixtureRows).withIndex()) {
            val line = "${shortSubstanceId(entry.substanceId)} ${formatNumber(entry.concentrationMolPerBucket)}"
            graphics.drawString(font, line.take(MaxMixtureLineLength), x + 154, y + 83 + index * 11, CreateTextColor, false)
        }
    }

    override fun renderForeground(graphics: GuiGraphics, mouseX: Int, mouseY: Int, partialTicks: Float) {
        super.renderForeground(graphics, mouseX, mouseY, partialTicks)
        val state = controllerState()
        if (state.formationState != ReactorControllerFormationState.FORMED && mouseX in leftPos + 18 until leftPos + 75 && mouseY in topPos + 32 until topPos + 56) {
            graphics.renderTooltip(font, Component.literal(state.diagnostic ?: "Reactor structure is not formed"), mouseX, mouseY)
        }
    }

    private fun controllerState(): ReactorControllerViewState =
        clientControllerBlockEntity()?.controllerScreenState() ?: menu.state

    private fun clientControllerBlockEntity(): ReactorMultiblockBlockEntity? {
        val position = menu.controllerPos ?: return null
        return minecraft?.level?.getBlockEntity(position) as? ReactorMultiblockBlockEntity
    }

    private fun ReactorControllerFormationState.label(): String =
        when (this) {
            ReactorControllerFormationState.FORMED -> "formed"
            ReactorControllerFormationState.NOT_FORMED -> "not formed"
            ReactorControllerFormationState.UNKNOWN -> "unknown"
        }

    private fun shortSubstanceId(substanceId: String): String =
        substanceId.substringAfter(':').replace('_', ' ').take(11)

    private fun formatTemperature(value: Double?): String =
        value?.let { "${formatNumber(it)} K" } ?: "temperature n/a"

    private fun formatPressure(value: Double?): String =
        value?.let { "${formatNumber(it / 1000.0)} kPa" } ?: "pressure n/a"

    private fun formatNumber(value: Double): String =
        when {
            value == 0.0 -> "0"
            kotlin.math.abs(value) >= 1000.0 -> String.format(Locale.ROOT, "%.0f", value)
            kotlin.math.abs(value) >= 10.0 -> String.format(Locale.ROOT, "%.1f", value)
            else -> String.format(Locale.ROOT, "%.3f", value)
        }

    companion object {
        private const val CreateTitleColor = 0x3D3C48
        private const val CreateTextColor = 0x3D3C48
        private const val CreateMutedColor = 0x6F6A75
        private const val MaxVisibleMixtureRows = 4
        private const val MaxMixtureLineLength = 13
        private const val BackgroundWidth = 232
        private const val BackgroundHeight = 140
        private val BackgroundTexture: ResourceLocation =
            ResourceLocation.fromNamespaceAndPath(CreateThermodynamicsMod.MOD_ID, "textures/gui/reactor_controller.png")
    }
}
