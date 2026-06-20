package dev.makargravanov.create_thermodynamics.neoforge.client

import com.simibubi.create.foundation.gui.menu.AbstractSimiContainerScreen
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerFormationState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerViewState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorZoneViewState
import dev.makargravanov.create_thermodynamics.neoforge.CreateThermodynamicsMod
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorControllerMenu
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlockEntity
import net.minecraft.client.gui.GuiGraphics
import net.minecraft.network.chat.Component
import net.minecraft.resources.ResourceLocation
import net.minecraft.world.entity.player.Inventory
import java.util.Locale

class ReactorControllerScreen(
    menu: ReactorControllerMenu,
    playerInventory: Inventory,
    title: Component,
) : AbstractSimiContainerScreen<ReactorControllerMenu>(menu, playerInventory, title) {
    private var selectedTab: ControllerTab = ControllerTab.OVERVIEW
    private var selectedZoneIndex: Int = 0

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

        renderTabs(graphics, x, y)
        when (selectedTab) {
            ControllerTab.OVERVIEW -> renderOverview(graphics, x, y, state)
            ControllerTab.ZONES -> renderZones(graphics, x, y, state)
            ControllerTab.MIXTURE -> renderMixture(graphics, x, y, state)
        }
    }

    private fun renderTabs(graphics: GuiGraphics, x: Int, y: Int) {
        for ((index, tab) in ControllerTab.entries.withIndex()) {
            val tabX = x + 16 + index * 66
            val tabY = y + 29
            val selected = tab == selectedTab
            graphics.fill(tabX, tabY, tabX + 58, tabY + 16, if (selected) CreatePanelColor else CreateTabColor)
            graphics.drawString(font, tab.label, tabX + 6, tabY + 4, if (selected) CreateTextColor else CreateMutedColor, false)
        }
    }

    private fun renderOverview(graphics: GuiGraphics, x: Int, y: Int, state: ReactorControllerViewState) {
        val zone = selectedZone(state)
        renderKeyValue(graphics, "State", state.formationState.label(), x + 18, y + 53)
        renderKeyValue(graphics, "Native", state.nativeBinding, x + 96, y + 53)
        renderKeyValue(graphics, "Zones", state.zoneCount.toString(), x + 174, y + 53)

        graphics.drawString(font, "Structure", x + 18, y + 80, CreateMutedColor, false)
        graphics.drawString(font, "blocks ${state.chamberBlockCount}", x + 18, y + 92, CreateTextColor, false)
        graphics.drawString(font, "ports ${state.portCount}", x + 18, y + 103, CreateTextColor, false)

        graphics.drawString(font, "Selected zone", x + 96, y + 80, CreateMutedColor, false)
        graphics.drawString(font, zone?.let { "zone ${it.index}" } ?: "no metrics", x + 96, y + 92, CreateTextColor, false)
        graphics.drawString(font, zone?.let { formatTemperature(it.temperatureKelvin) } ?: "temperature n/a", x + 96, y + 103, CreateTextColor, false)
        graphics.drawString(font, zone?.let { formatPressure(it.pressurePascal) } ?: "pressure n/a", x + 96, y + 114, CreateTextColor, false)

        graphics.drawString(font, "Mixture", x + 174, y + 80, CreateMutedColor, false)
        graphics.drawString(font, zone?.let { "${it.mixture.size} substances" } ?: "no data", x + 174, y + 92, CreateTextColor, false)
        graphics.drawString(font, zone?.let { "${formatNumber(it.totalConcentrationMolPerBucket)} mol/b" } ?: "", x + 174, y + 103, CreateTextColor, false)
    }

    private fun renderZones(graphics: GuiGraphics, x: Int, y: Int, state: ReactorControllerViewState) {
        graphics.drawString(font, "Zones", x + 18, y + 53, CreateMutedColor, false)
        val zones = state.zones.sortedBy { it.index }
        if (zones.isEmpty()) {
            graphics.drawString(font, "no native metrics yet", x + 18, y + 66, CreateTextColor, false)
            return
        }

        for ((row, zone) in zones.take(MaxVisibleZoneRows).withIndex()) {
            val rowY = y + 66 + row * 13
            val selected = zone.index == selectedZoneIndex
            if (selected) {
                graphics.fill(x + 16, rowY - 2, x + 68, rowY + 10, CreatePanelColor)
            }
            graphics.drawString(font, "zone ${zone.index}", x + 20, rowY, if (selected) CreateTextColor else CreateMutedColor, false)
        }

        val selected = selectedZone(state)
        if (selected == null) {
            graphics.drawString(font, "select a zone", x + 88, y + 66, CreateTextColor, false)
            return
        }

        graphics.drawString(font, "Zone ${selected.index}", x + 88, y + 53, CreateMutedColor, false)
        graphics.drawString(font, formatTemperature(selected.temperatureKelvin), x + 88, y + 66, CreateTextColor, false)
        graphics.drawString(font, formatPressure(selected.pressurePascal), x + 88, y + 78, CreateTextColor, false)
        graphics.drawString(font, "${selected.mixture.size} substances", x + 88, y + 90, CreateTextColor, false)
        graphics.drawString(font, "${formatNumber(selected.totalConcentrationMolPerBucket)} mol/b total", x + 88, y + 102, CreateTextColor, false)
    }

    private fun renderMixture(graphics: GuiGraphics, x: Int, y: Int, state: ReactorControllerViewState) {
        val zone = selectedZone(state)
        graphics.drawString(font, zone?.let { "Mixture: zone ${it.index}" } ?: "Mixture", x + 18, y + 53, CreateMutedColor, false)
        if (zone == null) {
            graphics.drawString(font, "no native metrics yet", x + 18, y + 66, CreateTextColor, false)
            return
        }
        if (zone.mixture.isEmpty()) {
            graphics.drawString(font, "empty", x + 18, y + 66, CreateTextColor, false)
            return
        }

        graphics.drawString(font, "Substance", x + 18, y + 66, CreateMutedColor, false)
        graphics.drawString(font, "mol/b", x + 160, y + 66, CreateMutedColor, false)

        val rows = zone.mixture.take(MaxVisibleMixtureRows)
        for ((index, entry) in rows.withIndex()) {
            val rowY = y + 79 + index * 10
            graphics.drawString(font, shortSubstanceId(entry.substanceId), x + 18, rowY, CreateTextColor, false)
            graphics.drawString(font, formatNumber(entry.concentrationMolPerBucket), x + 160, rowY, CreateTextColor, false)
        }

        val hidden = zone.mixture.size - rows.size
        if (hidden > 0) {
            graphics.drawString(font, "+$hidden more", x + 18, y + 122, CreateMutedColor, false)
        }
    }

    override fun renderForeground(graphics: GuiGraphics, mouseX: Int, mouseY: Int, partialTicks: Float) {
        super.renderForeground(graphics, mouseX, mouseY, partialTicks)
        val state = controllerState()
        if (state.formationState != ReactorControllerFormationState.FORMED && mouseX in leftPos + 18 until leftPos + 75 && mouseY in topPos + 32 until topPos + 56) {
            graphics.renderTooltip(font, Component.literal(state.diagnostic ?: "Reactor structure is not formed"), mouseX, mouseY)
        }
    }

    override fun mouseClicked(mouseX: Double, mouseY: Double, button: Int): Boolean {
        if (button == 0) {
            val tab = tabAt(mouseX.toInt(), mouseY.toInt())
            if (tab != null) {
                selectedTab = tab
                return true
            }
            if (selectedTab == ControllerTab.ZONES) {
                val zone = zoneAt(mouseX.toInt(), mouseY.toInt(), controllerState())
                if (zone != null) {
                    selectedZoneIndex = zone.index
                    return true
                }
            }
        }
        return super.mouseClicked(mouseX, mouseY, button)
    }

    private fun controllerState(): ReactorControllerViewState =
        clientControllerBlockEntity()?.controllerScreenState() ?: menu.state

    private fun selectedZone(state: ReactorControllerViewState): ReactorZoneViewState? =
        state.zones.firstOrNull { it.index == selectedZoneIndex }
            ?: state.zones.minByOrNull { it.index }?.also { selectedZoneIndex = it.index }

    private fun tabAt(mouseX: Int, mouseY: Int): ControllerTab? {
        for ((index, tab) in ControllerTab.entries.withIndex()) {
            val tabX = leftPos + 16 + index * 66
            val tabY = topPos + 29
            if (mouseX in tabX until tabX + 58 && mouseY in tabY until tabY + 16) {
                return tab
            }
        }
        return null
    }

    private fun zoneAt(mouseX: Int, mouseY: Int, state: ReactorControllerViewState): ReactorZoneViewState? {
        val zones = state.zones.sortedBy { it.index }.take(MaxVisibleZoneRows)
        for ((row, zone) in zones.withIndex()) {
            val rowY = topPos + 66 + row * 13
            if (mouseX in leftPos + 16 until leftPos + 68 && mouseY in rowY - 2 until rowY + 10) {
                return zone
            }
        }
        return null
    }

    private fun clientControllerBlockEntity(): ReactorMultiblockBlockEntity? {
        val position = menu.controllerPos ?: return null
        return minecraft?.level?.getBlockEntity(position) as? ReactorMultiblockBlockEntity
    }

    private fun renderKeyValue(graphics: GuiGraphics, key: String, value: String, x: Int, y: Int) {
        graphics.drawString(font, key, x, y, CreateMutedColor, false)
        graphics.drawString(font, value, x, y + 11, CreateTextColor, false)
    }

    private fun ReactorControllerFormationState.label(): String =
        when (this) {
            ReactorControllerFormationState.FORMED -> "formed"
            ReactorControllerFormationState.NOT_FORMED -> "not formed"
            ReactorControllerFormationState.UNKNOWN -> "unknown"
        }

    private fun shortSubstanceId(substanceId: String): String =
        substanceId.substringAfter(':').replace('_', ' ').take(22)

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
        private const val CreatePanelColor = 0xFFC9B28D.toInt()
        private const val CreateTabColor = 0xFFB69673.toInt()
        private const val MaxVisibleZoneRows = 5
        private const val MaxVisibleMixtureRows = 5
        private const val BackgroundWidth = 232
        private const val BackgroundHeight = 140
        private val BackgroundTexture: ResourceLocation =
            ResourceLocation.fromNamespaceAndPath(CreateThermodynamicsMod.MOD_ID, "textures/gui/reactor_controller.png")
    }

    private enum class ControllerTab(val label: String) {
        OVERVIEW("Overview"),
        ZONES("Zones"),
        MIXTURE("Mixture"),
    }
}
