package dev.makargravanov.create_thermodynamics.neoforge.client

import com.simibubi.create.foundation.gui.menu.AbstractSimiContainerScreen
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerFormationState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerViewState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorZoneViewState
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorControllerMenu
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlockEntity
import dev.makargravanov.create_thermodynamics.neoforge.client.ui.MinecraftCommandUiRenderer
import dev.makargravanov.create_thermodynamics.neoforge.client.ui.MinecraftUiTextMeasurer
import dev.makargravanov.create_thermodynamics.ui.layout.UiDrawCommand
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerCommandUi
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerTab
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerUiSnapshot
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorMixtureUiLine
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorZoneUiSnapshot
import net.minecraft.client.gui.GuiGraphics
import net.minecraft.network.chat.Component
import net.minecraft.world.entity.player.Inventory
import java.util.Locale

class ReactorControllerScreen(
    menu: ReactorControllerMenu,
    playerInventory: Inventory,
    title: Component,
) : AbstractSimiContainerScreen<ReactorControllerMenu>(menu, playerInventory, title) {
    private var selectedTab: ReactorControllerTab = ReactorControllerTab.Overview
    private var selectedZoneIndex: Int = 0
    private var lastCommands: List<UiDrawCommand> = emptyList()

    override fun init() {
        setWindowSize(ReactorControllerCommandUi.Width, ReactorControllerCommandUi.Height)
        setWindowOffset(-16, 0)
        super.init()
        clearWidgets()
    }

    override fun renderBg(
        graphics: GuiGraphics,
        partialTick: Float,
        mouseX: Int,
        mouseY: Int,
    ) {
        val snapshot = controllerSnapshot()
        val layout =
            ReactorControllerCommandUi.layout(
                state = snapshot,
                selectedTab = selectedTab,
                selectedZoneIndex = selectedZoneIndex,
                textMeasurer = MinecraftUiTextMeasurer(font),
            )
        lastCommands = layout.commands
        MinecraftCommandUiRenderer.render(graphics, font, leftPos, topPos, layout.commands)
    }

    override fun renderForeground(
        graphics: GuiGraphics,
        mouseX: Int,
        mouseY: Int,
        partialTicks: Float,
    ) {
        super.renderForeground(graphics, mouseX, mouseY, partialTicks)
        MinecraftCommandUiRenderer.tooltipAt(lastCommands, leftPos, topPos, mouseX, mouseY)
            ?.let { graphics.renderTooltip(font, it, mouseX, mouseY) }
    }

    override fun mouseClicked(
        mouseX: Double,
        mouseY: Double,
        button: Int,
    ): Boolean {
        if (button == 0) {
            val localX = mouseX.toInt() - leftPos
            val localY = mouseY.toInt() - topPos
            ReactorControllerCommandUi.tabAt(localX, localY)?.let { tab ->
                selectedTab = tab
                return true
            }
            if (selectedTab == ReactorControllerTab.Zones) {
                ReactorControllerCommandUi.zoneAt(localX, localY, controllerSnapshot())?.let { zoneIndex ->
                    selectedZoneIndex = zoneIndex
                    return true
                }
            }
        }
        return super.mouseClicked(mouseX, mouseY, button)
    }

    private fun controllerSnapshot(): ReactorControllerUiSnapshot {
        val state = controllerState()
        return ReactorControllerUiSnapshot(
            title = title.string,
            status = state.formationState.label(),
            active = state.formationState == ReactorControllerFormationState.FORMED,
            nativeBinding = state.nativeBinding,
            zoneCount = state.zoneCount,
            chamberBlocks = state.chamberBlockCount,
            portCount = state.portCount,
            zones = state.zones.map { it.toUiSnapshot() },
        )
    }

    private fun controllerState(): ReactorControllerViewState =
        clientControllerBlockEntity()?.controllerScreenState() ?: menu.state

    private fun clientControllerBlockEntity(): ReactorMultiblockBlockEntity? {
        val position = menu.controllerPos ?: return null
        return minecraft?.level?.getBlockEntity(position) as? ReactorMultiblockBlockEntity
    }

    private fun ReactorZoneViewState.toUiSnapshot(): ReactorZoneUiSnapshot =
        ReactorZoneUiSnapshot(
            index = index,
            temperature = formatTemperature(temperatureKelvin),
            pressure = formatPressure(pressurePascal),
            mixture =
                mixture.map { entry ->
                    ReactorMixtureUiLine(
                        substanceId = entry.substanceId,
                        concentration = formatNumber(entry.concentrationMolPerBucket),
                    )
                },
        )

    private fun ReactorControllerFormationState.label(): String =
        when (this) {
            ReactorControllerFormationState.FORMED -> "formed"
            ReactorControllerFormationState.NOT_FORMED -> "not formed"
            ReactorControllerFormationState.UNKNOWN -> "unknown"
        }

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
}
