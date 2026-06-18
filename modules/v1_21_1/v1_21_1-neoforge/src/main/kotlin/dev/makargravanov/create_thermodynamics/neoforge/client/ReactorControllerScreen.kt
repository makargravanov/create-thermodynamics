package dev.makargravanov.create_thermodynamics.neoforge.client

import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorControllerMenu
import dev.makargravanov.create_thermodynamics.neoforge.client.generated.GeneratedReactorControllerProgram
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerUiState
import net.minecraft.client.gui.GuiGraphics
import net.minecraft.client.gui.screens.inventory.AbstractContainerScreen
import net.minecraft.network.chat.Component
import net.minecraft.world.entity.player.Inventory
import ru.lazyhat.kraftui.program.ScreenRuntimeExecutor

class ReactorControllerScreen(
    menu: ReactorControllerMenu,
    playerInventory: Inventory,
    title: Component,
) : AbstractContainerScreen<ReactorControllerMenu>(menu, playerInventory, title) {
    private var executor: ScreenRuntimeExecutor? = null

    init {
        imageWidth = GeneratedReactorControllerProgram.Width
        imageHeight = GeneratedReactorControllerProgram.Height
    }

    override fun init() {
        super.init()
        rebuildUi()
    }

    override fun render(graphics: GuiGraphics, mouseX: Int, mouseY: Int, partialTick: Float) {
        renderBackground(graphics, mouseX, mouseY, partialTick)
        val currentExecutor = executor ?: rebuildUi()
        currentExecutor.updateMouse(mouseX, mouseY)
        currentExecutor.render(MinecraftUiRenderBackend(graphics, font))
        currentExecutor.activeTooltip?.let { tooltip ->
            graphics.renderTooltip(font, Component.literal(tooltip), mouseX, mouseY)
        }
    }

    override fun renderBg(graphics: GuiGraphics, partialTick: Float, mouseX: Int, mouseY: Int) {
    }

    override fun renderLabels(graphics: GuiGraphics, mouseX: Int, mouseY: Int) {
    }

    override fun mouseClicked(mouseX: Double, mouseY: Double, button: Int): Boolean =
        executor?.mouseClicked(mouseX.toInt(), mouseY.toInt()) == true || super.mouseClicked(mouseX, mouseY, button)

    override fun mouseDragged(
        mouseX: Double,
        mouseY: Double,
        button: Int,
        dragX: Double,
        dragY: Double,
    ): Boolean =
        executor?.mouseDragged(mouseX.toInt(), mouseY.toInt()) == true ||
            super.mouseDragged(mouseX, mouseY, button, dragX, dragY)

    override fun mouseReleased(mouseX: Double, mouseY: Double, button: Int): Boolean =
        executor?.mouseReleased(mouseX.toInt(), mouseY.toInt()) == true || super.mouseReleased(mouseX, mouseY, button)

    override fun mouseScrolled(mouseX: Double, mouseY: Double, scrollX: Double, scrollY: Double): Boolean =
        executor?.mouseScrolled(mouseX.toInt(), mouseY.toInt(), scrollY) == true ||
            super.mouseScrolled(mouseX, mouseY, scrollX, scrollY)

    override fun keyPressed(keyCode: Int, scanCode: Int, modifiers: Int): Boolean =
        executor?.keyPressed(keyCode, modifiers) == true || super.keyPressed(keyCode, scanCode, modifiers)

    override fun keyReleased(keyCode: Int, scanCode: Int, modifiers: Int): Boolean =
        executor?.keyReleased(keyCode) == true || super.keyReleased(keyCode, scanCode, modifiers)

    override fun charTyped(codePoint: Char, modifiers: Int): Boolean =
        executor?.charTyped(codePoint) == true || super.charTyped(codePoint, modifiers)

    private fun rebuildUi(): ScreenRuntimeExecutor {
        val program =
            GeneratedReactorControllerProgram.create(
                rootX = leftPos,
                rootY = topPos,
                state = ::uiState,
            )
        val previousFocus = executor?.focusedNodeId
        return ScreenRuntimeExecutor(program).also { newExecutor ->
            newExecutor.restoreFocus(previousFocus)
            executor = newExecutor
        }
    }

    private fun uiState(): ReactorControllerUiState {
        val state = menu.state
        return ReactorControllerUiState(
            title = title.string,
            status = if (state.formed) "formed" else "not formed",
            structureId = state.structureId ?: if (state.formed) "formed" else null,
            active = state.formed,
            zoneCount = state.zoneCount,
            chamberBlocks = state.chamberBlockCount,
            portCount = state.portCount,
        )
    }
}
