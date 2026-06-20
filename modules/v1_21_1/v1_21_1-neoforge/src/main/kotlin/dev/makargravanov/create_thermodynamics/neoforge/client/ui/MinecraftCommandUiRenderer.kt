package dev.makargravanov.create_thermodynamics.neoforge.client.ui

import dev.makargravanov.create_thermodynamics.ui.layout.UiAlignment
import dev.makargravanov.create_thermodynamics.ui.layout.UiDrawCommand
import dev.makargravanov.create_thermodynamics.ui.layout.UiTextMeasurer
import net.minecraft.client.gui.Font
import net.minecraft.client.gui.GuiGraphics
import net.minecraft.network.chat.Component

class MinecraftUiTextMeasurer(
    private val font: Font,
) : UiTextMeasurer {
    override fun width(text: String): Int =
        font.width(text)

    override fun lineHeight(): Int =
        font.lineHeight
}

object MinecraftCommandUiRenderer {
    fun render(
        graphics: GuiGraphics,
        font: Font,
        originX: Int,
        originY: Int,
        commands: List<UiDrawCommand>,
    ) {
        for (command in commands) {
            when (command) {
                is UiDrawCommand.DrawPanel -> {
                    val bounds = command.bounds
                    graphics.fill(
                        originX + bounds.x,
                        originY + bounds.y,
                        originX + bounds.right,
                        originY + bounds.bottom,
                        command.color,
                    )
                }
                is UiDrawCommand.DrawText -> {
                    val bounds = command.bounds
                    val textX =
                        when (command.alignment) {
                            UiAlignment.Start -> bounds.x
                            UiAlignment.Center -> bounds.x + (bounds.width - font.width(command.text)) / 2
                            UiAlignment.End -> bounds.right - font.width(command.text)
                        }
                    graphics.drawString(
                        font,
                        command.text,
                        originX + textX,
                        originY + bounds.y,
                        command.color,
                        false,
                    )
                }
            }
        }
    }

    fun tooltipAt(
        commands: List<UiDrawCommand>,
        originX: Int,
        originY: Int,
        mouseX: Int,
        mouseY: Int,
    ): Component? {
        val command =
            commands
                .asReversed()
                .filterIsInstance<UiDrawCommand.DrawText>()
                .firstOrNull { drawText ->
                    val bounds = drawText.bounds
                    drawText.tooltip != null &&
                        mouseX >= originX + bounds.x &&
                        mouseX < originX + bounds.right &&
                        mouseY >= originY + bounds.y &&
                        mouseY < originY + bounds.bottom
                }
        return command?.tooltip?.let(Component::literal)
    }
}
