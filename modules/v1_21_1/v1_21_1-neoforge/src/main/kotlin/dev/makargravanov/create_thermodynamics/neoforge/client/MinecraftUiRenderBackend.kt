package dev.makargravanov.create_thermodynamics.neoforge.client

import net.minecraft.client.gui.Font
import net.minecraft.client.gui.GuiGraphics
import ru.lazyhat.kraftui.editor.EditorViewModel
import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.program.RenderBackend
import kotlin.math.max
import kotlin.math.min

class MinecraftUiRenderBackend(
    private val graphics: GuiGraphics,
    private val font: Font,
) : RenderBackend {
    private data class ClipRect(
        val x: Int,
        val y: Int,
        val right: Int,
        val bottom: Int,
    )

    private val clips = ArrayDeque<ClipRect>()

    override fun fillRect(x: Int, y: Int, width: Int, height: Int, color: Color) {
        graphics.fill(x, y, x + width, y + height, color.value.toInt())
    }

    override fun drawText(x: Int, y: Int, text: String, color: Color) {
        graphics.drawString(font, text, x, y, color.value.toInt(), false)
    }

    override fun drawTerminalSurface(x: Int, y: Int, snapshot: Any) {
        throw UnsupportedOperationException("Terminal surfaces are not supported by reactor screens")
    }

    override fun pushClip(x: Int, y: Int, width: Int, height: Int) {
        val requested = ClipRect(x, y, x + width, y + height)
        val next =
            clips.lastOrNull()?.let { previous ->
                ClipRect(
                    x = max(previous.x, requested.x),
                    y = max(previous.y, requested.y),
                    right = min(previous.right, requested.right),
                    bottom = min(previous.bottom, requested.bottom),
                )
            } ?: requested
        clips.addLast(next)
        graphics.enableScissor(next.x, next.y, next.right, next.bottom)
    }

    override fun popClip() {
        check(clips.isNotEmpty()) { "Cannot pop UI clip: clip stack is empty" }
        clips.removeLast()
        val previous = clips.lastOrNull()
        if (previous == null) {
            graphics.disableScissor()
        } else {
            graphics.enableScissor(previous.x, previous.y, previous.right, previous.bottom)
        }
    }

    override fun drawCodeEditor(
        x: Int,
        y: Int,
        width: Int,
        height: Int,
        viewModel: EditorViewModel,
        fontWidth: Int,
        fontHeight: Int,
    ) {
        throw UnsupportedOperationException("Code editors are not supported by reactor screens")
    }

    override fun measureText(text: String): Int =
        font.width(text)
}
