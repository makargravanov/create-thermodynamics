package dev.makargravanov.create_thermodynamics.neoforge.client.ui

import net.minecraft.client.gui.Font
import net.minecraft.client.gui.GuiGraphics
import ru.lazyhat.kraftui.editor.EditorViewModel
import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.program.RenderBackend
import java.util.ArrayDeque

class MinecraftUiDslRenderBackend(
    private val graphics: GuiGraphics,
    private val font: Font,
    private val originX: Int,
    private val originY: Int,
) : RenderBackend {
    private val clipStack = ArrayDeque<ClipRect>()

    override fun fillRect(
        x: Int,
        y: Int,
        width: Int,
        height: Int,
        color: Color,
    ) {
        graphics.fill(
            originX + x,
            originY + y,
            originX + x + width,
            originY + y + height,
            color.value.toInt(),
        )
    }

    override fun drawText(
        x: Int,
        y: Int,
        text: String,
        color: Color,
    ) {
        graphics.drawString(font, text, originX + x, originY + y, color.value.toInt(), false)
    }

    override fun drawTerminalSurface(
        x: Int,
        y: Int,
        snapshot: Any,
    ) {
        throw UnsupportedOperationException("Terminal surfaces are not supported by Minecraft reactor screens")
    }

    override fun pushClip(
        x: Int,
        y: Int,
        width: Int,
        height: Int,
    ) {
        val requested = ClipRect(originX + x, originY + y, width, height)
        val next = clipStack.peekLast()?.intersect(requested) ?: requested
        clipStack.addLast(next)
        applyCurrentClip()
    }

    override fun popClip() {
        check(clipStack.isNotEmpty()) { "Cannot pop Minecraft UI clip: clip stack is empty" }
        clipStack.removeLast()
        applyCurrentClip()
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
        throw UnsupportedOperationException("Code editors are not supported by Minecraft reactor screens")
    }

    override fun measureText(text: String): Int =
        font.width(text)

    private fun applyCurrentClip() {
        graphics.disableScissor()
        val clip = clipStack.peekLast() ?: return
        graphics.enableScissor(clip.x, clip.y, clip.x + clip.width, clip.y + clip.height)
    }

    private data class ClipRect(
        val x: Int,
        val y: Int,
        val width: Int,
        val height: Int,
    ) {
        fun intersect(other: ClipRect): ClipRect {
            val left = maxOf(x, other.x)
            val top = maxOf(y, other.y)
            val right = minOf(x + width, other.x + other.width)
            val bottom = minOf(y + height, other.y + other.height)
            return ClipRect(left, top, maxOf(0, right - left), maxOf(0, bottom - top))
        }
    }
}
