package dev.makargravanov.create_thermodynamics.ui.preview

import dev.makargravanov.create_thermodynamics.ui.layout.TextOverflowPolicy
import dev.makargravanov.create_thermodynamics.ui.layout.UiAlignment
import dev.makargravanov.create_thermodynamics.ui.layout.UiDrawCommand
import dev.makargravanov.create_thermodynamics.ui.layout.UiLayoutDiagnostic
import dev.makargravanov.create_thermodynamics.ui.layout.UiRect
import java.awt.Font
import java.awt.RenderingHints
import java.awt.image.BufferedImage
import java.nio.file.Files
import java.nio.file.Path
import javax.imageio.ImageIO

data class CommandPreviewSpec(
    val id: String,
    val width: Int,
    val height: Int,
    val commands: List<UiDrawCommand>,
    val diagnostics: List<UiLayoutDiagnostic>,
)

data class CommandPreviewReport(
    val id: String,
    val diagnostics: List<UiLayoutDiagnostic>,
)

object CommandPreviewRenderer {
    private val previewFont = Font(Font.MONOSPACED, Font.PLAIN, 11)

    fun render(
        width: Int,
        height: Int,
        commands: List<UiDrawCommand>,
    ): BufferedImage {
        require(width > 0) { "UI preview width must be positive" }
        require(height > 0) { "UI preview height must be positive" }
        val image = BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB)
        val graphics = image.createGraphics()
        try {
            graphics.font = previewFont
            graphics.setRenderingHint(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_OFF)
            graphics.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_OFF)
            for (command in commands) {
                when (command) {
                    is UiDrawCommand.DrawPanel -> {
                        graphics.color = command.color.toAwtColor()
                        graphics.fillRect(command.bounds.x, command.bounds.y, command.bounds.width, command.bounds.height)
                    }
                    is UiDrawCommand.DrawText -> {
                        graphics.color = command.color.toAwtColor()
                        val x = command.alignedTextX(graphics.fontMetrics.stringWidth(command.text))
                        graphics.drawString(command.text, x, command.bounds.y + graphics.fontMetrics.ascent)
                    }
                }
            }
        } finally {
            graphics.dispose()
        }
        return image
    }

    fun renderAll(
        previews: Iterable<CommandPreviewSpec>,
        outputDirectory: Path,
    ): List<Path> {
        Files.createDirectories(outputDirectory)
        val outputs = mutableListOf<Path>()
        val reports = mutableListOf<CommandPreviewReport>()
        for (preview in previews) {
            val output = outputDirectory.resolve("${preview.id}.png")
            ImageIO.write(render(preview.width, preview.height, preview.commands), "png", output.toFile())
            outputs.add(output)
            reports.add(CommandPreviewReport(preview.id, preview.diagnostics))
        }
        writeReport(outputDirectory, reports)
        return outputs
    }

    fun writeReport(
        outputDirectory: Path,
        reports: List<CommandPreviewReport>,
    ): Path {
        Files.createDirectories(outputDirectory)
        val output = outputDirectory.resolve("layout-report.json")
        Files.writeString(output, reports.reportsToJson())
        return output
    }

    fun countDistinctColors(image: BufferedImage): Int {
        val colors = HashSet<Int>()
        for (y in 0 until image.height) {
            for (x in 0 until image.width) {
                colors += image.getRGB(x, y)
            }
        }
        return colors.size
    }

    private fun UiDrawCommand.DrawText.alignedTextX(textWidth: Int): Int =
        when (alignment) {
            UiAlignment.Start -> bounds.x
            UiAlignment.Center -> bounds.x + (bounds.width - textWidth) / 2
            UiAlignment.End -> bounds.right - textWidth
        }

    private fun Int.toAwtColor(): java.awt.Color =
        java.awt.Color(this, true)

    private fun List<CommandPreviewReport>.reportsToJson(): String =
        joinToString(prefix = "[", postfix = "]", separator = ",") { report ->
            """{"id":${report.id.json()},"diagnostics":${report.diagnostics.diagnosticsToJson()}}"""
        }

    private fun List<UiLayoutDiagnostic>.diagnosticsToJson(): String =
        joinToString(prefix = "[", postfix = "]", separator = ",") { diagnostic ->
            when (diagnostic) {
                is UiLayoutDiagnostic.TextWouldOverflow ->
                    """{"type":"TextWouldOverflow","nodeId":${diagnostic.nodeId.json()},"text":${diagnostic.text.json()},"rect":${diagnostic.rect.toJson()},"textWidth":${diagnostic.textWidth},"policy":${diagnostic.policy.toJson()}}"""
                is UiLayoutDiagnostic.NodeOutsideParent ->
                    """{"type":"NodeOutsideParent","nodeId":${diagnostic.nodeId.json()},"parent":${diagnostic.parent.toJson()},"child":${diagnostic.child.toJson()}}"""
            }
        }

    private fun TextOverflowPolicy.toJson(): String =
        when (this) {
            TextOverflowPolicy.FailInValidation -> "FailInValidation"
            TextOverflowPolicy.Clip -> "Clip"
            TextOverflowPolicy.Ellipsize -> "Ellipsize"
            TextOverflowPolicy.EllipsizeWithTooltip -> "EllipsizeWithTooltip"
            TextOverflowPolicy.HorizontalScroll -> "HorizontalScroll"
            is TextOverflowPolicy.WrapLines -> "WrapLines($maxLines)"
        }.json()

    private fun UiRect.toJson(): String =
        """{"x":$x,"y":$y,"width":$width,"height":$height}"""

    private fun String.json(): String =
        buildString {
            append('"')
            for (char in this@json) {
                when (char) {
                    '\\' -> append("\\\\")
                    '"' -> append("\\\"")
                    '\n' -> append("\\n")
                    '\r' -> append("\\r")
                    '\t' -> append("\\t")
                    else -> append(char)
                }
            }
            append('"')
        }
}
