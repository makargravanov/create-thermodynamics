package dev.makargravanov.create_thermodynamics.ui.preview

import ru.lazyhat.kraftui.foundation.UiElement
import ru.lazyhat.kraftui.program.FontMetrics
import ru.lazyhat.kraftui.program.ScreenProgramCompiler
import ru.lazyhat.kraftui.program.ScreenRuntimeExecutor
import java.awt.Font
import java.awt.image.BufferedImage
import java.nio.file.Files
import java.nio.file.Path
import javax.imageio.ImageIO

data class UiPreviewSpec(
    val id: String,
    val width: Int,
    val height: Int,
    val root: UiElement,
)

object UiPreviewRenderer {
    private val previewFont = Font(Font.MONOSPACED, Font.PLAIN, 11)

    fun render(spec: UiPreviewSpec): BufferedImage {
        require(spec.width > 0) { "UI preview width must be positive: ${spec.id}" }
        require(spec.height > 0) { "UI preview height must be positive: ${spec.id}" }

        val image = BufferedImage(spec.width, spec.height, BufferedImage.TYPE_INT_ARGB)
        val backend = ImageRenderBackend(image, font = previewFont)

        val program =
            ScreenProgramCompiler(fontMetrics = FontMetrics { text -> estimateTextWidth(text) })
                .compile(
                    root = spec.root,
                    rootWidth = spec.width,
                    rootHeight = spec.height,
                )
        ScreenRuntimeExecutor(program).render(backend)
        backend.close()
        return image
    }

    fun renderAll(previews: Iterable<UiPreviewSpec>, outputDirectory: Path): List<Path> {
        Files.createDirectories(outputDirectory)
        return previews.map { preview ->
            val output = outputDirectory.resolve("${preview.id}.png")
            ImageIO.write(render(preview), "png", output.toFile())
            output
        }
    }

    private fun estimateTextWidth(text: String): Int {
        val image = BufferedImage(1, 1, BufferedImage.TYPE_INT_ARGB)
        val graphics = image.createGraphics()
        try {
            graphics.font = previewFont
            return graphics.fontMetrics.stringWidth(text)
        } finally {
            graphics.dispose()
        }
    }
}
