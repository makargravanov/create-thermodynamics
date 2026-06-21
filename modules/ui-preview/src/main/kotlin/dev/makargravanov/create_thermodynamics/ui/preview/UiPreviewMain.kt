package dev.makargravanov.create_thermodynamics.ui.preview

import ru.lazyhat.kraftui.preview.MinecraftBitmapFont
import ru.lazyhat.kraftui.preview.UiPreviewRenderer
import java.nio.file.Path

fun main(args: Array<String>) {
    val outputDirectory =
        if (args.isNotEmpty()) {
            Path.of(args[0])
        } else {
            Path.of("build", "reports", "ui")
        }
    val renderer = UiPreviewRenderer(MinecraftBitmapFont.loadFromMinecraftClientJar())
    val outputs = renderer.renderAll(ReactorPortPreviews.all(), outputDirectory.resolve("ui-dsl"))
    outputs.forEach { output ->
        println("Rendered UI preview: $output")
    }
}
