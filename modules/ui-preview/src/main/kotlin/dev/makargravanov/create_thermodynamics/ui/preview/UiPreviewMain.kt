package dev.makargravanov.create_thermodynamics.ui.preview

import java.nio.file.Path

fun main(args: Array<String>) {
    val outputDirectory =
        if (args.isNotEmpty()) {
            Path.of(args[0])
        } else {
            Path.of("build", "reports", "ui")
        }
    val legacyOutputs = UiPreviewRenderer.renderAll(ReactorPortPreviews.all(), outputDirectory.resolve("legacy"))
    val commandOutputs = CommandPreviewRenderer.renderAll(ReactorControllerCommandPreviews.all(), outputDirectory.resolve("commands"))
    (legacyOutputs + commandOutputs).forEach { output ->
        println("Rendered UI preview: $output")
    }
}
