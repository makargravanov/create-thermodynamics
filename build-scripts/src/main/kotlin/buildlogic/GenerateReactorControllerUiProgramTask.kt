package buildlogic

import org.gradle.api.DefaultTask
import org.gradle.api.file.DirectoryProperty
import org.gradle.api.tasks.OutputDirectory
import org.gradle.api.tasks.TaskAction

abstract class GenerateReactorControllerUiProgramTask : DefaultTask() {
    @get:OutputDirectory
    abstract val outputDirectory: DirectoryProperty

    @TaskAction
    fun generate() {
        val outputFile =
            outputDirectory
                .file(
                    "dev/makargravanov/create_thermodynamics/neoforge/client/generated/" +
                        "GeneratedReactorControllerProgram.kt",
                ).get()
                .asFile
        outputFile.parentFile.mkdirs()
        outputFile.writeText(generatedSource())
    }

    private fun generatedSource(): String =
        """
        package dev.makargravanov.create_thermodynamics.neoforge.client.generated

        import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerUiState
        import ru.lazyhat.kraftui.foundation.Color
        import ru.lazyhat.kraftui.foundation.value
        import ru.lazyhat.kraftui.program.RenderFrame
        import ru.lazyhat.kraftui.program.RenderOp
        import ru.lazyhat.kraftui.program.ScreenProgram
        import ru.lazyhat.kraftui.foundation.modifier.TextAlignment

        /**
         * Generated from the reactor controller UI description.
         *
         * The Minecraft screen must execute this program directly. It must not
         * invoke ScreenProgramCompiler while the game is running.
         */
        object GeneratedReactorControllerProgram {
            const val Width: Int = 216
            const val Height: Int = 142

            fun create(
                rootX: Int,
                rootY: Int,
                state: () -> ReactorControllerUiState,
            ): ScreenProgram =
                ScreenProgram(
                    frames =
                        listOf(
                            RenderFrame(
                                ops =
                                    listOf(
                                        fill(rootX + 0, rootY + 0, Width, Height, 0xFF111317u),
                                        fill(rootX + 0, rootY + 0, Width, 18, 0xFF20252Bu),
                                        text(rootX + 8, rootY + 5, 160, 0xFFECF0F4u) { state().title },

                                        fill(rootX + 8, rootY + 28, 200, 28, 0xFF191D22u),
                                        text(rootX + 15, rootY + 33, 48, 0xFF84919Eu) { "State" },
                                        text(
                                            x = rootX + 64,
                                            y = rootY + 33,
                                            width = 120,
                                            color = { if (state().active) Color.hex(0xFF5ADC96u) else Color.hex(0xFFE8B45Au) },
                                            textValue = { state().status },
                                        ),

                                        fill(rootX + 8, rootY + 64, 96, 66, 0xFF191D22u),
                                        text(rootX + 15, rootY + 69, 72, 0xFF84919Eu) { "Structure" },
                                        text(rootX + 15, rootY + 86, 72, 0xFFD6DCE2u) { "Zones: ${'$'}{state().zoneCount}" },
                                        text(rootX + 15, rootY + 100, 72, 0xFFD6DCE2u) { "Blocks: ${'$'}{state().chamberBlocks}" },
                                        text(rootX + 15, rootY + 114, 72, 0xFFD6DCE2u) { "Ports: ${'$'}{state().portCount}" },

                                        fill(rootX + 112, rootY + 64, 96, 66, 0xFF191D22u),
                                        text(rootX + 119, rootY + 69, 72, 0xFF84919Eu) { "Binding" },
                                        text(rootX + 119, rootY + 86, 72, 0xFFD6DCE2u) {
                                            state().structureId?.take(8) ?: "not formed"
                                        },
                                        text(rootX + 119, rootY + 114, 72, 0xFF84919Eu) { "Native: pending" },
                                    ),
                            ),
                        ),
                    hitRegions = emptyList(),
                )

            private fun fill(
                x: Int,
                y: Int,
                width: Int,
                height: Int,
                color: UInt,
            ): RenderOp.FillRect =
                RenderOp.FillRect(x, y, width, height, Color.hex(color))

            private fun text(
                x: Int,
                y: Int,
                width: Int,
                color: UInt,
                textValue: () -> String,
            ): RenderOp.DrawText =
                text(x, y, width, { Color.hex(color) }, textValue)

            private fun text(
                x: Int,
                y: Int,
                width: Int,
                color: () -> Color,
                textValue: () -> String,
            ): RenderOp.DrawText =
                RenderOp.DrawText(
                    x = x,
                    y = y,
                    width = width,
                    value = value(textValue),
                    color = value(color),
                    alignment = TextAlignment.Start,
                )
        }
        """.trimIndent()
}
