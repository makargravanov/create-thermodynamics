package dev.makargravanov.create_thermodynamics.ui.preview

import dev.makargravanov.create_thermodynamics.ui.layout.TextOverflowPolicy
import dev.makargravanov.create_thermodynamics.ui.layout.UiLayoutDiagnostic
import dev.makargravanov.create_thermodynamics.ui.layout.UiLayoutEngine
import dev.makargravanov.create_thermodynamics.ui.layout.UiNode
import dev.makargravanov.create_thermodynamics.ui.layout.UiRect
import dev.makargravanov.create_thermodynamics.ui.layout.UiTextMeasurer
import kotlin.io.path.createTempDirectory
import kotlin.io.path.exists
import kotlin.io.path.readText
import kotlin.test.Test
import kotlin.test.assertContains
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CommandPreviewRendererTest {
    private val measurer = UiTextMeasurer { text -> text.length * 6 }

    @Test
    fun `renders layout commands to a png image`() {
        val result =
            UiLayoutEngine(measurer).layout(
                root =
                    UiNode.Panel(
                        id = "root",
                        color = 0xFF112233.toInt(),
                        children =
                            listOf(
                                UiNode.Text(
                                    id = "label",
                                    text = "hello",
                                    color = 0xFFFFFFFF.toInt(),
                                    overflow = TextOverflowPolicy.FailInValidation,
                                ),
                            ),
                    ),
                bounds = UiRect(0, 0, 64, 32),
            )

        val image = CommandPreviewRenderer.render(64, 32, result.commands)

        assertEquals(0xFF112233.toInt(), image.getRGB(4, 4))
        assertTrue(CommandPreviewRenderer.countDistinctColors(image) > 1)
    }

    @Test
    fun `writes machine readable diagnostics report`() {
        val diagnostic =
            UiLayoutDiagnostic.TextWouldOverflow(
                nodeId = "title",
                text = "very long title",
                rect = UiRect(0, 0, 10, 10),
                textWidth = 90,
                policy = TextOverflowPolicy.FailInValidation,
            )
        val output = createTempDirectory("ui-preview-report")

        CommandPreviewRenderer.writeReport(
            outputDirectory = output,
            reports =
                listOf(
                    CommandPreviewReport(
                        id = "overflow",
                        diagnostics = listOf(diagnostic),
                    ),
                ),
        )

        val report = output.resolve("layout-report.json")
        assertTrue(report.exists())
        assertContains(report.readText(), "\"nodeId\":\"title\"")
        assertContains(report.readText(), "\"type\":\"TextWouldOverflow\"")
    }
}
