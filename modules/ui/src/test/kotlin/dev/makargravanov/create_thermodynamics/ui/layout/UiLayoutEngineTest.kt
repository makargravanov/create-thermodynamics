package dev.makargravanov.create_thermodynamics.ui.layout

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertTrue

class UiLayoutEngineTest {
    private val measurer = UiTextMeasurer { text -> text.length * 6 }
    private val engine = UiLayoutEngine(measurer)

    @Test
    fun `fixed text with validation policy reports overflow before runtime`() {
        val result =
            engine.layout(
                root =
                    UiNode.Text(
                        id = "title",
                        text = "too long",
                        color = 0xFFFFFFFF.toInt(),
                        overflow = TextOverflowPolicy.FailInValidation,
                    ),
                bounds = UiRect(0, 0, 20, 10),
            )

        val diagnostic = result.diagnostics.single()
        assertIs<UiLayoutDiagnostic.TextWouldOverflow>(diagnostic)
        assertEquals("title", diagnostic.nodeId)
        assertEquals(48, diagnostic.textWidth)
    }

    @Test
    fun `dynamic text with ellipsis policy is clipped to the assigned rectangle without runtime diagnostic`() {
        val result =
            engine.layout(
                root =
                    UiNode.Text(
                        id = "substance",
                        text = "very long substance name",
                        color = 0xFFFFFFFF.toInt(),
                        overflow = TextOverflowPolicy.EllipsizeWithTooltip,
                    ),
                bounds = UiRect(0, 0, 42, 10),
            )

        assertEquals(emptyList(), result.diagnostics)
        val command = result.commands.single()
        assertIs<UiDrawCommand.DrawText>(command)
        assertEquals("very...", command.text)
        assertEquals(UiRect(0, 0, 42, 10), command.bounds)
        assertEquals("very long substance name", command.tooltip)
    }

    @Test
    fun `row layout keeps children inside their parent`() {
        val result =
            engine.layout(
                root =
                    UiNode.Row(
                        id = "row",
                        gap = 2,
                        children =
                            listOf(
                                UiNode.Panel(id = "left", width = UiLength.Weight(1), children = emptyList()),
                                UiNode.Panel(id = "right", width = UiLength.Weight(2), children = emptyList()),
                            ),
                    ),
                bounds = UiRect(10, 20, 62, 12),
            )

        assertEquals(emptyList(), result.diagnostics)
        val panels = result.commands.filterIsInstance<UiDrawCommand.DrawPanel>()
        assertEquals(UiRect(10, 20, 20, 12), panels.single { it.id == "left" }.bounds)
        assertEquals(UiRect(32, 20, 40, 12), panels.single { it.id == "right" }.bounds)
    }

    @Test
    fun `placed node is offset relative to its parent and still checked against parent bounds`() {
        val result =
            engine.layout(
                root =
                    UiNode.Place(
                        id = "placed",
                        rect = UiRect(4, 5, 12, 6),
                        child = UiNode.Panel(id = "panel", color = 0xFF112233.toInt(), children = emptyList()),
                    ),
                bounds = UiRect(10, 20, 64, 32),
            )

        assertEquals(emptyList(), result.diagnostics)
        val panel = result.commands.single()
        assertIs<UiDrawCommand.DrawPanel>(panel)
        assertEquals(UiRect(14, 25, 12, 6), panel.bounds)
    }

    @Test
    fun `table lays out only visible rows and gives every dynamic cell an overflow policy`() {
        val result =
            engine.layout(
                root =
                    UiNode.Table(
                        id = "mixture",
                        columns =
                            listOf(
                                UiTableColumn(id = "substance", title = "Substance", width = UiLength.Weight(1)),
                                UiTableColumn(id = "amount", title = "mol/b", width = UiLength.Fixed(38)),
                            ),
                        rows =
                            listOf(
                                listOf("water", "64.000"),
                                listOf("ethanol", "2.000"),
                                listOf("tetramethylsilane with intentionally long label", "0.125"),
                            ),
                        visibleRows = 2,
                    ),
                bounds = UiRect(0, 0, 120, 42),
            )

        assertEquals(emptyList(), result.diagnostics)
        val textCommands = result.commands.filterIsInstance<UiDrawCommand.DrawText>()
        assertTrue(textCommands.any { it.text == "Substance" })
        assertTrue(textCommands.any { it.text == "water" })
        assertTrue(textCommands.none { it.text.startsWith("tetramethyl") })
    }
}
