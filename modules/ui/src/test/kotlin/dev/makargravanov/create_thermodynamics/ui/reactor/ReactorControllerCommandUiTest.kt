package dev.makargravanov.create_thermodynamics.ui.reactor

import dev.makargravanov.create_thermodynamics.ui.layout.UiDrawCommand
import dev.makargravanov.create_thermodynamics.ui.layout.UiLayoutResult
import dev.makargravanov.create_thermodynamics.ui.layout.UiTextMeasurer
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotNull
import kotlin.test.assertTrue

class ReactorControllerCommandUiTest {
    private val measurer = UiTextMeasurer { text -> text.length * 6 }

    @Test
    fun `overview page has no layout diagnostics`() {
        val result =
            ReactorControllerCommandUi.layout(
                state = sampleState(),
                selectedTab = ReactorControllerTab.Overview,
                selectedZoneIndex = 0,
                textMeasurer = measurer,
            )

        assertEquals(emptyList(), result.diagnostics)
        val text = result.commands.filterIsInstance<UiDrawCommand.DrawText>().map { it.text }
        assertTrue("formed" in text)
        assertTrue("active" in text)
        assertTrue("zones 1" in text)
    }

    @Test
    fun `mixture page ellipsizes long dynamic substance names with tooltip`() {
        val result =
            ReactorControllerCommandUi.layout(
                state =
                    sampleState(
                        mixture =
                            listOf(
                                ReactorMixtureUiLine(
                                    substanceId = "destroy:tetramethylsilane_with_a_very_long_debug_name",
                                    concentration = "1.250",
                                ),
                            ),
                    ),
                selectedTab = ReactorControllerTab.Mixture,
                selectedZoneIndex = 0,
                textMeasurer = measurer,
            )

        assertEquals(emptyList(), result.diagnostics)
        val clipped =
            result.commands
                .filterIsInstance<UiDrawCommand.DrawText>()
                .singleOrNull { it.tooltip == "tetramethylsilane with a very long debug name" }
        assertNotNull(clipped)
        assertTrue(clipped.text.endsWith("..."))
    }

    @Test
    fun `tab hit testing is derived from the shared controller layout`() {
        assertEquals(ReactorControllerTab.Overview, ReactorControllerCommandUi.tabAt(18, 34))
        assertEquals(ReactorControllerTab.Zones, ReactorControllerCommandUi.tabAt(88, 34))
        assertEquals(ReactorControllerTab.Mixture, ReactorControllerCommandUi.tabAt(158, 34))
    }

    @Test
    fun `overview tabs metrics and cards share the same column grid`() {
        val result =
            ReactorControllerCommandUi.layout(
                state = sampleState(),
                selectedTab = ReactorControllerTab.Overview,
                selectedZoneIndex = 0,
                textMeasurer = measurer,
            )

        assertEquals(emptyList(), result.diagnostics)
        assertSameColumn(result, "tab.Overview.panel", "overview.state.panel", "overview.structure.panel")
        assertSameColumn(result, "tab.Zones.panel", "overview.native.panel", "overview.zone.panel")
        assertSameColumn(result, "tab.Mixture.panel", "overview.zones.panel", "overview.mixture.panel")
    }

    private fun assertSameColumn(
        result: UiLayoutResult,
        vararg panelIds: String,
    ) {
        val rects =
            panelIds.map { id ->
                result.commands
                    .filterIsInstance<UiDrawCommand.DrawPanel>()
                    .single { it.id == id }
                    .bounds
            }
        val first = rects.first()
        rects.drop(1).forEach { rect ->
            assertEquals(first.x, rect.x, "column x mismatch for ${panelIds.joinToString()}")
            assertEquals(first.right, rect.right, "column right edge mismatch for ${panelIds.joinToString()}")
        }
    }

    private fun sampleState(
        mixture: List<ReactorMixtureUiLine> =
            listOf(
                ReactorMixtureUiLine("destroy:water", "64.000"),
                ReactorMixtureUiLine("destroy:ethanol", "2.000"),
            ),
    ): ReactorControllerUiSnapshot =
        ReactorControllerUiSnapshot(
            title = "Reactor Controller",
            status = "formed",
            active = true,
            nativeBinding = "active",
            zoneCount = 1,
            chamberBlocks = 27,
            portCount = 2,
            zones =
                listOf(
                    ReactorZoneUiSnapshot(
                        index = 0,
                        temperature = "298.0 K",
                        pressure = "101.3 kPa",
                        mixture = mixture,
                    ),
                ),
        )
}
