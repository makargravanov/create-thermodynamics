package dev.makargravanov.create_thermodynamics.ui.reactor

import dev.makargravanov.create_thermodynamics.ui.layout.UiDrawCommand
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
        assertEquals(ReactorControllerTab.Zones, ReactorControllerCommandUi.tabAt(84, 34))
        assertEquals(ReactorControllerTab.Mixture, ReactorControllerCommandUi.tabAt(150, 34))
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
