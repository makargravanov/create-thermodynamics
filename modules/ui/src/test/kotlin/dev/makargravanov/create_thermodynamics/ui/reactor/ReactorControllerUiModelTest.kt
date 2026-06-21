package dev.makargravanov.create_thermodynamics.ui.reactor

import kotlin.test.Test
import kotlin.test.assertEquals

class ReactorControllerUiModelTest {
    @Test
    fun `snapshot creates flat generated state for selected zone`() {
        val state =
            activeState().toGeneratedState(
                selectedTab = ReactorControllerTab.Mixture,
                selectedZoneIndex = 0,
            )

        assertEquals("Reactor Controller", state.title)
        assertEquals("formed", state.status)
        assertEquals(false, state.overviewVisible)
        assertEquals(false, state.zonesVisible)
        assertEquals(true, state.mixtureVisible)
        assertEquals("blocks 27", state.chamberBlocksText)
        assertEquals("ports 2", state.portCountText)
        assertEquals("Mixture: zone 0", state.mixtureTitle)
        assertEquals("water", state.mixtureRow0Name)
        assertEquals("64.000", state.mixtureRow0Amount)
        assertEquals(ReactorControllerAction.SelectTab(ReactorControllerTab.Zones), state.selectZonesAction)
        assertEquals(ReactorControllerAction.SelectZone(0), state.selectZone0Action)
    }

    private fun activeState(): ReactorControllerUiSnapshot =
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
                        mixture =
                            listOf(
                                ReactorMixtureUiLine("destroy:water", "64.000"),
                                ReactorMixtureUiLine("destroy:ethanol", "2.000"),
                            ),
                    ),
                ),
        )
}
