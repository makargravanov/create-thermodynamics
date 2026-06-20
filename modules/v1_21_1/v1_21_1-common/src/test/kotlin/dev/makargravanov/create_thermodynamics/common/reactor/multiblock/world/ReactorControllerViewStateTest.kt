package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import java.util.UUID
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith

class ReactorControllerViewStateTest {
    @Test
    fun `formed controller state carries zone snapshots instead of one flat mixture`() {
        val zone = ReactorZoneViewState(
            index = 0,
            temperatureKelvin = 298.15,
            pressurePascal = 101_325.0,
            mixture = listOf(
                ReactorMixtureViewEntry("destroy:water", 64.0),
                ReactorMixtureViewEntry("destroy:ethanol", 2.0),
            ),
        )

        val state = ReactorControllerViewState(
            formationState = ReactorControllerFormationState.FORMED,
            structureId = ReactorStructureId(UUID.fromString("9bf1d7b3-08ff-492b-a199-031c7cb038a7")),
            zoneCount = 1,
            chamberBlockCount = 27,
            portCount = 2,
            diagnostic = null,
            nativeBinding = "active",
            zones = listOf(zone),
        )

        assertEquals(1, state.zones.size)
        assertEquals(0, state.zones.single().index)
        assertEquals(64.0, state.zones.single().mixture.single { it.substanceId == "destroy:water" }.concentrationMolPerBucket)
    }

    @Test
    fun `formed controller state rejects zone count that disagrees with snapshots`() {
        assertFailsWith<IllegalArgumentException> {
            ReactorControllerViewState(
                formationState = ReactorControllerFormationState.FORMED,
                structureId = ReactorStructureId(UUID.fromString("d3b8eb9e-a0e6-48ab-aad4-29b53aa9d731")),
                zoneCount = 2,
                chamberBlockCount = 27,
                portCount = 2,
                diagnostic = null,
                nativeBinding = "active",
                zones = listOf(
                    ReactorZoneViewState(
                        index = 0,
                        temperatureKelvin = 300.0,
                        pressurePascal = 100_000.0,
                        mixture = emptyList(),
                    ),
                ),
            )
        }
    }
}
