package dev.makargravanov.create_thermodynamics.ui.preview

import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerCommandUi
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerTab
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerUiSnapshot
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorMixtureUiLine
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorZoneUiSnapshot

object ReactorControllerCommandPreviews {
    fun all(): List<CommandPreviewSpec> =
        listOf(
            preview("reactor_controller_overview_active", activeState(), ReactorControllerTab.Overview),
            preview("reactor_controller_zones_many", manyZonesState(), ReactorControllerTab.Zones),
            preview("reactor_controller_mixture_long_names", longMixtureState(), ReactorControllerTab.Mixture),
        )

    private fun preview(
        id: String,
        state: ReactorControllerUiSnapshot,
        tab: ReactorControllerTab,
    ): CommandPreviewSpec {
        val result =
            ReactorControllerCommandUi.layout(
                state = state,
                selectedTab = tab,
                selectedZoneIndex = 0,
                textMeasurer = CommandPreviewRenderer.textMeasurer,
            )
        return CommandPreviewSpec(
            id = id,
            width = ReactorControllerCommandUi.Width,
            height = ReactorControllerCommandUi.Height,
            commands = result.commands,
            diagnostics = result.diagnostics,
        )
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

    private fun manyZonesState(): ReactorControllerUiSnapshot =
        activeState().copy(
            zoneCount = 6,
            zones =
                (0 until 6).map { index ->
                    ReactorZoneUiSnapshot(
                        index = index,
                        temperature = "${298 + index * 12}.0 K",
                        pressure = "${101 + index * 5}.0 kPa",
                        mixture =
                            listOf(
                                ReactorMixtureUiLine("destroy:water", "${64 - index}.000"),
                            ),
                    )
                },
        )

    private fun longMixtureState(): ReactorControllerUiSnapshot =
        activeState().copy(
            zones =
                listOf(
                    ReactorZoneUiSnapshot(
                        index = 0,
                        temperature = "711.5 K",
                        pressure = "4166.0 kPa",
                        mixture =
                            listOf(
                                ReactorMixtureUiLine("destroy:tetramethylsilane_with_a_very_long_debug_name", "1.250"),
                                ReactorMixtureUiLine("destroy:hexamethylenetetramine_intermediate", "0.750"),
                                ReactorMixtureUiLine("destroy:water", "64.000"),
                                ReactorMixtureUiLine("destroy:oxygen", "8.000"),
                                ReactorMixtureUiLine("destroy:carbon_dioxide", "2.125"),
                            ),
                    ),
                ),
        )
}
