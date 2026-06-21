package dev.makargravanov.create_thermodynamics.ui.reactor

import ru.lazyhat.kraftui.program.FontMetrics
import ru.lazyhat.kraftui.program.ScreenProgramCompiler
import ru.lazyhat.kraftui.program.ScreenRuntimeExecutor
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class ReactorControllerUiTest {
    @Test
    fun `controller ui compiles without layout diagnostics`() {
        val program =
            ScreenProgramCompiler(fontMetrics = FontMetrics { text -> text.length * 6 })
                .compile(
                    root =
                        ReactorControllerUi.build(
                            state = { activeState() },
                            selectedTab = { ReactorControllerTab.Overview },
                            selectedZoneIndex = { 0 },
                            onSelectTab = {},
                            onSelectZone = {},
                        ),
                    rootWidth = ReactorControllerUi.Width,
                    rootHeight = ReactorControllerUi.Height,
                )

        assertTrue(program.diagnostics.isEmpty(), program.diagnostics.joinToString())
    }

    @Test
    fun `controller tabs are real ui dsl click targets`() {
        var selectedTab = ReactorControllerTab.Overview
        val program =
            ScreenProgramCompiler(fontMetrics = FontMetrics { text -> text.length * 6 })
                .compile(
                    root =
                        ReactorControllerUi.build(
                            state = { activeState() },
                            selectedTab = { selectedTab },
                            selectedZoneIndex = { 0 },
                            onSelectTab = { selectedTab = it },
                            onSelectZone = {},
                        ),
                    rootWidth = ReactorControllerUi.Width,
                    rootHeight = ReactorControllerUi.Height,
                )

        assertTrue(ScreenRuntimeExecutor(program).mouseClicked(104, 37))
        assertEquals(ReactorControllerTab.Zones, selectedTab)
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
