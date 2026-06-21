package dev.makargravanov.create_thermodynamics.ui.reactor

import ru.lazyhat.kraftui.program.FontMetrics
import ru.lazyhat.kraftui.program.ScreenProgramCompiler
import ru.lazyhat.kraftui.program.ScreenRuntimeExecutor
import ru.lazyhat.kraftui.program.UiInputResult
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
                        ),
                    rootWidth = ReactorControllerUi.Width,
                    rootHeight = ReactorControllerUi.Height,
                )

        val result = ScreenRuntimeExecutor(program).mouseClicked(104, 37)
        assertEquals(UiInputResult.Action(ReactorControllerAction.SelectTab(ReactorControllerTab.Zones)), result)
        if (result is UiInputResult.Action) {
            when (val action = result.action) {
                is ReactorControllerAction.SelectTab -> selectedTab = action.tab
                is ReactorControllerAction.SelectZone -> error("unexpected zone action")
            }
        }
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
