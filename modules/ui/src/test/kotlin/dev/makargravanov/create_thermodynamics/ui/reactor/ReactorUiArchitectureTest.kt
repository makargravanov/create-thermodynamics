package dev.makargravanov.create_thermodynamics.ui.reactor

import kotlin.io.path.Path
import kotlin.io.path.readText
import kotlin.test.Test
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class ReactorUiArchitectureTest {
    @Test
    fun reactorControllerDoesNotDefineLocalStyleSystem() {
        val source =
            Path("src/uiSpec/kotlin/dev/makargravanov/create_thermodynamics/ui/reactor/ReactorControllerUi.kt")
                .readText()

        assertFalse("private object Colors" in source)
        assertFalse("style = \"" in source)
        assertTrue("ThermodynamicsUiTheme" in source)
        assertTrue("styledTab(" in source)
        assertTrue("metricCard(" in source)
    }
}
