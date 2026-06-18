package dev.makargravanov.create_thermodynamics.ui.reactor

import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.foundation.UiElement
import ru.lazyhat.kraftui.foundation.modifier.Modifier
import ru.lazyhat.kraftui.foundation.modifier.background
import ru.lazyhat.kraftui.foundation.modifier.offset
import ru.lazyhat.kraftui.foundation.modifier.size
import ru.lazyhat.kraftui.foundation.ui
import ru.lazyhat.kraftui.foundation.value

data class ReactorControllerUiState(
    val title: String,
    val status: String,
    val structureId: String?,
    val active: Boolean,
    val zoneCount: Int,
    val chamberBlocks: Int,
    val portCount: Int,
)

object ReactorControllerUi {
    const val Width: Int = 216
    const val Height: Int = 142

    fun build(state: () -> ReactorControllerUiState): UiElement =
        ui(Modifier.size(Width, Height).background(Colors.Background)) {
            box(Modifier.offset(0, 0).size(Width, 18).background(Colors.Header))
            text(
                text = value { state().title },
                modifier = Modifier.offset(8, 5),
                color = Colors.Title,
            )

            statusPanel(state)
            structurePanel(state)
            portPanel(state)
        }

    private fun ru.lazyhat.kraftui.foundation.UiScope.statusPanel(state: () -> ReactorControllerUiState) {
        box(Modifier.offset(8, 28).size(200, 28).background(Colors.Panel)) {
            text("State", modifier = Modifier.offset(7, 5), color = Colors.Muted)
            text(
                text = value { state().status },
                modifier = Modifier.offset(56, 5),
                color = value { if (state().active) Colors.Good else Colors.Warning },
            )
        }
    }

    private fun ru.lazyhat.kraftui.foundation.UiScope.structurePanel(state: () -> ReactorControllerUiState) {
        box(Modifier.offset(8, 64).size(96, 66).background(Colors.Panel)) {
            text("Structure", modifier = Modifier.offset(7, 5), color = Colors.Muted)
            text(
                text = value { "Zones: ${state().zoneCount}" },
                modifier = Modifier.offset(7, 22),
                color = Colors.Text,
            )
            text(
                text = value { "Blocks: ${state().chamberBlocks}" },
                modifier = Modifier.offset(7, 36),
                color = Colors.Text,
            )
            text(
                text = value { "Ports: ${state().portCount}" },
                modifier = Modifier.offset(7, 50),
                color = Colors.Text,
            )
        }
    }

    private fun ru.lazyhat.kraftui.foundation.UiScope.portPanel(state: () -> ReactorControllerUiState) {
        box(Modifier.offset(112, 64).size(96, 66).background(Colors.Panel)) {
            text("Binding", modifier = Modifier.offset(7, 5), color = Colors.Muted)
            text(
                text = value { state().structureId?.take(8) ?: "not formed" },
                modifier = Modifier.offset(7, 22),
                color = Colors.Text,
            )
            text("Native: pending", modifier = Modifier.offset(7, 50), color = Colors.Muted)
        }
    }

    private object Colors {
        val Background = Color.rgb(17, 19, 23)
        val Header = Color.rgb(32, 37, 43)
        val Panel = Color.rgb(25, 29, 34)
        val Title = Color.rgb(236, 240, 244)
        val Text = Color.rgb(214, 220, 226)
        val Muted = Color.rgb(132, 145, 158)
        val Good = Color.rgb(90, 220, 150)
        val Warning = Color.rgb(232, 180, 90)
    }
}
