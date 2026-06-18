package dev.makargravanov.create_thermodynamics.ui.preview

import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerUi
import dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerUiState
import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.foundation.modifier.Modifier
import ru.lazyhat.kraftui.foundation.modifier.background
import ru.lazyhat.kraftui.foundation.modifier.offset
import ru.lazyhat.kraftui.foundation.modifier.size
import ru.lazyhat.kraftui.foundation.ui

object ReactorPortPreviews {
    fun all(): List<UiPreviewSpec> =
        listOf(
            UiPreviewSpec(
                id = "reactor_controller",
                width = ReactorControllerUi.Width,
                height = ReactorControllerUi.Height,
                root =
                    ReactorControllerUi.build {
                        ReactorControllerUiState(
                            title = "Reactor Controller",
                            status = "formed",
                            structureId = "9f1a4c71-8b55-48df-9fb5-5f8ad7e1a334",
                            active = true,
                            zoneCount = 1,
                            chamberBlocks = 12,
                            portCount = 4,
                        )
                    },
            ),
            UiPreviewSpec(
                id = "reactor_port_inventory",
                width = 176,
                height = 90,
                root =
                    ui(Modifier.size(176, 90).background(Color.rgb(22, 24, 28))) {
                        text("Reactor Port", modifier = Modifier.offset(8, 6), color = Color.rgb(232, 236, 240))
                        box(Modifier.offset(7, 21).size(162, 1).background(Color.rgb(55, 62, 70)))
                        row(modifier = Modifier.offset(8, 30), gap = 2) {
                            repeat(9) {
                                box(Modifier.size(16, 16).background(Color.rgb(45, 51, 58))) {
                                    box(Modifier.offset(1, 1).size(14, 14).background(Color.rgb(30, 34, 39)))
                                }
                            }
                        }
                        text("Debug item buffer", modifier = Modifier.offset(8, 56), color = Color.rgb(148, 159, 170))
                    },
            ),
        )
}
