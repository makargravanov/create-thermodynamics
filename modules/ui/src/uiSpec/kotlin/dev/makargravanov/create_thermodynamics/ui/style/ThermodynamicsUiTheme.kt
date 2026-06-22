package dev.makargravanov.create_thermodynamics.ui.style

import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.foundation.TextureStyle
import ru.lazyhat.kraftui.foundation.modifier.TextOverflowPolicy
import ru.lazyhat.kraftui.style.CreateLikeTheme
import ru.lazyhat.kraftui.style.Insets
import ru.lazyhat.kraftui.style.MetricCardStyle
import ru.lazyhat.kraftui.style.PanelStyle
import ru.lazyhat.kraftui.style.SlotGridStyle
import ru.lazyhat.kraftui.style.StyleColor
import ru.lazyhat.kraftui.style.SurfaceStyle
import ru.lazyhat.kraftui.style.TextStyle
import ru.lazyhat.kraftui.style.TooltipStyle
import ru.lazyhat.kraftui.style.UiTheme

object ThermodynamicsUiTheme {
    private val windowPattern: TextureStyle =
        TextureStyle.Checkerboard(
            first = Color.rgb(184, 137, 101),
            second = Color.rgb(180, 132, 96),
            cellSize = 1,
        )

    private val windowFrame: TextureStyle =
        TextureStyle.BrassFrame(
            base = Color.rgb(187, 134, 58),
            borderWidth = 5,
            noiseStrength = 7,
            ornament = true,
            ornamentSpacing = 14,
            seed = 17,
        )

    val windowTexture: TextureStyle =
        TextureStyle.Layered(
            listOf(
                windowPattern,
                windowFrame,
            ),
        )

    private val controllerPanelSurface =
        SurfaceStyle(
            fill = StyleColor.Constant(CreateLikeTheme.theme.tokens.colors.panel),
            padding = Insets(3, 3, 3, 3),
        )

    val theme: UiTheme =
        CreateLikeTheme.theme.let { base ->
            val panel = base.styles.panel.copy(surface = controllerPanelSurface)
            val window =
                base.styles.window.copy(
                    surface =
                        base.styles.window.surface.copy(
                            fill = StyleColor.Constant(base.tokens.colors.background),
                        ),
                )
            base.copy(
                styles =
                    base.styles.copy(
                        window = window,
                        panel = panel,
                        button = base.styles.button,
                        tab = base.styles.tab,
                        metricCard = MetricCardStyle(panel),
                        listRow = base.styles.button,
                        slotGrid =
                            SlotGridStyle(
                                slot = controllerPanelSurface,
                                hoveredSlot = controllerPanelSurface.copy(fill = StyleColor.Constant(base.tokens.colors.selected)),
                                blockedSlot = base.styles.slotGrid.blockedSlot,
                                gap = base.styles.slotGrid.gap,
                            ),
                        tooltip = TooltipStyle(panel),
                    ),
            )
        }

    val goodText: TextStyle =
        theme.styles.panel.body.copy(
            color = StyleColor.Constant(Color.rgb(35, 134, 78)),
            overflow = TextOverflowPolicy.Ellipsize,
        )
}
