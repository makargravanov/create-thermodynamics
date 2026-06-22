package dev.makargravanov.create_thermodynamics.ui.reactor

import dev.makargravanov.create_thermodynamics.ui.style.ThermodynamicsUiTheme
import ru.lazyhat.kraftui.foundation.FixedGridTrack
import ru.lazyhat.kraftui.foundation.UiScope
import ru.lazyhat.kraftui.foundation.Value
import ru.lazyhat.kraftui.foundation.WeightedGridTrack
import ru.lazyhat.kraftui.foundation.modifier.Modifier
import ru.lazyhat.kraftui.foundation.modifier.TextAlignment
import ru.lazyhat.kraftui.foundation.modifier.UiAlignment
import ru.lazyhat.kraftui.foundation.modifier.align
import ru.lazyhat.kraftui.foundation.modifier.background
import ru.lazyhat.kraftui.foundation.modifier.fillMaxSize
import ru.lazyhat.kraftui.foundation.modifier.fillMaxWidth
import ru.lazyhat.kraftui.foundation.modifier.height
import ru.lazyhat.kraftui.foundation.modifier.padding
import ru.lazyhat.kraftui.foundation.modifier.size
import ru.lazyhat.kraftui.foundation.modifier.textAlign
import ru.lazyhat.kraftui.foundation.modifier.texture
import ru.lazyhat.kraftui.foundation.stateValue
import ru.lazyhat.kraftui.foundation.uiActions
import ru.lazyhat.kraftui.foundation.value
import ru.lazyhat.kraftui.style.StyleColor
import ru.lazyhat.kraftui.style.TextStyle
import ru.lazyhat.kraftui.style.asValue
import ru.lazyhat.kraftui.style.styledPanel
import ru.lazyhat.kraftui.style.styledTab
import ru.lazyhat.kraftui.style.styledText

object ReactorControllerUi {
    const val Width: Int = ReactorControllerUiSize.Width
    const val Height: Int = ReactorControllerUiSize.Height

    private const val OuterPadding = 8
    private const val SectionGap = 8
    private const val ColumnGap = 12
    private const val HeaderHeight = 18
    private const val TabHeight = 16
    private const val MetricHeight = 28
    private const val CardHeight = 36
    private const val PageHeight = 72
    private const val TextHeight = 9

    private val threeColumns = listOf(WeightedGridTrack(1f), WeightedGridTrack(1f), WeightedGridTrack(1f))
    private val theme = ThermodynamicsUiTheme.theme

    fun build(state: () -> ReactorControllerGeneratedState) =
        uiActions(
            Modifier
                .size(Width, Height)
                .background(theme.styles.window.surface.fill.asValue())
                .texture(ThermodynamicsUiTheme.windowTexture)
                .padding(OuterPadding),
        ) {
            column(
                modifier = Modifier.fillMaxSize(),
                gap = SectionGap,
                horizontalAlignment = UiAlignment.Stretch,
            ) {
                header(state)
                tabs(state)

                If(stateValue(ReactorControllerGeneratedState::overviewVisible, state)) {
                    overviewPage(state)
                }
                If(stateValue(ReactorControllerGeneratedState::zonesVisible, state)) {
                    zonesPage(state)
                }
                If(stateValue(ReactorControllerGeneratedState::mixtureVisible, state)) {
                    mixturePage(state)
                }
            }
        }

    private fun UiScope<ReactorControllerAction>.header(state: () -> ReactorControllerGeneratedState) {
        box(
            modifier =
                Modifier
                    .fillMaxWidth()
                    .height(HeaderHeight)
                    .background(StyleColor.Constant(theme.tokens.colors.panelAccent).asValue()),
        ) {
            label(
                modifier =
                    Modifier
                        .fillMaxWidth()
                        .height(TextHeight)
                        .align(UiAlignment.Center),
                text = stateValue(ReactorControllerGeneratedState::title, state),
                style = theme.styles.window.title,
                alignment = TextAlignment.Center,
            )
        }
    }

    private fun UiScope<ReactorControllerAction>.tabs(state: () -> ReactorControllerGeneratedState) {
        grid(
            modifier = Modifier.fillMaxWidth().height(TabHeight),
            columns = threeColumns,
            rows = listOf(FixedGridTrack(TabHeight)),
            columnGap = ColumnGap,
        ) {
            cell(column = 0, row = 0) {
                tab(
                    label = ReactorControllerTab.Overview.label,
                    selected = stateValue(ReactorControllerGeneratedState::overviewSelected, state),
                    action = stateValue(ReactorControllerGeneratedState::selectOverviewAction, state),
                )
            }
            cell(column = 1, row = 0) {
                tab(
                    label = ReactorControllerTab.Zones.label,
                    selected = stateValue(ReactorControllerGeneratedState::zonesSelected, state),
                    action = stateValue(ReactorControllerGeneratedState::selectZonesAction, state),
                )
            }
            cell(column = 2, row = 0) {
                tab(
                    label = ReactorControllerTab.Mixture.label,
                    selected = stateValue(ReactorControllerGeneratedState::mixtureSelected, state),
                    action = stateValue(ReactorControllerGeneratedState::selectMixtureAction, state),
                )
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.tab(
        label: String,
        selected: Value<Boolean>,
        action: Value<ReactorControllerAction?>,
    ) {
        styledTab(
            label = label,
            selected = selected,
            action = action,
            modifier = Modifier.fillMaxSize(),
            style = theme.styles.tab,
        )
    }

    private fun UiScope<ReactorControllerAction>.overviewPage(state: () -> ReactorControllerGeneratedState) {
        column(
            modifier = Modifier.fillMaxWidth().height(PageHeight),
            gap = SectionGap,
            horizontalAlignment = UiAlignment.Stretch,
        ) {
            grid(
                modifier = Modifier.fillMaxWidth().height(MetricHeight),
                columns = threeColumns,
                rows = listOf(FixedGridTrack(MetricHeight)),
                columnGap = ColumnGap,
            ) {
                cell(column = 0, row = 0) {
                    metric("State", stateValue(ReactorControllerGeneratedState::status, state))
                }
                cell(column = 1, row = 0) {
                    metric("Native", stateValue(ReactorControllerGeneratedState::nativeBinding, state))
                }
                cell(column = 2, row = 0) {
                    metric("Zones", stateValue(ReactorControllerGeneratedState::zoneCountText, state))
                }
            }
            grid(
                modifier = Modifier.fillMaxWidth().height(CardHeight),
                columns = threeColumns,
                rows = listOf(FixedGridTrack(CardHeight)),
                columnGap = ColumnGap,
            ) {
                cell(column = 0, row = 0) {
                    card(
                        title = "Structure",
                        lines =
                            listOf(
                                stateValue(ReactorControllerGeneratedState::chamberBlocksText, state),
                                stateValue(ReactorControllerGeneratedState::portCountText, state),
                            ),
                    )
                }
                cell(column = 1, row = 0) {
                    card(
                        title = "Zone",
                        lines =
                            listOf(
                                stateValue(ReactorControllerGeneratedState::selectedZoneTemperature, state),
                                stateValue(ReactorControllerGeneratedState::selectedZonePressure, state),
                            ),
                    )
                }
                cell(column = 2, row = 0) {
                    card(
                        title = "Mixture",
                        lines =
                            listOf(
                                stateValue(ReactorControllerGeneratedState::selectedZoneMixtureCompactCount, state),
                            ),
                    )
                }
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.zonesPage(state: () -> ReactorControllerGeneratedState) {
        grid(
            modifier = Modifier.fillMaxWidth().height(PageHeight),
            columns = threeColumns,
            rows = listOf(FixedGridTrack(PageHeight)),
            columnGap = ColumnGap,
        ) {
            cell(column = 0, row = 0) {
                zoneList(state)
            }
            cell(column = 1, row = 0, columnSpan = 2) {
                card(
                    title = "Zone",
                    lines =
                        listOf(
                            stateValue(ReactorControllerGeneratedState::selectedZoneTitle, state),
                            stateValue(ReactorControllerGeneratedState::selectedZoneTemperature, state),
                            stateValue(ReactorControllerGeneratedState::selectedZonePressure, state),
                            stateValue(ReactorControllerGeneratedState::selectedZoneMixtureCount, state),
                        ),
                )
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.zoneList(state: () -> ReactorControllerGeneratedState) {
        column(
            modifier = Modifier.fillMaxSize().background(theme.styles.panel.surface.fill.asValue()).padding(3),
            gap = 3,
            horizontalAlignment = UiAlignment.Stretch,
        ) {
            zoneRow(stateValue(ReactorControllerGeneratedState::zoneRow0Text, state), stateValue(ReactorControllerGeneratedState::selectZone0Action, state))
            zoneRow(stateValue(ReactorControllerGeneratedState::zoneRow1Text, state), stateValue(ReactorControllerGeneratedState::selectZone1Action, state))
            zoneRow(stateValue(ReactorControllerGeneratedState::zoneRow2Text, state), stateValue(ReactorControllerGeneratedState::selectZone2Action, state))
            zoneRow(stateValue(ReactorControllerGeneratedState::zoneRow3Text, state), stateValue(ReactorControllerGeneratedState::selectZone3Action, state))
            zoneRow(stateValue(ReactorControllerGeneratedState::zoneRow4Text, state), stateValue(ReactorControllerGeneratedState::selectZone4Action, state))
        }
    }

    private fun UiScope<ReactorControllerAction>.zoneRow(
        text: Value<String>,
        action: Value<ReactorControllerAction?>,
    ) {
        button(
            modifier = Modifier.fillMaxWidth().height(TextHeight),
            action = action,
        ) {
            label(Modifier.fillMaxSize(), text, theme.styles.panel.body)
        }
    }

    private fun UiScope<ReactorControllerAction>.mixturePage(state: () -> ReactorControllerGeneratedState) {
        column(
            modifier = Modifier.fillMaxWidth().height(PageHeight).background(theme.styles.panel.surface.fill.asValue()).padding(4),
            gap = 2,
            horizontalAlignment = UiAlignment.Stretch,
        ) {
            label(Modifier.fillMaxWidth().height(TextHeight), stateValue(ReactorControllerGeneratedState::mixtureTitle, state), theme.styles.panel.title)
            grid(
                modifier = Modifier.fillMaxWidth().height(TextHeight),
                columns = listOf(WeightedGridTrack(1f), FixedGridTrack(44)),
                rows = listOf(FixedGridTrack(TextHeight)),
            ) {
                cell(column = 0, row = 0) {
                    label(Modifier.fillMaxSize(), "Substance", theme.styles.panel.title)
                }
                cell(column = 1, row = 0) {
                    label(Modifier.fillMaxSize(), "mol/b", theme.styles.panel.title, TextAlignment.End)
                }
            }
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow0Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow0Amount, state))
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow1Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow1Amount, state))
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow2Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow2Amount, state))
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow3Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow3Amount, state))
        }
    }

    private fun UiScope<ReactorControllerAction>.mixtureRow(
        name: Value<String>,
        amount: Value<String>,
    ) {
        grid(
            modifier = Modifier.fillMaxWidth().height(TextHeight),
            columns = listOf(WeightedGridTrack(1f), FixedGridTrack(44)),
            rows = listOf(FixedGridTrack(TextHeight)),
        ) {
            cell(column = 0, row = 0) {
                label(Modifier.fillMaxSize(), name, theme.styles.panel.body)
            }
            cell(column = 1, row = 0) {
                label(Modifier.fillMaxSize(), amount, theme.styles.panel.body, TextAlignment.End)
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.metric(
        title: String,
        valueText: Value<String>,
    ) {
        styledPanel(
            modifier = Modifier.fillMaxSize(),
            style = theme.styles.panel,
        ) {
            column(
                modifier = Modifier.fillMaxSize(),
                gap = 1,
                horizontalAlignment = UiAlignment.Stretch,
            ) {
                label(Modifier.fillMaxWidth().height(TextHeight), title, theme.styles.panel.title)
                label(Modifier.fillMaxWidth().height(TextHeight), valueText, ThermodynamicsUiTheme.goodText)
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.card(
        title: String,
        lines: List<Value<String>>,
    ) {
        styledPanel(
            modifier = Modifier.fillMaxSize(),
            style = theme.styles.panel,
        ) {
            column(
                modifier = Modifier.fillMaxSize(),
                gap = 1,
                horizontalAlignment = UiAlignment.Stretch,
            ) {
                label(Modifier.fillMaxWidth().height(TextHeight), title, theme.styles.panel.title)
                lines.forEach { line ->
                    label(Modifier.fillMaxWidth().height(TextHeight), line, theme.styles.panel.body)
                }
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.label(
        modifier: Modifier,
        text: String,
        style: TextStyle,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        label(modifier, value(text), style, alignment)
    }

    private fun UiScope<ReactorControllerAction>.label(
        modifier: Modifier,
        text: Value<String>,
        style: TextStyle,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        styledText(
            modifier = modifier.textAlign(alignment),
            text = text,
            style = style,
        )
    }
}
