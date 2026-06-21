package dev.makargravanov.create_thermodynamics.ui.reactor

import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.foundation.FixedGridTrack
import ru.lazyhat.kraftui.foundation.UiScope
import ru.lazyhat.kraftui.foundation.WeightedGridTrack
import ru.lazyhat.kraftui.foundation.modifier.Modifier
import ru.lazyhat.kraftui.foundation.modifier.TextAlignment
import ru.lazyhat.kraftui.foundation.modifier.TextOverflowPolicy
import ru.lazyhat.kraftui.foundation.modifier.UiAlignment
import ru.lazyhat.kraftui.foundation.modifier.align
import ru.lazyhat.kraftui.foundation.modifier.background
import ru.lazyhat.kraftui.foundation.modifier.fillMaxHeight
import ru.lazyhat.kraftui.foundation.modifier.fillMaxSize
import ru.lazyhat.kraftui.foundation.modifier.fillMaxWidth
import ru.lazyhat.kraftui.foundation.modifier.height
import ru.lazyhat.kraftui.foundation.modifier.padding
import ru.lazyhat.kraftui.foundation.modifier.size
import ru.lazyhat.kraftui.foundation.modifier.textAlign
import ru.lazyhat.kraftui.foundation.modifier.textOverflow
import ru.lazyhat.kraftui.foundation.stateValue
import ru.lazyhat.kraftui.foundation.uiActions
import ru.lazyhat.kraftui.foundation.value

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

    fun build(state: () -> ReactorControllerGeneratedState) =
        uiActions(
            Modifier
                .size(Width, Height)
                .background(Colors.Background)
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
                    .background(Colors.Header),
        ) {
            label(
                modifier =
                    Modifier
                        .fillMaxWidth()
                        .height(TextHeight)
                        .align(UiAlignment.Center),
                text = stateValue(ReactorControllerGeneratedState::title, state),
                color = Colors.Title,
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
        selected: ru.lazyhat.kraftui.foundation.Value<Boolean>,
        action: ru.lazyhat.kraftui.foundation.Value<ReactorControllerAction?>,
    ) {
        button(
            modifier = Modifier.fillMaxSize(),
            action = action,
        ) {
            box(Modifier.fillMaxSize().background(Colors.Tab))
            If(selected) {
                box(Modifier.fillMaxSize().background(Colors.TabSelected))
                box(
                    Modifier
                        .fillMaxWidth()
                        .height(2)
                        .align(UiAlignment.End)
                        .background(Colors.TabAccent),
                )
                label(
                    modifier =
                        Modifier
                            .fillMaxWidth()
                            .height(TextHeight)
                            .align(UiAlignment.Center),
                    text = label,
                    color = Colors.Text,
                    alignment = TextAlignment.Center,
                )
            }
            label(
                modifier =
                    Modifier
                        .fillMaxWidth()
                        .height(TextHeight)
                        .align(UiAlignment.Center),
                text = label,
                color = Colors.Text,
                alignment = TextAlignment.Center,
            )
        }
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
            modifier = Modifier.fillMaxSize().background(Colors.Panel).padding(3),
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
        text: ru.lazyhat.kraftui.foundation.Value<String>,
        action: ru.lazyhat.kraftui.foundation.Value<ReactorControllerAction?>,
    ) {
        button(
            modifier = Modifier.fillMaxWidth().height(TextHeight),
            action = action,
        ) {
            label(Modifier.fillMaxSize(), text, Colors.Text)
        }
    }

    private fun UiScope<ReactorControllerAction>.mixturePage(state: () -> ReactorControllerGeneratedState) {
        column(
            modifier = Modifier.fillMaxWidth().height(PageHeight).background(Colors.Panel).padding(4),
            gap = 2,
            horizontalAlignment = UiAlignment.Stretch,
        ) {
            label(Modifier.fillMaxWidth().height(TextHeight), stateValue(ReactorControllerGeneratedState::mixtureTitle, state), Colors.Muted)
            grid(
                modifier = Modifier.fillMaxWidth().height(TextHeight),
                columns = listOf(WeightedGridTrack(1f), FixedGridTrack(44)),
                rows = listOf(FixedGridTrack(TextHeight)),
            ) {
                cell(column = 0, row = 0) {
                    label(Modifier.fillMaxSize(), "Substance", Colors.Muted)
                }
                cell(column = 1, row = 0) {
                    label(Modifier.fillMaxSize(), "mol/b", Colors.Muted, TextAlignment.End)
                }
            }
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow0Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow0Amount, state))
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow1Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow1Amount, state))
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow2Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow2Amount, state))
            mixtureRow(stateValue(ReactorControllerGeneratedState::mixtureRow3Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow3Amount, state))
        }
    }

    private fun UiScope<ReactorControllerAction>.mixtureRow(
        name: ru.lazyhat.kraftui.foundation.Value<String>,
        amount: ru.lazyhat.kraftui.foundation.Value<String>,
    ) {
        grid(
            modifier = Modifier.fillMaxWidth().height(TextHeight),
            columns = listOf(WeightedGridTrack(1f), FixedGridTrack(44)),
            rows = listOf(FixedGridTrack(TextHeight)),
        ) {
            cell(column = 0, row = 0) {
                label(Modifier.fillMaxSize(), name, Colors.Text)
            }
            cell(column = 1, row = 0) {
                label(Modifier.fillMaxSize(), amount, Colors.Text, TextAlignment.End)
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.metric(
        title: String,
        valueText: ru.lazyhat.kraftui.foundation.Value<String>,
    ) {
        column(
            modifier = Modifier.fillMaxSize().background(Colors.Panel).padding(4),
            gap = 1,
            horizontalAlignment = UiAlignment.Stretch,
        ) {
            label(Modifier.fillMaxWidth().height(TextHeight), title, Colors.Muted)
            label(Modifier.fillMaxWidth().height(TextHeight), valueText, Colors.Good)
        }
    }

    private fun UiScope<ReactorControllerAction>.card(
        title: String,
        lines: List<ru.lazyhat.kraftui.foundation.Value<String>>,
    ) {
        column(
            modifier = Modifier.fillMaxSize().background(Colors.Panel).padding(3),
            gap = 0,
            horizontalAlignment = UiAlignment.Stretch,
        ) {
            label(Modifier.fillMaxWidth().height(TextHeight), title, Colors.Muted)
            lines.forEach { line ->
                label(Modifier.fillMaxWidth().height(TextHeight), line, Colors.Text)
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.label(
        modifier: Modifier,
        text: String,
        color: Color,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        label(modifier, value(text), color, alignment)
    }

    private fun UiScope<ReactorControllerAction>.label(
        modifier: Modifier,
        text: ru.lazyhat.kraftui.foundation.Value<String>,
        color: Color,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        text(
            modifier =
                modifier
                    .textAlign(alignment)
                    .textOverflow(TextOverflowPolicy.Ellipsize),
            color = color,
            text = text,
        )
    }

    private object Colors {
        val Background: Color = Color.rgb(182, 135, 103)
        val Header: Color = Color.rgb(191, 209, 226)
        val Panel: Color = Color.rgb(205, 185, 148)
        val Tab: Color = Color.rgb(214, 195, 158)
        val TabSelected: Color = Color.rgb(229, 211, 176)
        val TabAccent: Color = Color.rgb(130, 91, 65)
        val Title: Color = Color.rgb(61, 60, 72)
        val Text: Color = Color.rgb(61, 60, 72)
        val Muted: Color = Color.rgb(111, 106, 117)
        val Good: Color = Color.rgb(35, 134, 78)
    }
}
