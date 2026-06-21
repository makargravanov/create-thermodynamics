package dev.makargravanov.create_thermodynamics.ui.reactor

import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.foundation.UiElement
import ru.lazyhat.kraftui.foundation.UiScope
import ru.lazyhat.kraftui.foundation.modifier.Modifier
import ru.lazyhat.kraftui.foundation.modifier.TextAlignment
import ru.lazyhat.kraftui.foundation.modifier.TextOverflowPolicy
import ru.lazyhat.kraftui.foundation.modifier.background
import ru.lazyhat.kraftui.foundation.modifier.offset
import ru.lazyhat.kraftui.foundation.modifier.size
import ru.lazyhat.kraftui.foundation.modifier.textAlign
import ru.lazyhat.kraftui.foundation.modifier.textOverflow
import ru.lazyhat.kraftui.foundation.uiActions
import ru.lazyhat.kraftui.foundation.value

data class ReactorControllerUiSnapshot(
    val title: String,
    val status: String,
    val active: Boolean,
    val nativeBinding: String,
    val zoneCount: Int,
    val chamberBlocks: Int,
    val portCount: Int,
    val zones: List<ReactorZoneUiSnapshot>,
)

data class ReactorZoneUiSnapshot(
    val index: Int,
    val temperature: String,
    val pressure: String,
    val mixture: List<ReactorMixtureUiLine>,
)

data class ReactorMixtureUiLine(
    val substanceId: String,
    val concentration: String,
) {
    val displayName: String =
        substanceId.substringAfter(':').replace('_', ' ')
}

enum class ReactorControllerTab(
    val label: String,
) {
    Overview("Overview"),
    Zones("Zones"),
    Mixture("Mixture"),
}

sealed interface ReactorControllerAction {
    data class SelectTab(
        val tab: ReactorControllerTab,
    ) : ReactorControllerAction

    data class SelectZone(
        val zoneIndex: Int,
    ) : ReactorControllerAction
}

object ReactorControllerUi {
    const val Width: Int = 232
    const val Height: Int = 140
    const val MaxVisibleZoneRows: Int = 5

    private val mainColumns = ColumnGrid(x = 16, width = 200, gap = 12, count = 3)
    private val tabRects =
        mapOf(
            ReactorControllerTab.Overview to mainColumns.rect(column = 0, y = 29, height = 16),
            ReactorControllerTab.Zones to mainColumns.rect(column = 1, y = 29, height = 16),
            ReactorControllerTab.Mixture to mainColumns.rect(column = 2, y = 29, height = 16),
        )

    fun build(
        state: () -> ReactorControllerUiSnapshot,
        selectedTab: () -> ReactorControllerTab,
        selectedZoneIndex: () -> Int,
    ): UiElement<ReactorControllerAction> =
        uiActions(Modifier.size(Width, Height).background(Colors.Background)) {
            panel(8, 6, 216, 18, Colors.Header)
            label(
                x = 8,
                y = 10,
                width = 216,
                text = { state().title },
                color = Colors.Title,
                alignment = TextAlignment.Center,
            )
            tabs(selectedTab)

            If(value { selectedTab() == ReactorControllerTab.Overview }) {
                overviewPage(state, selectedZoneIndex)
            }
            If(value { selectedTab() == ReactorControllerTab.Zones }) {
                zonesPage(state, selectedZoneIndex)
            }
            If(value { selectedTab() == ReactorControllerTab.Mixture }) {
                mixturePage(state, selectedZoneIndex)
            }
        }

    private fun UiScope<ReactorControllerAction>.tabs(
        selectedTab: () -> ReactorControllerTab,
    ) {
        for (tab in ReactorControllerTab.entries) {
            val rect = tabRects.getValue(tab)
            button(
                modifier = Modifier.offset(rect.x, rect.y).size(rect.width, rect.height),
                action = ReactorControllerAction.SelectTab(tab),
            ) {
                panel(0, 0, rect.width, rect.height, if (selectedTab() == tab) Colors.TabSelected else Colors.Tab)
                label(
                    x = 0,
                    y = 4,
                    width = rect.width,
                    text = { tab.label },
                    color = { if (selectedTab() == tab) Colors.Text else Colors.Muted },
                    alignment = TextAlignment.Center,
                )
            }
        }
    }

    private fun UiScope<ReactorControllerAction>.overviewPage(
        state: () -> ReactorControllerUiSnapshot,
        selectedZoneIndex: () -> Int,
    ) {
        metric(mainColumns.rect(0, 52, 28), "State", { state().status }, { state().active })
        metric(mainColumns.rect(1, 52, 28), "Native", { state().nativeBinding }, { state().nativeBinding == "active" })
        metric(mainColumns.rect(2, 52, 28), "Zones", { "zones ${state().zoneCount}" }, { true })
        card(
            rect = mainColumns.rect(0, 88, 36),
            title = "Structure",
            lines = { listOf("blocks ${state().chamberBlocks}", "ports ${state().portCount}") },
        )
        card(
            rect = mainColumns.rect(1, 88, 36),
            title = "Zone",
            lines = {
                val zone = state().selectedZone(selectedZoneIndex())
                listOf(zone?.temperature ?: "no data", zone?.pressure ?: "")
            },
        )
        card(
            rect = mainColumns.rect(2, 88, 36),
            title = "Mixture",
            lines = {
                val zone = state().selectedZone(selectedZoneIndex())
                listOf(zone?.let { "${it.mixture.size} subs" } ?: "no data", zone?.let { "${it.mixture.size} rows" } ?: "")
            },
        )
    }

    private fun UiScope<ReactorControllerAction>.zonesPage(
        state: () -> ReactorControllerUiSnapshot,
        selectedZoneIndex: () -> Int,
    ) {
        panel(mainColumns.rect(0, 52, 72), Colors.Panel)
        for (row in 0 until MaxVisibleZoneRows) {
            val rect = zoneRowRect(row)
            button(
                modifier = Modifier.offset(rect.x, rect.y).size(rect.width, rect.height),
                action =
                    value {
                        state()
                            .zones
                            .sortedBy { it.index }
                            .getOrNull(row)
                            ?.let { ReactorControllerAction.SelectZone(it.index) }
                    },
            ) {
                label(
                    x = 0,
                    y = 0,
                    width = rect.width,
                    text = { state().zones.sortedBy { it.index }.getOrNull(row)?.let { "zone ${it.index}" } ?: "" },
                    color = {
                        val zone = state().zones.sortedBy { it.index }.getOrNull(row)
                        if (zone?.index == state().selectedZone(selectedZoneIndex())?.index) Colors.Text else Colors.Muted
                    },
                )
            }
        }
        card(
            rect = mainColumns.span(1, 3, 52, 72),
            title = "Zone",
            lines = {
                val zone = state().selectedZone(selectedZoneIndex())
                zone?.let { listOf("zone ${it.index}", it.temperature, it.pressure, "${it.mixture.size} substances") }
                    ?: listOf("no native metrics yet")
            },
        )
    }

    private fun UiScope<ReactorControllerAction>.mixturePage(
        state: () -> ReactorControllerUiSnapshot,
        selectedZoneIndex: () -> Int,
    ) {
        val table = mainColumns.span(0, 3, 52, 72)
        panel(table, Colors.Panel)
        label(
            x = table.x + 4,
            y = table.y + 4,
            width = table.width - 8,
            text = { state().selectedZone(selectedZoneIndex())?.let { "Mixture: zone ${it.index}" } ?: "Mixture" },
            color = Colors.Muted,
        )
        label(table.x + 4, table.y + 18, table.width - 56, { "Substance" }, Colors.Muted)
        label(table.x + table.width - 48, table.y + 18, 44, { "mol/b" }, Colors.Muted, TextAlignment.End)
        for (row in 0 until 4) {
            label(
                x = table.x + 4,
                y = table.y + 30 + row * 10,
                width = table.width - 58,
                text = { state().selectedZone(selectedZoneIndex())?.mixture?.getOrNull(row)?.displayName ?: "" },
                color = Colors.Text,
            )
            label(
                x = table.x + table.width - 48,
                y = table.y + 30 + row * 10,
                width = 44,
                text = { state().selectedZone(selectedZoneIndex())?.mixture?.getOrNull(row)?.concentration ?: "" },
                color = Colors.Text,
                alignment = TextAlignment.End,
            )
        }
    }

    private fun UiScope<ReactorControllerAction>.metric(
        rect: UiRect,
        title: String,
        valueText: () -> String,
        good: () -> Boolean,
    ) {
        panel(rect, Colors.Panel)
        label(rect.x + 4, rect.y + 4, rect.width - 8, { title }, Colors.Muted)
        label(
            x = rect.x + 4,
            y = rect.y + 16,
            width = rect.width - 8,
            text = valueText,
            color = { if (good()) Colors.Good else Colors.Warning },
        )
    }

    private fun UiScope<ReactorControllerAction>.card(
        rect: UiRect,
        title: String,
        lines: () -> List<String>,
    ) {
        panel(rect, Colors.Panel)
        label(rect.x + 4, rect.y + 4, rect.width - 8, { title }, Colors.Muted)
        for (index in 0 until 3) {
            label(
                x = rect.x + 4,
                y = rect.y + 16 + index * 10,
                width = rect.width - 8,
                text = { lines().getOrNull(index) ?: "" },
                color = Colors.Text,
            )
        }
    }

    private fun UiScope<ReactorControllerAction>.panel(
        rect: UiRect,
        color: Color,
    ) {
        panel(rect.x, rect.y, rect.width, rect.height, color)
    }

    private fun UiScope<ReactorControllerAction>.panel(
        x: Int,
        y: Int,
        width: Int,
        height: Int,
        color: Color,
    ) {
        box(Modifier.offset(x, y).size(width, height).background(color))
    }

    private fun UiScope<ReactorControllerAction>.label(
        x: Int,
        y: Int,
        width: Int,
        text: () -> String,
        color: Color,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        label(x, y, width, text, { color }, alignment)
    }

    private fun UiScope<ReactorControllerAction>.label(
        x: Int,
        y: Int,
        width: Int,
        text: () -> String,
        color: () -> Color,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        text(
            modifier =
                Modifier
                    .offset(x, y)
                    .size(width, 9)
                    .textAlign(alignment)
                    .textOverflow(TextOverflowPolicy.Ellipsize),
            color = value(color),
            text = value(text),
        )
    }

    private fun ReactorControllerUiSnapshot.selectedZone(selectedZoneIndex: Int): ReactorZoneUiSnapshot? =
        zones.firstOrNull { it.index == selectedZoneIndex }
            ?: zones.minByOrNull { it.index }

    private fun zoneRowRect(row: Int): UiRect {
        val listRect = mainColumns.rect(column = 0, y = 52, height = 72)
        return UiRect(listRect.x + 2, 56 + row * 12, listRect.width - 4, 10)
    }

    private data class UiRect(
        val x: Int,
        val y: Int,
        val width: Int,
        val height: Int,
    ) {
        val right: Int = x + width
    }

    private class ColumnGrid(
        x: Int,
        width: Int,
        private val gap: Int,
        count: Int,
    ) {
        private val columns: List<UiRect>

        init {
            require(count > 0) { "column grid must have at least one column" }
            require(width >= gap * (count - 1)) { "column grid width is smaller than its gaps" }
            val available = width - gap * (count - 1)
            val baseWidth = available / count
            var spare = available - baseWidth * count
            var nextX = x
            columns =
                (0 until count).map { index ->
                    val columnWidth = baseWidth + if (index >= count - spare) 1 else 0
                    val rect = UiRect(nextX, 0, columnWidth, 0)
                    nextX += columnWidth + gap
                    rect
                }
        }

        fun rect(
            column: Int,
            y: Int,
            height: Int,
        ): UiRect {
            val source = columns[column]
            return UiRect(source.x, y, source.width, height)
        }

        fun span(
            startColumn: Int,
            endColumnExclusive: Int,
            y: Int,
            height: Int,
        ): UiRect {
            require(startColumn < endColumnExclusive) { "column span must include at least one column" }
            val start = columns[startColumn]
            val end = columns[endColumnExclusive - 1]
            return UiRect(start.x, y, end.right - start.x, height)
        }
    }

    private object Colors {
        val Background: Color = Color.rgb(182, 135, 103)
        val Header: Color = Color.rgb(191, 209, 226)
        val Panel: Color = Color.rgb(205, 185, 148)
        val Tab: Color = Color.rgb(192, 157, 119)
        val TabSelected: Color = Color.rgb(214, 195, 158)
        val Title: Color = Color.rgb(61, 60, 72)
        val Text: Color = Color.rgb(61, 60, 72)
        val Muted: Color = Color.rgb(111, 106, 117)
        val Good: Color = Color.rgb(35, 134, 78)
        val Warning: Color = Color.rgb(135, 85, 27)
    }
}
