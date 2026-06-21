package dev.makargravanov.create_thermodynamics.ui.reactor

import dev.makargravanov.create_thermodynamics.ui.layout.TextOverflowPolicy
import dev.makargravanov.create_thermodynamics.ui.layout.UiAlignment
import dev.makargravanov.create_thermodynamics.ui.layout.UiDrawCommand
import dev.makargravanov.create_thermodynamics.ui.layout.UiLayoutEngine
import dev.makargravanov.create_thermodynamics.ui.layout.UiLayoutResult
import dev.makargravanov.create_thermodynamics.ui.layout.UiLength
import dev.makargravanov.create_thermodynamics.ui.layout.UiNode
import dev.makargravanov.create_thermodynamics.ui.layout.UiRect
import dev.makargravanov.create_thermodynamics.ui.layout.UiTableColumn
import dev.makargravanov.create_thermodynamics.ui.layout.UiTextMeasurer

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

object ReactorControllerCommandUi {
    const val Width: Int = 232
    const val Height: Int = 140
    const val MaxVisibleZoneRows: Int = 5

    private val Bounds = UiRect(0, 0, Width, Height)
    private val MainColumns = ColumnGrid(x = 16, width = 200, gap = 12, count = 3)
    private val TabRects =
        mapOf(
            ReactorControllerTab.Overview to MainColumns.rect(column = 0, y = 29, height = 16),
            ReactorControllerTab.Zones to MainColumns.rect(column = 1, y = 29, height = 16),
            ReactorControllerTab.Mixture to MainColumns.rect(column = 2, y = 29, height = 16),
        )

    fun layout(
        state: ReactorControllerUiSnapshot,
        selectedTab: ReactorControllerTab,
        selectedZoneIndex: Int,
        textMeasurer: UiTextMeasurer,
    ): UiLayoutResult =
        UiLayoutEngine(textMeasurer).layout(
            root = buildNode(state, selectedTab, selectedZoneIndex),
            bounds = Bounds,
        )

    fun tabAt(
        x: Int,
        y: Int,
    ): ReactorControllerTab? =
        TabRects.entries.firstOrNull { (_, rect) -> x >= rect.x && x < rect.right && y >= rect.y && y < rect.bottom }?.key

    fun zoneAt(
        x: Int,
        y: Int,
        state: ReactorControllerUiSnapshot,
    ): Int? {
        val zones = state.zones.sortedBy { it.index }.take(MaxVisibleZoneRows)
        for ((row, zone) in zones.withIndex()) {
            val rect = zoneRowRect(row)
            if (x >= rect.x && x < rect.right && y >= rect.y && y < rect.bottom) {
                return zone.index
            }
        }
        return null
    }

    private fun buildNode(
        state: ReactorControllerUiSnapshot,
        selectedTab: ReactorControllerTab,
        selectedZoneIndex: Int,
    ): UiNode {
        val selectedZone = selectedZone(state, selectedZoneIndex)
        return UiNode.Panel(
            id = "reactor.root",
            color = Colors.Background,
            children =
                listOf(
                    place("header.panel", UiRect(8, 6, 216, 18), panel("header.panel.inner", Colors.Header)),
                    place(
                        "header.title",
                        UiRect(8, 9, 216, 10),
                        text("header.title.text", state.title, Colors.Title, TextOverflowPolicy.FailInValidation, UiAlignment.Center),
                    ),
                    tabs(selectedTab),
                    when (selectedTab) {
                        ReactorControllerTab.Overview -> overviewPage(state, selectedZone)
                        ReactorControllerTab.Zones -> zonesPage(state, selectedZone)
                        ReactorControllerTab.Mixture -> mixturePage(selectedZone)
                    },
                ),
        )
    }

    private fun tabs(selectedTab: ReactorControllerTab): UiNode.Panel =
        UiNode.Panel(
            id = "tabs",
            children =
                ReactorControllerTab.entries.map { tab ->
                    val rect = requireNotNull(TabRects[tab])
                    place(
                        id = "tab.${tab.name}",
                        rect = rect,
                        child =
                            UiNode.Panel(
                                id = "tab.${tab.name}.panel",
                                color = if (tab == selectedTab) Colors.TabSelected else Colors.Tab,
                                children =
                                    listOf(
                                        place(
                                            "tab.${tab.name}.text.place",
                                            UiRect(0, 4, rect.width, 10),
                                            text(
                                                id = "tab.${tab.name}.text",
                                                value = tab.label,
                                                color = if (tab == selectedTab) Colors.Text else Colors.Muted,
                                                overflow = TextOverflowPolicy.FailInValidation,
                                                alignment = UiAlignment.Center,
                                            ),
                                        ),
                                    ),
                            ),
                    )
                },
        )

    private fun overviewPage(
        state: ReactorControllerUiSnapshot,
        selectedZone: ReactorZoneUiSnapshot?,
    ): UiNode.Panel =
        UiNode.Panel(
            id = "overview",
            children =
                listOf(
                    metric("overview.state", MainColumns.rect(column = 0, y = 52, height = 28), "State", state.status, state.active),
                    metric("overview.native", MainColumns.rect(column = 1, y = 52, height = 28), "Native", state.nativeBinding, state.nativeBinding == "active"),
                    metric("overview.zones", MainColumns.rect(column = 2, y = 52, height = 28), "Zones", "zones ${state.zoneCount}", true),
                    card(
                        id = "overview.structure",
                        rect = MainColumns.rect(column = 0, y = 88, height = 36),
                        title = "Structure",
                        lines = listOf("blocks ${state.chamberBlocks}", "ports ${state.portCount}"),
                    ),
                    card(
                        id = "overview.zone",
                        rect = MainColumns.rect(column = 1, y = 88, height = 36),
                        title = "Zone",
                        lines = listOf(selectedZone?.temperature ?: "no data", selectedZone?.pressure ?: ""),
                    ),
                    card(
                        id = "overview.mixture",
                        rect = MainColumns.rect(column = 2, y = 88, height = 36),
                        title = "Mixture",
                        lines =
                            listOf(
                                selectedZone?.let { "${it.mixture.size} subs" } ?: "no data",
                                selectedZone?.let { "${it.mixture.size} rows" } ?: "",
                            ),
                    ),
                ),
        )

    private fun zonesPage(
        state: ReactorControllerUiSnapshot,
        selectedZone: ReactorZoneUiSnapshot?,
    ): UiNode.Panel =
        UiNode.Panel(
            id = "zones",
            children =
                listOf(
                    place("zones.list.panel", MainColumns.rect(column = 0, y = 52, height = 72), panel("zones.list.panel.inner", Colors.Panel)),
                    *state.zones.sortedBy { it.index }.take(MaxVisibleZoneRows).mapIndexed { row, zone ->
                        place(
                            id = "zones.row.${zone.index}",
                            rect = zoneRowRect(row),
                            child =
                                text(
                                    id = "zones.row.${zone.index}.text",
                                    value = "zone ${zone.index}",
                                    color = if (zone.index == selectedZone?.index) Colors.Text else Colors.Muted,
                                    overflow = TextOverflowPolicy.FailInValidation,
                                ),
                        )
                    }.toTypedArray(),
                    card(
                        id = "zones.detail",
                        rect = MainColumns.span(startColumn = 1, endColumnExclusive = 3, y = 52, height = 72),
                        title = selectedZone?.let { "Zone ${it.index}" } ?: "Zone",
                        lines =
                            selectedZone?.let {
                                listOf(it.temperature, it.pressure, "${it.mixture.size} substances")
                            } ?: listOf("no native metrics yet"),
                    ),
                ),
        )

    private fun mixturePage(selectedZone: ReactorZoneUiSnapshot?): UiNode.Panel =
        UiNode.Panel(
            id = "mixture",
            children =
                listOf(
                    place(
                        "mixture.title",
                        MainColumns.span(startColumn = 0, endColumnExclusive = 3, y = 52, height = 10),
                        text(
                            id = "mixture.title.text",
                            value = selectedZone?.let { "Mixture: zone ${it.index}" } ?: "Mixture",
                            color = Colors.Muted,
                            overflow = TextOverflowPolicy.FailInValidation,
                        ),
                    ),
                    place(
                        "mixture.table.place",
                        MainColumns.span(startColumn = 0, endColumnExclusive = 3, y = 66, height = 58),
                        UiNode.Table(
                            id = "mixture.table",
                            columns =
                                listOf(
                                    UiTableColumn("substance", "Substance", UiLength.Weight(1)),
                                    UiTableColumn("amount", "mol/b", UiLength.Fixed(44)),
                                ),
                            rows = selectedZone?.mixture?.map { listOf(it.displayName, it.concentration) } ?: emptyList(),
                            visibleRows = 4,
                            textColor = Colors.Text,
                            mutedTextColor = Colors.Muted,
                        ),
                    ),
                ),
        )

    private fun metric(
        id: String,
        rect: UiRect,
        title: String,
        value: String,
        good: Boolean,
    ): UiNode =
        UiNode.Place(
            id = "$id.place",
            rect = rect,
            child =
                UiNode.Panel(
                    id = "$id.panel",
                    color = Colors.Panel,
                    children =
                        listOf(
                            place("$id.title.place", UiRect(4, 4, rect.width - 8, 10), text("$id.title", title, Colors.Muted, TextOverflowPolicy.FailInValidation)),
                            place("$id.value.place", UiRect(4, 16, rect.width - 8, 10), text("$id.value", value, if (good) Colors.Good else Colors.Warning, TextOverflowPolicy.EllipsizeWithTooltip)),
                        ),
                ),
        )

    private fun card(
        id: String,
        rect: UiRect,
        title: String,
        lines: List<String>,
    ): UiNode =
        UiNode.Place(
            id = "$id.place",
            rect = rect,
            child =
                UiNode.Panel(
                    id = "$id.panel",
                    color = Colors.Panel,
                    children =
                        listOf(
                            place("$id.title.place", UiRect(4, 4, rect.width - 8, 10), text("$id.title", title, Colors.Muted, TextOverflowPolicy.EllipsizeWithTooltip)),
                            *lines.take(3).mapIndexed { index, line ->
                                place(
                                    "$id.line.$index.place",
                                    UiRect(4, 16 + index * 10, rect.width - 8, 10),
                                    text("$id.line.$index", line, Colors.Text, TextOverflowPolicy.EllipsizeWithTooltip),
                                )
                            }.toTypedArray(),
                        ),
                ),
        )

    private fun selectedZone(
        state: ReactorControllerUiSnapshot,
        selectedZoneIndex: Int,
    ): ReactorZoneUiSnapshot? =
        state.zones.firstOrNull { it.index == selectedZoneIndex }
            ?: state.zones.minByOrNull { it.index }

    private fun panel(
        id: String,
        color: Int,
    ): UiNode.Panel =
        UiNode.Panel(id = id, color = color, children = emptyList())

    private fun text(
        id: String,
        value: String,
        color: Int,
        overflow: TextOverflowPolicy,
        alignment: UiAlignment = UiAlignment.Start,
    ): UiNode.Text =
        UiNode.Text(id = id, text = value, color = color, overflow = overflow, alignment = alignment)

    private fun place(
        id: String,
        rect: UiRect,
        child: UiNode,
    ): UiNode.Place =
        UiNode.Place(id = id, rect = rect, child = child)

    private fun zoneRowRect(row: Int): UiRect {
        val listRect = MainColumns.rect(column = 0, y = 52, height = 72)
        return UiRect(listRect.x + 2, 56 + row * 12, listRect.width - 4, 10)
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
        const val Background: Int = 0xFFB68767.toInt()
        const val Header: Int = 0xFFBFD1E2.toInt()
        const val Panel: Int = 0xFFCDB994.toInt()
        const val Tab: Int = 0xFFC09D77.toInt()
        const val TabSelected: Int = 0xFFD6C39E.toInt()
        const val Title: Int = 0xFF3D3C48.toInt()
        const val Text: Int = 0xFF3D3C48.toInt()
        const val Muted: Int = 0xFF6F6A75.toInt()
        const val Good: Int = 0xFF23864E.toInt()
        const val Warning: Int = 0xFF87551B.toInt()
    }
}
