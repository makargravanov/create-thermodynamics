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
import ru.lazyhat.kraftui.foundation.stateValue
import ru.lazyhat.kraftui.foundation.uiActions
import ru.lazyhat.kraftui.foundation.value

object ReactorControllerUi {
    const val Width: Int = ReactorControllerUiSize.Width
    const val Height: Int = ReactorControllerUiSize.Height

    private val mainColumns = ColumnGrid(x = 16, width = 200, gap = 12, count = 3)
    private val tabRects =
        mapOf(
            ReactorControllerTab.Overview to mainColumns.rect(column = 0, y = 29, height = 16),
            ReactorControllerTab.Zones to mainColumns.rect(column = 1, y = 29, height = 16),
            ReactorControllerTab.Mixture to mainColumns.rect(column = 2, y = 29, height = 16),
        )

    fun build(state: () -> ReactorControllerGeneratedState): UiElement<ReactorControllerAction> =
        uiActions(Modifier.size(Width, Height).background(Colors.Background)) {
            panel(8, 6, 216, 18, Colors.Header)
            label(
                x = 8,
                y = 10,
                width = 216,
                text = stateValue(ReactorControllerGeneratedState::title, state),
                color = Colors.Title,
                alignment = TextAlignment.Center,
            )
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

    private fun UiScope<ReactorControllerAction>.tabs(state: () -> ReactorControllerGeneratedState) {
        tab(
            rect = tabRects.getValue(ReactorControllerTab.Overview),
            label = ReactorControllerTab.Overview.label,
            action = stateValue(ReactorControllerGeneratedState::selectOverviewAction, state),
        )
        tab(
            rect = tabRects.getValue(ReactorControllerTab.Zones),
            label = ReactorControllerTab.Zones.label,
            action = stateValue(ReactorControllerGeneratedState::selectZonesAction, state),
        )
        tab(
            rect = tabRects.getValue(ReactorControllerTab.Mixture),
            label = ReactorControllerTab.Mixture.label,
            action = stateValue(ReactorControllerGeneratedState::selectMixtureAction, state),
        )
    }

    private fun UiScope<ReactorControllerAction>.tab(
        rect: UiRect,
        label: String,
        action: ru.lazyhat.kraftui.foundation.Value<ReactorControllerAction?>,
    ) {
        button(
            modifier = Modifier.offset(rect.x, rect.y).size(rect.width, rect.height),
            action = action,
        ) {
            panel(0, 0, rect.width, rect.height, Colors.Tab)
            label(
                x = 0,
                y = 4,
                width = rect.width,
                text = label,
                color = Colors.Text,
                alignment = TextAlignment.Center,
            )
        }
    }

    private fun UiScope<ReactorControllerAction>.overviewPage(state: () -> ReactorControllerGeneratedState) {
        metric(mainColumns.rect(0, 52, 28), "State", stateValue(ReactorControllerGeneratedState::status, state))
        metric(mainColumns.rect(1, 52, 28), "Native", stateValue(ReactorControllerGeneratedState::nativeBinding, state))
        metric(mainColumns.rect(2, 52, 28), "Zones", stateValue(ReactorControllerGeneratedState::zoneCountText, state))
        card(
            rect = mainColumns.rect(0, 88, 36),
            title = "Structure",
            line0 = stateValue(ReactorControllerGeneratedState::chamberBlocksText, state),
            line1 = stateValue(ReactorControllerGeneratedState::portCountText, state),
        )
        card(
            rect = mainColumns.rect(1, 88, 36),
            title = "Zone",
            line0 = stateValue(ReactorControllerGeneratedState::selectedZoneTemperature, state),
            line1 = stateValue(ReactorControllerGeneratedState::selectedZonePressure, state),
        )
        card(
            rect = mainColumns.rect(2, 88, 36),
            title = "Mixture",
            line0 = stateValue(ReactorControllerGeneratedState::selectedZoneMixtureCount, state),
            line1 = value(""),
        )
    }

    private fun UiScope<ReactorControllerAction>.zonesPage(state: () -> ReactorControllerGeneratedState) {
        panel(mainColumns.rect(0, 52, 72), Colors.Panel)
        zoneRow(0, stateValue(ReactorControllerGeneratedState::zoneRow0Text, state), stateValue(ReactorControllerGeneratedState::selectZone0Action, state))
        zoneRow(1, stateValue(ReactorControllerGeneratedState::zoneRow1Text, state), stateValue(ReactorControllerGeneratedState::selectZone1Action, state))
        zoneRow(2, stateValue(ReactorControllerGeneratedState::zoneRow2Text, state), stateValue(ReactorControllerGeneratedState::selectZone2Action, state))
        zoneRow(3, stateValue(ReactorControllerGeneratedState::zoneRow3Text, state), stateValue(ReactorControllerGeneratedState::selectZone3Action, state))
        zoneRow(4, stateValue(ReactorControllerGeneratedState::zoneRow4Text, state), stateValue(ReactorControllerGeneratedState::selectZone4Action, state))
        card(
            rect = mainColumns.span(1, 3, 52, 72),
            title = "Zone",
            line0 = stateValue(ReactorControllerGeneratedState::selectedZoneTitle, state),
            line1 = stateValue(ReactorControllerGeneratedState::selectedZoneTemperature, state),
            line2 = stateValue(ReactorControllerGeneratedState::selectedZonePressure, state),
            line3 = stateValue(ReactorControllerGeneratedState::selectedZoneMixtureCount, state),
        )
    }

    private fun UiScope<ReactorControllerAction>.mixturePage(state: () -> ReactorControllerGeneratedState) {
        val table = mainColumns.span(0, 3, 52, 72)
        panel(table, Colors.Panel)
        label(table.x + 4, table.y + 4, table.width - 8, stateValue(ReactorControllerGeneratedState::mixtureTitle, state), Colors.Muted)
        label(table.x + 4, table.y + 18, table.width - 56, "Substance", Colors.Muted)
        label(table.x + table.width - 48, table.y + 18, 44, "mol/b", Colors.Muted, TextAlignment.End)
        mixtureRow(table, 0, stateValue(ReactorControllerGeneratedState::mixtureRow0Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow0Amount, state))
        mixtureRow(table, 1, stateValue(ReactorControllerGeneratedState::mixtureRow1Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow1Amount, state))
        mixtureRow(table, 2, stateValue(ReactorControllerGeneratedState::mixtureRow2Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow2Amount, state))
        mixtureRow(table, 3, stateValue(ReactorControllerGeneratedState::mixtureRow3Name, state), stateValue(ReactorControllerGeneratedState::mixtureRow3Amount, state))
    }

    private fun UiScope<ReactorControllerAction>.zoneRow(
        row: Int,
        text: ru.lazyhat.kraftui.foundation.Value<String>,
        action: ru.lazyhat.kraftui.foundation.Value<ReactorControllerAction?>,
    ) {
        val rect = zoneRowRect(row)
        button(
            modifier = Modifier.offset(rect.x, rect.y).size(rect.width, rect.height),
            action = action,
        ) {
            label(0, 0, rect.width, text, Colors.Text)
        }
    }

    private fun UiScope<ReactorControllerAction>.mixtureRow(
        table: UiRect,
        row: Int,
        name: ru.lazyhat.kraftui.foundation.Value<String>,
        amount: ru.lazyhat.kraftui.foundation.Value<String>,
    ) {
        label(table.x + 4, table.y + 30 + row * 10, table.width - 58, name, Colors.Text)
        label(table.x + table.width - 48, table.y + 30 + row * 10, 44, amount, Colors.Text, TextAlignment.End)
    }

    private fun UiScope<ReactorControllerAction>.metric(
        rect: UiRect,
        title: String,
        valueText: ru.lazyhat.kraftui.foundation.Value<String>,
    ) {
        panel(rect, Colors.Panel)
        label(rect.x + 4, rect.y + 4, rect.width - 8, title, Colors.Muted)
        label(rect.x + 4, rect.y + 16, rect.width - 8, valueText, Colors.Good)
    }

    private fun UiScope<ReactorControllerAction>.card(
        rect: UiRect,
        title: String,
        line0: ru.lazyhat.kraftui.foundation.Value<String>,
        line1: ru.lazyhat.kraftui.foundation.Value<String>,
        line2: ru.lazyhat.kraftui.foundation.Value<String> = value(""),
        line3: ru.lazyhat.kraftui.foundation.Value<String> = value(""),
    ) {
        panel(rect, Colors.Panel)
        label(rect.x + 4, rect.y + 4, rect.width - 8, title, Colors.Muted)
        label(rect.x + 4, rect.y + 16, rect.width - 8, line0, Colors.Text)
        label(rect.x + 4, rect.y + 26, rect.width - 8, line1, Colors.Text)
        label(rect.x + 4, rect.y + 36, rect.width - 8, line2, Colors.Text)
        label(rect.x + 4, rect.y + 46, rect.width - 8, line3, Colors.Text)
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
        text: String,
        color: Color,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        label(x, y, width, value(text), color, alignment)
    }

    private fun UiScope<ReactorControllerAction>.label(
        x: Int,
        y: Int,
        width: Int,
        text: ru.lazyhat.kraftui.foundation.Value<String>,
        color: Color,
        alignment: TextAlignment = TextAlignment.Start,
    ) {
        text(
            modifier =
                Modifier
                    .offset(x, y)
                    .size(width, 9)
                    .textAlign(alignment)
                    .textOverflow(TextOverflowPolicy.Ellipsize),
            color = color,
            text = text,
        )
    }

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
            val spare = available - baseWidth * count
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
        val Tab: Color = Color.rgb(214, 195, 158)
        val Title: Color = Color.rgb(61, 60, 72)
        val Text: Color = Color.rgb(61, 60, 72)
        val Muted: Color = Color.rgb(111, 106, 117)
        val Good: Color = Color.rgb(35, 134, 78)
    }
}
