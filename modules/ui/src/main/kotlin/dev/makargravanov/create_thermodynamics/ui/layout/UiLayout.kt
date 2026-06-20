package dev.makargravanov.create_thermodynamics.ui.layout

fun interface UiTextMeasurer {
    fun width(text: String): Int

    fun lineHeight(): Int = 10
}

data class UiRect(
    val x: Int,
    val y: Int,
    val width: Int,
    val height: Int,
) {
    init {
        require(width >= 0) { "UI rectangle width must be non-negative" }
        require(height >= 0) { "UI rectangle height must be non-negative" }
    }

    val right: Int get() = x + width
    val bottom: Int get() = y + height

    fun contains(child: UiRect): Boolean =
        child.x >= x && child.y >= y && child.right <= right && child.bottom <= bottom

    fun intersects(other: UiRect): Boolean =
        x < other.right && right > other.x && y < other.bottom && bottom > other.y
}

data class UiInsets(
    val left: Int,
    val top: Int,
    val right: Int,
    val bottom: Int,
) {
    init {
        require(left >= 0 && top >= 0 && right >= 0 && bottom >= 0) {
            "UI insets must be non-negative"
        }
    }

    companion object {
        val Zero = UiInsets(0, 0, 0, 0)
    }
}

enum class UiAlignment {
    Start,
    Center,
    End,
}

sealed interface UiLength {
    data class Fixed(val pixels: Int) : UiLength {
        init {
            require(pixels >= 0) { "fixed UI length must be non-negative" }
        }
    }

    data class Weight(val weight: Int) : UiLength {
        init {
            require(weight > 0) { "weighted UI length must be positive" }
        }
    }

    data object Fill : UiLength
}

sealed interface TextOverflowPolicy {
    data object FailInValidation : TextOverflowPolicy
    data object Clip : TextOverflowPolicy
    data object Ellipsize : TextOverflowPolicy
    data object EllipsizeWithTooltip : TextOverflowPolicy
    data class WrapLines(val maxLines: Int) : TextOverflowPolicy {
        init {
            require(maxLines > 0) { "wrapped text must allow at least one line" }
        }
    }

    data object HorizontalScroll : TextOverflowPolicy
}

sealed interface UiNode {
    val id: String

    data class Text(
        override val id: String,
        val text: String,
        val color: Int,
        val overflow: TextOverflowPolicy,
        val alignment: UiAlignment = UiAlignment.Start,
        val width: UiLength = UiLength.Fill,
        val height: UiLength = UiLength.Fill,
    ) : UiNode

    data class Panel(
        override val id: String,
        val color: Int = 0x00000000,
        val width: UiLength = UiLength.Fill,
        val height: UiLength = UiLength.Fill,
        val children: List<UiNode>,
    ) : UiNode

    data class Row(
        override val id: String,
        val gap: Int = 0,
        val children: List<UiNode>,
    ) : UiNode {
        init {
            require(gap >= 0) { "UI row gap must be non-negative" }
        }
    }

    data class Column(
        override val id: String,
        val gap: Int = 0,
        val children: List<UiNode>,
    ) : UiNode {
        init {
            require(gap >= 0) { "UI column gap must be non-negative" }
        }
    }

    data class Place(
        override val id: String,
        val rect: UiRect,
        val child: UiNode,
    ) : UiNode

    data class Table(
        override val id: String,
        val columns: List<UiTableColumn>,
        val rows: List<List<String>>,
        val visibleRows: Int,
        val rowHeight: Int = 10,
        val columnGap: Int = 2,
        val textColor: Int = 0xFFFFFFFF.toInt(),
        val mutedTextColor: Int = 0xFF808080.toInt(),
    ) : UiNode {
        init {
            require(columns.isNotEmpty()) { "UI table must contain at least one column" }
            require(visibleRows >= 0) { "visible row count must be non-negative" }
            require(rowHeight > 0) { "table row height must be positive" }
            require(columnGap >= 0) { "table column gap must be non-negative" }
            require(rows.all { it.size == columns.size }) {
                "every UI table row must have ${columns.size} cells"
            }
        }
    }
}

data class UiTableColumn(
    val id: String,
    val title: String,
    val width: UiLength,
)

sealed interface UiDrawCommand {
    val id: String
    val bounds: UiRect

    data class DrawPanel(
        override val id: String,
        override val bounds: UiRect,
        val color: Int,
    ) : UiDrawCommand

    data class DrawText(
        override val id: String,
        override val bounds: UiRect,
        val text: String,
        val color: Int,
        val alignment: UiAlignment,
        val tooltip: String? = null,
    ) : UiDrawCommand
}

sealed interface UiLayoutDiagnostic {
    val nodeId: String

    data class TextWouldOverflow(
        override val nodeId: String,
        val text: String,
        val rect: UiRect,
        val textWidth: Int,
        val policy: TextOverflowPolicy,
    ) : UiLayoutDiagnostic

    data class NodeOutsideParent(
        override val nodeId: String,
        val parent: UiRect,
        val child: UiRect,
    ) : UiLayoutDiagnostic
}

data class UiLayoutResult(
    val commands: List<UiDrawCommand>,
    val diagnostics: List<UiLayoutDiagnostic>,
) {
    fun requireValidationClean() {
        require(diagnostics.isEmpty()) {
            diagnostics.joinToString(separator = "\n")
        }
    }
}

class UiLayoutEngine(
    private val textMeasurer: UiTextMeasurer,
) {
    fun layout(
        root: UiNode,
        bounds: UiRect,
    ): UiLayoutResult {
        val commands = mutableListOf<UiDrawCommand>()
        val diagnostics = mutableListOf<UiLayoutDiagnostic>()
        layoutNode(root, bounds, parent = bounds, commands, diagnostics)
        return UiLayoutResult(commands = commands.toList(), diagnostics = diagnostics.toList())
    }

    private fun layoutNode(
        node: UiNode,
        bounds: UiRect,
        parent: UiRect,
        commands: MutableList<UiDrawCommand>,
        diagnostics: MutableList<UiLayoutDiagnostic>,
    ) {
        if (!parent.contains(bounds)) {
            diagnostics += UiLayoutDiagnostic.NodeOutsideParent(node.id, parent, bounds)
        }
        when (node) {
            is UiNode.Text -> layoutText(node, bounds, commands, diagnostics)
            is UiNode.Panel -> {
                commands += UiDrawCommand.DrawPanel(node.id, bounds, node.color)
                node.children.forEach { child -> layoutNode(child, bounds, bounds, commands, diagnostics) }
            }
            is UiNode.Row -> layoutRow(node, bounds, commands, diagnostics)
            is UiNode.Column -> layoutColumn(node, bounds, commands, diagnostics)
            is UiNode.Place -> {
                val childBounds = UiRect(
                    x = bounds.x + node.rect.x,
                    y = bounds.y + node.rect.y,
                    width = node.rect.width,
                    height = node.rect.height,
                )
                layoutNode(node.child, childBounds, bounds, commands, diagnostics)
            }
            is UiNode.Table -> layoutTable(node, bounds, commands, diagnostics)
        }
    }

    private fun layoutText(
        node: UiNode.Text,
        bounds: UiRect,
        commands: MutableList<UiDrawCommand>,
        diagnostics: MutableList<UiLayoutDiagnostic>,
    ) {
        val textWidth = textMeasurer.width(node.text)
        val rendered =
            if (textWidth <= bounds.width) {
                RenderedText(node.text, tooltip = null)
            } else {
                overflowText(node, bounds, textWidth, diagnostics)
            }

        commands +=
            UiDrawCommand.DrawText(
                id = node.id,
                bounds = bounds,
                text = rendered.text,
                color = node.color,
                alignment = node.alignment,
                tooltip = rendered.tooltip,
            )
    }

    private fun overflowText(
        node: UiNode.Text,
        bounds: UiRect,
        textWidth: Int,
        diagnostics: MutableList<UiLayoutDiagnostic>,
    ): RenderedText =
        when (node.overflow) {
            TextOverflowPolicy.FailInValidation -> {
                diagnostics +=
                    UiLayoutDiagnostic.TextWouldOverflow(
                        nodeId = node.id,
                        text = node.text,
                        rect = bounds,
                        textWidth = textWidth,
                        policy = node.overflow,
                    )
                RenderedText(node.text, tooltip = null)
            }
            TextOverflowPolicy.Clip,
            TextOverflowPolicy.HorizontalScroll,
            is TextOverflowPolicy.WrapLines,
            -> RenderedText(node.text, tooltip = null)
            TextOverflowPolicy.Ellipsize -> RenderedText(ellipsize(node.text, bounds.width), tooltip = null)
            TextOverflowPolicy.EllipsizeWithTooltip ->
                RenderedText(ellipsize(node.text, bounds.width), tooltip = node.text)
        }

    private fun ellipsize(
        text: String,
        width: Int,
    ): String {
        val marker = "..."
        val markerWidth = textMeasurer.width(marker)
        if (width <= 0) {
            return ""
        }
        if (markerWidth >= width) {
            return marker.take(width / maxOf(1, textMeasurer.width(".")))
        }
        var end = text.length
        while (end > 0 && textMeasurer.width(text.take(end)) + markerWidth > width) {
            end--
        }
        return text.take(end) + marker
    }

    private fun layoutRow(
        node: UiNode.Row,
        bounds: UiRect,
        commands: MutableList<UiDrawCommand>,
        diagnostics: MutableList<UiLayoutDiagnostic>,
    ) {
        val widths = distribute(bounds.width, node.gap, node.children.map { it.widthLength() })
        var x = bounds.x
        for ((index, child) in node.children.withIndex()) {
            val childBounds = UiRect(x, bounds.y, widths[index], bounds.height)
            layoutNode(child, childBounds, bounds, commands, diagnostics)
            x += widths[index] + node.gap
        }
    }

    private fun layoutColumn(
        node: UiNode.Column,
        bounds: UiRect,
        commands: MutableList<UiDrawCommand>,
        diagnostics: MutableList<UiLayoutDiagnostic>,
    ) {
        val heights = distribute(bounds.height, node.gap, node.children.map { it.heightLength() })
        var y = bounds.y
        for ((index, child) in node.children.withIndex()) {
            val childBounds = UiRect(bounds.x, y, bounds.width, heights[index])
            layoutNode(child, childBounds, bounds, commands, diagnostics)
            y += heights[index] + node.gap
        }
    }

    private fun layoutTable(
        node: UiNode.Table,
        bounds: UiRect,
        commands: MutableList<UiDrawCommand>,
        diagnostics: MutableList<UiLayoutDiagnostic>,
    ) {
        val widths = distribute(bounds.width, node.columnGap, node.columns.map { it.width })
        var x = bounds.x
        for ((index, column) in node.columns.withIndex()) {
            val cellBounds = UiRect(x, bounds.y, widths[index], node.rowHeight)
            layoutText(
                UiNode.Text(
                    id = "${node.id}.header.${column.id}",
                    text = column.title,
                    color = node.mutedTextColor,
                    overflow = TextOverflowPolicy.FailInValidation,
                ),
                cellBounds,
                commands,
                diagnostics,
            )
            x += widths[index] + node.columnGap
        }

        for ((rowIndex, row) in node.rows.take(node.visibleRows).withIndex()) {
            x = bounds.x
            val y = bounds.y + node.rowHeight * (rowIndex + 1)
            for ((columnIndex, value) in row.withIndex()) {
                val column = node.columns[columnIndex]
                val cellBounds = UiRect(x, y, widths[columnIndex], node.rowHeight)
                layoutText(
                    UiNode.Text(
                        id = "${node.id}.row.$rowIndex.${column.id}",
                        text = value,
                        color = node.textColor,
                        overflow = TextOverflowPolicy.EllipsizeWithTooltip,
                    ),
                    cellBounds,
                    commands,
                    diagnostics,
                )
                x += widths[columnIndex] + node.columnGap
            }
        }
    }

    private fun distribute(
        total: Int,
        gap: Int,
        lengths: List<UiLength>,
    ): List<Int> {
        if (lengths.isEmpty()) {
            return emptyList()
        }
        val available = maxOf(0, total - gap * (lengths.size - 1))
        val fixed = lengths.sumOf { if (it is UiLength.Fixed) it.pixels else 0 }
        val weighted = lengths.sumOf {
            when (it) {
                is UiLength.Weight -> it.weight
                UiLength.Fill -> 1
                is UiLength.Fixed -> 0
            }
        }
        val remaining = maxOf(0, available - fixed)
        var used = 0
        val widths =
            lengths.map { length ->
                val width =
                    when (length) {
                        is UiLength.Fixed -> length.pixels
                        is UiLength.Weight -> remaining * length.weight / weighted
                        UiLength.Fill -> remaining / weighted
                    }
                used += width
                width
            }.toMutableList()
        var spare = available - used
        var index = widths.lastIndex
        while (spare > 0 && index >= 0) {
            if (lengths[index] !is UiLength.Fixed) {
                widths[index] += 1
                spare--
            }
            index--
            if (index < 0) {
                index = widths.lastIndex
            }
        }
        return widths
    }

    private fun UiNode.widthLength(): UiLength =
        when (this) {
            is UiNode.Text -> width
            is UiNode.Panel -> width
            else -> UiLength.Fill
        }

    private fun UiNode.heightLength(): UiLength =
        when (this) {
            is UiNode.Text -> height
            is UiNode.Panel -> height
            else -> UiLength.Fill
        }

    private data class RenderedText(
        val text: String,
        val tooltip: String?,
    )
}
