package dev.makargravanov.create_thermodynamics.ui.preview

import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.foundation.modifier.Modifier
import ru.lazyhat.kraftui.foundation.modifier.background
import ru.lazyhat.kraftui.foundation.modifier.offset
import ru.lazyhat.kraftui.foundation.modifier.size
import ru.lazyhat.kraftui.foundation.ui
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertNotEquals

class UiPreviewRendererTest {
    @Test
    fun rendersBackgroundIntoImage() {
        val image =
            UiPreviewRenderer.render(
                UiPreviewSpec(
                    id = "background",
                    width = 8,
                    height = 8,
                    root = ui(Modifier.size(8, 8).background(Color.rgb(10, 20, 30))) {},
                ),
            )

        assertEquals(0xFF0A141E.toInt(), image.getRGB(4, 4))
    }

    @Test
    fun rendersChildOffsetInsideParent() {
        val image =
            UiPreviewRenderer.render(
                UiPreviewSpec(
                    id = "offset",
                    width = 16,
                    height = 16,
                    root =
                        ui(Modifier.size(16, 16).background(Color.Black)) {
                            box(Modifier.offset(4, 5).size(3, 2).background(Color.Red))
                        },
                ),
            )

        assertEquals(0xFFFF0000.toInt(), image.getRGB(4, 5))
        assertEquals(0xFFFF0000.toInt(), image.getRGB(6, 6))
        assertNotEquals(0xFFFF0000.toInt(), image.getRGB(7, 6))
    }

    @Test
    fun reactorPortPreviewIsNotBlank() {
        val image = UiPreviewRenderer.render(ReactorPortPreviews.all().single())
        val pixels =
            sequence {
                for (y in 0 until image.height) {
                    for (x in 0 until image.width) {
                        yield(image.getRGB(x, y))
                    }
                }
            }.toSet()

        assert(pixels.size > 4) { "preview should contain several colors" }
    }
}
