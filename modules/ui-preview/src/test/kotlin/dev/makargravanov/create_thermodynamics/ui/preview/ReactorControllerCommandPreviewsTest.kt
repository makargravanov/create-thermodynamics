package dev.makargravanov.create_thermodynamics.ui.preview

import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class ReactorControllerCommandPreviewsTest {
    @Test
    fun `reactor controller command previews are clean and non blank`() {
        val previews = ReactorControllerCommandPreviews.all()

        assertTrue(previews.size >= 3)
        for (preview in previews) {
            assertEquals(emptyList(), preview.diagnostics, "${preview.id} should not have layout diagnostics")
            val image = CommandPreviewRenderer.render(preview.width, preview.height, preview.commands)
            assertTrue(CommandPreviewRenderer.countDistinctColors(image) > 1, "${preview.id} should not be blank")
        }
    }
}
