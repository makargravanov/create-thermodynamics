package dev.makargravanov.create_thermodynamics.common.chemistry.binding

import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative
import kotlin.test.Test
import kotlin.test.assertEquals
import kotlin.test.assertFailsWith
import kotlin.test.assertTrue

class DefaultItemChemicalBindingsTest {
    @Test
    fun `default item bindings can be loaded into native registry`() {
        val staticSubstances = ThermodynamicsNative.staticSubstanceIds().toSet()
        assertTrue(DefaultItemChemicalBindings.bindings.all { it.substanceId in staticSubstances })

        ThermodynamicsNative.configureItemChemicalBindings(DefaultItemChemicalBindings.bindings)

        assertEquals(DefaultItemChemicalBindings.bindings.size, ThermodynamicsNative.itemChemicalBindingCount())
        assertTrue(ThermodynamicsNative.hasItemChemicalBinding("minecraft:water_bucket"))
    }

    @Test
    fun `item binding rejects invalid amounts`() {
        assertFailsWith<IllegalArgumentException> {
            ItemChemicalBinding(
                itemId = "minecraft:test",
                substanceId = "destroy:water",
                molPerItem = Double.NaN,
            )
        }
        assertFailsWith<IllegalArgumentException> {
            ItemChemicalBinding(
                itemId = "minecraft:test",
                substanceId = "destroy:water",
                molPerItem = 0.0,
            )
        }
    }
}
