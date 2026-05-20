package dev.makargravanov.create_thermodynamics.common.rust

import kotlin.test.Test
import kotlin.test.assertEquals

class NativeLibraryLoaderTest {
    @Test
    fun `maps windows amd64 to resource directory`() {
        assertEquals("windows-x86_64", NativeLibraryLoader.platformDirectory("Windows 11", "amd64"))
    }

    @Test
    fun `maps linux arm64 to resource directory`() {
        assertEquals("linux-aarch64", NativeLibraryLoader.platformDirectory("Linux", "arm64"))
    }

    @Test
    fun `maps macos arm64 to resource directory`() {
        assertEquals("macos-aarch64", NativeLibraryLoader.platformDirectory("Mac OS X", "aarch64"))
    }
}
