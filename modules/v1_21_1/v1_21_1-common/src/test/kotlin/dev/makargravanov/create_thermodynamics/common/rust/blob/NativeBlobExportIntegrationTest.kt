package dev.makargravanov.create_thermodynamics.common.rust.blob

import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative
import java.nio.file.Files
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertIs
import kotlin.test.assertTrue

class NativeBlobExportIntegrationTest {
    @Test
    fun `catalog checkpoint exported by native layer is stored as opaque blob`() {
        val bytes = ThermodynamicsNative.exportCatalogCheckpoint(contentVersion = 1)
        val storage = NativeBlobStorage(Files.createTempDirectory("ct-native-blob-export-test"))

        val stored = assertIs<NativeBlobStorageResult.Stored>(
            storage.store(
                storageKey = "catalog/checkpoint_000001.bin.zst",
                bytes = bytes,
            ),
        )
        val loaded = assertIs<NativeBlobStorageResult.Loaded>(storage.load(stored.ref))

        assertTrue(bytes.isNotEmpty())
        assertEquals(NativeBlobKind.DynamicCatalogCheckpoint, stored.ref.kind)
        assertEquals("create-thermodynamics:dynamic-catalog:1", stored.ref.modelVersion)
        assertContentEquals(bytes, loaded.bytes)
    }

    @Test
    fun `reactor checkpoint exported by native layer is stored as opaque blob`() {
        val reactorId = ThermodynamicsNative.createSingleZoneReactor(
            volumeCubicMeters = 0.001,
            itemInputCount = 1,
            itemOutputCount = 1,
            fluidInputCount = 0,
            fluidOutputCount = 0,
        )
        try {
            val bytes = ThermodynamicsNative.exportReactorCheckpoint(reactorId, contentVersion = 1)
            val storage = NativeBlobStorage(Files.createTempDirectory("ct-native-blob-export-test"))

            val stored = assertIs<NativeBlobStorageResult.Stored>(
                storage.store(
                    storageKey = "reactors/reactor_000001_snapshot_000001.bin.zst",
                    bytes = bytes,
                ),
            )
            val loaded = assertIs<NativeBlobStorageResult.Loaded>(storage.load(stored.ref))

            assertTrue(bytes.isNotEmpty())
            assertEquals(NativeBlobKind.ReactorSnapshot, stored.ref.kind)
            assertEquals("create-thermodynamics:reactor:1", stored.ref.modelVersion)
            assertContentEquals(bytes, loaded.bytes)
        } finally {
            ThermodynamicsNative.removeReactor(reactorId)
        }
    }
}
