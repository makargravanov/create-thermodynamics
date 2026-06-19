package dev.makargravanov.create_thermodynamics.common.rust.blob

import java.nio.file.Files
import kotlin.test.Test
import kotlin.test.assertContentEquals
import kotlin.test.assertEquals
import kotlin.test.assertIs

class NativeBlobStorageTest {
    @Test
    fun `stores and loads native blob without unpacking domain data`() {
        val storage = NativeBlobStorage(Files.createTempDirectory("ct-native-blob-test"))
        val bytes = nativeBlobBytes(
            kind = NativeBlobKind.ReactorSnapshot,
            contentVersion = 7,
            modelVersion = "test-model",
            uncompressedByteSize = 128,
            compressedPayload = byteArrayOf(1, 2, 3, 4, 5),
        )

        val stored = assertIs<NativeBlobStorageResult.Stored>(
            storage.store(
                storageKey = "reactors/reactor_1_snapshot_000007.bin.zst",
                bytes = bytes,
            ),
        )
        val loaded = assertIs<NativeBlobStorageResult.Loaded>(storage.load(stored.ref))

        assertEquals(NativeBlobKind.ReactorSnapshot, stored.ref.kind)
        assertEquals(bytes.size.toLong(), stored.ref.encodedByteSize)
        assertEquals(5, stored.ref.compressedPayloadByteSize)
        assertEquals(128, stored.ref.uncompressedByteSize)
        assertContentEquals(bytes, loaded.bytes)
    }

    @Test
    fun `rejects storage key escaping root directory`() {
        val storage = NativeBlobStorage(Files.createTempDirectory("ct-native-blob-test"))

        val result = storage.store(
            storageKey = "../catalog.bin.zst",
            bytes = nativeBlobBytes(NativeBlobKind.DynamicCatalogCheckpoint),
        )

        val rejected = assertIs<NativeBlobStorageResult.Rejected>(result)
        assertEquals(NativeBlobStorageRejection.INVALID_STORAGE_KEY, rejected.reason)
    }

    @Test
    fun `rejects blob larger than configured limit`() {
        val storage = NativeBlobStorage(
            rootDirectory = Files.createTempDirectory("ct-native-blob-test"),
            limits = NativeBlobStorageLimits(maxEncodedBytes = 2),
        )

        val result = storage.store(
            storageKey = "catalog/delta_000002.bin.zst",
            bytes = nativeBlobBytes(NativeBlobKind.DynamicCatalogDelta),
        )

        val rejected = assertIs<NativeBlobStorageResult.Rejected>(result)
        assertEquals(NativeBlobStorageRejection.BLOB_TOO_LARGE, rejected.reason)
    }

    @Test
    fun `detects stored blob mutation by hash`() {
        val root = Files.createTempDirectory("ct-native-blob-test")
        val storage = NativeBlobStorage(root)
        val stored = assertIs<NativeBlobStorageResult.Stored>(
            storage.store(
                storageKey = "catalog/checkpoint_000003.bin.zst",
                bytes = nativeBlobBytes(
                    kind = NativeBlobKind.DynamicCatalogCheckpoint,
                    contentVersion = 3,
                    uncompressedByteSize = 64,
                    compressedPayload = byteArrayOf(9, 8, 7),
                ),
            ),
        )
        Files.write(
            root.resolve("catalog/checkpoint_000003.bin.zst"),
            nativeBlobBytes(
                kind = NativeBlobKind.DynamicCatalogCheckpoint,
                contentVersion = 3,
                uncompressedByteSize = 64,
                compressedPayload = byteArrayOf(9, 8, 6),
            ),
        )

        val result = storage.load(stored.ref)

        val rejected = assertIs<NativeBlobStorageResult.Rejected>(result)
        assertEquals(NativeBlobStorageRejection.HASH_MISMATCH, rejected.reason)
    }

    private fun nativeBlobBytes(
        kind: NativeBlobKind,
        contentVersion: Long = 1,
        modelVersion: String = "test-model",
        uncompressedByteSize: Long = 16,
        compressedPayload: ByteArray = byteArrayOf(1),
    ): ByteArray =
        buildList<Int> {
            addAll(byteArrayOf(0x43, 0x54, 0x4e, 0x42, 0x4c, 0x42, 0x31, 0x00).map { it.toInt() and 0xff })
            addU16(1)
            addU16(kind.wireId)
            addU16(modelVersion.encodeToByteArray().size)
            addU64(contentVersion)
            addU64(uncompressedByteSize)
            addU64(compressedPayload.size.toLong())
            repeat(32) { add(0) }
            addAll(modelVersion.encodeToByteArray().map { it.toInt() and 0xff })
            addAll(compressedPayload.map { it.toInt() and 0xff })
        }.map { it.toByte() }.toByteArray()

    private fun MutableList<Int>.addU16(value: Int) {
        add(value and 0xff)
        add((value ushr 8) and 0xff)
    }

    private fun MutableList<Int>.addU64(value: Long) {
        for (index in 0 until 8) {
            add(((value ushr (index * 8)) and 0xff).toInt())
        }
    }
}
