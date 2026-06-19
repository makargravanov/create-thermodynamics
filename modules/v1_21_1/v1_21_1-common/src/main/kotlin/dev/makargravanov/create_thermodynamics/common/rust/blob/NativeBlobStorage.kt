package dev.makargravanov.create_thermodynamics.common.rust.blob

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption.ATOMIC_MOVE
import java.nio.file.StandardOpenOption.CREATE_NEW
import java.nio.file.StandardOpenOption.WRITE
import java.security.MessageDigest
import java.util.UUID

data class NativeBlobStorageLimits(
    val maxEncodedBytes: Long = 64L * 1024L * 1024L,
) {
    init {
        require(maxEncodedBytes > 0) { "maxEncodedBytes must be positive" }
    }
}

sealed interface NativeBlobStorageResult {
    data class Stored(val ref: NativeBlobRef) : NativeBlobStorageResult
    data class Loaded(val bytes: ByteArray) : NativeBlobStorageResult {
        override fun equals(other: Any?): Boolean =
            other is Loaded && bytes.contentEquals(other.bytes)

        override fun hashCode(): Int = bytes.contentHashCode()
    }
    data class Rejected(val reason: NativeBlobStorageRejection, val message: String) : NativeBlobStorageResult
}

enum class NativeBlobStorageRejection {
    INVALID_STORAGE_KEY,
    INVALID_BLOB_HEADER,
    BLOB_TOO_LARGE,
    HASH_MISMATCH,
    SIZE_MISMATCH,
    FILE_NOT_FOUND,
    IO_FAILURE,
}

class NativeBlobStorage(
    rootDirectory: Path,
    private val limits: NativeBlobStorageLimits = NativeBlobStorageLimits(),
) {
    private val rootDirectory: Path = rootDirectory.toAbsolutePath().normalize()

    fun store(
        storageKey: String,
        bytes: ByteArray,
    ): NativeBlobStorageResult {
        if (bytes.isEmpty()) {
            return rejected(NativeBlobStorageRejection.SIZE_MISMATCH, "native blob bytes must not be empty")
        }
        if (bytes.size.toLong() > limits.maxEncodedBytes) {
            return rejected(
                NativeBlobStorageRejection.BLOB_TOO_LARGE,
                "native blob has ${bytes.size} bytes, limit is ${limits.maxEncodedBytes}",
            )
        }
        val target = resolveStorageKey(storageKey)
            ?: return rejected(NativeBlobStorageRejection.INVALID_STORAGE_KEY, "invalid native blob storage key: $storageKey")
        val header = try {
            NativeBlobHeader.read(bytes)
        } catch (error: IllegalArgumentException) {
            return rejected(NativeBlobStorageRejection.INVALID_BLOB_HEADER, error.message ?: "invalid native blob header")
        }

        val ref = try {
            NativeBlobRef(
                kind = header.kind,
                schemaVersion = header.schemaVersion,
                modelVersion = header.modelVersion,
                contentVersion = header.contentVersion,
                payloadHash = header.payloadHash,
                encodedHash = sha256(bytes),
                encodedByteSize = bytes.size.toLong(),
                compressedPayloadByteSize = header.compressedPayloadByteSize,
                uncompressedByteSize = header.uncompressedByteSize,
                storageKey = storageKey,
            )
        } catch (error: IllegalArgumentException) {
            return rejected(NativeBlobStorageRejection.INVALID_STORAGE_KEY, error.message ?: "invalid native blob reference")
        }

        return try {
            Files.createDirectories(target.parent)
            val temporary = target.resolveSibling(".${target.fileName}.tmp-${UUID.randomUUID()}")
            try {
                Files.write(temporary, bytes, CREATE_NEW, WRITE)
                Files.move(temporary, target, ATOMIC_MOVE)
            } finally {
                Files.deleteIfExists(temporary)
            }
            NativeBlobStorageResult.Stored(ref)
        } catch (error: Exception) {
            rejected(NativeBlobStorageRejection.IO_FAILURE, "failed to write native blob $storageKey: ${error.message}")
        }
    }

    fun load(ref: NativeBlobRef): NativeBlobStorageResult {
        if (ref.encodedByteSize > limits.maxEncodedBytes) {
            return rejected(
                NativeBlobStorageRejection.BLOB_TOO_LARGE,
                "native blob ${ref.storageKey} declares ${ref.encodedByteSize} bytes, limit is ${limits.maxEncodedBytes}",
            )
        }
        val target = resolveStorageKey(ref.storageKey)
            ?: return rejected(NativeBlobStorageRejection.INVALID_STORAGE_KEY, "invalid native blob storage key: ${ref.storageKey}")
        if (!Files.exists(target)) {
            return rejected(NativeBlobStorageRejection.FILE_NOT_FOUND, "native blob ${ref.storageKey} does not exist")
        }

        return try {
            val actualSize = Files.size(target)
            if (actualSize != ref.encodedByteSize) {
                return rejected(
                    NativeBlobStorageRejection.SIZE_MISMATCH,
                    "native blob ${ref.storageKey} has $actualSize bytes, expected ${ref.encodedByteSize}",
                )
            }
            val bytes = Files.readAllBytes(target)
            val header = try {
                NativeBlobHeader.read(bytes)
            } catch (error: IllegalArgumentException) {
                return rejected(NativeBlobStorageRejection.INVALID_BLOB_HEADER, error.message ?: "invalid native blob header")
            }
            if (header.kind != ref.kind ||
                header.schemaVersion != ref.schemaVersion ||
                header.modelVersion != ref.modelVersion ||
                header.contentVersion != ref.contentVersion ||
                header.payloadHash != ref.payloadHash ||
                header.compressedPayloadByteSize != ref.compressedPayloadByteSize ||
                header.uncompressedByteSize != ref.uncompressedByteSize
            ) {
                return rejected(NativeBlobStorageRejection.INVALID_BLOB_HEADER, "native blob ${ref.storageKey} header does not match reference")
            }
            val actualHash = sha256(bytes)
            if (actualHash != ref.encodedHash) {
                return rejected(NativeBlobStorageRejection.HASH_MISMATCH, "native blob ${ref.storageKey} hash mismatch")
            }
            NativeBlobStorageResult.Loaded(bytes)
        } catch (error: Exception) {
            rejected(NativeBlobStorageRejection.IO_FAILURE, "failed to read native blob ${ref.storageKey}: ${error.message}")
        }
    }

    private fun resolveStorageKey(storageKey: String): Path? {
        if (storageKey.isBlank() || storageKey.contains('\\')) {
            return null
        }
        val parts = storageKey.split('/')
        if (parts.any { it.isBlank() || it == "." || it == ".." }) {
            return null
        }
        val relative = Path.of(storageKey)
        if (relative.isAbsolute) {
            return null
        }
        val resolved = rootDirectory.resolve(relative).normalize()
        return resolved.takeIf { it.startsWith(rootDirectory) }
    }

    private fun rejected(
        reason: NativeBlobStorageRejection,
        message: String,
    ): NativeBlobStorageResult.Rejected =
        NativeBlobStorageResult.Rejected(reason, message)

    private fun sha256(bytes: ByteArray): NativeBlobHash {
        val digest = MessageDigest.getInstance("SHA-256").digest(bytes)
        return NativeBlobHash(digest.toHex())
    }
}
