package dev.makargravanov.create_thermodynamics.common.rust.blob

private val NATIVE_BLOB_MAGIC = byteArrayOf(0x43, 0x54, 0x4e, 0x42, 0x4c, 0x42, 0x31, 0x00)
private const val FIXED_HEADER_BYTES = 8 + 2 + 2 + 2 + 8 + 8 + 8 + 32
private const val MAX_MODEL_VERSION_BYTES = 256

data class NativeBlobHeader(
    val kind: NativeBlobKind,
    val schemaVersion: Int,
    val modelVersion: String,
    val contentVersion: Long,
    val payloadHash: NativeBlobHash,
    val compressedPayloadByteSize: Long,
    val uncompressedByteSize: Long,
) {
    companion object {
        fun read(bytes: ByteArray): NativeBlobHeader {
            require(bytes.size >= FIXED_HEADER_BYTES) {
                "native blob is too short: ${bytes.size} bytes"
            }
            require(bytes.copyOfRange(0, NATIVE_BLOB_MAGIC.size).contentEquals(NATIVE_BLOB_MAGIC)) {
                "native blob has invalid magic"
            }
            val cursor = Cursor(bytes)
            cursor.skip(NATIVE_BLOB_MAGIC.size)
            val schemaVersion = cursor.readU16()
            val kind = kindFromWire(cursor.readU16())
            val modelVersionByteSize = cursor.readU16()
            require(modelVersionByteSize <= MAX_MODEL_VERSION_BYTES) {
                "native blob model version is too long: $modelVersionByteSize bytes"
            }
            val contentVersion = cursor.readU64()
            val uncompressedByteSize = cursor.readU64()
            val compressedPayloadByteSize = cursor.readU64()
            val payloadHash = NativeBlobHash(cursor.readBytes(32).toHex())
            val modelVersion = cursor.readBytes(modelVersionByteSize).decodeToString()
            val expectedSize = FIXED_HEADER_BYTES + modelVersionByteSize + compressedPayloadByteSize
            require(expectedSize <= Int.MAX_VALUE) { "native blob encoded size does not fit into Int" }
            require(bytes.size == expectedSize.toInt()) {
                "native blob has ${bytes.size} encoded bytes, header declares $expectedSize"
            }
            require(modelVersion.isNotBlank()) { "native blob model version must not be blank" }
            return NativeBlobHeader(
                kind = kind,
                schemaVersion = schemaVersion,
                modelVersion = modelVersion,
                contentVersion = contentVersion,
                payloadHash = payloadHash,
                compressedPayloadByteSize = compressedPayloadByteSize,
                uncompressedByteSize = uncompressedByteSize,
            )
        }

        private fun kindFromWire(value: Int): NativeBlobKind =
            NativeBlobKind.entries.firstOrNull { it.wireId == value }
                ?: throw IllegalArgumentException("unknown native blob kind $value")
    }
}

private class Cursor(
    private val bytes: ByteArray,
) {
    private var offset = 0

    fun skip(count: Int) {
        offset += count
    }

    fun readU16(): Int {
        val result = (bytes[offset].toInt() and 0xff) or
            ((bytes[offset + 1].toInt() and 0xff) shl 8)
        offset += 2
        return result
    }

    fun readU64(): Long {
        var result = 0L
        for (index in 0 until 8) {
            result = result or ((bytes[offset + index].toLong() and 0xffL) shl (index * 8))
        }
        offset += 8
        return result
    }

    fun readBytes(count: Int): ByteArray {
        val result = bytes.copyOfRange(offset, offset + count)
        offset += count
        return result
    }
}

internal fun ByteArray.toHex(): String =
    joinToString(separator = "") { "%02x".format(it.toInt() and 0xff) }
