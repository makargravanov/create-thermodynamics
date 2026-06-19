package dev.makargravanov.create_thermodynamics.common.rust.blob

private val HEX_64 = Regex("[0-9a-f]{64}")

enum class NativeBlobKind(val wireId: Int) {
    DynamicCatalogCheckpoint(1),
    DynamicCatalogDelta(2),
    ReactorSnapshot(3),
    ReactionSpaceSnapshot(4),
    PlannerSnapshot(5),
}

@JvmInline
value class NativeBlobHash(val value: String) {
    init {
        require(HEX_64.matches(value)) {
            "native blob hash must be a lowercase 64-character hexadecimal SHA-256 string"
        }
    }
}

data class NativeBlobRef(
    val kind: NativeBlobKind,
    val schemaVersion: Int,
    val modelVersion: String,
    val contentVersion: Long,
    val payloadHash: NativeBlobHash,
    val encodedHash: NativeBlobHash,
    val encodedByteSize: Long,
    val compressedPayloadByteSize: Long,
    val uncompressedByteSize: Long,
    val storageKey: String,
) {
    init {
        require(schemaVersion > 0) { "native blob schemaVersion must be positive" }
        require(modelVersion.isNotBlank()) { "native blob modelVersion must not be blank" }
        require(contentVersion >= 0) { "native blob contentVersion must be non-negative" }
        require(encodedByteSize > 0) { "native blob encodedByteSize must be positive" }
        require(compressedPayloadByteSize > 0) { "native blob compressedPayloadByteSize must be positive" }
        require(uncompressedByteSize >= 0) { "native blob uncompressedByteSize must be non-negative" }
        require(storageKey.isNotBlank()) { "native blob storageKey must not be blank" }
    }
}
