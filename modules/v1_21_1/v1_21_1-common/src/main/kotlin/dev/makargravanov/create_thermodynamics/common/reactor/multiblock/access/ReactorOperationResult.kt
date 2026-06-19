package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.access

import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobRef

enum class ReactorOperationRejection {
    STRUCTURE_NOT_FOUND,
    STRUCTURE_NOT_ACTIVE,
    STRUCTURE_NOT_SUSPENDED,
    PORT_NOT_FOUND,
    WRONG_PORT_KIND,
    INVALID_ITEM_ID,
    INVALID_ITEM_COUNT,
    ITEM_BUFFER_FULL,
    INVALID_CONTENT_VERSION,
    STALE_REPORT,
    SNAPSHOT_STORAGE_REJECTED,
    OPERATION_NOT_SUPPORTED,
}

sealed interface ReactorOperationResult {
    data object Completed : ReactorOperationResult

    data class ItemInserted(
        val molInserted: Double,
    ) : ReactorOperationResult

    data class ItemBuffered(
        val itemId: String,
        val itemCount: Int,
        val message: String,
    ) : ReactorOperationResult

    data class ReactorSuspended(
        val message: String,
    ) : ReactorOperationResult

    data class ReactorResumed(
        val message: String,
    ) : ReactorOperationResult

    data class ReactorCheckpointExported(
        val checkpoint: NativeBlobRef,
    ) : ReactorOperationResult

    data class Rejected(
        val reason: ReactorOperationRejection,
        val message: String,
    ) : ReactorOperationResult

    data class Failed(
        val message: String,
    ) : ReactorOperationResult
}
