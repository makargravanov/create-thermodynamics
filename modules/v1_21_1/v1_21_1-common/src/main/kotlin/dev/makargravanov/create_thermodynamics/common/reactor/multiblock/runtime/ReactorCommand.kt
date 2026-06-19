package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobRef

sealed interface ReactorCommand {
    val commandId: ReactorCommandId
    val structureId: ReactorStructureId
    val expectedSnapshotVersion: ReactorSnapshotVersion?

    data class EnsureLoaded(
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val expectedSnapshotVersion: ReactorSnapshotVersion?,
        val checkpoint: NativeBlobRef?,
    ) : ReactorCommand

    data class UnloadAndExportSnapshot(
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val expectedSnapshotVersion: ReactorSnapshotVersion,
        val reason: String,
    ) : ReactorCommand {
        init {
            require(reason.isNotBlank()) { "snapshot export reason must not be blank" }
        }
    }

    data class RemoveReactor(
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val expectedSnapshotVersion: ReactorSnapshotVersion?,
        val reason: String,
    ) : ReactorCommand {
        init {
            require(reason.isNotBlank()) { "reactor removal reason must not be blank" }
        }
    }

    data class Tick(
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val expectedSnapshotVersion: ReactorSnapshotVersion,
        val dtSeconds: Double,
    ) : ReactorCommand {
        init {
            require(dtSeconds.isFinite() && dtSeconds >= 0.0) {
                "reactor tick duration must be non-negative and finite"
            }
        }
    }

    data class InsertItem(
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val expectedSnapshotVersion: ReactorSnapshotVersion,
        val portPosition: ReactorBlockPosition,
        val itemId: String,
        val itemCount: Int,
    ) : ReactorCommand {
        init {
            require(itemId.isNotBlank()) { "reactor item id must not be blank" }
            require(itemCount > 0) { "reactor item count must be positive" }
        }
    }

    data class RequestSnapshot(
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val expectedSnapshotVersion: ReactorSnapshotVersion?,
        val reason: String,
    ) : ReactorCommand {
        init {
            require(reason.isNotBlank()) { "snapshot request reason must not be blank" }
        }
    }

    data class RequestMetrics(
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val expectedSnapshotVersion: ReactorSnapshotVersion?,
        val mask: Int,
    ) : ReactorCommand {
        init {
            require(mask >= 0) { "reactor metrics mask must be non-negative" }
        }
    }
}
