package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobRef

sealed interface ReactorReport {
    val reportId: ReactorReportId
    val commandId: ReactorCommandId?
    val structureId: ReactorStructureId
    val snapshotVersion: ReactorSnapshotVersion
    val requiredDelivery: Boolean

    data class CommandAccepted(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true
    }

    data class CommandRejected(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val reason: String,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true

        init {
            require(reason.isNotBlank()) { "reactor command rejection reason must not be blank" }
        }
    }

    data class TickCompleted(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId?,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val metrics: ReactorTickMetrics,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = false
    }

    data class SnapshotReady(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId?,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val checkpoint: NativeBlobRef,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true
    }

    data class PortInputAccepted(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val portPosition: ReactorBlockPosition,
        val itemId: String,
        val acceptedCount: Int,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true

        init {
            require(itemId.isNotBlank()) { "accepted item id must not be blank" }
            require(acceptedCount > 0) { "accepted item count must be positive" }
        }
    }

    data class PortInputRejected(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val portPosition: ReactorBlockPosition,
        val reason: String,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true

        init {
            require(reason.isNotBlank()) { "port input rejection reason must not be blank" }
        }
    }

    data class PortOutputAccepted(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val portPosition: ReactorBlockPosition,
        val itemId: String,
        val extractedCount: Int,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true

        init {
            require(itemId.isNotBlank()) { "output item id must not be blank" }
            require(extractedCount >= 0) { "extracted output item count must be non-negative" }
        }
    }

    data class PortOutputRejected(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val portPosition: ReactorBlockPosition,
        val reason: String,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true

        init {
            require(reason.isNotBlank()) { "port output rejection reason must not be blank" }
        }
    }

    data class DiagnosticsUpdated(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId?,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val entries: List<String>,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = false

        init {
            require(entries.all { it.isNotBlank() }) { "reactor diagnostic entry must not be blank" }
        }
    }

    data class ReactorFailed(
        override val reportId: ReactorReportId,
        override val commandId: ReactorCommandId?,
        override val structureId: ReactorStructureId,
        override val snapshotVersion: ReactorSnapshotVersion,
        val message: String,
    ) : ReactorReport {
        override val requiredDelivery: Boolean = true

        init {
            require(message.isNotBlank()) { "reactor failure message must not be blank" }
        }
    }
}

data class ReactorTickMetrics(
    val simulatedSeconds: Double,
    val temperatureKelvin: Double?,
    val pressurePascal: Double?,
    val substances: List<ReactorMixtureSubstanceMetric> = emptyList(),
) {
    init {
        require(simulatedSeconds.isFinite() && simulatedSeconds >= 0.0) {
            "simulatedSeconds must be non-negative and finite"
        }
        require(temperatureKelvin == null || temperatureKelvin.isFinite()) {
            "temperatureKelvin must be finite when present"
        }
        require(pressurePascal == null || pressurePascal.isFinite()) {
            "pressurePascal must be finite when present"
        }
    }
}

data class ReactorMixtureSubstanceMetric(
    val substanceId: String,
    val concentrationMolPerBucket: Double,
) {
    init {
        require(substanceId.isNotBlank()) { "substanceId must not be blank" }
        require(concentrationMolPerBucket.isFinite() && concentrationMolPerBucket >= 0.0) {
            "concentrationMolPerBucket must be non-negative and finite"
        }
    }
}
