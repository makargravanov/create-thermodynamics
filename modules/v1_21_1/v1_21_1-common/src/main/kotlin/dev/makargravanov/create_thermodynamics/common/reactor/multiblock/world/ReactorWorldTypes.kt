package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockDirection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId

data class ReactorWorldBlockSnapshot(
    val position: ReactorBlockPosition,
    val kind: ReactorMultiblockBlockKind?,
    val facing: ReactorBlockDirection?,
    val loaded: Boolean,
) {
    init {
        require(loaded || kind == null) { "unloaded reactor world snapshot must not contain a block kind" }
    }
}

data class ReactorAssemblyScan(
    val blocks: Set<ReactorWorldBlockSnapshot>,
    val hasUnknownBoundary: Boolean,
) {
    val loadedReactorPositions: Set<ReactorBlockPosition> =
        blocks.filterTo(linkedSetOf()) { it.loaded && it.kind != null }.mapTo(linkedSetOf()) { it.position }
}

data class ReactorAssemblyDiagnostic(
    val controllerPosition: ReactorBlockPosition?,
    val errors: List<String>,
) {
    init {
        require(errors.isNotEmpty()) { "reactor assembly diagnostic must contain at least one error" }
    }
}

data class ReactorAssemblyPlan(
    val definitions: List<ReactorMultiblockDefinition>,
    val memberships: Map<ReactorBlockPosition, ReactorBlockMembership>,
    val diagnostics: List<ReactorAssemblyDiagnostic>,
    val hasUnknownBoundary: Boolean,
)

enum class ReactorBlockRole {
    CONTROLLER,
    CHAMBER,
    PORT,
}

data class ReactorStructureSummary(
    val structureId: ReactorStructureId,
    val zoneCount: Int,
    val chamberBlockCount: Int,
    val portCount: Int,
) {
    init {
        require(zoneCount > 0) { "formed reactor structure must contain at least one zone" }
        require(chamberBlockCount > 0) { "formed reactor structure must contain at least one volume block" }
        require(portCount >= 0) { "reactor port count must be non-negative" }
    }
}

data class ReactorBlockMembership(
    val structureId: ReactorStructureId,
    val role: ReactorBlockRole,
    val activeVolumeBlock: Boolean,
    val summary: ReactorStructureSummary,
) {
    init {
        require(structureId == summary.structureId) { "reactor block membership and summary structure ids differ" }
    }
}

enum class ReactorControllerFormationState {
    FORMED,
    NOT_FORMED,
    UNKNOWN,
}

data class ReactorControllerViewState(
    val formationState: ReactorControllerFormationState,
    val structureId: ReactorStructureId?,
    val zoneCount: Int,
    val chamberBlockCount: Int,
    val portCount: Int,
    val diagnostic: String?,
    val nativeBinding: String = "pending",
    val temperatureKelvin: Double? = null,
    val pressurePascal: Double? = null,
    val mixture: List<ReactorMixtureViewEntry> = emptyList(),
) {
    init {
        if (formationState == ReactorControllerFormationState.FORMED) {
            require(zoneCount > 0) { "formed reactor controller view state must contain a zone" }
            require(chamberBlockCount > 0) { "formed reactor controller view state must contain chamber volume" }
        }
        require(temperatureKelvin == null || temperatureKelvin.isFinite()) {
            "temperatureKelvin must be finite when present"
        }
        require(pressurePascal == null || pressurePascal.isFinite()) {
            "pressurePascal must be finite when present"
        }
    }

    companion object {
        val NotFormed = ReactorControllerViewState(
            formationState = ReactorControllerFormationState.NOT_FORMED,
            structureId = null,
            zoneCount = 0,
            chamberBlockCount = 0,
            portCount = 0,
            diagnostic = null,
        )
    }
}

data class ReactorMixtureViewEntry(
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
