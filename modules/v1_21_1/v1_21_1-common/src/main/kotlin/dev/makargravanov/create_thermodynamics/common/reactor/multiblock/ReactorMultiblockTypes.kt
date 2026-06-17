package dev.makargravanov.create_thermodynamics.common.reactor.multiblock

import java.util.UUID

enum class ReactorMultiblockBlockKind {
    CONTROLLER,
    CHAMBER,
    ITEM_INPUT_PORT,
    ITEM_OUTPUT_PORT,
    FLUID_INPUT_PORT,
    FLUID_OUTPUT_PORT,
}

enum class ReactorPortKind {
    ITEM_INPUT,
    ITEM_OUTPUT,
    FLUID_INPUT,
    FLUID_OUTPUT,
}

data class ReactorMultiblockBlock(
    val position: ReactorBlockPosition,
    val kind: ReactorMultiblockBlockKind,
)

data class ReactorMultiblockRules(
    val chamberVolumeCubicMeters: Double,
    val minimumChamberBlocks: Int = 1,
    val maximumChamberBlocks: Int? = null,
) {
    init {
        require(chamberVolumeCubicMeters.isFinite() && chamberVolumeCubicMeters > 0.0) {
            "chamberVolumeCubicMeters must be positive and finite"
        }
        require(minimumChamberBlocks > 0) {
            "minimumChamberBlocks must be positive"
        }
        require(maximumChamberBlocks == null || maximumChamberBlocks >= minimumChamberBlocks) {
            "maximumChamberBlocks must be absent or not smaller than minimumChamberBlocks"
        }
    }
}

data class ReactorZoneDescriptor(
    val zoneIndex: Int,
    val chamberPositions: Set<ReactorBlockPosition>,
    val volumeCubicMeters: Double,
)

data class ReactorPortDescriptor(
    val portIndex: Int,
    val kind: ReactorPortKind,
    val position: ReactorBlockPosition,
    val zoneIndex: Int,
    val attachedChamberPosition: ReactorBlockPosition,
)

data class ReactorMultiblockDefinition(
    val structureId: UUID,
    val controllerPosition: ReactorBlockPosition,
    val zone: ReactorZoneDescriptor,
    val ports: List<ReactorPortDescriptor>,
) {
    val totalVolumeCubicMeters: Double
        get() = zone.volumeCubicMeters

    fun portsOfKind(kind: ReactorPortKind): List<ReactorPortDescriptor> =
        ports.filter { it.kind == kind }
}

class ReactorMultiblockValidationException(
    val validationErrors: List<String>,
) : IllegalArgumentException(validationErrors.joinToString(separator = "; "))
