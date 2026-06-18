package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model

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
    val facing: ReactorBlockDirection? = null,
)

data class ReactorZoneDescriptor(
    val zoneIndex: Int,
    val volumePositions: Set<ReactorBlockPosition>,
    val plainChamberPositions: Set<ReactorBlockPosition>,
    val volumeCubicMeters: Double,
)

data class ReactorPortDescriptor(
    val portIndex: Int,
    val kind: ReactorPortKind,
    val position: ReactorBlockPosition,
    val zoneIndex: Int,
    val attachedChamberPosition: ReactorBlockPosition,
    val contactDirection: ReactorBlockDirection,
)

@JvmInline
value class ReactorStructureId(val value: UUID)

data class ReactorMultiblockDefinition(
    val structureId: ReactorStructureId,
    val controllerPosition: ReactorBlockPosition,
    val controllerContactDirection: ReactorBlockDirection?,
    val zone: ReactorZoneDescriptor,
    val ports: List<ReactorPortDescriptor>,
    val inactiveChamberPositions: Set<ReactorBlockPosition> = emptySet(),
) {
    val totalVolumeCubicMeters: Double
        get() = zone.volumeCubicMeters

    fun portsOfKind(kind: ReactorPortKind): List<ReactorPortDescriptor> =
        ports.filter { it.kind == kind }
}

class ReactorMultiblockValidationException(
    val validationErrors: List<String>,
) : IllegalArgumentException(validationErrors.joinToString(separator = "; "))
