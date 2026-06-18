package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId

class ReactorStructureWorldStore {
    private val definitions = linkedMapOf<ReactorStructureId, ReactorMultiblockDefinition>()
    private val membershipsByPosition = linkedMapOf<ReactorBlockPosition, ReactorBlockMembership>()
    private val controllerByPosition = linkedMapOf<ReactorBlockPosition, ReactorStructureId>()
    private val diagnosticsByController = linkedMapOf<ReactorBlockPosition, ReactorAssemblyDiagnostic>()
    private val unknownControllers = linkedSetOf<ReactorBlockPosition>()

    fun applyPlan(
        scannedPositions: Set<ReactorBlockPosition>,
        plan: ReactorAssemblyPlan,
        removeMissingStructures: Boolean,
    ): Set<ReactorBlockPosition> {
        val changed = linkedSetOf<ReactorBlockPosition>()
        val newIds = plan.definitions.mapTo(linkedSetOf()) { it.structureId }
        val intersectingIds = scannedPositions
            .mapNotNullTo(linkedSetOf()) { membershipsByPosition[it]?.structureId }

        if (removeMissingStructures) {
            for (structureId in intersectingIds - newIds) {
                changed += removeStructure(structureId)
            }
            for (position in scannedPositions) {
                if (position !in plan.memberships) {
                    changed += clearPosition(position)
                }
            }
        }

        for (definition in plan.definitions) {
            changed += removeStructure(definition.structureId)
            definitions[definition.structureId] = definition
            controllerByPosition[definition.controllerPosition] = definition.structureId
        }
        for ((position, membership) in plan.memberships) {
            if (membershipsByPosition[position] != membership) {
                membershipsByPosition[position] = membership
                changed += position
            }
        }

        for (position in scannedPositions) {
            diagnosticsByController.remove(position)
            unknownControllers.remove(position)
        }
        for (diagnostic in plan.diagnostics) {
            val controller = diagnostic.controllerPosition ?: continue
            diagnosticsByController[controller] = diagnostic
            changed += controller
        }
        if (plan.hasUnknownBoundary) {
            unknownControllers += scannedPositions.filter { position ->
                membershipsByPosition[position]?.role == ReactorBlockRole.CONTROLLER ||
                    diagnosticsByController.containsKey(position)
            }
            changed += unknownControllers
        }

        return changed
    }

    fun removePosition(position: ReactorBlockPosition): Set<ReactorBlockPosition> {
        val structureId = membershipsByPosition[position]?.structureId ?: return setOf(position)
        return removeStructure(structureId) + position
    }

    fun membershipAt(position: ReactorBlockPosition): ReactorBlockMembership? =
        membershipsByPosition[position]

    fun controllerViewState(position: ReactorBlockPosition): ReactorControllerViewState {
        if (position in unknownControllers) {
            val membership = membershipsByPosition[position]
            return ReactorControllerViewState(
                formationState = ReactorControllerFormationState.UNKNOWN,
                structureId = membership?.structureId,
                zoneCount = membership?.summary?.zoneCount ?: 0,
                chamberBlockCount = membership?.summary?.chamberBlockCount ?: 0,
                portCount = membership?.summary?.portCount ?: 0,
                diagnostic = "reactor structure touches an unloaded chunk; state is unknown",
            )
        }

        val membership = membershipsByPosition[position]
        if (membership?.role == ReactorBlockRole.CONTROLLER) {
            return ReactorControllerViewState(
                formationState = ReactorControllerFormationState.FORMED,
                structureId = membership.structureId,
                zoneCount = membership.summary.zoneCount,
                chamberBlockCount = membership.summary.chamberBlockCount,
                portCount = membership.summary.portCount,
                diagnostic = null,
            )
        }

        val diagnostic = diagnosticsByController[position]
        if (diagnostic != null) {
            return ReactorControllerViewState(
                formationState = ReactorControllerFormationState.NOT_FORMED,
                structureId = null,
                zoneCount = 0,
                chamberBlockCount = 0,
                portCount = 0,
                diagnostic = diagnostic.errors.joinToString(separator = "; "),
            )
        }

        return ReactorControllerViewState.NotFormed
    }

    private fun removeStructure(structureId: ReactorStructureId): Set<ReactorBlockPosition> {
        definitions.remove(structureId)
        controllerByPosition.entries.removeAll { it.value == structureId }
        val removed = membershipsByPosition
            .filterValues { it.structureId == structureId }
            .keys
            .toList()
        for (position in removed) {
            membershipsByPosition.remove(position)
            diagnosticsByController.remove(position)
            unknownControllers.remove(position)
        }
        return removed.toSet()
    }

    private fun clearPosition(position: ReactorBlockPosition): Set<ReactorBlockPosition> {
        val changed = linkedSetOf<ReactorBlockPosition>()
        if (membershipsByPosition.remove(position) != null) {
            changed += position
        }
        if (diagnosticsByController.remove(position) != null) {
            changed += position
        }
        if (unknownControllers.remove(position)) {
            changed += position
        }
        controllerByPosition.remove(position)
        return changed
    }
}
