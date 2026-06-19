package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorAssemblyResult
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockAssembler
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly.ReactorMultiblockRules
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlock
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import java.util.UUID

class ReactorAssemblyPlanner(
    rules: ReactorMultiblockRules,
) {
    private val assembler = ReactorMultiblockAssembler(rules)

    fun buildPlan(scan: ReactorAssemblyScan): ReactorAssemblyPlan {
        val blocksByPosition = scan.blocks
            .filter { it.loaded && it.kind != null }
            .associateBy { it.position }
        val componentPositions = blocksByPosition.keys.toSortedSet()
        val used = linkedSetOf<ReactorBlockPosition>()
        val definitions = mutableListOf<ReactorMultiblockDefinition>()
        val diagnostics = mutableListOf<ReactorAssemblyDiagnostic>()

        val controllers = blocksByPosition.values
            .filter { it.kind == ReactorMultiblockBlockKind.CONTROLLER }
            .map { it.position }
            .sorted()

        for (controller in controllers) {
            if (controller in used) {
                continue
            }
            val candidatePositions = collectCandidateForController(
                blocksByPosition = blocksByPosition,
                componentPositions = componentPositions,
                controller = controller,
                alreadyUsed = used,
            )
            val blocks = candidatePositions.map { position ->
                val snapshot = requireNotNull(blocksByPosition[position]) {
                    "missing reactor block snapshot for candidate position $position"
                }
                ReactorMultiblockBlock(
                    position = position,
                    kind = requireNotNull(snapshot.kind),
                    facing = snapshot.facing,
                )
            }
            val structureId = structureIdForCandidate(controller, blocks)
            when (val result = assembler.tryAssemble(structureId, blocks)) {
                is ReactorAssemblyResult.Formed -> {
                    definitions += result.definition
                    used += result.definition.structurePositions()
                }

                is ReactorAssemblyResult.Rejected -> {
                    diagnostics += ReactorAssemblyDiagnostic(
                        controllerPosition = controller,
                        errors = result.errors,
                    )
                }
            }
        }

        return ReactorAssemblyPlan(
            definitions = definitions,
            memberships = definitions
                .flatMap { definition -> definition.memberships().entries }
                .associate { it.toPair() },
            diagnostics = diagnostics,
            hasUnknownBoundary = scan.hasUnknownBoundary,
        )
    }

    private fun collectCandidateForController(
        blocksByPosition: Map<ReactorBlockPosition, ReactorWorldBlockSnapshot>,
        componentPositions: Set<ReactorBlockPosition>,
        controller: ReactorBlockPosition,
        alreadyUsed: Set<ReactorBlockPosition>,
    ): Set<ReactorBlockPosition> {
        val visited = linkedSetOf<ReactorBlockPosition>()
        val queue = ArrayDeque<ReactorBlockPosition>()
        queue += controller

        while (queue.isNotEmpty()) {
            val position = queue.removeFirst()
            if (!visited.add(position)) {
                continue
            }
            for (next in position.faceNeighbours()) {
                val nextKind = blocksByPosition[next]?.kind
                if (
                    next in componentPositions &&
                    next !in alreadyUsed &&
                    (next == controller || nextKind != ReactorMultiblockBlockKind.CONTROLLER)
                ) {
                    queue += next
                }
            }
        }
        return visited
    }

    private fun ReactorMultiblockDefinition.memberships(): Map<ReactorBlockPosition, ReactorBlockMembership> {
        val summary = ReactorStructureSummary(
            structureId = structureId,
            zoneCount = 1,
            chamberBlockCount = zone.volumePositions.size,
            portCount = ports.size,
        )
        val result = linkedMapOf<ReactorBlockPosition, ReactorBlockMembership>()
        result[controllerPosition] = ReactorBlockMembership(
            structureId = structureId,
            role = ReactorBlockRole.CONTROLLER,
            activeVolumeBlock = controllerPosition in zone.volumePositions,
            summary = summary,
        )
        for (position in zone.plainChamberPositions) {
            result[position] = ReactorBlockMembership(
                structureId = structureId,
                role = ReactorBlockRole.CHAMBER,
                activeVolumeBlock = position in zone.volumePositions,
                summary = summary,
            )
        }
        for (port in ports) {
            result[port.position] = ReactorBlockMembership(
                structureId = structureId,
                role = ReactorBlockRole.PORT,
                activeVolumeBlock = port.position in zone.volumePositions,
                summary = summary,
            )
        }
        return result
    }

    private fun ReactorMultiblockDefinition.structurePositions(): Set<ReactorBlockPosition> =
        linkedSetOf<ReactorBlockPosition>().also { positions ->
            positions += zone.volumePositions
            positions += controllerPosition
            positions += ports.map { it.position }
        }

    private fun structureIdForCandidate(
        controller: ReactorBlockPosition,
        blocks: List<ReactorMultiblockBlock>,
    ): ReactorStructureId {
        val structureFingerprint = buildString {
            append("reactor-structure:")
            append(controller.x)
            append(',')
            append(controller.y)
            append(',')
            append(controller.z)
            for (block in blocks.sortedBy { it.position }) {
                append(';')
                append(block.position.x)
                append(',')
                append(block.position.y)
                append(',')
                append(block.position.z)
                append('=')
                append(block.kind.name)
                append('@')
                append(block.facing?.name ?: "none")
            }
        }
        return ReactorStructureId(
            UUID.nameUUIDFromBytes(
                structureFingerprint.toByteArray(Charsets.UTF_8),
            ),
        )
    }
}
