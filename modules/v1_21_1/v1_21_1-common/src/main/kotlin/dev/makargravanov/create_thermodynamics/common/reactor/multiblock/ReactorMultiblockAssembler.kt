package dev.makargravanov.create_thermodynamics.common.reactor.multiblock

import java.util.UUID

class ReactorMultiblockAssembler(
    private val rules: ReactorMultiblockRules,
) {
    fun assemble(
        structureId: UUID,
        blocks: Iterable<ReactorMultiblockBlock>,
    ): ReactorMultiblockDefinition {
        val errors = mutableListOf<String>()
        val blocksByPosition = mutableMapOf<ReactorBlockPosition, ReactorMultiblockBlockKind>()

        for (block in blocks) {
            val previous = blocksByPosition.put(block.position, block.kind)
            if (previous != null) {
                errors += "duplicate reactor multiblock block at ${block.position}: $previous and ${block.kind}"
            }
        }

        val controllers = blocksByPosition.entries
            .filter { it.value == ReactorMultiblockBlockKind.CONTROLLER }
            .map { it.key }
            .sorted()
        val chambers = blocksByPosition.entries
            .filter { it.value == ReactorMultiblockBlockKind.CHAMBER }
            .map { it.key }
            .toSet()

        if (controllers.size != 1) {
            errors += "reactor multiblock must contain exactly one controller, got ${controllers.size}"
        }
        if (chambers.size < rules.minimumChamberBlocks) {
            errors += "reactor multiblock must contain at least ${rules.minimumChamberBlocks} chamber blocks, got ${chambers.size}"
        }
        if (rules.maximumChamberBlocks != null && chambers.size > rules.maximumChamberBlocks) {
            errors += "reactor multiblock may contain at most ${rules.maximumChamberBlocks} chamber blocks, got ${chambers.size}"
        }

        if (chambers.isNotEmpty()) {
            val connectedChambers = connectedChamberComponent(chambers)
            if (connectedChambers.size != chambers.size) {
                val disconnected = (chambers - connectedChambers).sorted()
                errors += "reactor chamber blocks must form one face-connected zone; disconnected chambers: $disconnected"
            }
        }

        val controller = controllers.firstOrNull()
        if (controller != null && controller.faceNeighbours().none { it in chambers }) {
            errors += "reactor controller at $controller must touch a chamber block by a face"
        }

        val portDescriptors = buildPortDescriptors(blocksByPosition, chambers, errors)

        if (errors.isNotEmpty()) {
            throw ReactorMultiblockValidationException(errors)
        }

        val chamberPositions = chambers.toSortedSet()
        val volume = chamberPositions.size * rules.chamberVolumeCubicMeters
        return ReactorMultiblockDefinition(
            structureId = structureId,
            controllerPosition = requireNotNull(controller),
            zone = ReactorZoneDescriptor(
                zoneIndex = 0,
                chamberPositions = chamberPositions,
                volumeCubicMeters = volume,
            ),
            ports = portDescriptors,
        )
    }

    private fun buildPortDescriptors(
        blocksByPosition: Map<ReactorBlockPosition, ReactorMultiblockBlockKind>,
        chambers: Set<ReactorBlockPosition>,
        errors: MutableList<String>,
    ): List<ReactorPortDescriptor> {
        val portEntries = blocksByPosition.entries
            .mapNotNull { (position, kind) -> kind.toPortKind()?.let { position to it } }
            .sortedWith(compareBy<Pair<ReactorBlockPosition, ReactorPortKind>> { it.second }.thenBy { it.first })

        val nextIndexByKind = mutableMapOf<ReactorPortKind, Int>()
        val descriptors = mutableListOf<ReactorPortDescriptor>()
        for ((position, portKind) in portEntries) {
            val attachedChambers = position.faceNeighbours().filter { it in chambers }.sorted()
            if (attachedChambers.isEmpty()) {
                errors += "reactor port $portKind at $position must touch a chamber block by a face"
                continue
            }
            val portIndex = nextIndexByKind.getOrDefault(portKind, 0)
            nextIndexByKind[portKind] = portIndex + 1
            descriptors += ReactorPortDescriptor(
                portIndex = portIndex,
                kind = portKind,
                position = position,
                zoneIndex = 0,
                attachedChamberPosition = attachedChambers.first(),
            )
        }
        return descriptors
    }

    private fun connectedChamberComponent(chambers: Set<ReactorBlockPosition>): Set<ReactorBlockPosition> {
        val start = chambers.minOrNull() ?: return emptySet()
        val visited = mutableSetOf<ReactorBlockPosition>()
        val queue = ArrayDeque<ReactorBlockPosition>()
        queue += start
        while (queue.isNotEmpty()) {
            val current = queue.removeFirst()
            if (!visited.add(current)) {
                continue
            }
            for (neighbour in current.faceNeighbours()) {
                if (neighbour in chambers && neighbour !in visited) {
                    queue += neighbour
                }
            }
        }
        return visited
    }

    private fun ReactorMultiblockBlockKind.toPortKind(): ReactorPortKind? =
        when (this) {
            ReactorMultiblockBlockKind.CONTROLLER,
            ReactorMultiblockBlockKind.CHAMBER,
            -> null

            ReactorMultiblockBlockKind.ITEM_INPUT_PORT -> ReactorPortKind.ITEM_INPUT
            ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT -> ReactorPortKind.ITEM_OUTPUT
            ReactorMultiblockBlockKind.FLUID_INPUT_PORT -> ReactorPortKind.FLUID_INPUT
            ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT -> ReactorPortKind.FLUID_OUTPUT
        }
}
