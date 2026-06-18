package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockDirection
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlock
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockValidationException
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import java.util.UUID

class ReactorMultiblockAssembler(
    private val rules: ReactorMultiblockRules,
) {
    fun assemble(
        structureId: UUID,
        blocks: Iterable<ReactorMultiblockBlock>,
    ): ReactorMultiblockDefinition =
        assemble(ReactorStructureId(structureId), blocks)

    fun assemble(
        structureId: ReactorStructureId,
        blocks: Iterable<ReactorMultiblockBlock>,
    ): ReactorMultiblockDefinition =
        when (val result = tryAssemble(structureId, blocks)) {
            is ReactorAssemblyResult.Formed -> result.definition
            is ReactorAssemblyResult.Rejected -> throw ReactorMultiblockValidationException(result.errors)
        }

    fun tryAssemble(
        structureId: ReactorStructureId,
        blocks: Iterable<ReactorMultiblockBlock>,
    ): ReactorAssemblyResult {
        val errors = mutableListOf<String>()
        val blocksByPosition = mutableMapOf<ReactorBlockPosition, ReactorMultiblockBlockKind>()
        val facingByPosition = mutableMapOf<ReactorBlockPosition, ReactorBlockDirection?>()

        for (block in blocks) {
            val previous = blocksByPosition.put(block.position, block.kind)
            if (previous != null) {
                errors += "duplicate reactor multiblock block at ${block.position}: $previous and ${block.kind}"
            }
            facingByPosition[block.position] = block.facing
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
        val controller = controllers.firstOrNull()
        val shapeResult = rules.chamberShapeStrategy.buildZone(
            chambers = chambers,
            volumeCapableBlocks = blocksByPosition.filterValues { it.isVolumeCapable() },
            controller = controller,
            chamberVolumeCubicMeters = rules.chamberVolumeCubicMeters,
            maximumVolumeBlocks = rules.maximumVolumeBlocks,
        )
        errors += shapeResult.errors
        val zone = shapeResult.zone

        if (zone != null) {
            val selectedVolume = zone.volumePositions
            if (selectedVolume.size < rules.minimumVolumeBlocks) {
                errors += "reactor multiblock must contain at least ${rules.minimumVolumeBlocks} volume blocks, got ${selectedVolume.size}"
            }
            if (controller != null && controller !in zone.volumePositions && contactDirections(controller, zone.plainChamberPositions).isEmpty()) {
                errors += "reactor controller at $controller must touch a chamber block by a face"
            }
        }

        val portDescriptors = buildPortDescriptors(
            blocksByPosition = blocksByPosition,
            facingByPosition = facingByPosition,
            zoneVolume = zone?.volumePositions.orEmpty(),
            plainChambers = zone?.plainChamberPositions.orEmpty(),
            errors = errors,
        )

        if (errors.isNotEmpty()) {
            return ReactorAssemblyResult.Rejected(errors)
        }

        return ReactorAssemblyResult.Formed(
            ReactorMultiblockDefinition(
                structureId = structureId,
                controllerPosition = requireNotNull(controller),
                controllerContactDirection = contactDirections(controller, requireNotNull(zone).plainChamberPositions).singleOrNull(),
                zone = requireNotNull(zone),
                ports = portDescriptors,
                inactiveChamberPositions = shapeResult.inactiveChamberPositions,
            ),
        )
    }

    private fun buildPortDescriptors(
        blocksByPosition: Map<ReactorBlockPosition, ReactorMultiblockBlockKind>,
        facingByPosition: Map<ReactorBlockPosition, ReactorBlockDirection?>,
        zoneVolume: Set<ReactorBlockPosition>,
        plainChambers: Set<ReactorBlockPosition>,
        errors: MutableList<String>,
    ): List<ReactorPortDescriptor> {
        val portEntries = blocksByPosition.entries
            .mapNotNull { (position, kind) -> kind.toPortKind()?.let { position to it } }
            .sortedWith(compareBy<Pair<ReactorBlockPosition, ReactorPortKind>> { it.second }.thenBy { it.first })

        val nextIndexByKind = mutableMapOf<ReactorPortKind, Int>()
        val descriptors = mutableListOf<ReactorPortDescriptor>()
        for ((position, portKind) in portEntries) {
            val facing = facingByPosition[position]
            val contact = portContact(position, portKind, facing, zoneVolume, plainChambers, errors) ?: continue
            val portIndex = nextIndexByKind.getOrDefault(portKind, 0)
            nextIndexByKind[portKind] = portIndex + 1
            descriptors += ReactorPortDescriptor(
                portIndex = portIndex,
                kind = portKind,
                position = position,
                zoneIndex = 0,
                attachedChamberPosition = contact.attachedPosition,
                contactDirection = contact.contactDirection,
            )
        }
        return descriptors
    }

    private fun portContact(
        position: ReactorBlockPosition,
        portKind: ReactorPortKind,
        facing: ReactorBlockDirection?,
        zoneVolume: Set<ReactorBlockPosition>,
        plainChambers: Set<ReactorBlockPosition>,
        errors: MutableList<String>,
    ): PortContact? {
        if (position in zoneVolume) {
            if (facing == null) {
                errors += "embedded reactor port $portKind at $position must have explicit facing"
                return null
            }
            return PortContact(attachedPosition = position, contactDirection = facing)
        }

        val contactDirections = contactDirections(position, plainChambers)
        if (contactDirections.isEmpty()) {
            errors += "external reactor port $portKind at $position must touch a chamber block by a face"
            return null
        }
        if (contactDirections.size > 1) {
            errors += "external reactor port $portKind at $position must touch exactly one chamber face, got ${contactDirections.size}"
            return null
        }
        val contactDirection = contactDirections.single()
        return PortContact(attachedPosition = position.neighbour(contactDirection), contactDirection = contactDirection)
    }

    private fun contactDirections(
        position: ReactorBlockPosition,
        chambers: Set<ReactorBlockPosition>,
    ) = position.faceNeighbours()
        .mapNotNull { neighbour -> position.directionTo(neighbour)?.takeIf { neighbour in chambers } }
        .sorted()

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

    private fun ReactorMultiblockBlockKind.isVolumeCapable(): Boolean =
        when (this) {
            ReactorMultiblockBlockKind.CONTROLLER,
            ReactorMultiblockBlockKind.CHAMBER,
            ReactorMultiblockBlockKind.ITEM_INPUT_PORT,
            ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT,
            ReactorMultiblockBlockKind.FLUID_INPUT_PORT,
            ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT,
            -> true
        }

    private data class PortContact(
        val attachedPosition: ReactorBlockPosition,
        val contactDirection: ReactorBlockDirection,
    )
}
