package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorBlockPosition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortDescriptor
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId

enum class ReactorStructureState {
    ACTIVE,
    INVALID,
    REMOVED,
}

data class ReactorStructureRecord(
    val structureId: ReactorStructureId,
    val definition: ReactorMultiblockDefinition,
    val nativeBinding: NativeReactorMultiblockBinding,
    val state: ReactorStructureState,
) {
    private val portsByPosition: Map<ReactorBlockPosition, ReactorPortDescriptor> =
        definition.ports.associateBy { it.position }

    fun portAt(position: ReactorBlockPosition): ReactorPortDescriptor? =
        portsByPosition[position]

    fun portOfKind(position: ReactorBlockPosition, kind: ReactorPortKind): ReactorPortDescriptor? =
        portAt(position)?.takeIf { it.kind == kind }
}
