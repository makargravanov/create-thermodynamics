package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorPortKind

enum class ReactorMultiblockKind(val modelKind: ReactorMultiblockBlockKind) {
    CHAMBER(ReactorMultiblockBlockKind.CHAMBER),
    CONTROLLER(ReactorMultiblockBlockKind.CONTROLLER),
    ITEM_INPUT_PORT(ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
    ITEM_OUTPUT_PORT(ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT),
    FLUID_INPUT_PORT(ReactorMultiblockBlockKind.FLUID_INPUT_PORT),
    FLUID_OUTPUT_PORT(ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT),
    ;

    val hasFacing: Boolean
        get() = when (this) {
            CONTROLLER,
            ITEM_INPUT_PORT,
            ITEM_OUTPUT_PORT,
            FLUID_INPUT_PORT,
            FLUID_OUTPUT_PORT,
            -> true

            CHAMBER,
            -> false
        }

    val isPort: Boolean
        get() = portKind != null

    val opensMenu: Boolean
        get() = this == CONTROLLER || isPort

    val portKind: ReactorPortKind?
        get() = when (this) {
            ITEM_INPUT_PORT -> ReactorPortKind.ITEM_INPUT
            ITEM_OUTPUT_PORT -> ReactorPortKind.ITEM_OUTPUT
            FLUID_INPUT_PORT -> ReactorPortKind.FLUID_INPUT
            FLUID_OUTPUT_PORT -> ReactorPortKind.FLUID_OUTPUT
            CONTROLLER,
            CHAMBER,
            -> null
        }

}
