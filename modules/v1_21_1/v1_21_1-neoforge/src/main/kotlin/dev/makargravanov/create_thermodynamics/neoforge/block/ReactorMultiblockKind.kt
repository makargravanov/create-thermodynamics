package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockBlockKind

enum class ReactorMultiblockKind(val modelKind: ReactorMultiblockBlockKind) {
    CHAMBER(ReactorMultiblockBlockKind.CHAMBER),
    CONTROLLER(ReactorMultiblockBlockKind.CONTROLLER),
    ITEM_INPUT_PORT(ReactorMultiblockBlockKind.ITEM_INPUT_PORT),
    ITEM_OUTPUT_PORT(ReactorMultiblockBlockKind.ITEM_OUTPUT_PORT),
    FLUID_INPUT_PORT(ReactorMultiblockBlockKind.FLUID_INPUT_PORT),
    FLUID_OUTPUT_PORT(ReactorMultiblockBlockKind.FLUID_OUTPUT_PORT),
}
