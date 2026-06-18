package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorMultiblockDefinition

sealed interface ReactorAssemblyResult {
    data class Formed(
        val definition: ReactorMultiblockDefinition,
    ) : ReactorAssemblyResult

    data class Rejected(
        val errors: List<String>,
    ) : ReactorAssemblyResult {
        init {
            require(errors.isNotEmpty()) { "rejected reactor assembly must contain diagnostics" }
        }
    }
}
