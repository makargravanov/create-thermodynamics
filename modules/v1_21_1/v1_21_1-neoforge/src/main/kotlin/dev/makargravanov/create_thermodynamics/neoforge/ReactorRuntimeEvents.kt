package dev.makargravanov.create_thermodynamics.neoforge

import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorRuntimeWorlds
import net.neoforged.neoforge.event.tick.ServerTickEvent

object ReactorRuntimeEvents {
    private const val SECONDS_PER_SERVER_TICK = 1.0 / 20.0

    fun tickReactors(event: ServerTickEvent.Post) {
        ReactorRuntimeWorlds.tickExisting(SECONDS_PER_SERVER_TICK)
    }
}
