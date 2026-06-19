package dev.makargravanov.create_thermodynamics.neoforge

import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorRuntimeWorlds
import net.neoforged.neoforge.event.tick.ServerTickEvent

object ReactorRuntimeEvents {
    fun tickReactors(event: ServerTickEvent.Post) {
        for (level in event.server.allLevels) {
            ReactorRuntimeWorlds.tickLevel(level)
        }
    }
}
