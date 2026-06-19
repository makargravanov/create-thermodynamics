package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorWorldRuntime
import dev.makargravanov.create_thermodynamics.common.rust.blob.NativeBlobStorage
import net.minecraft.server.level.ServerLevel
import net.minecraft.world.level.storage.LevelResource
import java.util.WeakHashMap

object ReactorRuntimeWorlds {
    private val runtimesByLevel = WeakHashMap<ServerLevel, ReactorWorldRuntime>()

    fun runtime(level: ServerLevel): ReactorWorldRuntime =
        runtimesByLevel.getOrPut(level) {
            ReactorWorldRuntime(
                blobStorage = NativeBlobStorage(
                    level.server.getWorldPath(LevelResource.ROOT)
                        .resolve("create_thermodynamics")
                        .resolve("native"),
                ),
            )
        }

    fun tickExisting(dtSeconds: Double) {
        for (runtime in runtimesByLevel.values) {
            runtime.tickAll(dtSeconds)
        }
    }
}
