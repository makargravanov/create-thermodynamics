package com.example.examplemod.neoforge

import com.example.examplemod.common.CommonModule
import net.neoforged.bus.api.IEventBus
import net.neoforged.fml.common.Mod
import net.neoforged.fml.event.lifecycle.FMLCommonSetupEvent
import org.slf4j.LoggerFactory

@Mod(ExampleMod.MOD_ID)
class ExampleMod(modEventBus: IEventBus) {
    init {
        modEventBus.addListener(::commonSetup)
    }

    private fun commonSetup(event: FMLCommonSetupEvent) {
        LOGGER.info(
            "Loaded {} with Rust JNI ABI {} and demo pressure {} kPa",
            MOD_ID,
            CommonModule.nativeAbiVersion(),
            CommonModule.demoPressureKilopascals(),
        )
    }

    companion object {
        const val MOD_ID = "examplemod"
        private val LOGGER = LoggerFactory.getLogger(ExampleMod::class.java)
    }
}
