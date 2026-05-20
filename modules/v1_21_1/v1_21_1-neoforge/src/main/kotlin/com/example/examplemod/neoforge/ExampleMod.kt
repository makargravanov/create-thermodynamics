package com.example.examplemod.neoforge

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
        LOGGER.info("Loaded {}", MOD_ID)
    }

    companion object {
        const val MOD_ID = "examplemod"
        private val LOGGER = LoggerFactory.getLogger(ExampleMod::class.java)
    }
}
