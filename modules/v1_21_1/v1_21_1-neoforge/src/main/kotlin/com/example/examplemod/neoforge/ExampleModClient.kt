package com.example.examplemod.neoforge

import net.neoforged.api.distmarker.Dist
import net.neoforged.bus.api.SubscribeEvent
import net.neoforged.fml.common.EventBusSubscriber
import net.neoforged.fml.event.lifecycle.FMLClientSetupEvent

@EventBusSubscriber(modid = ExampleMod.MOD_ID, value = [Dist.CLIENT])
object ExampleModClient {
    @SubscribeEvent
    fun onClientSetup(event: FMLClientSetupEvent) {
    }
}
