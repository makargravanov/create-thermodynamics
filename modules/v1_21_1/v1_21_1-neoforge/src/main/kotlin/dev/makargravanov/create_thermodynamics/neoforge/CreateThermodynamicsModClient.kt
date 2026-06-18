package dev.makargravanov.create_thermodynamics.neoforge

import net.neoforged.api.distmarker.Dist
import net.neoforged.bus.api.SubscribeEvent
import net.neoforged.fml.common.EventBusSubscriber
import net.neoforged.fml.event.lifecycle.FMLClientSetupEvent

@EventBusSubscriber(modid = CreateThermodynamicsMod.MOD_ID, value = [Dist.CLIENT])
object CreateThermodynamicsModClient {
    @SubscribeEvent
    @JvmStatic
    fun onClientSetup(event: FMLClientSetupEvent) {
    }
}
