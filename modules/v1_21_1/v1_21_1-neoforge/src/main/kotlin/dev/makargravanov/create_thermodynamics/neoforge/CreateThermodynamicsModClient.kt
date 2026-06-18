package dev.makargravanov.create_thermodynamics.neoforge

import dev.makargravanov.create_thermodynamics.neoforge.client.ReactorChamberConnectedTexture
import net.neoforged.api.distmarker.Dist
import net.neoforged.bus.api.SubscribeEvent
import net.neoforged.fml.common.EventBusSubscriber.Bus
import net.neoforged.fml.common.EventBusSubscriber
import net.neoforged.fml.event.lifecycle.FMLClientSetupEvent
import net.neoforged.neoforge.client.event.ModelEvent

@Suppress("DEPRECATION")
@EventBusSubscriber(modid = CreateThermodynamicsMod.MOD_ID, value = [Dist.CLIENT], bus = Bus.MOD)
object CreateThermodynamicsModClient {
    @SubscribeEvent
    @JvmStatic
    fun onClientSetup(event: FMLClientSetupEvent) {
    }

    @SubscribeEvent
    @JvmStatic
    fun onModifyBakingResult(event: ModelEvent.ModifyBakingResult) {
        ReactorChamberConnectedTexture.onModifyBakingResult(event)
    }
}
