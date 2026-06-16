package dev.makargravanov.create_thermodynamics.neoforge

import dev.makargravanov.create_thermodynamics.common.CommonModule
import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative
import net.neoforged.bus.api.IEventBus
import net.neoforged.fml.common.Mod
import net.neoforged.fml.event.lifecycle.FMLCommonSetupEvent
import org.slf4j.LoggerFactory

@Mod(CreateThermodynamicsMod.MOD_ID)
class CreateThermodynamicsMod(modEventBus: IEventBus) {
    init {
        modEventBus.addListener(::commonSetup)
    }

    private fun commonSetup(event: FMLCommonSetupEvent) {
        CommonModule.initializeNativeChemistry()
        LOGGER.info(
            "Loaded {} with Rust JNI ABI {}, {} default item bindings and {} static substances; demo pressure {} kPa",
            MOD_ID,
            CommonModule.nativeAbiVersion(),
            ThermodynamicsNative.itemChemicalBindingCount(),
            CommonModule.nativeStaticSubstanceIds().size,
            CommonModule.demoPressureKilopascals(),
        )
    }

    companion object {
        const val MOD_ID = "create_thermodynamics"
        private val LOGGER = LoggerFactory.getLogger(CreateThermodynamicsMod::class.java)
    }
}
