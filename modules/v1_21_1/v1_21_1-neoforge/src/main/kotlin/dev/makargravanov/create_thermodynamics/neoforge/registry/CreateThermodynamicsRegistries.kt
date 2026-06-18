package dev.makargravanov.create_thermodynamics.neoforge.registry

import dev.makargravanov.create_thermodynamics.neoforge.CreateThermodynamicsMod
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorChamberBlock
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlock
import net.minecraft.core.registries.Registries
import net.minecraft.network.chat.Component
import net.minecraft.world.item.BlockItem
import net.minecraft.world.item.CreativeModeTab
import net.minecraft.world.item.Item
import net.minecraft.world.item.ItemStack
import net.minecraft.world.level.block.Block
import net.minecraft.world.level.block.SoundType
import net.minecraft.world.level.block.state.BlockBehaviour
import net.minecraft.world.level.material.MapColor
import net.neoforged.bus.api.IEventBus
import net.neoforged.neoforge.registries.DeferredHolder
import net.neoforged.neoforge.registries.DeferredRegister
import java.util.function.Supplier

object CreateThermodynamicsRegistries {
    private val blocks = DeferredRegister.create(Registries.BLOCK, CreateThermodynamicsMod.MOD_ID)
    private val items = DeferredRegister.create(Registries.ITEM, CreateThermodynamicsMod.MOD_ID)
    private val creativeModeTabs = DeferredRegister.create(Registries.CREATIVE_MODE_TAB, CreateThermodynamicsMod.MOD_ID)

    val reactorChamber: DeferredHolder<Block, ReactorChamberBlock> =
        blocks.register("reactor_chamber", Supplier { ReactorChamberBlock(reactorBlockProperties()) })
    val reactorController: DeferredHolder<Block, ReactorMultiblockBlock> = registerReactorBlock("reactor_controller")
    val reactorItemInputPort: DeferredHolder<Block, ReactorMultiblockBlock> = registerReactorBlock("reactor_item_input_port")
    val reactorItemOutputPort: DeferredHolder<Block, ReactorMultiblockBlock> = registerReactorBlock("reactor_item_output_port")
    val reactorFluidInputPort: DeferredHolder<Block, ReactorMultiblockBlock> = registerReactorBlock("reactor_fluid_input_port")
    val reactorFluidOutputPort: DeferredHolder<Block, ReactorMultiblockBlock> = registerReactorBlock("reactor_fluid_output_port")

    val reactorChamberItem: DeferredHolder<Item, BlockItem> = registerBlockItem("reactor_chamber", reactorChamber)
    val reactorControllerItem: DeferredHolder<Item, BlockItem> = registerBlockItem("reactor_controller", reactorController)
    val reactorItemInputPortItem: DeferredHolder<Item, BlockItem> = registerBlockItem("reactor_item_input_port", reactorItemInputPort)
    val reactorItemOutputPortItem: DeferredHolder<Item, BlockItem> = registerBlockItem("reactor_item_output_port", reactorItemOutputPort)
    val reactorFluidInputPortItem: DeferredHolder<Item, BlockItem> = registerBlockItem("reactor_fluid_input_port", reactorFluidInputPort)
    val reactorFluidOutputPortItem: DeferredHolder<Item, BlockItem> = registerBlockItem("reactor_fluid_output_port", reactorFluidOutputPort)

    val mainCreativeTab: DeferredHolder<CreativeModeTab, CreativeModeTab> =
        creativeModeTabs.register(
            "main",
            Supplier {
                CreativeModeTab.builder()
                    .title(Component.translatable("itemGroup.create_thermodynamics.main"))
                    .icon { ItemStack(reactorControllerItem.get()) }
                    .displayItems { _, output ->
                        output.accept(reactorChamberItem.get())
                        output.accept(reactorControllerItem.get())
                        output.accept(reactorItemInputPortItem.get())
                        output.accept(reactorItemOutputPortItem.get())
                        output.accept(reactorFluidInputPortItem.get())
                        output.accept(reactorFluidOutputPortItem.get())
                    }
                    .build()
            },
        )

    fun register(eventBus: IEventBus) {
        blocks.register(eventBus)
        items.register(eventBus)
        creativeModeTabs.register(eventBus)
    }

    private fun registerReactorBlock(id: String): DeferredHolder<Block, ReactorMultiblockBlock> =
        blocks.register(id, Supplier { ReactorMultiblockBlock(reactorBlockProperties()) })

    private fun registerBlockItem(id: String, block: DeferredHolder<Block, out Block>): DeferredHolder<Item, BlockItem> =
        items.register(id, Supplier { BlockItem(block.get(), Item.Properties()) })

    private fun reactorBlockProperties(): BlockBehaviour.Properties =
        BlockBehaviour.Properties.of()
            .mapColor(MapColor.METAL)
            .strength(3.0f, 6.0f)
            .sound(SoundType.METAL)
            .requiresCorrectToolForDrops()
}
