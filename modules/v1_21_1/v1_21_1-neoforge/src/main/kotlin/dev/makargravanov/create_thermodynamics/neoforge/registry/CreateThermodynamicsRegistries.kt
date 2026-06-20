package dev.makargravanov.create_thermodynamics.neoforge.registry

import dev.makargravanov.create_thermodynamics.neoforge.CreateThermodynamicsMod
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlock
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockBlockEntity
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorMultiblockKind
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorControllerMenu
import dev.makargravanov.create_thermodynamics.neoforge.block.ReactorPortMenu
import net.minecraft.core.registries.Registries
import net.minecraft.network.chat.Component
import net.minecraft.world.flag.FeatureFlags
import net.minecraft.world.item.BlockItem
import net.minecraft.world.item.CreativeModeTab
import net.minecraft.world.item.Item
import net.minecraft.world.item.ItemStack
import net.minecraft.world.inventory.MenuType
import net.minecraft.world.level.block.Block
import net.minecraft.world.level.block.entity.BlockEntityType
import net.minecraft.world.level.block.SoundType
import net.minecraft.world.level.block.state.BlockBehaviour
import net.minecraft.world.level.material.MapColor
import net.neoforged.neoforge.common.extensions.IMenuTypeExtension
import net.neoforged.bus.api.IEventBus
import net.neoforged.neoforge.registries.DeferredHolder
import net.neoforged.neoforge.registries.DeferredRegister
import java.util.function.Supplier

object CreateThermodynamicsRegistries {
    private val blocks = DeferredRegister.create(Registries.BLOCK, CreateThermodynamicsMod.MOD_ID)
    private val items = DeferredRegister.create(Registries.ITEM, CreateThermodynamicsMod.MOD_ID)
    private val creativeModeTabs = DeferredRegister.create(Registries.CREATIVE_MODE_TAB, CreateThermodynamicsMod.MOD_ID)
    private val blockEntityTypes = DeferredRegister.create(Registries.BLOCK_ENTITY_TYPE, CreateThermodynamicsMod.MOD_ID)
    private val menuTypes = DeferredRegister.create(Registries.MENU, CreateThermodynamicsMod.MOD_ID)

    val reactorChamber: DeferredHolder<Block, ReactorMultiblockBlock> =
        registerReactorBlock("reactor_chamber", ReactorMultiblockKind.CHAMBER)
    val reactorController: DeferredHolder<Block, ReactorMultiblockBlock> =
        registerReactorBlock("reactor_controller", ReactorMultiblockKind.CONTROLLER)
    val reactorItemInputPort: DeferredHolder<Block, ReactorMultiblockBlock> =
        registerReactorBlock("reactor_item_input_port", ReactorMultiblockKind.ITEM_INPUT_PORT)
    val reactorItemOutputPort: DeferredHolder<Block, ReactorMultiblockBlock> =
        registerReactorBlock("reactor_item_output_port", ReactorMultiblockKind.ITEM_OUTPUT_PORT)
    val reactorFluidInputPort: DeferredHolder<Block, ReactorMultiblockBlock> =
        registerReactorBlock("reactor_fluid_input_port", ReactorMultiblockKind.FLUID_INPUT_PORT)
    val reactorFluidOutputPort: DeferredHolder<Block, ReactorMultiblockBlock> =
        registerReactorBlock("reactor_fluid_output_port", ReactorMultiblockKind.FLUID_OUTPUT_PORT)

    @Suppress("NULLABILITY_MISMATCH_BASED_ON_JAVA_ANNOTATIONS")
    val reactorMultiblockBlockEntity: DeferredHolder<BlockEntityType<*>, BlockEntityType<ReactorMultiblockBlockEntity>> =
        blockEntityTypes.register(
            "reactor_multiblock",
            Supplier {
                BlockEntityType.Builder.of(
                    ::ReactorMultiblockBlockEntity,
                    reactorChamber.get(),
                    reactorController.get(),
                    reactorItemInputPort.get(),
                    reactorItemOutputPort.get(),
                    reactorFluidInputPort.get(),
                    reactorFluidOutputPort.get(),
                ).build(null)
            },
        )

    val reactorControllerMenu: DeferredHolder<MenuType<*>, MenuType<ReactorControllerMenu>> =
        menuTypes.register(
            "reactor_controller",
            Supplier {
                IMenuTypeExtension.create { containerId, playerInventory, extraData ->
                    ReactorControllerMenu(containerId, playerInventory, extraData)
                }
            },
        )

    val reactorPortMenu: DeferredHolder<MenuType<*>, MenuType<ReactorPortMenu>> =
        menuTypes.register(
            "reactor_port",
            Supplier {
                MenuType({ containerId, playerInventory -> ReactorPortMenu(containerId, playerInventory) }, FeatureFlags.DEFAULT_FLAGS)
            },
        )

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
        blockEntityTypes.register(eventBus)
        menuTypes.register(eventBus)
    }

    private fun registerReactorBlock(
        id: String,
        kind: ReactorMultiblockKind,
    ): DeferredHolder<Block, ReactorMultiblockBlock> =
        blocks.register(id, Supplier { ReactorMultiblockBlock(reactorBlockProperties(), kind) })

    private fun registerBlockItem(id: String, block: DeferredHolder<Block, out Block>): DeferredHolder<Item, BlockItem> =
        items.register(id, Supplier { BlockItem(block.get(), Item.Properties()) })

    private fun reactorBlockProperties(): BlockBehaviour.Properties =
        BlockBehaviour.Properties.of()
            .mapColor(MapColor.METAL)
            .strength(3.0f, 6.0f)
            .sound(SoundType.METAL)
            .requiresCorrectToolForDrops()
}
