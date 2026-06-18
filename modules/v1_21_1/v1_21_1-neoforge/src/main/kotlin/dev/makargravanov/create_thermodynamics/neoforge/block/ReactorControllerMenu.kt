package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerFormationState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerViewState
import dev.makargravanov.create_thermodynamics.neoforge.registry.CreateThermodynamicsRegistries
import net.minecraft.world.entity.player.Player
import net.minecraft.world.inventory.AbstractContainerMenu
import net.minecraft.world.inventory.ContainerData
import net.minecraft.world.inventory.SimpleContainerData
import net.minecraft.world.item.ItemStack

class ReactorControllerMenu(
    containerId: Int,
    private val blockEntity: ReactorMultiblockBlockEntity? = null,
    private val data: ContainerData = blockEntity?.controllerMenuData() ?: SimpleContainerData(DataSlotCount),
) : AbstractContainerMenu(CreateThermodynamicsRegistries.reactorControllerMenu.get(), containerId) {
    val state: ReactorControllerViewState
        get() {
            val localState = blockEntity?.controllerScreenState()
            return ReactorControllerViewState(
                formationState = ReactorControllerFormationState.entries.getOrElse(data.get(FormationStateSlot)) {
                    ReactorControllerFormationState.NOT_FORMED
                },
                structureId = localState?.structureId,
                zoneCount = data.get(ZoneCountSlot),
                chamberBlockCount = data.get(ChamberBlockCountSlot),
                portCount = data.get(PortCountSlot),
                diagnostic = localState?.diagnostic,
            )
        }

    init {
        check(data.count == DataSlotCount) {
            "reactor controller menu data must expose $DataSlotCount slots, got ${data.count}"
        }
        addDataSlots(data)
    }

    override fun stillValid(player: Player): Boolean =
        blockEntity?.stillValid(player) ?: true

    override fun quickMoveStack(player: Player, index: Int): ItemStack =
        ItemStack.EMPTY

    companion object {
        private const val FormationStateSlot = 0
        private const val ZoneCountSlot = 1
        private const val ChamberBlockCountSlot = 2
        private const val PortCountSlot = 3
        private const val DataSlotCount = 4
    }
}
