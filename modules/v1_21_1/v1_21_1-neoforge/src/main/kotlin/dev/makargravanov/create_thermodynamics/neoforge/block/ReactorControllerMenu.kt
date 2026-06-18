package dev.makargravanov.create_thermodynamics.neoforge.block

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
    val state: ReactorControllerScreenState
        get() {
            val localState = blockEntity?.controllerScreenState()
            return ReactorControllerScreenState(
                structureId = localState?.structureId,
                formed = data.get(FormedSlot) != 0,
                zoneCount = data.get(ZoneCountSlot),
                chamberBlockCount = data.get(ChamberBlockCountSlot),
                portCount = data.get(PortCountSlot),
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
        private const val FormedSlot = 0
        private const val ZoneCountSlot = 1
        private const val ChamberBlockCountSlot = 2
        private const val PortCountSlot = 3
        private const val DataSlotCount = 4
    }
}

data class ReactorControllerScreenState(
    val structureId: String?,
    val formed: Boolean,
    val zoneCount: Int,
    val chamberBlockCount: Int,
    val portCount: Int,
) {
    companion object {
        val Empty = ReactorControllerScreenState(
            structureId = null,
            formed = false,
            zoneCount = 0,
            chamberBlockCount = 0,
            portCount = 0,
        )
    }
}
