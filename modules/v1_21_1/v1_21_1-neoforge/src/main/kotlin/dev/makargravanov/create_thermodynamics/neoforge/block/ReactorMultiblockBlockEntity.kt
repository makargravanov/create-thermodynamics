package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.neoforge.registry.CreateThermodynamicsRegistries
import net.minecraft.core.BlockPos
import net.minecraft.core.HolderLookup
import net.minecraft.core.NonNullList
import net.minecraft.nbt.CompoundTag
import net.minecraft.network.chat.Component
import net.minecraft.network.protocol.Packet
import net.minecraft.network.protocol.game.ClientGamePacketListener
import net.minecraft.network.protocol.game.ClientboundBlockEntityDataPacket
import net.minecraft.world.Container
import net.minecraft.world.ContainerHelper
import net.minecraft.world.entity.player.Inventory
import net.minecraft.world.entity.player.Player
import net.minecraft.world.inventory.AbstractContainerMenu
import net.minecraft.world.inventory.ChestMenu
import net.minecraft.world.inventory.MenuType
import net.minecraft.world.item.ItemStack
import net.minecraft.world.level.block.entity.BlockEntity
import net.minecraft.world.level.block.state.BlockState
import net.minecraft.world.level.block.Block
import net.minecraft.world.MenuProvider
import java.util.UUID

class ReactorMultiblockBlockEntity(pos: BlockPos, state: BlockState) :
    BlockEntity(CreateThermodynamicsRegistries.reactorMultiblockBlockEntity.get(), pos, state),
    Container,
    MenuProvider {
    private val items: NonNullList<ItemStack> = NonNullList.withSize(CONTAINER_SIZE, ItemStack.EMPTY)

    var structureId: UUID? = null
        private set
    var activeVolumeBlock: Boolean = false
        private set

    fun visualGroupKey(): UUID? =
        structureId?.takeIf { activeVolumeBlock }

    fun setStructureMembership(newStructureId: UUID?, newActiveVolumeBlock: Boolean): Boolean {
        val normalizedActive = newStructureId != null && newActiveVolumeBlock
        if (structureId == newStructureId && activeVolumeBlock == normalizedActive) {
            return false
        }
        structureId = newStructureId
        activeVolumeBlock = normalizedActive
        setChanged()
        refreshVisualModel()
        return true
    }

    override fun loadAdditional(tag: CompoundTag, registries: HolderLookup.Provider) {
        val oldStructureId = structureId
        val oldActiveVolumeBlock = activeVolumeBlock
        super.loadAdditional(tag, registries)
        structureId = if (tag.hasUUID(STRUCTURE_ID_TAG)) tag.getUUID(STRUCTURE_ID_TAG) else null
        activeVolumeBlock = structureId != null && tag.getBoolean(ACTIVE_VOLUME_TAG)
        ContainerHelper.loadAllItems(tag, items, registries)
        if (structureId != oldStructureId || activeVolumeBlock != oldActiveVolumeBlock) {
            refreshVisualModel()
        }
    }

    override fun saveAdditional(tag: CompoundTag, registries: HolderLookup.Provider) {
        super.saveAdditional(tag, registries)
        structureId?.let { tag.putUUID(STRUCTURE_ID_TAG, it) }
        tag.putBoolean(ACTIVE_VOLUME_TAG, activeVolumeBlock)
        ContainerHelper.saveAllItems(tag, items, registries)
    }

    override fun getUpdatePacket(): Packet<ClientGamePacketListener> =
        ClientboundBlockEntityDataPacket.create(this)

    override fun getUpdateTag(registries: HolderLookup.Provider): CompoundTag =
        saveWithoutMetadata(registries)

    private fun refreshVisualModel() {
        requestModelDataUpdate()
        level?.sendBlockUpdated(blockPos, blockState, blockState, Block.UPDATE_CLIENTS)
    }

    override fun getContainerSize(): Int =
        items.size

    override fun isEmpty(): Boolean =
        items.all(ItemStack::isEmpty)

    override fun getItem(slot: Int): ItemStack =
        items[slot]

    override fun removeItem(slot: Int, amount: Int): ItemStack =
        ContainerHelper.removeItem(items, slot, amount).also { removed ->
            if (!removed.isEmpty) {
                setChanged()
            }
        }

    override fun removeItemNoUpdate(slot: Int): ItemStack =
        ContainerHelper.takeItem(items, slot)

    override fun setItem(slot: Int, stack: ItemStack) {
        items[slot] = stack
        if (stack.count > maxStackSize) {
            stack.count = maxStackSize
        }
        setChanged()
    }

    override fun stillValid(player: Player): Boolean =
        Container.stillValidBlockEntity(this, player)

    override fun clearContent() {
        items.clear()
        setChanged()
    }

    override fun getDisplayName(): Component =
        Component.translatable("container.create_thermodynamics.reactor_port")

    override fun createMenu(containerId: Int, playerInventory: Inventory, player: Player): AbstractContainerMenu =
        ChestMenu(MenuType.GENERIC_9x1, containerId, playerInventory, this, 1)

    companion object {
        private const val CONTAINER_SIZE = 9
        private const val STRUCTURE_ID_TAG = "structure_id"
        private const val ACTIVE_VOLUME_TAG = "active_volume"
    }
}
