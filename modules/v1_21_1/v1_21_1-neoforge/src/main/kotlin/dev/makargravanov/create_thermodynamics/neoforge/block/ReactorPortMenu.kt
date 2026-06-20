package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.neoforge.registry.CreateThermodynamicsRegistries
import net.minecraft.world.Container
import net.minecraft.world.SimpleContainer
import net.minecraft.world.entity.player.Inventory
import net.minecraft.world.entity.player.Player
import net.minecraft.world.inventory.AbstractContainerMenu
import net.minecraft.world.inventory.Slot
import net.minecraft.world.item.ItemStack

class ReactorPortMenu(
    containerId: Int,
    private val playerInventory: Inventory,
    private val portInventory: Container = SimpleContainer(PortSlotCount),
    private val blockEntity: ReactorMultiblockBlockEntity? = null,
) : AbstractContainerMenu(CreateThermodynamicsRegistries.reactorPortMenu.get(), containerId) {
    init {
        check(portInventory.containerSize == PortSlotCount) {
            "reactor port menu requires $PortSlotCount buffer slots, got ${portInventory.containerSize}"
        }
        addPortSlots()
        addPlayerSlots()
    }

    private fun addPortSlots() {
        for (row in 0 until PortRows) {
            for (column in 0 until SlotsPerRow) {
                val slot = column + row * SlotsPerRow
                addSlot(Slot(portInventory, slot, PortSlotX + column * SlotStep, PortSlotY + row * SlotStep))
            }
        }
    }

    private fun addPlayerSlots() {
        for (row in 0 until PlayerInventoryRows) {
            for (column in 0 until SlotsPerRow) {
                val slot = column + row * SlotsPerRow + SlotsPerRow
                addSlot(Slot(playerInventory, slot, PlayerSlotX + column * SlotStep, PlayerSlotY + row * SlotStep))
            }
        }
        for (column in 0 until SlotsPerRow) {
            addSlot(Slot(playerInventory, column, PlayerSlotX + column * SlotStep, PlayerHotbarY))
        }
    }

    override fun stillValid(player: Player): Boolean =
        blockEntity?.stillValid(player) ?: true

    override fun quickMoveStack(player: Player, index: Int): ItemStack {
        val slot = slots.getOrNull(index) ?: return ItemStack.EMPTY
        if (!slot.hasItem()) {
            return ItemStack.EMPTY
        }
        val original = slot.item
        val copy = original.copy()
        val moved = if (index < PortSlotCount) {
            moveItemStackTo(original, PortSlotCount, slots.size, true)
        } else {
            moveItemStackTo(original, 0, PortSlotCount, false)
        }
        if (!moved) {
            return ItemStack.EMPTY
        }
        if (original.isEmpty) {
            slot.setByPlayer(ItemStack.EMPTY)
        } else {
            slot.setChanged()
        }
        return copy
    }

    companion object {
        const val PortSlotCount = 27
        private const val PortRows = 3
        private const val SlotsPerRow = 9
        private const val SlotStep = 18
        private const val PortSlotX = 8
        private const val PortSlotY = 18
        private const val PlayerSlotX = 8
        private const val PlayerSlotY = 85
        private const val PlayerHotbarY = 143
        private const val PlayerInventoryRows = 3
    }
}
