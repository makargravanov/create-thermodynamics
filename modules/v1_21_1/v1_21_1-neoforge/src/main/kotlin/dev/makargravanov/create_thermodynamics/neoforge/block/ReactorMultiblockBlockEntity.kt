package dev.makargravanov.create_thermodynamics.neoforge.block

import com.simibubi.create.foundation.blockEntity.SmartBlockEntity
import com.simibubi.create.foundation.blockEntity.behaviour.BlockEntityBehaviour
import com.simibubi.create.foundation.blockEntity.behaviour.CenteredSideValueBoxTransform
import com.simibubi.create.foundation.blockEntity.behaviour.filtering.FilteringBehaviour
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime.ReactorTickMetrics
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model.ReactorStructureId
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorBlockMembership
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerFormationState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorControllerViewState
import dev.makargravanov.create_thermodynamics.common.reactor.multiblock.world.ReactorMixtureViewEntry
import dev.makargravanov.create_thermodynamics.neoforge.registry.CreateThermodynamicsRegistries
import net.minecraft.core.BlockPos
import net.minecraft.core.HolderLookup
import net.minecraft.core.NonNullList
import net.minecraft.core.registries.BuiltInRegistries
import net.minecraft.nbt.CompoundTag
import net.minecraft.nbt.ListTag
import net.minecraft.nbt.Tag
import net.minecraft.network.chat.Component
import net.minecraft.network.protocol.game.ClientboundBlockEntityDataPacket
import net.minecraft.resources.ResourceLocation
import net.minecraft.world.Container
import net.minecraft.world.ContainerHelper
import net.minecraft.world.entity.player.Inventory
import net.minecraft.world.entity.player.Player
import net.minecraft.world.inventory.AbstractContainerMenu
import net.minecraft.world.inventory.ContainerData
import net.minecraft.world.item.ItemStack
import net.minecraft.world.item.Items
import net.minecraft.world.level.block.state.BlockState
import net.minecraft.world.level.block.Block
import net.minecraft.world.MenuProvider
import java.util.UUID

class ReactorMultiblockBlockEntity(pos: BlockPos, state: BlockState) :
    SmartBlockEntity(CreateThermodynamicsRegistries.reactorMultiblockBlockEntity.get(), pos, state),
    Container,
    MenuProvider {
    private val items: NonNullList<ItemStack> = NonNullList.withSize(CONTAINER_SIZE, ItemStack.EMPTY)

    var structureId: UUID? = null
        private set
    var activeVolumeBlock: Boolean = false
        private set
    var zoneCount: Int = 0
        private set
    var chamberBlockCount: Int = 0
        private set
    var portCount: Int = 0
        private set
    var formationState: ReactorControllerFormationState = ReactorControllerFormationState.NOT_FORMED
        private set
    var diagnostic: String? = null
        private set
    private var nativeBinding: String = "pending"
    private var temperatureKelvin: Double? = null
    private var pressurePascal: Double? = null
    private var mixtureEntries: List<ReactorMixtureViewEntry> = emptyList()

    fun visualGroupKey(): UUID? =
        structureId?.takeIf { activeVolumeBlock }

    fun applyWorldProjection(
        membership: ReactorBlockMembership?,
        controllerViewState: ReactorControllerViewState? = null,
    ): Boolean {
        val newStructureId = membership?.structureId?.value ?: controllerViewState?.structureId?.value
        val normalizedActive = membership?.activeVolumeBlock == true
        val normalizedZoneCount = membership?.summary?.zoneCount ?: controllerViewState?.zoneCount ?: 0
        val normalizedChamberBlockCount = membership?.summary?.chamberBlockCount ?: controllerViewState?.chamberBlockCount ?: 0
        val normalizedPortCount = membership?.summary?.portCount ?: controllerViewState?.portCount ?: 0
        val normalizedFormationState = controllerViewState?.formationState
            ?: if (membership != null) ReactorControllerFormationState.FORMED else ReactorControllerFormationState.NOT_FORMED
        val normalizedDiagnostic = controllerViewState?.diagnostic
        if (
            structureId == newStructureId &&
            activeVolumeBlock == normalizedActive &&
            zoneCount == normalizedZoneCount &&
            chamberBlockCount == normalizedChamberBlockCount &&
            portCount == normalizedPortCount &&
            formationState == normalizedFormationState &&
            diagnostic == normalizedDiagnostic
        ) {
            return false
        }
        structureId = newStructureId
        activeVolumeBlock = normalizedActive
        zoneCount = normalizedZoneCount
        chamberBlockCount = normalizedChamberBlockCount
        portCount = normalizedPortCount
        formationState = normalizedFormationState
        diagnostic = normalizedDiagnostic
        setChanged()
        refreshVisualModel()
        return true
    }

    override fun read(tag: CompoundTag, registries: HolderLookup.Provider, clientPacket: Boolean) {
        val oldStructureId = structureId
        val oldActiveVolumeBlock = activeVolumeBlock
        super.read(tag, registries, clientPacket)
        structureId = if (tag.hasUUID(STRUCTURE_ID_TAG)) tag.getUUID(STRUCTURE_ID_TAG) else null
        activeVolumeBlock = structureId != null && tag.getBoolean(ACTIVE_VOLUME_TAG)
        zoneCount = if (structureId != null) tag.getInt(ZONE_COUNT_TAG) else 0
        chamberBlockCount = if (structureId != null) tag.getInt(CHAMBER_BLOCK_COUNT_TAG) else 0
        portCount = if (structureId != null) tag.getInt(PORT_COUNT_TAG) else 0
        formationState = tag.getString(FORMATION_STATE_TAG)
            .takeIf { it.isNotBlank() }
            ?.let(ReactorControllerFormationState::valueOf)
            ?: if (structureId != null) ReactorControllerFormationState.FORMED else ReactorControllerFormationState.NOT_FORMED
        diagnostic = tag.getString(DIAGNOSTIC_TAG).takeIf { it.isNotBlank() }
        nativeBinding = tag.getString(NATIVE_BINDING_TAG).takeIf { it.isNotBlank() } ?: "pending"
        temperatureKelvin = tag.getDoubleOrNull(TEMPERATURE_KELVIN_TAG)
        pressurePascal = tag.getDoubleOrNull(PRESSURE_PASCAL_TAG)
        mixtureEntries = tag.getList(MIXTURE_TAG, Tag.TAG_COMPOUND.toInt())
            .mapNotNull { entryTag ->
                val entry = entryTag as? CompoundTag ?: return@mapNotNull null
                val substanceId = entry.getString(MIXTURE_SUBSTANCE_ID_TAG)
                if (substanceId.isBlank()) {
                    return@mapNotNull null
                }
                ReactorMixtureViewEntry(
                    substanceId = substanceId,
                    concentrationMolPerBucket = entry.getDouble(MIXTURE_CONCENTRATION_TAG),
                )
            }
        ContainerHelper.loadAllItems(tag, items, registries)
        if (structureId != oldStructureId || activeVolumeBlock != oldActiveVolumeBlock) {
            refreshVisualModel()
        }
    }

    override fun write(tag: CompoundTag, registries: HolderLookup.Provider, clientPacket: Boolean) {
        super.write(tag, registries, clientPacket)
        structureId?.let { tag.putUUID(STRUCTURE_ID_TAG, it) }
        tag.putBoolean(ACTIVE_VOLUME_TAG, activeVolumeBlock)
        tag.putInt(ZONE_COUNT_TAG, zoneCount)
        tag.putInt(CHAMBER_BLOCK_COUNT_TAG, chamberBlockCount)
        tag.putInt(PORT_COUNT_TAG, portCount)
        tag.putString(FORMATION_STATE_TAG, formationState.name)
        diagnostic?.let { tag.putString(DIAGNOSTIC_TAG, it) }
        tag.putString(NATIVE_BINDING_TAG, nativeBinding)
        temperatureKelvin?.let { tag.putDouble(TEMPERATURE_KELVIN_TAG, it) }
        pressurePascal?.let { tag.putDouble(PRESSURE_PASCAL_TAG, it) }
        val mixtureTag = ListTag()
        for (entry in mixtureEntries) {
            mixtureTag.add(
                CompoundTag().also { entryTag ->
                    entryTag.putString(MIXTURE_SUBSTANCE_ID_TAG, entry.substanceId)
                    entryTag.putDouble(MIXTURE_CONCENTRATION_TAG, entry.concentrationMolPerBucket)
                },
            )
        }
        tag.put(MIXTURE_TAG, mixtureTag)
        ContainerHelper.saveAllItems(tag, items, registries)
    }

    override fun getUpdatePacket(): ClientboundBlockEntityDataPacket =
        ClientboundBlockEntityDataPacket.create(this)

    override fun getUpdateTag(registries: HolderLookup.Provider): CompoundTag =
        saveWithoutMetadata(registries)

    private fun refreshVisualModel() {
        requestModelDataUpdate()
        level?.sendBlockUpdated(blockPos, blockState, blockState, Block.UPDATE_CLIENTS)
    }

    fun applyNativeMetrics(metrics: ReactorTickMetrics): Boolean {
        val nextMixture = metrics.substances
            .take(MAX_CONTROLLER_MIXTURE_ENTRIES)
            .map { substance ->
                ReactorMixtureViewEntry(
                    substanceId = substance.substanceId,
                    concentrationMolPerBucket = substance.concentrationMolPerBucket,
                )
            }
        if (
            nativeBinding == "active" &&
            temperatureKelvin == metrics.temperatureKelvin &&
            pressurePascal == metrics.pressurePascal &&
            mixtureEntries == nextMixture
        ) {
            return false
        }
        nativeBinding = "active"
        temperatureKelvin = metrics.temperatureKelvin
        pressurePascal = metrics.pressurePascal
        mixtureEntries = nextMixture
        setChanged()
        level?.sendBlockUpdated(blockPos, blockState, blockState, Block.UPDATE_CLIENTS)
        return true
    }

    override fun getContainerSize(): Int =
        items.size

    override fun isEmpty(): Boolean =
        items.all(ItemStack::isEmpty)

    override fun addBehaviours(behaviours: MutableList<BlockEntityBehaviour>) {
        if (!reactorKind().isItemPort) {
            return
        }
        val filter = FilteringBehaviour(
            this,
            CenteredSideValueBoxTransform { state, direction ->
                val kind = (state.block as? ReactorMultiblockBlock)?.kind
                kind.isItemPort && direction == state.getValue(ReactorMultiblockBlock.FACING)
            },
        )
            .withPredicate { stack -> stack.isEmpty || stack.item != Items.AIR }
            .onlyActiveWhen { reactorKind().isItemPort }
            .showCount()
        filter.setLabel(Component.translatable("container.create_thermodynamics.reactor_port.filter"))
        behaviours += filter
    }

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

    fun firstPortInputStack(): PortItemStack? {
        check(reactorKind() == ReactorMultiblockKind.ITEM_INPUT_PORT) {
            "reactor block entity at $blockPos is not an item input port"
        }
        for (slot in BUFFER_SLOT_RANGE) {
            val stack = items[slot]
            val itemId = if (stack.isEmpty) null else BuiltInRegistries.ITEM.getKey(stack.item).toString()
            if (itemId != null && portFilterAllows(stack)) {
                return PortItemStack(
                    slot = slot,
                    itemId = itemId,
                    count = stack.count,
                )
            }
        }
        return null
    }

    fun removeConfirmedPortInput(itemId: String, count: Int): Int {
        check(reactorKind() == ReactorMultiblockKind.ITEM_INPUT_PORT) {
            "reactor block entity at $blockPos is not an item input port"
        }
        require(itemId.isNotBlank()) { "itemId must not be blank" }
        require(count >= 0) { "count must be non-negative" }
        if (count == 0) {
            return 0
        }
        val available = items
            .asSequence()
            .filter { !it.isEmpty && BuiltInRegistries.ITEM.getKey(it.item).toString() == itemId }
            .sumOf { it.count }
        check(available >= count) {
            "reactor input port at $blockPos cannot remove $count of $itemId after native acceptance; only $available remain"
        }
        var remaining = count
        var removed = 0
        for (slot in BUFFER_SLOT_RANGE) {
            if (remaining == 0) {
                break
            }
            val stack = items[slot]
            if (stack.isEmpty || BuiltInRegistries.ITEM.getKey(stack.item).toString() != itemId) {
                continue
            }
            val taken = minOf(stack.count, remaining)
            stack.shrink(taken)
            if (stack.isEmpty) {
                items[slot] = ItemStack.EMPTY
            }
            remaining -= taken
            removed += taken
        }
        if (removed > 0) {
            setChanged()
        }
        return removed
    }

    fun insertablePortOutputCount(itemId: String, maxCount: Int): Int {
        check(reactorKind() == ReactorMultiblockKind.ITEM_OUTPUT_PORT) {
            "reactor block entity at $blockPos is not an item output port"
        }
        require(maxCount >= 0) { "maxCount must be non-negative" }
        if (maxCount == 0) {
            return 0
        }
        val template = stackForItemId(itemId)
        var remaining = maxCount
        for (slot in BUFFER_SLOT_RANGE) {
            val stack = items[slot]
            if (remaining == 0) {
                break
            }
            if (stack.isEmpty) {
                remaining -= minOf(template.maxStackSize, remaining)
            } else if (ItemStack.isSameItemSameComponents(stack, template)) {
                remaining -= minOf(stack.maxStackSize - stack.count, remaining)
            }
        }
        return maxCount - remaining
    }

    fun insertConfirmedPortOutput(itemId: String, count: Int): Int {
        check(reactorKind() == ReactorMultiblockKind.ITEM_OUTPUT_PORT) {
            "reactor block entity at $blockPos is not an item output port"
        }
        require(count >= 0) { "count must be non-negative" }
        if (count == 0) {
            return 0
        }
        val template = stackForItemId(itemId)
        val insertable = insertablePortOutputCount(itemId, count)
        check(insertable >= count) {
            "reactor output port at $blockPos cannot accept confirmed output $count of $itemId; only $insertable items fit"
        }
        var remaining = count
        var inserted = 0
        for (slot in BUFFER_SLOT_RANGE) {
            if (remaining == 0) {
                break
            }
            val stack = items[slot]
            if (!stack.isEmpty && ItemStack.isSameItemSameComponents(stack, template)) {
                val added = minOf(stack.maxStackSize - stack.count, remaining)
                if (added > 0) {
                    stack.grow(added)
                    remaining -= added
                    inserted += added
                }
            }
        }
        for (slot in BUFFER_SLOT_RANGE) {
            if (remaining == 0) {
                break
            }
            if (items[slot].isEmpty) {
                val added = minOf(template.maxStackSize, remaining)
                val insertedStack = template.copy()
                insertedStack.count = added
                items[slot] = insertedStack
                remaining -= added
                inserted += added
            }
        }
        if (inserted > 0) {
            setChanged()
        }
        return inserted
    }

    fun portFilterAllowsItemId(itemId: String): Boolean {
        check(reactorKind()?.isPort == true) {
            "reactor block entity at $blockPos is not a reactor port"
        }
        return portFilterAllows(stackForItemId(itemId))
    }

    private fun portFilterAllows(stack: ItemStack): Boolean {
        check(reactorKind()?.isPort == true) {
            "reactor block entity at $blockPos is not a reactor port"
        }
        if (stack.isEmpty) {
            return false
        }
        return getBehaviour(FilteringBehaviour.TYPE)?.test(stack) ?: true
    }

    override fun getDisplayName(): Component =
        if (reactorKind() == ReactorMultiblockKind.CONTROLLER) {
            Component.translatable("container.create_thermodynamics.reactor_controller")
        } else {
            Component.translatable("container.create_thermodynamics.reactor_port")
        }

    override fun createMenu(containerId: Int, playerInventory: Inventory, player: Player): AbstractContainerMenu =
        when (reactorKind()) {
            ReactorMultiblockKind.CONTROLLER -> ReactorControllerMenu(containerId, playerInventory, this)
            ReactorMultiblockKind.ITEM_INPUT_PORT,
            ReactorMultiblockKind.ITEM_OUTPUT_PORT,
            -> ReactorPortMenu(containerId, playerInventory, this, this)

            ReactorMultiblockKind.CHAMBER,
            ReactorMultiblockKind.FLUID_INPUT_PORT,
            ReactorMultiblockKind.FLUID_OUTPUT_PORT,
            null,
            -> error("reactor block entity at $blockPos cannot create a menu for ${blockState.block}")
        }

    fun controllerScreenState(): ReactorControllerViewState =
        ReactorControllerViewState(
            formationState = formationState,
            structureId = structureId?.let(::ReactorStructureId),
            zoneCount = zoneCount,
            chamberBlockCount = chamberBlockCount,
            portCount = portCount,
            diagnostic = diagnostic,
            nativeBinding = nativeBinding,
            temperatureKelvin = temperatureKelvin,
            pressurePascal = pressurePascal,
            mixture = mixtureEntries,
        )

    fun controllerMenuData(): ContainerData =
        object : ContainerData {
            override fun get(index: Int): Int =
                when (index) {
                    CONTROLLER_FORMATION_STATE_DATA_SLOT -> formationState.ordinal
                    CONTROLLER_ZONE_COUNT_DATA_SLOT -> zoneCount
                    CONTROLLER_CHAMBER_BLOCK_COUNT_DATA_SLOT -> chamberBlockCount
                    CONTROLLER_PORT_COUNT_DATA_SLOT -> portCount
                    else -> error("unknown reactor controller data slot $index")
                }

            override fun set(index: Int, value: Int) {
                check(index in 0 until CONTROLLER_DATA_SLOT_COUNT) {
                    "unknown reactor controller data slot $index"
                }
            }

            override fun getCount(): Int =
                CONTROLLER_DATA_SLOT_COUNT
        }

    private fun reactorKind(): ReactorMultiblockKind? =
        (blockState.block as? ReactorMultiblockBlock)?.kind

    private fun stackForItemId(itemId: String): ItemStack {
        require(itemId.isNotBlank()) { "itemId must not be blank" }
        val item = BuiltInRegistries.ITEM.get(ResourceLocation.parse(itemId))
        require(item != Items.AIR) { "unknown item id $itemId" }
        return ItemStack(item)
    }

    companion object {
        private const val CONTAINER_SIZE = 27
        private val BUFFER_SLOT_RANGE = 0 until CONTAINER_SIZE
        private const val STRUCTURE_ID_TAG = "structure_id"
        private const val ACTIVE_VOLUME_TAG = "active_volume"
        private const val ZONE_COUNT_TAG = "zone_count"
        private const val CHAMBER_BLOCK_COUNT_TAG = "chamber_block_count"
        private const val PORT_COUNT_TAG = "port_count"
        private const val FORMATION_STATE_TAG = "formation_state"
        private const val DIAGNOSTIC_TAG = "diagnostic"
        private const val NATIVE_BINDING_TAG = "native_binding"
        private const val TEMPERATURE_KELVIN_TAG = "temperature_kelvin"
        private const val PRESSURE_PASCAL_TAG = "pressure_pascal"
        private const val MIXTURE_TAG = "mixture"
        private const val MIXTURE_SUBSTANCE_ID_TAG = "substance_id"
        private const val MIXTURE_CONCENTRATION_TAG = "concentration_mol_per_bucket"
        private const val MAX_CONTROLLER_MIXTURE_ENTRIES = 6
        private const val CONTROLLER_FORMATION_STATE_DATA_SLOT = 0
        private const val CONTROLLER_ZONE_COUNT_DATA_SLOT = 1
        private const val CONTROLLER_CHAMBER_BLOCK_COUNT_DATA_SLOT = 2
        private const val CONTROLLER_PORT_COUNT_DATA_SLOT = 3
        private const val CONTROLLER_DATA_SLOT_COUNT = 4
    }
}

private fun CompoundTag.getDoubleOrNull(key: String): Double? =
    if (contains(key, Tag.TAG_DOUBLE.toInt())) getDouble(key) else null

private val ReactorMultiblockKind?.isItemPort: Boolean
    get() = this == ReactorMultiblockKind.ITEM_INPUT_PORT || this == ReactorMultiblockKind.ITEM_OUTPUT_PORT

data class PortItemStack(
    val slot: Int,
    val itemId: String,
    val count: Int,
) {
    init {
        require(slot >= 0) { "slot must be non-negative" }
        require(itemId.isNotBlank()) { "itemId must not be blank" }
        require(count > 0) { "count must be positive" }
    }
}
