package dev.makargravanov.create_thermodynamics.neoforge.block

import dev.makargravanov.create_thermodynamics.neoforge.registry.CreateThermodynamicsRegistries
import net.minecraft.core.BlockPos
import net.minecraft.core.HolderLookup
import net.minecraft.nbt.CompoundTag
import net.minecraft.network.protocol.Packet
import net.minecraft.network.protocol.game.ClientGamePacketListener
import net.minecraft.network.protocol.game.ClientboundBlockEntityDataPacket
import net.minecraft.world.level.block.entity.BlockEntity
import net.minecraft.world.level.block.state.BlockState
import net.minecraft.world.level.block.Block
import java.util.UUID

class ReactorMultiblockBlockEntity(pos: BlockPos, state: BlockState) :
    BlockEntity(CreateThermodynamicsRegistries.reactorMultiblockBlockEntity.get(), pos, state) {
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
        if (structureId != oldStructureId || activeVolumeBlock != oldActiveVolumeBlock) {
            refreshVisualModel()
        }
    }

    override fun saveAdditional(tag: CompoundTag, registries: HolderLookup.Provider) {
        super.saveAdditional(tag, registries)
        structureId?.let { tag.putUUID(STRUCTURE_ID_TAG, it) }
        tag.putBoolean(ACTIVE_VOLUME_TAG, activeVolumeBlock)
    }

    override fun getUpdatePacket(): Packet<ClientGamePacketListener> =
        ClientboundBlockEntityDataPacket.create(this)

    override fun getUpdateTag(registries: HolderLookup.Provider): CompoundTag =
        saveWithoutMetadata(registries)

    private fun refreshVisualModel() {
        requestModelDataUpdate()
        level?.sendBlockUpdated(blockPos, blockState, blockState, Block.UPDATE_CLIENTS)
    }

    companion object {
        private const val STRUCTURE_ID_TAG = "structure_id"
        private const val ACTIVE_VOLUME_TAG = "active_volume"
    }
}
