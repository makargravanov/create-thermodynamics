package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative

interface NativeCatalogBridge {
    fun exportCatalogCheckpoint(contentVersion: Long): ByteArray
}

object ThermodynamicsNativeCatalogBridge : NativeCatalogBridge {
    override fun exportCatalogCheckpoint(contentVersion: Long): ByteArray =
        ThermodynamicsNative.exportCatalogCheckpoint(contentVersion)
}
