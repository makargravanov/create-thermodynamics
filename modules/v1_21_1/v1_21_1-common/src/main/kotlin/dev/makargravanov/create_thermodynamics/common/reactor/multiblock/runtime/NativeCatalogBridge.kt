package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

import dev.makargravanov.create_thermodynamics.common.rust.ThermodynamicsNative

interface NativeCatalogBridge {
    fun exportCatalogCheckpoint(contentVersion: Long): ByteArray
    fun importCatalogCheckpoint(encoded: ByteArray)
}

object ThermodynamicsNativeCatalogBridge : NativeCatalogBridge {
    override fun exportCatalogCheckpoint(contentVersion: Long): ByteArray =
        ThermodynamicsNative.exportCatalogCheckpoint(contentVersion)

    override fun importCatalogCheckpoint(encoded: ByteArray) {
        ThermodynamicsNative.importCatalogCheckpoint(encoded)
    }
}
