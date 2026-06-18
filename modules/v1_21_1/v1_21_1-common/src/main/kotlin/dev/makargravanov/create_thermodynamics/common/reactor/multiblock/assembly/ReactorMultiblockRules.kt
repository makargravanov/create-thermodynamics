package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

data class ReactorMultiblockRules(
    val chamberVolumeCubicMeters: Double,
    val chamberShapeStrategy: ReactorChamberShapeStrategy = VerticalTankChamberShapeStrategy,
    val minimumVolumeBlocks: Int = 1,
    val maximumVolumeBlocks: Int? = null,
) {
    init {
        require(chamberVolumeCubicMeters.isFinite() && chamberVolumeCubicMeters > 0.0) {
            "chamberVolumeCubicMeters must be positive and finite"
        }
        require(minimumVolumeBlocks > 0) {
            "minimumVolumeBlocks must be positive"
        }
        require(maximumVolumeBlocks == null || maximumVolumeBlocks >= minimumVolumeBlocks) {
            "maximumVolumeBlocks must be absent or not smaller than minimumVolumeBlocks"
        }
    }
}
