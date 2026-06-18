package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.assembly

data class ReactorMultiblockRules(
    val chamberVolumeCubicMeters: Double,
    val chamberShapeStrategy: ReactorChamberShapeStrategy = VerticalTankChamberShapeStrategy,
    val minimumChamberBlocks: Int = 1,
    val maximumChamberBlocks: Int? = null,
) {
    init {
        require(chamberVolumeCubicMeters.isFinite() && chamberVolumeCubicMeters > 0.0) {
            "chamberVolumeCubicMeters must be positive and finite"
        }
        require(minimumChamberBlocks > 0) {
            "minimumChamberBlocks must be positive"
        }
        require(maximumChamberBlocks == null || maximumChamberBlocks >= minimumChamberBlocks) {
            "maximumChamberBlocks must be absent or not smaller than minimumChamberBlocks"
        }
    }
}
