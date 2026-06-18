package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model

data class ReactorBlockPosition(
    val x: Int,
    val y: Int,
    val z: Int,
) : Comparable<ReactorBlockPosition> {
    fun faceNeighbours(): List<ReactorBlockPosition> =
        listOf(
            copy(x = x + 1),
            copy(x = x - 1),
            copy(y = y + 1),
            copy(y = y - 1),
            copy(z = z + 1),
            copy(z = z - 1),
        )

    override fun compareTo(other: ReactorBlockPosition): Int =
        compareValuesBy(this, other, ReactorBlockPosition::x, ReactorBlockPosition::y, ReactorBlockPosition::z)
}
