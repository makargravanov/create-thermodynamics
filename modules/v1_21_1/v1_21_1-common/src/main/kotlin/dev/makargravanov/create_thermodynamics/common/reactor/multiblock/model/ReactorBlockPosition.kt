package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.model

data class ReactorBlockPosition(
    val x: Int,
    val y: Int,
    val z: Int,
) : Comparable<ReactorBlockPosition> {
    fun neighbour(direction: ReactorBlockDirection): ReactorBlockPosition =
        when (direction) {
            ReactorBlockDirection.EAST -> copy(x = x + 1)
            ReactorBlockDirection.WEST -> copy(x = x - 1)
            ReactorBlockDirection.UP -> copy(y = y + 1)
            ReactorBlockDirection.DOWN -> copy(y = y - 1)
            ReactorBlockDirection.SOUTH -> copy(z = z + 1)
            ReactorBlockDirection.NORTH -> copy(z = z - 1)
        }

    fun faceNeighbours(): List<ReactorBlockPosition> =
        ReactorBlockDirection.entries.map(::neighbour)

    fun directionTo(neighbour: ReactorBlockPosition): ReactorBlockDirection? =
        ReactorBlockDirection.entries.singleOrNull { this.neighbour(it) == neighbour }

    override fun compareTo(other: ReactorBlockPosition): Int =
        compareValuesBy(this, other, ReactorBlockPosition::x, ReactorBlockPosition::y, ReactorBlockPosition::z)
}

enum class ReactorBlockDirection {
    EAST,
    WEST,
    UP,
    DOWN,
    SOUTH,
    NORTH,
}
