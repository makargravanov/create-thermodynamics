package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

data class ReactorCommandOutboxLimits(
    val maxCommands: Int = 1024,
    val maxDrainBatch: Int = 256,
) {
    init {
        require(maxCommands > 0) { "maxCommands must be positive" }
        require(maxDrainBatch > 0) { "maxDrainBatch must be positive" }
    }
}

sealed interface ReactorCommandOutboxResult {
    data class Enqueued(val size: Int) : ReactorCommandOutboxResult
    data class Rejected(val reason: ReactorQueueRejection, val message: String) : ReactorCommandOutboxResult
}

enum class ReactorQueueRejection {
    QUEUE_FULL,
}

class ReactorCommandOutbox(
    private val limits: ReactorCommandOutboxLimits = ReactorCommandOutboxLimits(),
) {
    private val commands = ArrayDeque<ReactorCommand>()

    val size: Int
        get() = commands.size

    val remainingCapacity: Int
        get() = limits.maxCommands - commands.size

    fun isEmpty(): Boolean =
        commands.isEmpty()

    fun enqueue(command: ReactorCommand): ReactorCommandOutboxResult {
        if (commands.size >= limits.maxCommands) {
            return ReactorCommandOutboxResult.Rejected(
                ReactorQueueRejection.QUEUE_FULL,
                "reactor command queue is full: ${commands.size}/${limits.maxCommands}",
            )
        }
        commands.addLast(command)
        return ReactorCommandOutboxResult.Enqueued(commands.size)
    }

    fun drainableCount(maxCommands: Int = limits.maxDrainBatch): Int {
        require(maxCommands > 0) { "drain maxCommands must be positive" }
        return minOf(maxCommands, limits.maxDrainBatch, commands.size)
    }

    fun drain(maxCommands: Int = limits.maxDrainBatch): List<ReactorCommand> {
        val count = drainableCount(maxCommands)
        return buildList(count) {
            repeat(count) {
                add(commands.removeFirst())
            }
        }
    }
}
