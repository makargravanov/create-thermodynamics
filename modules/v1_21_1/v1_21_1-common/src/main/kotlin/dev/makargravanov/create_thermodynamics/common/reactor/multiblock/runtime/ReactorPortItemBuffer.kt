package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

class ReactorPortItemBuffer(
    val capacity: Int,
) {
    private val countsByItemId = linkedMapOf<String, Int>()

    val storedItemCount: Int
        get() = countsByItemId.values.sum()

    fun storedCount(itemId: String): Int =
        countsByItemId[itemId] ?: 0

    fun canAccept(itemCount: Int): Boolean =
        itemCount > 0 && storedItemCount + itemCount <= capacity

    fun add(itemId: String, itemCount: Int) {
        require(itemId.isNotBlank()) { "itemId must not be blank" }
        require(canAccept(itemCount)) {
            "cannot accept $itemCount items into reactor port buffer with $storedItemCount/$capacity occupied"
        }
        countsByItemId[itemId] = storedCount(itemId) + itemCount
    }

    fun remove(itemId: String, itemCount: Int) {
        require(itemCount > 0) { "itemCount must be positive" }
        val current = storedCount(itemId)
        require(current >= itemCount) {
            "cannot remove $itemCount items of $itemId from reactor port buffer with only $current stored"
        }
        if (current == itemCount) {
            countsByItemId.remove(itemId)
        } else {
            countsByItemId[itemId] = current - itemCount
        }
    }
}
