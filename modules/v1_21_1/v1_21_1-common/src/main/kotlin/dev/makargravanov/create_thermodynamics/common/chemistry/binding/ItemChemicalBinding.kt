package dev.makargravanov.create_thermodynamics.common.chemistry.binding

data class ItemChemicalBinding(
    val itemId: String,
    val substanceId: String,
    val molPerItem: Double,
) {
    init {
        require(itemId.isNotBlank()) {
            "itemId must not be blank"
        }
        require(substanceId.isNotBlank()) {
            "substanceId must not be blank"
        }
        require(molPerItem.isFinite() && molPerItem > 0.0) {
            "molPerItem must be positive and finite"
        }
    }
}
