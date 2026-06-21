package dev.makargravanov.create_thermodynamics.ui.reactor

data class ReactorControllerUiSnapshot(
    val title: String,
    val status: String,
    val active: Boolean,
    val nativeBinding: String,
    val zoneCount: Int,
    val chamberBlocks: Int,
    val portCount: Int,
    val zones: List<ReactorZoneUiSnapshot>,
) {
    fun toGeneratedState(
        selectedTab: ReactorControllerTab,
        selectedZoneIndex: Int,
    ): ReactorControllerGeneratedState {
        val zone = selectedZone(selectedZoneIndex)
        val mixture = zone?.mixture.orEmpty()
        return ReactorControllerGeneratedState(
            title = title,
            status = status,
            nativeBinding = nativeBinding,
            zoneCountText = "zones $zoneCount",
            chamberBlocksText = "blocks $chamberBlocks",
            portCountText = "ports $portCount",
            overviewVisible = selectedTab == ReactorControllerTab.Overview,
            zonesVisible = selectedTab == ReactorControllerTab.Zones,
            mixtureVisible = selectedTab == ReactorControllerTab.Mixture,
            selectedTabLabel = selectedTab.label,
            selectedZoneTitle = zone?.let { "zone ${it.index}" } ?: "no zone",
            selectedZoneTemperature = zone?.temperature ?: "no data",
            selectedZonePressure = zone?.pressure ?: "",
            selectedZoneMixtureCount = "${mixture.size} substances",
            mixtureTitle = zone?.let { "Mixture: zone ${it.index}" } ?: "Mixture",
            mixtureRow0Name = mixture.getOrNull(0)?.displayName.orEmpty(),
            mixtureRow0Amount = mixture.getOrNull(0)?.concentration.orEmpty(),
            mixtureRow1Name = mixture.getOrNull(1)?.displayName.orEmpty(),
            mixtureRow1Amount = mixture.getOrNull(1)?.concentration.orEmpty(),
            mixtureRow2Name = mixture.getOrNull(2)?.displayName.orEmpty(),
            mixtureRow2Amount = mixture.getOrNull(2)?.concentration.orEmpty(),
            mixtureRow3Name = mixture.getOrNull(3)?.displayName.orEmpty(),
            mixtureRow3Amount = mixture.getOrNull(3)?.concentration.orEmpty(),
            selectOverviewAction = ReactorControllerAction.SelectTab(ReactorControllerTab.Overview),
            selectZonesAction = ReactorControllerAction.SelectTab(ReactorControllerTab.Zones),
            selectMixtureAction = ReactorControllerAction.SelectTab(ReactorControllerTab.Mixture),
            selectZone0Action = zones.sortedBy { it.index }.getOrNull(0)?.let { ReactorControllerAction.SelectZone(it.index) },
            selectZone1Action = zones.sortedBy { it.index }.getOrNull(1)?.let { ReactorControllerAction.SelectZone(it.index) },
            selectZone2Action = zones.sortedBy { it.index }.getOrNull(2)?.let { ReactorControllerAction.SelectZone(it.index) },
            selectZone3Action = zones.sortedBy { it.index }.getOrNull(3)?.let { ReactorControllerAction.SelectZone(it.index) },
            selectZone4Action = zones.sortedBy { it.index }.getOrNull(4)?.let { ReactorControllerAction.SelectZone(it.index) },
            zoneRow0Text = zones.sortedBy { it.index }.getOrNull(0)?.let { "zone ${it.index}" }.orEmpty(),
            zoneRow1Text = zones.sortedBy { it.index }.getOrNull(1)?.let { "zone ${it.index}" }.orEmpty(),
            zoneRow2Text = zones.sortedBy { it.index }.getOrNull(2)?.let { "zone ${it.index}" }.orEmpty(),
            zoneRow3Text = zones.sortedBy { it.index }.getOrNull(3)?.let { "zone ${it.index}" }.orEmpty(),
            zoneRow4Text = zones.sortedBy { it.index }.getOrNull(4)?.let { "zone ${it.index}" }.orEmpty(),
        )
    }

    private fun selectedZone(selectedZoneIndex: Int): ReactorZoneUiSnapshot? =
        zones.firstOrNull { it.index == selectedZoneIndex }
            ?: zones.minByOrNull { it.index }
}

data class ReactorZoneUiSnapshot(
    val index: Int,
    val temperature: String,
    val pressure: String,
    val mixture: List<ReactorMixtureUiLine>,
)

data class ReactorMixtureUiLine(
    val substanceId: String,
    val concentration: String,
) {
    val displayName: String =
        substanceId.substringAfter(':').replace('_', ' ')
}

data class ReactorControllerGeneratedState(
    val title: String,
    val status: String,
    val nativeBinding: String,
    val zoneCountText: String,
    val chamberBlocksText: String,
    val portCountText: String,
    val overviewVisible: Boolean,
    val zonesVisible: Boolean,
    val mixtureVisible: Boolean,
    val selectedTabLabel: String,
    val selectedZoneTitle: String,
    val selectedZoneTemperature: String,
    val selectedZonePressure: String,
    val selectedZoneMixtureCount: String,
    val mixtureTitle: String,
    val mixtureRow0Name: String,
    val mixtureRow0Amount: String,
    val mixtureRow1Name: String,
    val mixtureRow1Amount: String,
    val mixtureRow2Name: String,
    val mixtureRow2Amount: String,
    val mixtureRow3Name: String,
    val mixtureRow3Amount: String,
    val selectOverviewAction: ReactorControllerAction?,
    val selectZonesAction: ReactorControllerAction?,
    val selectMixtureAction: ReactorControllerAction?,
    val selectZone0Action: ReactorControllerAction?,
    val selectZone1Action: ReactorControllerAction?,
    val selectZone2Action: ReactorControllerAction?,
    val selectZone3Action: ReactorControllerAction?,
    val selectZone4Action: ReactorControllerAction?,
    val zoneRow0Text: String,
    val zoneRow1Text: String,
    val zoneRow2Text: String,
    val zoneRow3Text: String,
    val zoneRow4Text: String,
)

enum class ReactorControllerTab(
    val label: String,
) {
    Overview("Overview"),
    Zones("Zones"),
    Mixture("Mixture"),
}

sealed interface ReactorControllerAction {
    data class SelectTab(
        val tab: ReactorControllerTab,
    ) : ReactorControllerAction

    data class SelectZone(
        val zoneIndex: Int,
    ) : ReactorControllerAction
}

object ReactorControllerUiSize {
    const val Width: Int = 232
    const val Height: Int = 140
}
