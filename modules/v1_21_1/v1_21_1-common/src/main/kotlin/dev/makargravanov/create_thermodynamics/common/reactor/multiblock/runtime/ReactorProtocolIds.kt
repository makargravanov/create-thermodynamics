package dev.makargravanov.create_thermodynamics.common.reactor.multiblock.runtime

@JvmInline
value class ReactorCommandId(val value: Long) {
    init {
        require(value >= 0) { "reactor command id must be non-negative" }
    }
}

@JvmInline
value class ReactorReportId(val value: Long) {
    init {
        require(value >= 0) { "reactor report id must be non-negative" }
    }
}

@JvmInline
value class ReactorSnapshotVersion(val value: Long) {
    init {
        require(value >= 0) { "reactor snapshot version must be non-negative" }
    }
}
