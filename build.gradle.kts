plugins {
    base
}

tasks.register("checkArchitecture") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Runs repository architecture boundary checks."
    dependsOn(
        ":modules:core:checkPureCoreBoundary",
        ":modules:v1_21_1:v1_21_1-common:checkCommonLoaderBoundary",
        ":modules:v1_21_1:v1_21_1-neoforge:checkThinLoaderBoundary",
    )
}

tasks.named("check") {
    dependsOn("checkArchitecture")
}
