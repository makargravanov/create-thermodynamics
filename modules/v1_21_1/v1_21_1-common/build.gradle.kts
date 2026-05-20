import buildlogic.CheckArchitectureBoundaryTask

plugins {
    id("mod.kotlin-convention")
}

dependencies {
    api(project(":modules:core"))
}

tasks.register<CheckArchitectureBoundaryTask>("checkCommonLoaderBoundary") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Ensures the common module uses Minecraft APIs only and stays loader-agnostic."

    sourceRoot.set(layout.projectDirectory.dir("src/main/kotlin"))
    forbiddenImports.set(listOf("net.neoforged."))
    forbiddenPathSegments.set(emptyList())
    failureMessage.set("Common module must not import NeoForge APIs")
}

tasks.named("check") {
    dependsOn("checkCommonLoaderBoundary")
}
