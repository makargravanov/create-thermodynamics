import buildlogic.CheckArchitectureBoundaryTask

plugins {
    id("mod.kotlin-convention")
    id("mod.rust-jni-convention")
}

dependencies {
    api(project(":modules:core"))
}

rustJni {
    crateDirectory.set(rootProject.layout.projectDirectory.dir("native/thermodynamics-jni"))
    libraryBaseName.set("create_thermodynamics_jni")
}

tasks.register<CheckArchitectureBoundaryTask>("checkCommonLoaderBoundary") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Ensures the common module uses Minecraft APIs only and stays loader-agnostic."

    sourceRoot.set(layout.projectDirectory.dir("src/main/kotlin"))
    forbiddenImports.set(listOf("net.neoforged."))
    forbiddenPathSegments.set(emptyList())
    forbiddenText.set(listOf("ReactorMultiblockDefinition?"))
    failureMessage.set("Common module must not import NeoForge APIs")
}

tasks.named("check") {
    dependsOn("checkCommonLoaderBoundary")
}
