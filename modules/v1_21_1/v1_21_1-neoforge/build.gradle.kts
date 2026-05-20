import buildlogic.CheckArchitectureBoundaryTask

plugins {
    id("mod.neoforge-convention")
}

dependencies {
    implementation(project(":modules:v1_21_1:v1_21_1-common"))
}

tasks.register<CheckArchitectureBoundaryTask>("checkThinLoaderBoundary") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Ensures the NeoForge module remains a thin loader leaf."

    sourceRoot.set(layout.projectDirectory.dir("src/main/kotlin"))
    forbiddenImports.set(emptyList())
    forbiddenPathSegments.set(listOf("/content/", "\\content\\"))
    failureMessage.set("NeoForge loader module must not contain content packages")
}

tasks.named("check") {
    dependsOn("checkThinLoaderBoundary")
}
