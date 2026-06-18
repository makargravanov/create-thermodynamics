import buildlogic.CheckArchitectureBoundaryTask

plugins {
    id("mod.neoforge-convention")
}

dependencies {
    implementation(project(":modules:v1_21_1:v1_21_1-common"))
}

val generatedResourcesRoot = layout.projectDirectory.dir("src/generated/resources")

val generateReactorChamberAssets by tasks.registering(Exec::class) {
    group = LifecycleBasePlugin.BUILD_GROUP
    description = "Generates reactor chamber blockstate, models and split textures from Blockbench sources."

    workingDir = rootProject.projectDir
    commandLine(
        "node",
        "blockbench/export_reactor_chamber_assets.mjs",
        generatedResourcesRoot.dir("assets/create_thermodynamics").asFile.absolutePath,
    )

    inputs.file(rootProject.layout.projectDirectory.file("blockbench/reactor_chamber.bbmodel"))
    inputs.file(rootProject.layout.projectDirectory.file("blockbench/reactor_chamber_connected_template.png"))
    inputs.file(rootProject.layout.projectDirectory.file("blockbench/export_reactor_chamber_assets.mjs"))
    outputs.dir(generatedResourcesRoot.dir("assets/create_thermodynamics/blockstates"))
    outputs.dir(generatedResourcesRoot.dir("assets/create_thermodynamics/models"))
    outputs.dir(generatedResourcesRoot.dir("assets/create_thermodynamics/textures"))
}

tasks.named("processResources") {
    dependsOn(generateReactorChamberAssets)
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
