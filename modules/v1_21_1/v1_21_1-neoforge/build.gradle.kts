import buildlogic.CheckArchitectureBoundaryTask
import buildlogic.GenerateReactorChamberAssetsTask

plugins {
    id("mod.neoforge-convention")
}

val modProperties = extensions.extraProperties["modProperties"] as Map<*, *>

fun modProperty(name: String): String =
    requireNotNull(modProperties[name]?.toString()) { "Missing '$name' in mod properties" }

dependencies {
    implementation(project(":modules:v1_21_1:v1_21_1-common"))

    implementation("com.simibubi.create:create-${modProperty("minecraft_version")}:${modProperty("create_version")}:slim") {
        isTransitive = false
    }
    implementation("net.createmod.ponder:ponder-neoforge:${modProperty("ponder_version")}+mc${modProperty("minecraft_version")}")
    compileOnly("dev.engine-room.flywheel:flywheel-neoforge-api-${modProperty("minecraft_version")}:${modProperty("flywheel_version")}")
    runtimeOnly("dev.engine-room.flywheel:flywheel-neoforge-${modProperty("minecraft_version")}:${modProperty("flywheel_version")}")
    implementation("com.tterrag.registrate:Registrate:${modProperty("registrate_version")}")
}

val generatedResourcesRoot = layout.projectDirectory.dir("src/generated/resources")

val generateReactorChamberAssets by tasks.registering(GenerateReactorChamberAssetsTask::class) {
    group = LifecycleBasePlugin.BUILD_GROUP
    description = "Generates reactor chamber resources and connection previews from Blockbench sources."

    namespace.set("create_thermodynamics")
    blockId.set("reactor_chamber")
    blockbenchModel.set(rootProject.layout.projectDirectory.file("blockbench/reactor_chamber.bbmodel"))
    connectedTemplate.set(rootProject.layout.projectDirectory.file("blockbench/reactor_chamber_connected_template.png"))
    outputAssetsRoot.set(generatedResourcesRoot.dir("assets/create_thermodynamics"))
    previewDirectory.set(layout.buildDirectory.dir("reports/reactor-assets"))
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
