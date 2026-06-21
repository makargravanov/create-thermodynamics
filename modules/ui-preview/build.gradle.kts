plugins {
    id("mod.kotlin-convention")
}

dependencies {
    implementation(project(":modules:ui"))
    implementation("ru.lazyhat:kraft-ui-dsl")
}

val minecraftVersion = providers.gradleProperty("minecraft_version")

val minecraftClientJar = providers.provider {
    val version = minecraftVersion.get()
    val home = providers.systemProperty("user.home").get()
    listOf(
        file("$home/.gradle/caches/neoformruntime/artifacts/minecraft_${version}_client.jar"),
        file("$home/.gradle/caches/fabric-loom/$version/minecraft-client.jar"),
    ).firstOrNull { it.isFile }
        ?: error("Minecraft client jar for UI previews was not found for version $version")
}

tasks.register<JavaExec>("renderUiPreviews") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Renders compile-time UI previews to PNG images."

    classpath = sourceSets.main.get().runtimeClasspath
    mainClass.set("dev.makargravanov.create_thermodynamics.ui.preview.UiPreviewMainKt")
    args(layout.buildDirectory.dir("reports/ui").get().asFile.absolutePath)
    systemProperty("kraftui.minecraftVersion", minecraftVersion.get())
    systemProperty("kraftui.minecraftClientJar", minecraftClientJar.get().absolutePath)
}

tasks.withType<Test>().configureEach {
    systemProperty("kraftui.minecraftVersion", minecraftVersion.get())
    systemProperty("kraftui.minecraftClientJar", minecraftClientJar.get().absolutePath)
}
