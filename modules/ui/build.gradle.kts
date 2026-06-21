import buildlogic.CheckArchitectureBoundaryTask

plugins {
    id("mod.kotlin-convention")
}

val uiSpec by sourceSets.creating {
    kotlin.srcDir("src/uiSpec/kotlin")
    compileClasspath += sourceSets.main.get().output
    runtimeClasspath += output + compileClasspath
}

dependencies {
    "uiSpecImplementation"("ru.lazyhat:kraft-ui-dsl")
}

val generatedMinecraftUiKotlin = layout.buildDirectory.dir("generated/kraftui/neoforge/kotlin")
val generatedMinecraftUiResources = layout.buildDirectory.dir("generated/kraftui/neoforge/resources")
val generatedMinecraftUiReports = layout.buildDirectory.dir("reports/kraftui")

val generateReactorControllerMinecraftUi by tasks.registering(JavaExec::class) {
    group = LifecycleBasePlugin.BUILD_GROUP
    description = "Generates Minecraft source and resources for the reactor controller UI."

    dependsOn(tasks.named(uiSpec.classesTaskName))
    classpath = uiSpec.runtimeClasspath
    mainClass.set("dev.makargravanov.create_thermodynamics.ui.reactor.GenerateReactorControllerMinecraftUi")
    args(
        generatedMinecraftUiKotlin.get().asFile.absolutePath,
        generatedMinecraftUiResources.get().asFile.absolutePath,
        generatedMinecraftUiReports.get().asFile.absolutePath,
    )
    outputs.dir(generatedMinecraftUiKotlin)
    outputs.dir(generatedMinecraftUiResources)
    outputs.dir(generatedMinecraftUiReports)
}

tasks.named("check") {
    dependsOn("checkUiModelBoundary")
    dependsOn(generateReactorControllerMinecraftUi)
}

tasks.register<CheckArchitectureBoundaryTask>("checkUiModelBoundary") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Ensures runtime UI model classes do not depend on the compile-time UI DSL."

    sourceRoot.set(layout.projectDirectory.dir("src/main/kotlin"))
    forbiddenImports.set(listOf("ru.lazyhat.kraftui."))
    forbiddenPathSegments.set(emptyList())
    forbiddenText.set(listOf("ru.lazyhat.kraftui"))
    failureMessage.set("UI model source set must not depend on the compile-time UI DSL")
}
