plugins {
    id("mod.kotlin-convention")
}

dependencies {
    implementation(project(":modules:ui"))
    implementation("ru.lazyhat:kraft-ui-dsl")
}

tasks.register<JavaExec>("renderUiPreviews") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Renders compile-time UI previews to PNG images."

    classpath = sourceSets.main.get().runtimeClasspath
    mainClass.set("dev.makargravanov.create_thermodynamics.ui.preview.UiPreviewMainKt")
    args(layout.buildDirectory.dir("reports/ui").get().asFile.absolutePath)
}
