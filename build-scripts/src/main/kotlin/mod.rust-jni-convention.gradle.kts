import buildlogic.BuildRustJniTask
import buildlogic.RustJniExtension
import buildlogic.RustJniPlatforms
import buildlogic.TestRustJniTask
import org.gradle.language.jvm.tasks.ProcessResources

plugins {
    base
}

val rustTargetsProvider = providers.gradleProperty("rustTargets")
    .map { property ->
        property.split(',')
            .map(String::trim)
            .filter(String::isNotEmpty)
    }
    .orElse(listOf(RustJniPlatforms.hostTargetTriple()))

val rustJni = extensions.create<RustJniExtension>("rustJni")
rustJni.crateDirectory.convention(rootProject.layout.projectDirectory.dir("native/thermodynamics-jni"))
rustJni.libraryBaseName.convention("create_thermodynamics_jni")
rustJni.resourceRoot.convention("natives")
rustJni.targets.convention(rustTargetsProvider)

val buildRustJni = tasks.register<BuildRustJniTask>("buildRustJni") {
    group = LifecycleBasePlugin.BUILD_GROUP
    description = "Builds Rust JNI libraries and stages them as JVM resources."

    crateDirectory.set(rustJni.crateDirectory)
    libraryBaseName.set(rustJni.libraryBaseName)
    resourceRoot.set(rustJni.resourceRoot)
    targets.set(rustJni.targets)
    outputDirectory.set(layout.buildDirectory.dir("generated/resources/rustJni"))
}

tasks.register<TestRustJniTask>("testRustJni") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Runs Rust crate tests for the JNI bridge."

    crateDirectory.set(rustJni.crateDirectory)
}

extensions.configure<SourceSetContainer>("sourceSets") {
    named("main") {
        resources.srcDir(buildRustJni)
    }
}

tasks.named<ProcessResources>("processResources") {
    dependsOn(buildRustJni)
}

tasks.named("check") {
    dependsOn("testRustJni")
}
