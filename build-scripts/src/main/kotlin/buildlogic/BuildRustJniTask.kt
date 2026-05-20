package buildlogic

import javax.inject.Inject
import org.gradle.api.DefaultTask
import org.gradle.api.file.DirectoryProperty
import org.gradle.api.provider.ListProperty
import org.gradle.api.provider.Property
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.InputDirectory
import org.gradle.api.tasks.OutputDirectory
import org.gradle.api.tasks.PathSensitive
import org.gradle.api.tasks.PathSensitivity
import org.gradle.api.tasks.TaskAction
import org.gradle.process.ExecOperations

abstract class BuildRustJniTask : DefaultTask() {
    @get:Inject
    abstract val execOperations: ExecOperations

    @get:InputDirectory
    @get:PathSensitive(PathSensitivity.RELATIVE)
    abstract val crateDirectory: DirectoryProperty

    @get:Input
    abstract val libraryBaseName: Property<String>

    @get:Input
    abstract val resourceRoot: Property<String>

    @get:Input
    abstract val targets: ListProperty<String>

    @get:OutputDirectory
    abstract val outputDirectory: DirectoryProperty

    @TaskAction
    fun buildLibraries() {
        val crateDir = crateDirectory.asFile.get()
        val outputDir = outputDirectory.asFile.get()
        project.delete(outputDir)

        targets.get().distinct().forEach { target ->
            execOperations.exec {
                workingDir(crateDir)
                commandLine(
                    "cargo",
                    "build",
                    "--manifest-path",
                    crateDir.resolve("Cargo.toml").absolutePath,
                    "--lib",
                    "--release",
                    "--target",
                    target,
                )
            }.assertNormalExitValue()

            val libraryName = RustJniPlatforms.libraryFileName(target, libraryBaseName.get())
            val compiledLibrary = crateDir
                .resolve("target")
                .resolve(target)
                .resolve("release")
                .resolve(libraryName)

            check(compiledLibrary.exists()) {
                "Rust JNI artifact not found after build: ${compiledLibrary.absolutePath}"
            }

            val resourceDir = outputDir
                .resolve(resourceRoot.get())
                .resolve(RustJniPlatforms.resourceDirectoryForTriple(target))
            resourceDir.mkdirs()
            compiledLibrary.copyTo(resourceDir.resolve(libraryName), overwrite = true)
        }
    }
}
