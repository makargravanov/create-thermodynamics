package buildlogic

import javax.inject.Inject
import org.gradle.api.DefaultTask
import org.gradle.api.file.DirectoryProperty
import org.gradle.api.tasks.InputDirectory
import org.gradle.api.tasks.PathSensitive
import org.gradle.api.tasks.PathSensitivity
import org.gradle.api.tasks.TaskAction
import org.gradle.process.ExecOperations

abstract class TestRustJniTask : DefaultTask() {
    @get:Inject
    abstract val execOperations: ExecOperations

    @get:InputDirectory
    @get:PathSensitive(PathSensitivity.RELATIVE)
    abstract val crateDirectory: DirectoryProperty

    @TaskAction
    fun runTests() {
        val crateDir = crateDirectory.asFile.get()
        execOperations.exec {
            workingDir(crateDir)
            commandLine(
                "cargo",
                "test",
                "--manifest-path",
                crateDir.resolve("Cargo.toml").absolutePath,
                "--lib",
            )
        }.assertNormalExitValue()
    }
}
