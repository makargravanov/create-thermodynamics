package buildlogic

import org.gradle.api.DefaultTask
import org.gradle.api.file.DirectoryProperty
import org.gradle.api.provider.ListProperty
import org.gradle.api.provider.Property
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.InputDirectory
import org.gradle.api.tasks.Optional
import org.gradle.api.tasks.PathSensitive
import org.gradle.api.tasks.PathSensitivity
import org.gradle.api.tasks.TaskAction

abstract class CheckArchitectureBoundaryTask : DefaultTask() {
    @get:InputDirectory
    @get:PathSensitive(PathSensitivity.RELATIVE)
    abstract val sourceRoot: DirectoryProperty

    @get:Input
    abstract val forbiddenImports: ListProperty<String>

    @get:Input
    abstract val forbiddenPathSegments: ListProperty<String>

    @get:Input
    abstract val forbiddenText: ListProperty<String>

    @get:Input
    @get:Optional
    abstract val failureMessage: Property<String>

    @TaskAction
    fun checkBoundary() {
        val root = sourceRoot.asFile.get()
        if (!root.exists()) {
            return
        }

        val violations = root.walkTopDown()
            .filter { it.isFile && it.extension == "kt" }
            .flatMap { file ->
                val pathViolation = forbiddenPathSegments.get().any { segment ->
                    file.path.contains(segment)
                }
                val importViolations = file.readLines().mapIndexedNotNull { index, line ->
                    val trimmed = line.trim()
                    val forbidden = forbiddenImports.get().firstOrNull { trimmed.startsWith("import $it") }
                    if (forbidden == null) null else "${file.relativeTo(root)}:${index + 1}: $trimmed"
                }
                val textViolations = file.readLines().mapIndexedNotNull { index, line ->
                    val forbidden = forbiddenText.getOrElse(emptyList()).firstOrNull { text -> line.contains(text) }
                    if (forbidden == null) null else "${file.relativeTo(root)}:${index + 1}: forbidden text '$forbidden'"
                }

                sequence {
                    if (pathViolation) yield(file.relativeTo(root).path)
                    yieldAll(importViolations)
                    yieldAll(textViolations)
                }
            }
            .toList()

        check(violations.isEmpty()) {
            "${failureMessage.getOrElse("Architecture boundary violation")}:\n${violations.joinToString("\n")}"
        }
    }
}
