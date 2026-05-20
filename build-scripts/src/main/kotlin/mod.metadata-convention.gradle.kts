import java.util.Properties
import org.gradle.language.jvm.tasks.ProcessResources

val modPropertiesFile = rootProject.file("config/mod.properties")
val modProperties = Properties().apply {
    modPropertiesFile.inputStream().use(::load)
}
val expandedModProperties = modProperties.entries.associate { (key, value) ->
    key.toString() to value.toString()
}

extensions.extraProperties["modProperties"] = expandedModProperties

val generateModMetadata = tasks.register<ProcessResources>("generateModMetadata") {
    inputs.file(modPropertiesFile)
    inputs.properties(expandedModProperties)

    expand(expandedModProperties)
    from("src/main/templates")
    into(layout.buildDirectory.dir("generated/sources/modMetadata"))
}

extensions.configure<SourceSetContainer>("sourceSets") {
    named("main") {
        resources.srcDir(generateModMetadata)
    }
}
