import java.util.Properties
import net.neoforged.moddevgradle.dsl.NeoForgeExtension
import org.slf4j.event.Level

plugins {
    id("mod.kotlin-convention")
    id("mod.metadata-convention")
    id("net.neoforged.moddev")
    id("maven-publish")
    idea
}

val modPropertiesFile = rootProject.file("config/mod.properties")
val modProperties = Properties().apply {
    modPropertiesFile.inputStream().use(::load)
}

fun prop(name: String): String = requireNotNull(modProperties.getProperty(name)) {
    "Missing '$name' in ${modPropertiesFile.path}"
}

version = prop("mod_version")
group = prop("mod_group_id")

base {
    archivesName.set(prop("mod_id"))
}

extensions.configure<NeoForgeExtension>("neoForge") {
    version = prop("neo_version")

    parchment {
        minecraftVersion.set(prop("parchment_minecraft_version"))
        mappingsVersion.set(prop("parchment_mappings_version"))
    }

    runs {
        register("client") {
            client()
            systemProperty("neoforge.enabledGameTestNamespaces", prop("mod_id"))
        }

        register("server") {
            server()
            programArgument("--nogui")
            systemProperty("neoforge.enabledGameTestNamespaces", prop("mod_id"))
        }

        register("gameTestServer") {
            type.set("gameTestServer")
            systemProperty("neoforge.enabledGameTestNamespaces", prop("mod_id"))
        }

        register("data") {
            data()
            programArguments.addAll(
                "--mod",
                prop("mod_id"),
                "--all",
                "--output",
                file("src/generated/resources/").absolutePath,
                "--existing",
                file("src/main/resources/").absolutePath,
            )
        }

        configureEach {
            systemProperty("forge.logging.markers", "REGISTRIES")
            logLevel.set(Level.DEBUG)
        }
    }

    mods {
        register(prop("mod_id")) {
            sourceSet(sourceSets.main.get())
        }
    }
}

sourceSets.main {
    resources {
        srcDir("src/generated/resources")
        exclude("**/*.bbmodel")
        exclude("src/generated/**/.cache")
    }
}

configurations {
    create("localRuntime")
    named("runtimeClasspath") {
        extendsFrom(configurations.named("localRuntime").get())
    }
}

tasks.named<Jar>("jar") {
    from(project(":modules:core").sourceSets.main.get().output)
    from(project(":modules:v1_21_1:v1_21_1-common").sourceSets.main.get().output)
}

tasks.named("neoForgeIdeSync") {
    dependsOn("generateModMetadata")
}

publishing {
    publications {
        register<MavenPublication>("mavenJava") {
            from(components["java"])
        }
    }
    repositories {
        maven {
            url = uri(layout.projectDirectory.dir("repo"))
        }
    }
}

idea {
    module {
        isDownloadSources = true
        isDownloadJavadoc = true
    }
}
