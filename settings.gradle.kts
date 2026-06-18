pluginManagement {
    repositories {
        gradlePluginPortal()
        mavenCentral()
    }
}

plugins {
    id("org.gradle.toolchains.foojay-resolver-convention") version "1.0.0"
}

dependencyResolutionManagement {
    repositories {
        mavenCentral()
        maven("https://maven.createmod.net")
        maven("https://maven.ithundxr.dev/snapshots")
    }
}

includeBuild("build-scripts")
includeBuild("external/ui-dsl")

rootProject.name = "create-thermodynamics"

include(":modules:core")
include(":modules:ui-preview")
include(":modules:v1_21_1:v1_21_1-common")
include(":modules:v1_21_1:v1_21_1-neoforge")
