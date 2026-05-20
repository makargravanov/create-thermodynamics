import buildlogic.CheckArchitectureBoundaryTask

plugins {
    id("mod.kotlin-convention")
}

tasks.register<CheckArchitectureBoundaryTask>("checkPureCoreBoundary") {
    group = LifecycleBasePlugin.VERIFICATION_GROUP
    description = "Ensures core stays independent from Minecraft and loader APIs."

    sourceRoot.set(layout.projectDirectory.dir("src/main/kotlin"))
    forbiddenImports.set(listOf("net.minecraft.", "net.neoforged."))
    forbiddenPathSegments.set(emptyList())
    failureMessage.set("Core module must not import Minecraft or loader APIs")
}

tasks.named("check") {
    dependsOn("checkPureCoreBoundary")
}
