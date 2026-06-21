package dev.makargravanov.create_thermodynamics.ui.reactor

import dev.makargravanov.create_thermodynamics.ui.style.ThermodynamicsUiTheme
import ru.lazyhat.kraftui.program.FontMetrics
import ru.lazyhat.kraftui.program.PrimitiveOptimizationOptions
import ru.lazyhat.kraftui.program.PrimitiveOptimizationPass
import ru.lazyhat.kraftui.program.PrimitiveSourceTargets
import ru.lazyhat.kraftui.program.PrimitiveStaticTextureBakingOptions
import ru.lazyhat.kraftui.program.PrimitiveTargetSourceRequest
import ru.lazyhat.kraftui.program.ScreenProgramCompiler
import ru.lazyhat.kraftui.program.generateTargetSource
import ru.lazyhat.kraftui.program.toPrimitiveScreenProgram
import ru.lazyhat.kraftui.style.BakeHint
import ru.lazyhat.kraftui.style.StyleOptimizationHint
import ru.lazyhat.kraftui.style.StyleReport
import ru.lazyhat.kraftui.style.StyleUsage
import java.nio.file.Files
import java.nio.file.Path
import kotlin.io.path.createDirectories
import kotlin.io.path.writeBytes
import kotlin.io.path.writeText

object GenerateReactorControllerMinecraftUi {
    @JvmStatic
    fun main(args: Array<String>) {
        require(args.size == 3) {
            "Expected arguments: <generated kotlin directory> <generated resources directory> <report directory>"
        }
        val kotlinDirectory = Path.of(args[0])
        val resourcesDirectory = Path.of(args[1])
        val reportDirectory = Path.of(args[2])
        clearGeneratedDirectory(kotlinDirectory)
        clearGeneratedDirectory(resourcesDirectory)
        clearGeneratedDirectory(reportDirectory)

        val sampleState = previewState()
        val primitive =
            ScreenProgramCompiler(fontMetrics = FontMetrics { text -> text.length * 6 })
                .compile(
                    root = ReactorControllerUi.build { sampleState },
                    rootWidth = ReactorControllerUi.Width,
                    rootHeight = ReactorControllerUi.Height,
                ).toPrimitiveScreenProgram()

        val result =
            primitive.generateTargetSource(
                target = PrimitiveSourceTargets.minecraftGuiGraphics,
                request =
                    PrimitiveTargetSourceRequest(
                        packageName = "dev.makargravanov.create_thermodynamics.neoforge.client.generated",
                        className = "GeneratedReactorControllerScreen",
                        stateType = "dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerGeneratedState",
                        actionType = "dev.makargravanov.create_thermodynamics.ui.reactor.ReactorControllerAction",
                        optimization =
                            PrimitiveOptimizationOptions(
                                passes = PrimitiveOptimizationPass.default + PrimitiveOptimizationPass.StaticTextureBaking,
                                staticTextureBaking =
                                    PrimitiveStaticTextureBakingOptions.Enabled(
                                        minInstructionCount = 4,
                                        textureNamespace = "create_thermodynamics",
                                        texturePathPrefix = "textures/gui/generated",
                                    ),
                            ),
                    ),
            )

        val sourcePath =
            kotlinDirectory.resolve(
                result.source.packageName.replace('.', '/') + "/" + result.source.className + ".kt",
            )
        sourcePath.parent.createDirectories()
        sourcePath.writeText(result.source.source)

        for (asset in result.source.assets) {
            val assetPath = resourcesDirectory.resolve(asset.path)
            assetPath.parent.createDirectories()
            assetPath.writeBytes(asset.bytes)
        }

        reportDirectory.createDirectories()
        reportDirectory.resolve("reactor-controller-analysis.txt").writeText(result.analysisReport.asText())
        reportDirectory.resolve("reactor-controller-optimization.txt").writeText(result.optimizationReport.asText())
        reportDirectory.resolve("reactor-controller-style-report.txt").writeText(styleReport().asText())
        reportDirectory.resolve("reactor-controller-generated-files.txt").writeText(
            buildString {
                appendLine(sourcePath)
                result.source.assets.forEach { asset ->
                    appendLine(resourcesDirectory.resolve(asset.path))
                }
            },
        )
    }

    private fun styleReport(): StyleReport =
        StyleReport(
            themeName = "ThermodynamicsUiTheme",
            usages =
                listOf(
                    StyleUsage("window", 1),
                    StyleUsage("tab", 3),
                    StyleUsage("metricCard", 6),
                    StyleUsage("panel", 2),
                    StyleUsage("text", 22),
                ),
            diagnostics = emptyList(),
            optimizationHints =
                listOf(
                    StyleOptimizationHint("/root/window", ThermodynamicsUiTheme.theme.styles.window.surface.bakeHint),
                    StyleOptimizationHint("/root/cards", BakeHint.PreferBakedTexture),
                ),
        )

    private fun clearGeneratedDirectory(root: Path) {
        if (!Files.exists(root)) return
        Files.walk(root).use { paths ->
            paths
                .sorted(Comparator.reverseOrder())
                .filter { it != root }
                .forEach(Files::delete)
        }
    }

    private fun previewState(): ReactorControllerGeneratedState =
        ReactorControllerUiSnapshot(
            title = "Reactor Controller",
            status = "formed",
            active = true,
            nativeBinding = "active",
            zoneCount = 1,
            chamberBlocks = 27,
            portCount = 2,
            zones =
                listOf(
                    ReactorZoneUiSnapshot(
                        index = 0,
                        temperature = "298.0 K",
                        pressure = "101.3 kPa",
                        mixture =
                            listOf(
                                ReactorMixtureUiLine("destroy:water", "64.000"),
                                ReactorMixtureUiLine("destroy:ethanol", "2.000"),
                            ),
                    ),
                ),
        ).toGeneratedState(
            selectedTab = ReactorControllerTab.Overview,
            selectedZoneIndex = 0,
        )
}
