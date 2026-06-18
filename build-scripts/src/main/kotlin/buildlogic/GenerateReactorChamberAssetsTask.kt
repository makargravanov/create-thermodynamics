package buildlogic

import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.awt.Color
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import java.io.ByteArrayInputStream
import java.io.File
import java.util.Base64
import javax.imageio.ImageIO
import org.gradle.api.DefaultTask
import org.gradle.api.file.DirectoryProperty
import org.gradle.api.file.RegularFileProperty
import org.gradle.api.provider.Property
import org.gradle.api.tasks.Input
import org.gradle.api.tasks.InputFile
import org.gradle.api.tasks.OutputDirectory
import org.gradle.api.tasks.PathSensitive
import org.gradle.api.tasks.PathSensitivity
import org.gradle.api.tasks.TaskAction

abstract class GenerateReactorChamberAssetsTask : DefaultTask() {
    @get:InputFile
    @get:PathSensitive(PathSensitivity.RELATIVE)
    abstract val blockbenchModel: RegularFileProperty

    @get:InputFile
    @get:PathSensitive(PathSensitivity.RELATIVE)
    abstract val connectedTemplate: RegularFileProperty

    @get:Input
    abstract val namespace: Property<String>

    @get:Input
    abstract val blockId: Property<String>

    @get:OutputDirectory
    abstract val outputAssetsRoot: DirectoryProperty

    @get:OutputDirectory
    abstract val previewDirectory: DirectoryProperty

    @TaskAction
    fun generate() {
        val generator = ReactorChamberAssetGenerator(
            namespace = namespace.get(),
            blockId = blockId.get(),
            blockbenchModel = blockbenchModel.asFile.get(),
            connectedTemplate = connectedTemplate.asFile.get(),
            outputAssetsRoot = outputAssetsRoot.asFile.get(),
            previewDirectory = previewDirectory.asFile.get(),
        )
        generator.generate()
    }
}

private class ReactorChamberAssetGenerator(
    private val namespace: String,
    private val blockId: String,
    private val blockbenchModel: File,
    private val connectedTemplate: File,
    private val outputAssetsRoot: File,
    private val previewDirectory: File,
) {
    private val faces = BlockFace.entries

    fun generate() {
        val model = loadModel()
        val baseImage = readEmbeddedBaseTexture(model)
        val template = readPng(connectedTemplate)

        check(baseImage.width == 64 && baseImage.height == 64) {
            "Expected ${blockbenchModel.path} base texture to be 64x64, got ${baseImage.width}x${baseImage.height}"
        }
        check(template.width == 64 && template.height == 64) {
            "Expected ${connectedTemplate.path} to be 64x64, got ${template.width}x${template.height}"
        }

        cleanupGeneratedFiles()
        val baseFaceTextures = writeBaseFaceTextures(model, baseImage)
        val connectedTexture = writeConnectedTexture(template)
        writeBlockModels(model)
        writeBlockState()
        writeAtlas()
        validateGeneratedResources(baseFaceTextures, connectedTexture)
        writePreview(connectedTexture)
    }

    private fun loadModel(): BlockbenchCubeModel {
        val parsed = JsonSlurper().parse(blockbenchModel) as? Map<*, *>
            ?: error("Blockbench model root must be a JSON object: ${blockbenchModel.path}")
        val elements = parsed["elements"] as? List<*>
            ?: error("Blockbench model does not contain an elements array: ${blockbenchModel.path}")
        val element = elements.singleOrNull() as? Map<*, *>
            ?: error("Reactor chamber model must contain exactly one cube element: ${blockbenchModel.path}")
        val textures = parsed["textures"] as? List<*>
            ?: error("Blockbench model does not contain a textures array: ${blockbenchModel.path}")
        val texture = textures
            .filterIsInstance<Map<*, *>>()
            .singleOrNull { it["id"] == "0" }
            ?: error("Blockbench model must contain exactly one base texture with id=0: ${blockbenchModel.path}")

        val from = numberList(element["from"], "element.from", 3)
        val to = numberList(element["to"], "element.to", 3)
        val rawFaces = element["faces"] as? Map<*, *>
            ?: error("Cube element does not contain faces: ${blockbenchModel.path}")
        val faceUvs = faces.associateWith { face ->
            val rawFace = rawFaces[face.id] as? Map<*, *>
                ?: error("Cube element does not contain face '${face.id}': ${blockbenchModel.path}")
            val uv = numberList(rawFace["uv"], "face.${face.id}.uv", 4)
            FaceUv(uv[0], uv[1], uv[2], uv[3])
        }
        val source = texture["source"] as? String
            ?: error("Base texture id=0 does not contain an embedded source: ${blockbenchModel.path}")

        return BlockbenchCubeModel(from = from, to = to, faceUvs = faceUvs, baseTextureSource = source)
    }

    private fun numberList(value: Any?, name: String, size: Int): List<Double> {
        val list = value as? List<*>
            ?: error("$name must be an array in ${blockbenchModel.path}")
        check(list.size == size) {
            "$name must contain $size numbers in ${blockbenchModel.path}"
        }
        return list.mapIndexed { index, item ->
            (item as? Number)?.toDouble()
                ?: error("$name[$index] must be a number in ${blockbenchModel.path}")
        }
    }

    private fun readEmbeddedBaseTexture(model: BlockbenchCubeModel): BufferedImage {
        val prefix = "data:image/png;base64,"
        check(model.baseTextureSource.startsWith(prefix)) {
            "Base texture id=0 must be an embedded PNG data URI in ${blockbenchModel.path}"
        }
        val bytes = Base64.getDecoder().decode(model.baseTextureSource.removePrefix(prefix))
        return readPng(bytes)
    }

    private fun writeBaseFaceTextures(model: BlockbenchCubeModel, baseImage: BufferedImage): Map<BlockFace, BufferedImage> =
        faces.associateWith { face ->
            val texture = cropUv(baseImage, model.faceUvs.getValue(face))
            writePng(texturePath("${blockId}_${face.id}.png"), texture)
            texture
        }

    private fun writeConnectedTexture(template: BufferedImage): BufferedImage {
        val connected = BufferedImage(64, 64, BufferedImage.TYPE_INT_ARGB)
        for (mask in ConnectionMask.all()) {
            val sourceTile = cropMaskTile(template, mask)
            val createIndex = rectangleCtIndex(mask)
            val targetX = (createIndex % 4) * 16
            val targetY = (createIndex / 4) * 16
            for (y in 0 until 16) {
                for (x in 0 until 16) {
                    connected.setRGB(targetX + x, targetY + y, sourceTile.getRGB(x, y))
                }
            }
        }
        writePng(texturePath("${blockId}_connected.png"), connected)
        return connected
    }

    private fun writeBlockModels(model: BlockbenchCubeModel) {
        writeJson(modelPath("$blockId.json"), baseBlockModel(model))
        writeJson(modelPath("../item/$blockId.json"), mapOf("parent" to "$namespace:block/$blockId"))
    }

    private fun baseBlockModel(model: BlockbenchCubeModel): Map<String, Any> =
        mapOf(
            "ambientocclusion" to false,
            "textures" to linkedMapOf<String, String>().apply {
                put("particle", "$namespace:block/${blockId}_north")
                faces.forEach { face -> put(face.id, "$namespace:block/${blockId}_${face.id}") }
            },
            "elements" to listOf(elementJson(model) { face -> "#${face.id}" }),
        )

    private fun elementJson(model: BlockbenchCubeModel, textureFor: (BlockFace) -> String): Map<String, Any> =
        mapOf(
            "from" to model.from,
            "to" to model.to,
            "shade" to false,
            "faces" to faces.associate { face ->
                face.id to mapOf(
                    "uv" to listOf(0, 0, 16, 16),
                    "texture" to textureFor(face),
                    "cullface" to face.id,
                )
            },
        )

    private fun writeBlockState() {
        writeJson(
            outputAssetsRoot.resolve("blockstates/$blockId.json"),
            mapOf(
                "variants" to listOf("down", "up", "north", "south", "west", "east").associate { facing ->
                    "facing=$facing" to mapOf("model" to "$namespace:block/$blockId")
                },
            ),
        )
    }

    private fun writeAtlas() {
        writeJson(
            outputAssetsRoot.resolve("atlases/blocks.json"),
            mapOf(
                "sources" to listOf(
                    mapOf(
                        "type" to "minecraft:single",
                        "resource" to "$namespace:block/${blockId}_connected",
                    ),
                ),
            ),
        )
    }

    private fun validateGeneratedResources(
        baseFaceTextures: Map<BlockFace, BufferedImage>,
        connectedTexture: BufferedImage,
    ) {
        check(baseFaceTextures.size == 6) {
            "Expected six base face textures, got ${baseFaceTextures.size}"
        }
        check(connectedTexture.width == 64 && connectedTexture.height == 64) {
            "Connected texture sheet must be 64x64, got ${connectedTexture.width}x${connectedTexture.height}"
        }
        for ((face, image) in baseFaceTextures) {
            check(image.width == 16 && image.height == 16) {
                "Base texture for ${face.id} must be 16x16, got ${image.width}x${image.height}"
            }
        }

        val sample = PreviewVolume.box(width = 3, height = 3, depth = 2)
        sample.assertVisibleFaceMasksAreConsistent()
    }

    private fun writePreview(connectedTexture: BufferedImage) {
        previewDirectory.mkdirs()
        val volume = PreviewVolume.box(width = 3, height = 3, depth = 2)
        val scale = 4
        val tile = 16 * scale
        val gap = 12
        val labelHeight = 18
        val panelWidth = tile * 3
        val panelHeight = tile * 3 + labelHeight
        val image = BufferedImage(panelWidth * 3 + gap * 4, panelHeight + gap * 2, BufferedImage.TYPE_INT_ARGB)
        val graphics = image.createGraphics()
        graphics.color = Color(24, 28, 30)
        graphics.fillRect(0, 0, image.width, image.height)

        drawPanel(graphics, "north wall", gap, gap, BlockFace.NORTH, volume, connectedTexture, scale)
        drawPanel(graphics, "top", gap * 2 + panelWidth, gap, BlockFace.UP, volume, connectedTexture, scale)
        drawPanel(graphics, "east wall", gap * 3 + panelWidth * 2, gap, BlockFace.EAST, volume, connectedTexture, scale)
        graphics.dispose()
        ImageIO.write(image, "png", previewDirectory.resolve("${blockId}_connections_preview.png"))
    }

    private fun drawPanel(
        graphics: Graphics2D,
        label: String,
        originX: Int,
        originY: Int,
        face: BlockFace,
        volume: PreviewVolume,
        connectedTexture: BufferedImage,
        scale: Int,
    ) {
        val tile = 16 * scale
        graphics.color = Color(235, 238, 238)
        graphics.drawString(label, originX, originY + 12)
        for (slot in volume.visibleSlots(face)) {
            val state = volume.stateAt(slot.position)
            val mask = face.connectionMask(state)
            val texture = cropCreateTile(connectedTexture, rectangleCtIndex(mask))
            graphics.drawImage(scaleNearest(texture, scale), originX + slot.x * tile, originY + 18 + slot.y * tile, null)
            if (slot.embeddedController) {
                graphics.color = Color(218, 154, 44, 210)
                graphics.drawRect(originX + slot.x * tile + tile / 4, originY + 18 + slot.y * tile + tile / 4, tile / 2, tile / 2)
            }
        }
    }

    private fun scaleNearest(source: BufferedImage, scale: Int): BufferedImage {
        val image = BufferedImage(source.width * scale, source.height * scale, BufferedImage.TYPE_INT_ARGB)
        for (y in 0 until image.height) {
            for (x in 0 until image.width) {
                image.setRGB(x, y, source.getRGB(x / scale, y / scale))
            }
        }
        return image
    }

    private fun cleanupGeneratedFiles() {
        outputAssetsRoot.resolve("blockstates/$blockId.json").delete()
        outputAssetsRoot.resolve("models/item/$blockId.json").delete()
        outputAssetsRoot.resolve("models/block/$blockId.json").delete()
        texturePath("${blockId}.png").delete()
        texturePath("${blockId}_connected_template.png").delete()
        texturePath("${blockId}_connected.png").delete()
        outputAssetsRoot.resolve("atlases/blocks.json").delete()
        faces.forEach { face -> texturePath("${blockId}_${face.id}.png").delete() }
        for (mask in 0 until 16) {
            texturePath("${blockId}_connected_${mask.toString(16)}.png").delete()
            faces.forEach { face ->
                texturePath("${blockId}_connected_${face.id}_${mask.toString(16)}.png").delete()
            }
        }
        for (mask in 0 until 64) {
            modelPath("${blockId}_${mask.toString(16).padStart(2, '0')}.json").delete()
        }
    }

    private fun cropMaskTile(template: BufferedImage, mask: ConnectionMask): BufferedImage {
        val x = (mask.value % 4) * 16
        val y = (mask.value / 4) * 16
        return cropRect(template, x, y, 16, 16)
    }

    private fun cropCreateTile(connectedTexture: BufferedImage, index: Int): BufferedImage {
        val x = (index % 4) * 16
        val y = (index / 4) * 16
        return cropRect(connectedTexture, x, y, 16, 16)
    }

    private fun rectangleCtIndex(mask: ConnectionMask): Int {
        val x = when {
            mask.left && mask.right -> 2
            mask.left -> 3
            mask.right -> 1
            else -> 0
        }
        val y = when {
            mask.up && mask.down -> 1
            mask.up -> 2
            mask.down -> 0
            else -> 3
        }
        return x + y * 4
    }

    private fun cropUv(image: BufferedImage, uv: FaceUv): BufferedImage {
        val output = BufferedImage(16, 16, BufferedImage.TYPE_INT_ARGB)
        for (y in 0 until 16) {
            for (x in 0 until 16) {
                val sourceX = (uv.u1 + ((uv.u2 - uv.u1) * (x + 0.5)) / 16.0)
                    .toInt()
                    .coerceIn(0, image.width - 1)
                val sourceY = (uv.v1 + ((uv.v2 - uv.v1) * (y + 0.5)) / 16.0)
                    .toInt()
                    .coerceIn(0, image.height - 1)
                output.setRGB(x, y, image.getRGB(sourceX, sourceY))
            }
        }
        return output
    }

    private fun cropRect(image: BufferedImage, x: Int, y: Int, width: Int, height: Int): BufferedImage {
        val output = BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB)
        for (yy in 0 until height) {
            for (xx in 0 until width) {
                output.setRGB(xx, yy, image.getRGB(x + xx, y + yy))
            }
        }
        return output
    }

    private fun readPng(file: File): BufferedImage =
        ImageIO.read(file)?.toArgb()
            ?: error("Unable to read PNG image: ${file.path}")

    private fun readPng(bytes: ByteArray): BufferedImage =
        ImageIO.read(ByteArrayInputStream(bytes))?.toArgb()
            ?: error("Unable to read embedded PNG image from ${blockbenchModel.path}")

    private fun BufferedImage.toArgb(): BufferedImage {
        if (type == BufferedImage.TYPE_INT_ARGB) {
            return this
        }
        val converted = BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB)
        val graphics = converted.createGraphics()
        graphics.drawImage(this, 0, 0, null)
        graphics.dispose()
        return converted
    }

    private fun writePng(file: File, image: BufferedImage) {
        file.parentFile.mkdirs()
        check(ImageIO.write(image, "png", file)) {
            "No PNG writer is available for ${file.path}"
        }
    }

    private fun writeJson(file: File, value: Any) {
        file.parentFile.mkdirs()
        file.writeText(JsonOutput.prettyPrint(JsonOutput.toJson(value)) + "\n")
    }

    private fun texturePath(name: String): File =
        outputAssetsRoot.resolve("textures/block/$name")

    private fun modelPath(name: String): File =
        outputAssetsRoot.resolve("models/block/$name")
}

private enum class BlockFace(val id: String) {
    NORTH("north"),
    EAST("east"),
    SOUTH("south"),
    WEST("west"),
    UP("up"),
    DOWN("down");

    fun connectionMask(state: ChamberState): ConnectionMask =
        when (this) {
            NORTH -> ConnectionMask(
                up = state.up,
                right = state.west,
                down = state.down,
                left = state.east,
            )
            EAST -> ConnectionMask(
                up = state.up,
                right = state.north,
                down = state.down,
                left = state.south,
            )
            SOUTH -> ConnectionMask(
                up = state.up,
                right = state.east,
                down = state.down,
                left = state.west,
            )
            WEST -> ConnectionMask(
                up = state.up,
                right = state.south,
                down = state.down,
                left = state.north,
            )
            UP -> ConnectionMask(
                up = state.north,
                right = state.east,
                down = state.south,
                left = state.west,
            )
            DOWN -> ConnectionMask(
                up = state.south,
                right = state.east,
                down = state.north,
                left = state.west,
            )
        }
}

private data class ConnectionMask(val value: Int) {
    val hex: String = value.toString(16)
    val up: Boolean = value and UP != 0
    val right: Boolean = value and RIGHT != 0
    val down: Boolean = value and DOWN != 0
    val left: Boolean = value and LEFT != 0

    constructor(up: Boolean, right: Boolean, down: Boolean, left: Boolean) : this(
        (if (up) UP else 0) or
            (if (right) RIGHT else 0) or
            (if (down) DOWN else 0) or
            (if (left) LEFT else 0),
    )

    companion object {
        private const val UP = 1
        private const val RIGHT = 2
        private const val DOWN = 4
        private const val LEFT = 8

        fun all(): List<ConnectionMask> = (0 until 16).map(::ConnectionMask)
    }
}

private data class ChamberState(
    val north: Boolean,
    val east: Boolean,
    val south: Boolean,
    val west: Boolean,
    val up: Boolean,
    val down: Boolean,
)

private data class FaceUv(val u1: Double, val v1: Double, val u2: Double, val v2: Double)

private data class BlockbenchCubeModel(
    val from: List<Double>,
    val to: List<Double>,
    val faceUvs: Map<BlockFace, FaceUv>,
    val baseTextureSource: String,
)

private data class PreviewPosition(val x: Int, val y: Int, val z: Int)

private data class PreviewSlot(
    val position: PreviewPosition,
    val x: Int,
    val y: Int,
    val embeddedController: Boolean,
)

private class PreviewVolume(private val blocks: Set<PreviewPosition>) {
    fun stateAt(position: PreviewPosition): ChamberState =
        ChamberState(
            north = blocks.contains(position.copy(z = position.z - 1)),
            east = blocks.contains(position.copy(x = position.x + 1)),
            south = blocks.contains(position.copy(z = position.z + 1)),
            west = blocks.contains(position.copy(x = position.x - 1)),
            up = blocks.contains(position.copy(y = position.y + 1)),
            down = blocks.contains(position.copy(y = position.y - 1)),
        )

    fun visibleSlots(face: BlockFace): List<PreviewSlot> =
        when (face) {
            BlockFace.NORTH -> blocks
                .filter { it.z == blocks.minOf(PreviewPosition::z) }
                .sortedWith(compareByDescending<PreviewPosition> { it.y }.thenBy { it.x })
                .map { PreviewSlot(it, it.x, blocks.maxOf(PreviewPosition::y) - it.y, it.x == 1 && it.y == 1) }
            BlockFace.EAST -> blocks
                .filter { it.x == blocks.maxOf(PreviewPosition::x) }
                .sortedWith(compareByDescending<PreviewPosition> { it.y }.thenBy { it.z })
                .map { PreviewSlot(it, it.z, blocks.maxOf(PreviewPosition::y) - it.y, false) }
            BlockFace.UP -> blocks
                .filter { it.y == blocks.maxOf(PreviewPosition::y) }
                .sortedWith(compareBy<PreviewPosition> { it.z }.thenBy { it.x })
                .map { PreviewSlot(it, it.x, it.z, false) }
            else -> emptyList()
        }

    fun assertVisibleFaceMasksAreConsistent() {
        for (block in blocks) {
            val state = stateAt(block)
            check(state.east == blocks.contains(block.copy(x = block.x + 1)))
            check(state.west == blocks.contains(block.copy(x = block.x - 1)))
            check(state.up == blocks.contains(block.copy(y = block.y + 1)))
            check(state.down == blocks.contains(block.copy(y = block.y - 1)))
            check(state.north == blocks.contains(block.copy(z = block.z - 1)))
            check(state.south == blocks.contains(block.copy(z = block.z + 1)))
        }
    }

    companion object {
        fun box(width: Int, height: Int, depth: Int): PreviewVolume {
            val blocks = buildSet {
                for (x in 0 until width) {
                    for (y in 0 until height) {
                        for (z in 0 until depth) {
                            add(PreviewPosition(x, y, z))
                        }
                    }
                }
            }
            return PreviewVolume(blocks)
        }
    }
}
