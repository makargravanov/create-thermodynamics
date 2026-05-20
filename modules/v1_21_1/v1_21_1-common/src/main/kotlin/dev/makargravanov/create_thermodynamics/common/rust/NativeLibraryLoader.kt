package dev.makargravanov.create_thermodynamics.common.rust

import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption

internal object NativeLibraryLoader {
    private const val RESOURCE_ROOT = "/natives"
    private const val LIBRARY_BASE_NAME = "create_thermodynamics_jni"

    @Volatile
    private var loaded = false

    fun load() {
        if (loaded) {
            return
        }

        synchronized(this) {
            if (loaded) {
                return
            }

            val fileName = libraryFileNameFor(currentOsName())
            val resourcePath = "$RESOURCE_ROOT/${platformDirectory(currentOsName(), currentArchitecture())}/$fileName"
            val extractedLibrary = extractLibrary(resourcePath, fileName)
            System.load(extractedLibrary.toAbsolutePath().toString())
            loaded = true
        }
    }

    internal fun platformDirectory(osName: String, architecture: String): String {
        val os = osName.lowercase()
        val arch = normalizeArchitecture(architecture)

        return when {
            os.contains("win") -> "windows-$arch"
            os.contains("mac") -> "macos-$arch"
            os.contains("linux") -> "linux-$arch"
            else -> error("Unsupported OS for native thermodynamics library: $osName")
        }
    }

    private fun currentOsName(): String = System.getProperty("os.name")

    private fun currentArchitecture(): String = System.getProperty("os.arch")

    private fun libraryFileNameFor(osName: String): String = when {
        osName.contains("win", ignoreCase = true) -> "$LIBRARY_BASE_NAME.dll"
        osName.contains("mac", ignoreCase = true) -> "lib$LIBRARY_BASE_NAME.dylib"
        osName.contains("linux", ignoreCase = true) -> "lib$LIBRARY_BASE_NAME.so"
        else -> error("Unsupported OS for native thermodynamics library: $osName")
    }

    private fun normalizeArchitecture(value: String): String = when (value.lowercase()) {
        "x86_64", "amd64" -> "x86_64"
        "aarch64", "arm64" -> "aarch64"
        else -> error("Unsupported architecture for native thermodynamics library: $value")
    }

    private fun extractLibrary(resourcePath: String, fileName: String): Path {
        val resource = requireNotNull(javaClass.getResourceAsStream(resourcePath)) {
            "Native library resource not found: $resourcePath. " +
                "Build with -PrustTargets to package the required platform binary."
        }

        Files.createTempDirectory("create-thermodynamics-native").also { tempDir ->
            tempDir.toFile().deleteOnExit()
            val extracted = tempDir.resolve(fileName)
            resource.use { input ->
                Files.copy(input, extracted, StandardCopyOption.REPLACE_EXISTING)
            }
            extracted.toFile().deleteOnExit()
            return extracted
        }
    }
}
