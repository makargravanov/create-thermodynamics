package buildlogic

import java.util.Locale

object RustJniPlatforms {
    fun hostTargetTriple(): String {
        val os = System.getProperty("os.name").lowercase(Locale.ROOT)
        val arch = normalizeArchitecture(System.getProperty("os.arch"))

        return when {
            os.contains("win") -> "$arch-pc-windows-msvc"
            os.contains("mac") -> "$arch-apple-darwin"
            os.contains("linux") -> "$arch-unknown-linux-gnu"
            else -> error("Unsupported host OS for Rust JNI build: ${System.getProperty("os.name")}")
        }
    }

    fun resourceDirectoryForTriple(targetTriple: String): String {
        val triple = targetTriple.lowercase(Locale.ROOT)
        val os = when {
            triple.contains("windows") -> "windows"
            triple.contains("apple-darwin") -> "macos"
            triple.contains("linux") -> "linux"
            else -> error("Unsupported Rust target triple: $targetTriple")
        }
        val arch = normalizeArchitecture(triple.substringBefore('-'))
        return "$os-$arch"
    }

    fun libraryFileName(targetTriple: String, libraryBaseName: String): String = when {
        targetTriple.contains("windows", ignoreCase = true) -> "$libraryBaseName.dll"
        targetTriple.contains("apple-darwin", ignoreCase = true) -> "lib$libraryBaseName.dylib"
        else -> "lib$libraryBaseName.so"
    }

    private fun normalizeArchitecture(value: String): String = when (value.lowercase(Locale.ROOT)) {
        "x86_64", "amd64" -> "x86_64"
        "aarch64", "arm64" -> "aarch64"
        else -> error("Unsupported architecture for Rust JNI build: $value")
    }
}
