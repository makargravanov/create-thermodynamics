package buildlogic

import javax.inject.Inject
import org.gradle.api.file.DirectoryProperty
import org.gradle.api.model.ObjectFactory
import org.gradle.api.provider.ListProperty
import org.gradle.api.provider.Property

abstract class RustJniExtension @Inject constructor(
    objects: ObjectFactory,
) {
    val crateDirectory: DirectoryProperty = objects.directoryProperty()
    val libraryBaseName: Property<String> = objects.property(String::class.java)
    val resourceRoot: Property<String> = objects.property(String::class.java)
    val targets: ListProperty<String> = objects.listProperty(String::class.java)
}
