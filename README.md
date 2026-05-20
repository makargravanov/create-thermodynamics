
Installation information
=======

This template repository can be directly cloned to get you started with a new
mod. Simply create a new repository cloned from this one, by following the
instructions provided by [GitHub](https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-repository-from-a-template).

Once you have your clone, simply open the repository in the IDE of your choice. The usual recommendation for an IDE is either IntelliJ IDEA or Eclipse.

If at any point you are missing libraries in your IDE, or you've run into problems you can
run `gradlew --refresh-dependencies` to refresh the local cache. `gradlew clean` to reset everything 
{this does not affect your code} and then start the process again.

Mapping Names:
============
By default, the MDK is configured to use the official mapping names from Mojang for methods and fields 
in the Minecraft codebase. These names are covered by a specific license. All modders should be aware of this
license. For the latest license text, refer to the mapping file itself, or the reference copy here:
https://github.com/NeoForged/NeoForm/blob/main/Mojang.md

Additional Resources: 
==========
Community Documentation: https://docs.neoforged.net/  
NeoForged Discord: https://discord.neoforged.net/

Rust JNI
==========
The project now contains a Rust crate at `native/thermodynamics-jni` and a JNI bridge in `modules/v1_21_1/v1_21_1-common`.

Default Gradle builds package the host platform native library into the mod jar. For example:

`./gradlew build`

To package several targets into one artifact, pass Rust target triples through `rustTargets`:

`./gradlew :modules:v1_21_1:v1_21_1-neoforge:build -PrustTargets=x86_64-pc-windows-msvc,x86_64-unknown-linux-gnu,aarch64-apple-darwin`

Each requested target must already be installed in the local Rust toolchain, for example with:

`rustup target add x86_64-unknown-linux-gnu aarch64-apple-darwin`
