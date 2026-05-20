# Architecture

This repository uses three Kotlin module layers:

- `modules/core`: pure Kotlin code with no Minecraft or loader imports.
- `modules/v1_21_1/v1_21_1-common`: Minecraft-facing code for Minecraft 1.21.1 that is not tied to a specific loader.
- `modules/v1_21_1/v1_21_1-neoforge`: NeoForge bootstrap, event wiring, metadata, run configs, and thin compatibility glue.

The root build is an aggregator. Reusable Gradle behavior lives in `build-scripts` as precompiled convention plugins.
