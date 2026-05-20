# Architecture

Create Thermodynamics is structured as a native-first simulation engine with a
thin JVM and Minecraft integration layer. The repository must grow around that
boundary: Rust owns the scientific model, state, numerical execution, and
performance-sensitive storage; Kotlin owns game integration, lifecycle wiring,
Minecraft/Create adaptation, resource registration, and user-facing glue.

## Module Tree

Current repository shape:

```text
build-scripts/
  src/main/kotlin/
    buildlogic/                         Gradle tasks and boundary checks
    mod.*-convention.gradle.kts         reusable Gradle conventions
config/
  mod.properties                        mod metadata source
docs/                                  architecture and design contracts
modules/
  core/                                pure Kotlin support; no Minecraft imports
  v1_21_1/
    v1_21_1-common/                    Minecraft-facing, loader-agnostic code
    v1_21_1-neoforge/                  NeoForge bootstrap and shims
native/
  thermodynamics-jni/                  Rust native library exposed through JNI
```

Expected long-term shape:

```text
native/
  thermodynamics-core/                 pure Rust domain model and solvers
  thermodynamics-data/                 Rust data schemas, loading, serialization
  thermodynamics-jni/                  ABI-stable JNI facade over native crates
  thermodynamics-bench/                native profiling and benchmark harnesses
modules/
  core/                                Kotlin-side shared contracts only
  v1_21_1/
    v1_21_1-common/                    Minecraft abstractions for one MC version
    v1_21_1-neoforge/                  loader-specific entrypoints and events
```

Only add a module when it owns a real runtime, dependency, or version boundary.
Do not split by vague technical category if the code still changes together.

## Responsibility Boundaries

### Rust Native Core

Rust is the source of truth for:

- chemical and physical domain entities;
- thermodynamic and kinetic state;
- numerical solvers and process simulation;
- native memory layout and persistent simulation state;
- deterministic stepping and rollback-compatible state transitions;
- validation of physical constraints;
- serialization formats for model state and scientific data;
- profiling, benchmarking, and numerical regression tests.

The native core must be usable without Minecraft. A command-line harness,
benchmark suite, or future CAD/education integration should be able to link it
without depending on JVM or game classes.

### JNI Facade

`native/thermodynamics-jni` exposes a stable, narrow ABI. It translates JVM calls
into native operations but must not contain the domain model itself.

The facade is responsible for:

- ABI version reporting;
- opaque native handles for long-lived state;
- conversion between JVM buffers and native data;
- error transfer across the JNI boundary;
- keeping JNI symbols stable for released mod versions;
- minimizing per-tick allocation and copying.

### Kotlin Core

`modules/core` is pure Kotlin. It may contain small JVM-side contracts, utility
types, or test helpers, but it must not own chemistry or physics rules.

Rules:

- no `net.minecraft.*` imports;
- no loader imports;
- no Create-specific imports;
- no attempt to duplicate native solver behavior in Kotlin.

### Minecraft Common Layer

`modules/v1_21_1/v1_21_1-common` is Minecraft-facing and loader-agnostic for
Minecraft 1.21.1. It may import `net.minecraft.*`, but must not import
`net.neoforged.*`.

This layer adapts native simulation concepts to game concepts:

- level and chunk lifecycle;
- block entities, contraptions, and Create kinetic networks;
- Minecraft registry and data-pack discovery;
- passing native snapshots through Minecraft save data;
- server tick scheduling;
- client sync payloads.

It should pass coarse commands and packed state to Rust, not run chemistry in
Kotlin collections on every tick.

### Loader Layer

`modules/v1_21_1/v1_21_1-neoforge` is a leaf module. Keep it limited to:

- mod entrypoints;
- event bus wiring;
- registration bootstrap;
- run configs and generated metadata;
- client bootstrap;
- tiny compatibility shims.

Gameplay, simulation rules, native calls, and domain translation must move
inward to common or native modules.

## Data Loading And Chemical Serialization

Scientific and gameplay data are separate concerns.

Scientific data includes elements, isotopes, species, phases, equations of
state, transport models, reaction mechanisms, empirical correlations, and
validated constants. Loading, validation, schema handling, and serialization of
chemical/scientific data belong to a dedicated Rust crate, expected to be
`thermodynamics-data` or an equivalent native workspace member.

Gameplay data includes Minecraft blocks, fluids, items, tags, Create components,
recipes, and server configuration. Kotlin may discover game-side records and
call native loading APIs, but it must not become the serializer or validator for
chemical data. When Kotlin needs chemical data, it should reference native IDs,
handles, compact descriptors, or already-validated data returned by Rust.

Every data source needs:

- schema version;
- provenance when data is physical or empirical;
- unit declarations;
- validation diagnostics;
- deterministic load order;
- stable IDs that are not tied to localized names.

## Serialization

Simulation state and chemical data must be versioned independently from
Minecraft save formats. Minecraft save data may contain a native snapshot, but
the snapshot format is a native contract produced and consumed by Rust.

Kotlin's role is transport and lifecycle integration: request a snapshot from
native code, store the resulting bytes or reference through Minecraft's save
pipeline, and pass it back to native code on load. Kotlin should not parse,
reinterpret, migrate, or author the chemical serialization format.

Required properties:

- explicit schema and ABI versioning;
- deterministic binary representation for hot state;
- human-readable diagnostic exports for debugging;
- forward migration strategy for released worlds;
- no reliance on JVM object serialization;
- no hidden dependence on map iteration order.

## Simulation Lifecycle

A normal server lifecycle should look like this:

1. Minecraft and loader bootstrap initialize Kotlin entrypoints.
2. Kotlin loads mod metadata, game registries, and server config.
3. Kotlin initializes the native library and verifies ABI compatibility.
4. Kotlin passes game data-pack/resource references to native loading APIs when
   chemical data is involved.
5. Native code loads scientific data, validates model schemas, and owns
   chemical serialization.
6. Kotlin creates or restores native simulation worlds by passing native
   snapshot bytes or handles back to Rust.
7. On each scheduled simulation step, Kotlin passes batched world events and
   receives compact state deltas.
8. Kotlin applies game-facing deltas to block entities, Create networks, client
   sync packets, particles, sounds, and UI state.
9. Save operations request native snapshots and write opaque snapshot payloads
   through Minecraft's save pipeline.
10. Shutdown releases native handles in a defined order.

## Boundary Enforcement

Current boundary checks should remain active:

- `modules/core` must not import Minecraft or loader packages.
- `v1_21_1-common` may import Minecraft packages but not NeoForge packages.
- `v1_21_1-neoforge` must stay thin and must not become the home of gameplay
  logic.

Future checks should verify:

- JNI facade code does not contain solver implementations;
- native core crates do not depend on JNI;
- Kotlin code does not reimplement native physical equations except in tests;
- public native ABI changes require an ABI version bump and migration note;
- generated resources do not encode scientific constants without provenance.

## Package Placement Rules

Use feature-first packages once real features exist. Until then, prefer neutral
infrastructure packages over speculative hierarchies.

Suggested Kotlin package roles:

- `common.rust`: native loading and JNI-facing wrappers;
- `common.platform`: game/platform abstractions shared by features;
- `common.registry`: Minecraft registry integration;
- `common.content`: block, item, fluid, and component definitions;
- `common.network`: packet and sync contracts;
- `common.ui`: screens, menus, overlays;
- `common.util`: small utilities only when no stronger owner exists.

Rust packages should follow domain boundaries, not JNI symbols. For example:

- `chemistry::species`;
- `chemistry::reaction`;
- `thermodynamics::phase`;
- `transport::diffusion`;
- `numerics::solver`;
- `state::snapshot`;
- `data::schema`;
- `data::serialization`;
- `ffi::jni`.

## Gradle Convention Responsibilities

`build-scripts` owns repeated build behavior:

- Kotlin/JVM compiler settings;
- NeoForge setup and run configuration;
- metadata generation from `config/mod.properties`;
- Rust JNI build and test tasks;
- architecture boundary checks;
- future generated-resource pipelines.

Build scripts should express version and loader context centrally. Individual
feature modules should not carry scattered literals for Minecraft versions,
loader versions, mod IDs, or Rust targets.

## Open Decisions

- Whether the Rust core should live in one workspace crate first or be split
  immediately into `core`, `data`, `ffi`, and `bench` crates.
- Which binary format should be used for native snapshots and chemical data
  packages.
- Whether JVM-to-native calls should use JNI arrays, direct byte buffers, Panama
  later, or a mixed strategy.
- Which parts of Create integration require per-tick coupling and which can be
  event-driven.
