# Rust Native Strategy

Rust is the main implementation language for the model because the project
needs predictable performance, explicit ownership, strong typing, safe
concurrency, and low-level control over memory layout.

Kotlin should call into Rust. Kotlin should not become the home of the
simulation.

## Crate Responsibilities

The current native crate is `native/thermodynamics-jni`. It is a starting JNI
bridge, not the final shape of the native model.

Expected direction:

```text
native/
  thermodynamics-core/       pure domain model, state, and solvers
  thermodynamics-data/       loading, validation, serialization, provenance
  thermodynamics-jni/        JNI facade and JVM interop
  thermodynamics-bench/      benchmarks, profiling scenarios, fixtures
```

The JNI crate should depend on the core and data crates. Core and data crates
should not depend on JNI, Minecraft, or JVM concepts.

## FFI Principles

The FFI boundary should be narrow, versioned, and batch-oriented.

Prefer:

- opaque native handles for long-lived worlds and systems;
- direct buffers or packed arrays for large transfers;
- explicit create/destroy lifecycle;
- ABI version checks on startup;
- structured error codes plus diagnostic retrieval;
- batched events and deltas instead of per-block calls.

Avoid:

- exposing deep Rust object graphs to JVM;
- returning large object trees through JNI;
- calling native code once per tiny field;
- making JVM allocation part of hot solver loops;
- letting JNI naming leak into domain module names.

## State And Data Ownership

Native code owns simulation state and chemical data serialization. Kotlin may
hold native handles, stable IDs, opaque snapshot bytes, and Minecraft-facing
cached views, but it should not be the authoritative owner, parser, validator,
or migrator of chemical state.

Long-lived native state should have:

- explicit allocation and release;
- generation or version checks for handles;
- deterministic snapshot export;
- clear thread-affinity rules;
- validation hooks after deserialization;
- panic containment at the FFI boundary.

## Memory Layout

Hot state should be designed for traversal, not object comfort.

Use data layouts that match access patterns:

- dense arrays for species amounts and intensive state;
- stable integer IDs for indexed lookup;
- structure-of-arrays where vectorized or cache-local traversal benefits;
- arena or slab allocation for stable handles;
- sparse structures only when sparsity is real and measured;
- compact graph representations for reactor and network topology.

Readable code still matters. Low-level layout should be wrapped in named domain
types and documented invariants.

## Deterministic Execution

Native execution should be deterministic for a fixed input set, data version,
configuration, and supported platform class.

Design requirements:

- stable ordering for IDs, reactions, phases, and topology traversal;
- no hidden dependence on hash map iteration order;
- controlled use of parallelism;
- explicit floating-point tolerance policies;
- deterministic serialization;
- regression fixtures for known scenarios.

## Serialization

Native snapshots and chemical data packages are part of the engine contract.
They should be versioned independently from Minecraft's save representation.
The Rust data crate should own schemas, validation, migrations, provenance, and
serialization for chemistry data.

Snapshot requirements:

- schema version;
- native ABI version;
- data registry version;
- endian and format rules if binary;
- migration hooks;
- validation on load;
- diagnostic export for debugging.

Kotlin's save/load role is deliberately small: ask Rust to serialize, store the
opaque payload through Minecraft's save system, then pass that payload back to
Rust. Do not use JVM object serialization for authoritative simulation state or
chemical data.

## Error Handling

Inside Rust, prefer typed errors that preserve context. Across JNI, transfer
errors through stable codes and retrievable diagnostics.

At the boundary:

- panics must not unwind into JVM;
- invalid handles must be detected;
- ABI mismatch must fail startup clearly;
- non-convergence must be a solver result, not a crash;
- missing data must identify the missing entity and source.

## Versioning

Version all externally visible native contracts:

- JNI ABI;
- snapshot format;
- scientific data schema;
- scientific data serialization format;
- solver configuration schema;
- domain registry IDs.

Any breaking native contract change requires a migration note. Released worlds
must not depend on undocumented binary accidents.

## Testing

Native tests should cover:

- unit conversion;
- domain invariants;
- reaction balancing;
- conservation of mass, charge, and energy;
- solver reference cases;
- serialization round trips;
- chemical data schema and migration fixtures;
- ABI smoke tests;
- deterministic replay fixtures.

Kotlin tests should verify integration boundaries, loading, packaging, and
Minecraft-facing adaptation. They should not be the primary proof of scientific
correctness.

## Profiling

Performance claims require measurement.

Maintain benchmark scenarios for:

- small machine systems;
- large factory networks;
- high species-count mixtures;
- stiff reaction systems;
- phase equilibrium cases;
- serialization snapshots;
- JVM/native transfer overhead.

Track allocation, cache behavior, step time distribution, and scaling behavior.
Optimization should target measured bottlenecks, but data structures should be
chosen from the start with hot-path access patterns in mind.

## Open Decisions

- Exact native workspace split and crate naming.
- JNI transfer mechanism for large data: direct buffers, primitive arrays, or
  generated bindings.
- Binary snapshot and chemical data package formats.
- Minimum supported Rust toolchain policy.
- Platform matrix for packaged native libraries.
