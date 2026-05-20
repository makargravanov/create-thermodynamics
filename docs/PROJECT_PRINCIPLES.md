# Project Principles

Create Thermodynamics is a hardcore Create addon whose main product is not a
set of recipes. The main product is a serious, high-performance physico-chemical
simulation engine that happens to be integrated into Minecraft.

The game integration must serve the model. The model must not be weakened just
because the first runtime is a game.

## Engineering Values

### Native-First Model

Rust is the primary implementation language for chemistry, thermodynamics,
transport, state storage, and numerical solving. Kotlin exists to integrate the
native engine with Minecraft, Create, Gradle, registries, UI, and platform
lifecycle.

This is an architectural rule, not a temporary optimization. If a behavior is
part of the scientific model, it belongs in Rust unless there is a deliberate
and documented exception.

### Scientific Seriousness

The project should start from engineering and scientific expectations:

- explicit quantities and units;
- physically meaningful state variables;
- traceable constants and empirical correlations;
- validation against reference cases;
- clear assumptions and known limits;
- numerical methods chosen for stability and performance, not convenience.

Minecraft balance is allowed at the adapter layer. It must not corrupt the
scientific model.

### Performance as a Design Constraint

Performance is a first-order requirement. The code should be designed around:

- predictable memory layout;
- cache-local state traversal;
- bounded allocation;
- batched JVM/native calls;
- deterministic stepping;
- profiling before broad optimization claims;
- algorithms that scale to large factories, networks, and worlds.

Slow but simple implementations are acceptable only as clearly marked prototype
steps. They must not become the public architecture.

### Programmable Clarity

The codebase should be readable by strong software engineers who are not domain
experts. Scientific code must not become a wall of unnamed symbols.

Mathematical code should expose:

- named intermediate values;
- documented equations;
- references or derivations for non-obvious formulas;
- explicit units in names or types;
- tests that tie formulas to known physical cases.

The goal is not to make the science fake-simple. The goal is to make real
complexity inspectable.

### Runtime Boundaries

The native core must be independently useful. It should be possible to reuse it
for education, engineering tools, offline simulation, or future CAD-style
integration without dragging Minecraft classes into the model.

The Minecraft mod is one runtime. It is not the whole architecture.

## Anti-Goals

The project explicitly rejects:

- recipe-only chemistry disguised as simulation;
- "game balance" constants embedded inside physical equations;
- ad hoc unit conversion through naked `Double` values;
- Kotlin-side reimplementation of native solvers;
- per-tick JNI chatter for each block entity when batching is possible;
- formulas written as unreadable expressions with no naming or explanation;
- hidden mutable global state that makes worlds non-deterministic;
- save formats that cannot be migrated;
- loader modules that accumulate gameplay logic;
- architecture that depends on one Minecraft version forever.

## Quality Bar

The code should be good enough that publishing it teaches useful software
engineering habits. A reader should be able to distinguish domain complexity
from accidental complexity.

Bad scientific code often hides weak architecture behind notation. This project
should do the opposite: make equations, assumptions, data, and performance
constraints explicit enough that both programmers and domain experts can review
them.

## Practical Rule

When there is tension between scientific fidelity, performance, and gameplay,
do not hide the tradeoff in code. Put the strict model in Rust, put adaptation
policy at the Minecraft boundary, and document the decision.
