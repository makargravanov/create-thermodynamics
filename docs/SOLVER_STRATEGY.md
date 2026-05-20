# Solver Strategy

The project should treat solving as a family of cooperating numerical systems,
not as one universal function. Chemistry, heat transfer, transport, phase
behavior, electrical effects, and metallurgy have different mathematical shapes
and failure modes.

The first design rule: separate the physical model, numerical method, and game
runtime policy.

## Solver Layers

### Physical Model

The physical model defines what is being solved:

- conserved quantities;
- phases and interfaces;
- species and reactions;
- energy, charge, and mass balance;
- equations of state;
- transport laws;
- kinetic and equilibrium constraints.

This layer should be as strict as practical. It is not where Minecraft balance
belongs.

### Numerical Method

The numerical method defines how the model is solved:

- nonlinear equation solving;
- ODE/DAE integration;
- constrained minimization;
- sparse linear solves;
- graph and network solving;
- operator splitting;
- event detection;
- convergence and stability control.

Methods must expose diagnostics. A solver that silently fails and returns
plausible-looking state is worse than a solver that refuses a bad step.

### Runtime Policy

Runtime policy decides how the solver is scheduled inside Minecraft:

- tick budgets;
- batching;
- level and chunk activation;
- degraded update frequency for distant systems;
- player-visible error handling;
- synchronization to clients;
- save and rollback points.

Runtime policy may approximate scheduling. It must not rewrite the physical
model to make scheduling easier.

## Required Solver Families

### Thermodynamic Equilibrium

Needed for phase splits, reaction equilibrium, vapor-liquid equilibrium,
aqueous speciation, metallurgical slag/alloy systems, and later nuclear or
plasma extensions.

Expected capabilities:

- mass and charge balance;
- phase stability checks;
- activity/fugacity models;
- temperature and pressure dependence;
- robust failure diagnostics.

### Reaction Kinetics

Needed for irreversible and finite-rate chemistry, combustion, catalysis,
polymerization, corrosion, ore processing, and industrial reactors.

Expected capabilities:

- stoichiometric mechanisms;
- rate laws;
- Arrhenius-style temperature dependence where applicable;
- catalyst and surface terms;
- stiffness-aware integration;
- event handling for depletion and phase changes.

### Heat Transfer

Needed across almost every feature: machines, pipes, tanks, furnaces, fluids,
contraptions, and environmental exchange.

Expected capabilities:

- conduction;
- convection abstractions;
- radiation where relevant;
- heat capacity and latent heat;
- coupling to phase changes and reactions;
- stable integration under large temperature gradients.

### Mass Transport

Needed for diffusion, flow, separation, membranes, distillation, extraction,
electrolysis, and gas handling.

Expected capabilities:

- advection and diffusion;
- interphase transfer;
- pressure-driven flow;
- concentration-driven transport;
- network coupling;
- boundary conditions from Minecraft machines and Create contraptions.

### Electrochemical Solvers

Needed for electrolysis, batteries, corrosion, plating, refining, and future
advanced machinery.

Expected capabilities:

- charge balance;
- electrode reactions;
- potentials and overpotentials;
- electrolyte composition;
- coupling to heat and mass transport.

### Metallurgy and Materials

Needed for ores, smelting, alloying, slag behavior, heat treatment, and
materials engineering.

Expected capabilities:

- multiphase high-temperature systems;
- phase diagrams and empirical models;
- solid solution and precipitate handling;
- impurity and inclusion tracking;
- mechanical/material property derivation where data permits.

### Network and Reactor Coupling

Needed to connect local reactors into factories. A pipe, tank, boiler,
distillation column, and Create contraption should not be solved as unrelated
toys when they exchange conserved quantities.

Expected capabilities:

- graph-based ownership of connected systems;
- batched stepping;
- conservation checks across boundaries;
- partitioning for cache locality and parallel execution;
- deterministic merge/split behavior when machines connect or disconnect.

## Coupling Strategy

The initial architecture should allow operator splitting: solve subsystems in a
controlled order with explicit exchange variables. Fully coupled solves can be
added for domains where splitting produces unacceptable error or instability.

Recommended default step shape:

1. Gather batched runtime events from Minecraft.
2. Update topology and boundary conditions.
3. Resolve fast algebraic constraints.
4. Advance transport and heat transfer.
5. Advance kinetic reactions.
6. Resolve phase and equilibrium constraints when required.
7. Validate conservation and bounds.
8. Emit compact deltas and diagnostics.

## Determinism

The simulation should be deterministic for the same inputs, platform class, data
versions, and configuration. Determinism matters for saves, debugging,
multiplayer, testing, and future replay tooling.

Rules:

- stable iteration order for domain IDs and graph traversal;
- explicit random seeds when stochastic behavior exists;
- no hidden dependence on thread scheduling;
- controlled floating-point assumptions;
- regression tests for reference scenarios.

## Failure Handling

Solver failure is a valid outcome and must be represented explicitly.

Failures should include:

- non-convergence;
- invalid physical input;
- missing data;
- impossible phase state;
- conservation violation;
- step rejected by stability criteria;
- native ABI or serialization mismatch.

Minecraft-facing code may translate failures into machine shutdown, warning UI,
log diagnostics, or safe state. It must not silently continue with corrupted
state.

## Game Runtime Boundary

The game runtime can decide when and how often a system is stepped. It can
choose display granularity, UI simplifications, and player-facing affordances.

The game runtime must not:

- mutate native state outside declared APIs;
- invent physical constants;
- bypass conservation checks for convenience;
- run independent Kotlin chemistry that disagrees with Rust;
- encode progression balance inside solver equations.

## Open Decisions

- Which solver family should be implemented first as the reference vertical
  slice.
- Whether native parallelism should be internal to Rust from the start or
  introduced after deterministic single-thread baselines.
- Which tolerance policies are global and which are process-specific.
- How aggressively Minecraft tick budgets may defer simulation without
  violating player expectations.
