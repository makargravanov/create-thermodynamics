# Domain Model

This document names the first-class concepts the simulation should use. It is a
conceptual contract, not a promise that every type already exists.

The model should separate scientific entities from Minecraft representations.
A Minecraft fluid, item, block, or Create contraption may refer to domain
objects, but it is not itself the chemistry model.

## Foundational Concepts

### Element

An element is a chemical element identified by atomic number. Element metadata
should include symbol, standard name, atomic number, and references for physical
data.

Element identity must not depend on localized names or Minecraft registry IDs.

### Nuclide and Isotope

A nuclide is a specific nucleus identified by proton count, neutron count, and
possibly nuclear excitation state. An isotope is a group of nuclides sharing the
same element.

The domain should reserve space for nuclear processes even if they are not part
of the first implementation.

### Species

A species is a chemically distinct entity: atom, molecule, ion, radical,
electron, defect, complex, surface site, or other identifiable participant in
chemical processes.

Species data should capture:

- composition;
- charge;
- phase applicability;
- thermodynamic reference data;
- transport data when available;
- identifiers for mechanisms and reactions.

### Phase

A phase describes a physically distinct form of matter with shared intensive
properties. Examples: gas, liquid, solid, aqueous phase, plasma, surface phase,
crystal phase, slag phase, alloy phase.

Do not reduce phase to a Minecraft fluid or block state. Minecraft states are
rendering and interaction artifacts; phase is a physical concept.

### Substance

A substance is a material identity under a specific domain description. It may
be a pure species, a defined compound, a mineral, a polymer family, an alloy
grade, or a named industrial material.

Substance is useful when the player or data pack needs a named material while
the solver still uses species, phases, and mixtures internally.

### Material

A material is an engineering-facing aggregate with properties relevant to use:
thermal conductivity, heat capacity, density, strength, corrosion behavior,
phase diagram references, and processing constraints.

Material may map to one or more substances and phases.

### Mixture

A mixture is a set of species amounts within one phase or a modeled effective
phase. It owns composition and intensive state such as temperature and pressure.

Mixture composition must use explicit amount units. Avoid percentage-only state
as the authoritative representation.

### Heterogeneous System

A heterogeneous system is a collection of phases and interfaces that may
exchange heat, mass, charge, and momentum. Examples: boiling solution, ore plus
slag plus gas, catalyst bed, electrolytic cell, condenser, distillation column
segment.

This is likely the central runtime object for many machines and networks.

## Process Concepts

### Reaction

A reaction is a transformation relation between species with stoichiometry,
directionality, thermodynamic constraints, kinetic law candidates, and catalyst
or surface dependencies.

Reaction data must not imply that all reactions are always active. Solver
context decides which reactions participate in a process step.

### Process

A process is a physical or physico-chemical operation over a system. Examples:

- heating and cooling;
- compression and expansion;
- distillation;
- extraction;
- electrolysis;
- precipitation;
- crystallization;
- roasting;
- smelting;
- alloying;
- corrosion;
- gas absorption;
- diffusion through a membrane.

Processes combine domain equations, data, and solver policies.

### Reactor

A reactor is a modeled runtime boundary where processes occur. It may represent
a machine, pipe segment, tank, heat exchanger, furnace, electrolyzer, Create
contraption volume, or an abstract control volume.

Reactors own geometry and coupling rules. They should not own global scientific
truths.

### Solver

A solver advances or resolves state under a specific mathematical contract:
equilibrium, kinetics, transport, heat transfer, phase split, network flow, or
coupled multiphysics.

Solvers should be composable and explicit about required inputs, convergence
criteria, failure modes, and approximation boundaries.

### State

State is the complete set of values required to continue a simulation
deterministically. It includes composition, phase state, temperature, pressure,
energy, charge, geometry-derived parameters, and solver-owned internal
variables where necessary.

State must be serializable, versioned, and suitable for native memory ownership.

## Measurement Concepts

### Quantity

A quantity is a numeric value with a physical dimension. The model should avoid
raw scalar values at API boundaries where units matter.

Internal hot loops may use optimized scalar storage, but the type and field
names must preserve dimensional meaning.

### Unit

A unit is a representation of a physical dimension and scale. The authoritative
model should use coherent internal units, preferably SI-derived, and convert at
the edges.

Unit conversion must be centralized and tested. Do not scatter conversion
factors through formulas.

## Identity Rules

Domain IDs must be stable, deterministic, and independent of display names.
Minecraft registry IDs may reference domain IDs, but domain IDs should survive
renames, localization changes, and loader/version changes.

## Open Decisions

- Exact split between `species`, `substance`, and `material` in data files.
- Whether isotope-aware composition is mandatory from day one or introduced as
  an extension of element composition.
- How to represent charged surfaces and catalyst active sites.
- How much dimensional typing should exist in hot Rust paths versus validation
  and API layers.
