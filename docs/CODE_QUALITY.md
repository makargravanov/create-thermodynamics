# Code Quality

The project must be readable as software, not only as applied chemistry or
physics. Scientific complexity is expected. Accidental complexity is not.

## Naming

Use domain names when the domain concept is real:

- `temperature_kelvin`, not `t`;
- `pressure_pascal`, not `p` outside a very small formula scope;
- `species_amount_mol`, not `n` in public or persistent state;
- `reaction_extent`, not `x`;
- `heat_capacity_joule_per_kelvin`, not `cp` unless the type and scope make it
  unambiguous.

Short mathematical symbols may appear in tiny, local scopes only when they
match a documented equation immediately above the code. Do not let notation
escape into APIs, persisted state, or large functions.

## Mathematical Code

Non-trivial formulas require structure:

1. Name intermediate values.
2. State the equation or reference.
3. Make units visible through types, field names, or comments.
4. Add a reference test or dimensional sanity test.
5. Keep empirical constants traceable to data files or cited sources.

Avoid code shaped like:

```rust
let y = x0 * (3.1415926 * dt / (i * q * dvxy).sqrt());
```

Prefer code shaped like:

```rust
let diffusion_time_scale_seconds = time_step_seconds / cell_width_meters.powi(2);
let normalized_flux = PI * diffusion_time_scale_seconds / charge_transport_factor;
let concentration_delta = initial_concentration * normalized_flux.sqrt();
```

This example is illustrative, not a domain recommendation.

## Units and Quantities

Every physical value crossing an API boundary must make its unit obvious.

Accepted strategies:

- strong quantity types where practical;
- unit suffixes in field and parameter names;
- schema-level unit declarations for data files;
- coherent internal units in hot paths;
- centralized conversion functions.

Rejected strategies:

- naked conversion factors inside formulas;
- comments that contradict field names;
- UI units leaking into solver state;
- Minecraft-specific units as authoritative physical units.

## Types and Invariants

Represent invariants in types where doing so does not harm hot-path layout.

Examples:

- distinct IDs for species, phase, reaction, reactor, and material;
- non-zero or positive wrappers for values that require them;
- bounded fractions when stored as fractions;
- explicit signed quantities for charge and energy flow;
- separate absolute temperature from temperature difference where useful.

When optimized storage uses raw arrays, wrap access behind APIs that preserve
domain meaning.

## Documentation

Document:

- why a model or equation was selected;
- valid input ranges;
- assumptions and approximation boundaries;
- data provenance;
- numerical tolerances;
- failure behavior;
- performance-sensitive layout decisions.

Do not document obvious assignments. Spend comments on decisions that would be
expensive to rediscover.

## Tests

Tests should be close to the risk.

Required categories:

- domain invariant tests;
- unit conversion tests;
- reference physical cases;
- conservation checks;
- serialization round trips;
- deterministic replay tests;
- FFI smoke tests;
- architecture boundary checks.

A solver without reference cases is not finished. A formula without at least one
sanity test should be treated as provisional.

## Performance Code

Optimized code must still explain itself.

For performance-sensitive areas, record:

- access pattern;
- chosen memory layout;
- expected complexity;
- reason for avoiding a simpler representation;
- benchmark that protects the decision.

Micro-optimizations without a benchmark should not make code obscure.

## Error Messages

Errors should identify:

- operation;
- entity ID or source record;
- relevant units and values;
- violated invariant;
- whether the failure is recoverable.

Native errors crossing into Kotlin should remain structured long enough for
Minecraft-facing code to produce useful logs and diagnostics.

## Review Checklist

Before merging model or solver code, check:

- Are units explicit?
- Are formulas readable and named?
- Are constants traceable?
- Are invariants represented or validated?
- Are failure modes explicit?
- Does this belong in Rust rather than Kotlin?
- Does this allocation or JNI call happen in a hot path?
- Is there a test that would fail if the equation or unit were wrong?
