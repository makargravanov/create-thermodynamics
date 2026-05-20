# Package Guidelines

Use feature-first packages only when real features exist.

Neutral shared packages in the common module:

- `content`
- `network`
- `platform`
- `registry`
- `ui`
- `util`

Do not put content or gameplay logic in the NeoForge module. Keep the loader module limited to entrypoints, event bus wiring, registration bootstrap, config registration, client bootstrap, and tiny loader compatibility shims.

The `core` module must not import Minecraft or loader packages.
