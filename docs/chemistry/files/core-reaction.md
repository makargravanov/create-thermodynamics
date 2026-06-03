# Reaction — модель химической реакции

Исходный код: `core/reaction.rs`

## Назначение

Определяет структуру данных `Reaction` и её построитель `ReactionBuilder`. Это центральный объект химического ядра: именно реакция несёт стехиометрию, кинетику (Аррениус), термодинамику, условия и связи с внешним контекстом мода.

## Ключевые типы

`ReactionId` — newtype-обёртка над `String`. Не может быть пустой; создаётся через `ReactionId::new(...)` или `From<&str>`.

`StoichiometricTerm` — пара (`SubstanceId`, `coefficient: u32`). Coefficient > 0 — инвариант, проверяется при `validate_shape`.

`ProductDistributionVariant` / `ProductDistribution` — вероятностное распределение продуктов: список вариантов с `fraction: f64`, сумма долей должна быть ровно 1.0 (±1e-9).

`ExternalRequirement` — реагент/катализатор/продукт вне смеси (игровой контекст). Несёт `description`, `moles_per_reaction`, опционально `molar_mass_grams` и `charge` для проверки баланса. Если масса/заряд неизвестны — требуется `unchecked_mass_reason`.

`ReactionResult` — побочный игровой результат (например, нагрев блока): `description` + `moles_per_reaction`.

`ReactionPhaseAccess` — список фаз `MixturePhase`, из которых читается концентрация реагента. По умолчанию — `[Aqueous, Organic]`.

`Reaction` — основная структура. Ключевые поля:

| Поле | Тип | Роль |
|---|---|---|
| `reactants` / `products` | `Vec<StoichiometricTerm>` | Прямая стехиометрия |
| `product_distribution` | `Option<ProductDistribution>` | Вероятностный выход — взаимоисключает `products` |
| `channels` | `Vec<ReactionChannel>` | Множественные пути — взаимоисключают `products` и `reverse_reaction_id` |
| `orders` | `BTreeMap<SubstanceId, u32>` | Порядки кинетики; может включать катализаторы |
| `pre_exponential_factor` | `f64` | A в формуле Аррениуса, по умолч. 10 000 |
| `activation_energy_kj_per_mol` | `f64` | Ea, по умолч. 25 кДж/моль |
| `enthalpy_change_kj_per_mol` | `f64` | ΔH; применяется при `mixture.heat()` |
| `thermodynamics` | `Option<ReactionThermodynamics>` | Gibbs / Keq для термодинамического торможения |
| `reverse_reaction_id` | `Option<ReactionId>` | Ссылка на обратную реакцию (пара) |
| `conditions` | `Vec<ReactionCondition>` | Условия среды — см. [[core-condition\|Condition]] |
| `surface_requirements` | `Vec<SurfaceRequirement>` | Катализаторные поверхности — см. [[core-catalysis\|Catalysis]] |
| `selectivity_profile` | `Option<SelectivityProfile>` | Профиль селективности — см. [[selectivity-engine\|SelectivityEngine]] |
| `phase_access` | `BTreeMap<SubstanceId, ReactionPhaseAccess>` | Переопределение фаз для реагентов |
| `product_phases` | `BTreeMap<SubstanceId, MixturePhase>` | Целевая фаза для продуктов |

## Публичные входы

`Reaction::builder(id)` → `ReactionBuilder` — единственный способ создания. Методы построителя:

- `.reactant(id, coeff, order)` — добавляет реагент и его кинетический порядок.
- `.product(id, coeff)` / `.product_distribution_variant(frac, [(id, coeff)])` — прямые или распределённые продукты.
- `.channel(ReactionChannel)` — добавляет канал пути реакции.
- `.condition(ReactionCondition)` — условие среды.
- `.external_reactant(desc, moles)` / `.chemical_external_reactant(desc, moles, mass, charge)` — внешние реагенты без/с балансом массы.
- `.electron_reactant(n)` / `.electron_product(n)` — электроны для полуреакций (заряд −1, масса 0).
- `.surface_requirement(id, sites)` / `.surface_adsorption(...)` / `.surface_desorption(...)` — каталитические поверхности.
- `.gibbs_free_energy_change_kj_per_mol(ΔG)` / `.equilibrium_constant(Keq)` — термодинамика.
- `.reverse_reaction_id(id)` — связывает с обратной реакцией.
- `.build()` → `Reaction` (без валидации).

`rate_constant_per_second(T)` — Аррениус: `A · exp(−Ea·1000 / (R·T))`.

`has_external_context()` — `true`, если реакция требует UV, внешних реагентов/катализаторов/поверхностей.

`requires_context_to_proceed()` — строже: только UV + реагенты + катализаторы + поверхности (без продуктов и результатов).

## Поток данных / Алгоритм

```mermaid
flowchart LR
    B[ReactionBuilder] -->|.build()| R[Reaction]
    R -->|validate_shape| V{OK?}
    V -- ошибка --> E[ChemistryError]
    V -- OK --> REG[ChemistryRegistry.build]
    REG -->|validate_reactions| RC{mass/charge/redox/reverse}
    RC --> SIM[Simulation]
```

`validate_shape` — структурная проверка формы (без реестра): непустые реагенты/продукты, суммы долей, коэффициенты > 0, поля Аррениуса конечны и ≥ 0, redox-аннотация электробалансирована.

Полная содержательная валидация (масса, заряд, ссылки) — в `ChemistryRegistry.validate_reactions`.

## Инварианты и ошибки

- `ReactionId` не может быть пустым.
- `coefficient > 0` для всех `StoichiometricTerm`.
- Если есть `channels` → нет `products`/`product_distribution`/`reverse_reaction_id`.
- Если есть `product_distribution` → нет `products`/`reverse_reaction_id`; сумма `fraction` = 1.0 ± 1e-9.
- `pre_exponential_factor > 0`, `activation_energy_kj_per_mol ≥ 0`, `enthalpy_change_kj_per_mol` конечен.
- `thermodynamics` допустим только при наличии прямого набора продуктов (не каналы, не распределение).
- `redox.transferred_electrons > 0`; `allow_charge_imbalance` несовместим с redox.
- `ExternalRequirement`: либо (mass, charge) оба заданы и корректны, либо `unchecked_mass_reason` не None.

## Связи

- [[core-condition|Condition]] — `Reaction.conditions` содержит `Vec<ReactionCondition>`.
- [[core-kinetics|Kinetics]] — `ReactionChannel` используется в `channels`.
- [[core-catalysis|Catalysis]] — `SurfaceRequirement`, `SurfaceStep` встроены в `Reaction`.
- [[core-registry|Registry]] — принимает `Reaction` через `ChemistryRegistryBuilder::reaction()`.
- [[core-simulation|Simulation]] — использует `rate_constant_per_second`, `enthalpy_change_kj_per_mol`.
- [[core-redox|Redox]] — `RedoxAnnotation` в `reaction.redox`.
- [[selectivity-engine|SelectivityEngine]] — `selectivity_profile` применяется в симуляции.
