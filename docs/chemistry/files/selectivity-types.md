# Типы системы селективности

Исходный код: `kinetics/selectivity/types.rs`, `kinetics/selectivity/mod.rs`

## Назначение

Базовые структуры данных для оценки кинетической селективности органических реакций. Описывают реакционный центр (стерика, электроника, степень замещения) и условия среды, в которых протекает реакция. Результат оценки — числовые поправки к скорости и энергии активации.

Публичный API модуля: `pub use engine::SelectivityEngine; pub use engine::SiteDescriptorBuilder; pub use types::*;`

## Ключевые типы

### SubstitutionDegree

Степень замещения у реакционного центра. Варианты: `Primary` (≤1 C), `Secondary` (2 C), `Tertiary` (3 C), `Benzylic` (Ph-CH₂-), `Allylic` (C=C-CH₂-).

Метод `base_steric_score() -> f64` возвращает базовый стерический балл (0.0–1.0):

| Вариант   | Балл |
|-----------|------|
| Primary   | 0.0  |
| Secondary | 0.3  |
| Tertiary  | 0.6  |
| Benzylic  | 0.2  |
| Allylic   | 0.25 |

### ElectronicEnvironment

Электронное окружение центра. Поля: `electron_donating_groups` (EDG, +I/+M эффекты), `electron_withdrawing_groups` (EWG, -I/-M), `resonance_stabilization`, `aromatic`.

`net_effect() -> f64`: каждый EDG даёт +0.1, каждый EWG −0.1, резонанс +0.15.

### SiteDescriptor

Полный дескриптор реакционного центра. Агрегирует `site_kind` (`ReactiveSiteKind`), `degree`, `electronics`, `steric_score` и флаг `has_beta_hydrogen`.

Конструктор `new(...)` вычисляет итоговый стерический балл:

```
steric_score = min(base_steric + bulky_substituents × 0.1, 1.0)
```

`steric_accessibility() = max(1.0 − steric_score, 0.1)` — нижняя граница 0.1 гарантирует, что полная блокировка не обнуляет реакцию при `never_suppress`.

### SelectivityContext

Условия оценки: температура, pH, тип растворителя, активности ключевых веществ (вода, O₂, H₂, фторид), redox-баланс, наличие металлов и поверхностей.

`from_mixture(registry, mixture, reaction_context)` — единственный «серьёзный» конструктор. Читает pH, парциальные давления газов, окислительно-восстановительные роли всех компонентов, концентрацию свободных и скомплексованных металлов, мощность УФ и суммарную светимость из `ReactionContext`. Если классификацию растворителя не удаётся определить однозначно, используется `SolventType::Neutral`.

Вспомогательные предикаты:

| Метод | Условие |
|---|---|
| `is_acidic()` | pH < 6 |
| `is_basic()` | pH > 8 |
| `is_high_temperature()` | T > 353 K (80 °C) |
| `is_very_high_temperature()` | T > 423 K (150 °C) |
| `is_water_rich()` | water_activity ≥ 0.35 |
| `is_water_poor()` | water_activity ≤ 0.20 |
| `is_oxidizing()` | oxidizing > reducing + trace |
| `is_reducing()` | reducing > oxidizing + trace |

### SolventType

`Protic` (вода, спирты — благоприятствует SN1/E1), `AproticPolar` (DMSO, DMF — SN2), `NonPolar`, `Basic`, `Acidic`, `Neutral`.

Тип определяется из смеси: если водная фаза доминирует → `Protic`; если органическая с апротонным растворителем → `AproticPolar`; иначе `NonPolar`. Кислотность/основность pH переопределяет этот классификатор.

### ReactionType

Перечисление механизмов реакций: `SN1`, `SN2`, `E1`, `E2`, `FischerEsterification`, `CarbonylAddition`, `CarbonylReduction`, олефинирование (`Wittig`, `HWE`, `Julia`), реакции альфа-углерода (`AlphaHalogenation`, `AldolAddition`, `MichaelAddition`, …), защитные группы (`SilylEther*`, `Acetal*`, `Carbamate*`, `Ester*`).

### NucleophileStrength

`VeryStrong` (Grignard, RLi), `Strong` (NaBH₄), `Moderate` (OH⁻, CN⁻, амины), `Weak` (вода, спирты). Определяет чувствительность к внешней среде: только `VeryStrong` нуклеофилы гасятся водой и металлами в `apply_inorganic_environment_to_carbonyl_score`.

### SelectivityProfile

Описание профиля конкретной реакции для передачи в движок. Содержит `mechanism`, `primary_site`, опциональный `secondary_site` (требуется для `FischerEsterification`, `Wittig`/`HWE`/`Julia`), `nucleophile_strength` и `suppression_policy`.

Billder-методы: `with_secondary_site`, `with_nucleophile_strength`, `never_suppress` (переключает политику на `NeverSuppress`).

### SelectivityRuntimeEffect

Итог оценки профиля: `rate_multiplier`, `activation_delta_kj_per_mol`, `pre_exp_multiplier`, флаг `suppressed`, текстовая причина. Потребляется генераторами органических реакций.

### ReactivityScore

Числовой результат: `value` (0.0…∞, >1.0 — быстрее опорной), `activation_delta` (кДж/моль, < 0 ускоряет), `pre_exp_multiplier`, список конкурирующих механизмов.

Преобразование между `value` и `activation_delta` через уравнение Аррениуса при T = 298 К:

```
ΔEa = −RT·ln(k/k₀)   // value → delta
k/k₀ = exp(−ΔEa/RT)  // delta → value
```

### SelectivityResult

Агрегирует `primary: ReactivityScore`, словарь всех механизмов и рекомендацию `SelectivityRecommendation`:

| Отношение primary/max_competitor | Рекомендация |
|---|---|
| ≥ 10 | `Exclusive` |
| ≥ 3 | `Preferred` |
| ≥ 0.5 | `Mixed` |
| ≥ 0.1 | `Suppressed` |
| < 0.1 | `None` |

## Инварианты и ошибки

- `steric_accessibility()` ограничена снизу 0.1 — движок не производит нулевых скоростей при `NeverSuppress`.
- `SelectivityContext::from_mixture` возвращает `ChemistryResult`; если вещество не найдено в реестре, его концентрация трактуется как 0.
- Классификация `palladium_available` проверяет как `external_catalysts` (подстрока "palladium"), так и свободные сайты поверхности.

## Связи

- [[selectivity-engine|Движок селективности]] — потребляет все типы этого модуля.
- [[selectivity-reactions|Реакции]] — строят `SelectivityProfile` и вызывают функции оценки.
- [[core-simulation|ReactionContext]] — источник данных для `from_mixture`.
- [[core-mixture|Mixture]] — pH, активности, фазовые концентрации.
- [[molecule-reactive-site|ReactiveSiteKind]] — используется в `SiteDescriptor.site_kind`.
