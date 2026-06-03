# Condition — условия протекания реакции

Исходный код: `core/condition.rs`

## Назначение

Описывает набор физико-химических условий, при которых реакция разрешена или усилена. Каждое условие независимо: нарушение любого блокирует реакцию и фиксирует причину. Выполненные условия перемножают свои `rate_multiplier`.

## Ключевые типы

`AcidityCondition` — перечисление `Acidic | Neutral | Basic`. Пороги pH: кислая < 6, нейтральная 6–8, щелочная > 8.

`AtmosphereCondition` — `Air | Inert`. `Inert` требует отсутствия кислорода выше трейс-уровня.

`ReactionCondition` — одно условие с набором опциональных ограничений:

| Поле | Тип | Смысл |
|---|---|---|
| `phases` | `Vec<MixturePhase>` | Хотя бы одна фаза должна содержать вещество выше трейс-уровня |
| `acidity` | `Option<AcidityCondition>` | Требование кислотности/щёлочности |
| `solvent` | `Option<SubstanceId>` | Растворитель должен присутствовать |
| `min_temperature_kelvin` | `Option<f64>` | Минимальная температура |
| `max_temperature_kelvin` | `Option<f64>` | Максимальная температура |
| `min_water_activity` | `Option<f64>` | Минимальная активность воды |
| `max_water_activity` | `Option<f64>` | Максимальная активность воды |
| `max_oxygen_activity` | `Option<f64>` | Максимальная активность кислорода |
| `gas_pressure_atm` | `Option<f64>` | Минимальное давление газа (атм) |
| `atmosphere` | `Option<AtmosphereCondition>` | Требование атмосферы |
| `rate_multiplier` | `f64` | Множитель скорости если условие выполнено, по умолч. 1.0 |
| `reason` | `String` | Обязательное описание — попадает в `blocked_reasons` |

`ReactionConditionEvaluation` — результат оценки: `allowed: bool`, суммарный `multiplier`, `blocked_reasons: Vec<String>`.

## Публичные входы

`ReactionCondition::new(reason)` → fluent builder:

```
.in_phase(MixturePhase)
.acidity(AcidityCondition)
.solvent(id)
.min_temperature_kelvin(T)
.max_temperature_kelvin(T)
.min_water_activity(a) / .max_water_activity(a)
.max_oxygen_activity(a)
.gas_pressure_atm(p)
.atmosphere(AtmosphereCondition)
.rate_multiplier(k)
```

`evaluate_reaction_conditions(registry, mixture, conditions)` → `ChemistryResult<ReactionConditionEvaluation>`.

Вызывается из [[core-simulation|Simulation]] внутри расчёта скорости. Если `!evaluation.allowed` → rate = 0; иначе `rate *= evaluation.multiplier`.

## Поток данных / Алгоритм

```mermaid
flowchart TD
    I[conditions: &[ReactionCondition]] --> L{для каждого}
    L --> P1{phases пусты или есть фаза выше трейс?}
    P1 -->|нет| BLK[blocked = true]
    L --> P2{acidity?}
    P2 -->|pH не совпадает| BLK
    L --> P3{solvent present?}
    P3 -->|ниже трейс| BLK
    L --> P4{temperature in range?}
    P4 -->|нет| BLK
    L --> P5{water_activity in range?}
    P5 -->|нет| BLK
    L --> P6{oxygen_activity <= max?}
    P6 -->|нет| BLK
    L --> P7{gas_pressure >= min?}
    P7 -->|нет| BLK
    BLK --> R1[blocked_reasons.push reason]
    L -->|не блокировано| R2[multiplier *= rate_multiplier]
    R1 & R2 --> E[ReactionConditionEvaluation]
```

Каждое условие проверяется независимо. Итог: `allowed = blocked_reasons.is_empty()`.

pH для `AcidityCondition` вычисляется через `mixture.ph(registry)` — может вернуть `None` при отсутствии кислотно-основной системы, что блокирует условие.

Активность кислорода проверяется только при наличии `max_oxygen_activity` или `atmosphere == Inert`; иначе пропускается (без поиска вещества).

## Инварианты и ошибки

- `reason` не может быть пустым (пробелы) — `InvalidReaction`.
- `min_temperature`, `max_temperature`, `gas_pressure_atm` — конечны и > 0.
- `min/max_water_activity`, `max_oxygen_activity` — конечны, в диапазоне [0.0, 1.0].
- `min_temperature ≤ max_temperature` (если оба заданы).
- `min_water_activity ≤ max_water_activity` (если оба заданы).
- `rate_multiplier ≥ 0` и конечен.

## Связи

- [[core-reaction|Reaction]] — `Reaction.conditions: Vec<ReactionCondition>`.
- [[core-simulation|Simulation]] — вызывает `evaluate_reaction_conditions` в каждом расчёте rate.
- [[core-mixture|Mixture]] — читает `temperature_kelvin()`, `ph()`, `activity_of()`, `gas_pressure_pascal()`.
