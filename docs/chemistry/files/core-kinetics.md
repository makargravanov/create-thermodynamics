# Каналы реакций и энергетическая модель (kinetics)

Исходный код: `core/kinetics.rs`

## Назначение

Описывает мультиканальные реакции: одни реагенты могут давать разные наборы продуктов в
зависимости от температуры, освещённости, состояния поверхности и профиля селективности.
Предоставляет `EnergyModel` для термодинамически равновесного распределения изомеров.

## Ключевые типы

| Тип | Роль |
|-----|------|
| `ReactionChannelId` | строковый идентификатор канала, не может быть пустым |
| `ReactionChannelMode` | `Kinetic / Thermodynamic / Mixed / Photochemical / Surface` |
| `LightBand` | `Ultraviolet / Visible / Infrared` — источник фотохимической активации |
| `ChannelConditionEffect` | модификатор веса канала: `Phase`, `Light`, `Surface` |
| `ReactionChannel` | канал: продукты, E_a, A, режим, список `ChannelConditionEffect`, опциональный `SelectivityProfile` |
| `IsomerEnergy` | относительная энергия Гиббса изомера (вещество, фаза, поверхность) |
| `TransitionStateEnergy` | энергия переходного состояния для конкретного канала |
| `EnergyModel` | контейнер `IsomerEnergy` + `TransitionStateEnergy`; вычисляет равновесное распределение |

## Публичные входы

- `ReactionChannel::rate_weight(registry, mixture, context)` — вес канала при текущих условиях
- `channel_product_distribution(registry, reaction_id, channels, mixture, context)` — нормированное распределение `ProductDistribution`
- `channel_rate_sum_per_second(...)` — суммарная скорость всех каналов реакции
- `EnergyModel::equilibrium_distribution(substances, T, phase)` — Больцман-распределение по изомерам
- `EnergyModel::isomer_energy / transition_state_energy` — поиск энергии с fallback по (фаза, поверхность)

## Поток данных / Алгоритм

### Вес канала (`channel_rate_weight`)

```
1. Если есть SelectivityProfile:
     эффект = SelectivityEngine::evaluate_profile(profile, context)
     если suppressed → 0
     A'  = A · effect.pre_exp_multiplier
     Ea' = max(0, Ea + effect.activation_delta)
     mul = effect.rate_multiplier
2. weight = A' · exp(−Ea'·1000 / (R·T))
3. weight *= mul
4. Для каждого ChannelConditionEffect:
     Phase   → 0 если фаза пуста, иначе multiplier
     Light   → 0 если мощность < minimum_power, иначе power · multiplier
     Surface → 0 если свободных сайтов ≤ TRACE, иначе multiplier
     первый нулевой множитель → ранний выход
```

### Распределение продуктов (`channel_product_distribution`)

```
weights[i] = channel[i].rate_weight(...)
fractions  = normalize_weights(weights)   # weight_i / Σ weights
→ ProductDistribution { variants: [{ fraction, products }] }
```

Варианты с нулевой долей отфильтровываются.

### Равновесное распределение изомеров (`EnergyModel::equilibrium_distribution`)

```
w_i = exp(−G_i·1000 / (R·T))   # G_i — relative_gibbs_kj_per_mol
fraction_i = w_i / Σ w_j
```

`EnergyModel` поддерживает специфику (фаза, поверхность): поиск ведётся от наиболее
специфичного ключа до самого общего `(substance, None, None)`.

## Инварианты и ошибки

| Условие | Ошибка |
|---------|--------|
| Пустой список продуктов в канале | `InvalidReaction` |
| `activation_gibbs < 0` или нечисловое | `InvalidReaction` |
| `pre_exponential_factor ≤ 0` | `InvalidReaction` |
| Мультипликатор условия < 0 или нечисловой | `InvalidReaction` |
| Суммарный вес каналов = 0 (normalize_weights) | `InvalidReaction` |
| Дублирующаяся запись в `EnergyModel` | `InvalidReaction` |
| Отсутствует энергия изомера при equilibrium_distribution | `InvalidReaction` |
| Температура ≤ 0 | `InvalidReaction` |

## Связи

- [[core-reaction|Reaction]] — `ReactionChannel` хранится в `Reaction.channels`; `ProductDistribution` описан там же
- [[core-catalysis|catalysis]] — `CatalystSurfaceId` используется в `ChannelConditionEffect::Surface` и `IsomerEnergy`
- [[core-simulation|simulation]] — вызывает `channel_product_distribution` и `channel_rate_sum_per_second`
- [[selectivity-engine|SelectivityEngine]] — применяется внутри `channel_rate_weight`
- [[core-mixture|Mixture]] — `total_concentration_in_phase`, `temperature_kelvin`
- [[core-error|ChemistryError]]
