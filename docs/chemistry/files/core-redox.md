# Окислительно-восстановительная химия (redox)

Исходный код: `core/redox.rs`

## Назначение

Полная модель редокс-химии: присвоение степеней окисления, описание полуреакций и парных
реакций, расчёт потенциала по уравнению Нернста, симуляция электролиза по закону Фарадея.
Является основой для валидации электронного баланса при регистрации реакций в
[[core-registry|ChemistryRegistry]].

## Ключевые типы

### Степени окисления

| Тип | Роль |
|-----|------|
| `OxidationStateRule` | `Electronegativity` — по разнице ЭО связей; `Explicit` — задана вручную |
| `AtomOxidationState` | степень окисления одного атома с индексом и элементом |
| `OxidationStateAssignment` | полный набор состояний молекулы + суммарный заряд |
| `ExplicitOxidationState` | явное задание (элемент, значение, кол-во атомов) для неорганики |

### Полуреакции

| Тип | Роль |
|-----|------|
| `RedoxHalfReaction` | одна полуреакция: реагенты, продукты, кол-во электронов, `ElectronSide`, среда, стандартный потенциал |
| `ElectronSide` | `Reactant` (восстановление) / `Product` (окисление) |
| `RedoxEnvironment` | `Acidic / Basic / Neutral / Any` |

### Парные реакции

| Тип | Роль |
|-----|------|
| `RedoxAnnotation` | метаданные редокс-реакции: окислитель/восстановитель, число e⁻, ссылки на полуреакции, флаг `electron_balance_checked` |
| `RedoxPair` | готовая закрытая реакция + ссылки на обе полуреакции |
| `RedoxPairSpec` | спецификация для автогенерации парной реакции из двух полуреакций |
| `RedoxRole` | `Oxidant / Reductant / OxidantAndReductant` |

### Потенциал и электролиз

| Тип | Роль |
|-----|------|
| `RedoxPotentialEvaluation` | результат расчёта по Нернсту: потенциалы обоих электродов, E_cell, K_eq, `thermodynamic_rate_factor` |
| `ElectrolysisCell` | пара электродов + приложенное напряжение |
| `ElectrodeProcess` | полуреакция на электроде с КПД и перенапряжением |
| `ElectrolysisReport` | результат: заряд, моль e⁻, степени превращения, напряжения |

## Публичные входы

- `assign_oxidation_states(structure)` — по структуре молекулы, алгоритм по ЭО
- `explicit_oxidation_assignment(id, charge, states)` — ручное задание для ионов без графа
- `evaluate_redox_potential(registry, mixture, reaction)` — уравнение Нернста → `RedoxPotentialEvaluation`
- `apply_electrolysis_cell(registry, mixture, cell, current_A, duration_s)` — закон Фарадея → мутирует `Mixture`
- `build_redox_pair_reaction(spec, halves)` — автогенерация сбалансированной реакции из двух полуреакций
- `electron_reactant / electron_product(builder, count)` — хелперы добавления свободных электронов в реакцию
- Внутренние валидаторы: `validate_half_reaction_shape`, `validate_half_reaction_conservation`, `validate_redox_annotation`, `validate_redox_pair`

## Поток данных / Алгоритм

### Присвоение степеней окисления

```
MolecularStructure → assign_oxidation_states
  для каждой связи: ΔОС = bond.order · sign(ΔЭО)
  → OxidationStateAssignment (сумма = total_charge)
```

Для неорганических ионов без молекулярного графа используется `explicit_oxidation_assignment`;
сумма `state × atom_count` обязана совпадать с зарядом вещества (допуск 1e-6).

### Расчёт потенциала (Нернст)

```
E_half = E°_half − (RT / nF) · ln Q_half
E_cell = E_oxidation + E_reduction
exponent = n · F · E_cell / (R · T)
K_eq = exp(exponent)
thermodynamic_rate_factor = sigmoid(exponent)  # обрезается в [0, 1]
```

`thermodynamic_rate_factor` умножается на кинетическую скорость в [[core-simulation|simulation]]:
при отрицательном потенциале клетки фактор → 0, реакция заглушается.

### Генерация парной реакции (`build_redox_pair_reaction`)

```
n_lcm = lcm(oxidation.electron_count, reduction.electron_count)
scaled_left  = oxidation.reactants × (n_lcm / ox.n) + reduction.reactants × (n_lcm / red.n)
scaled_right = oxidation.products  × (n_lcm / ox.n) + reduction.products  × (n_lcm / red.n)
cancel_common_terms(left, right)   # убирает виды с обеих сторон
→ Reaction с RedoxAnnotation.from_halves
```

### Электролиз (закон Фарадея)

```
Q  = I · t  [Кл]
n_e = Q / F
ξ_anode   = n_e · η_anode   / anode.electron_count
ξ_cathode = n_e · η_cathode / cathode.electron_count
Mixture мутируется атомарно: сначала staged-клон, затем swap
```

Проверяется: `V_applied ≥ V_reversible + η_anode + η_cathode`; при недостаточном напряжении
`Mixture` не мутируется.

## Инварианты и ошибки

| Условие | Ошибка |
|---------|--------|
| Сумма степеней окисления ≠ заряд (явное задание) | `InvalidSubstance` |
| Нарушение сохранения массы в полуреакции | `MassNotConserved` |
| Нарушение сохранения заряда в полуреакции | `ChargeNotConserved` |
| Закрытая редокс-реакция со свободными электронами | `InvalidReaction` |
| Редокс-реакция с `allow_charge_imbalance` | `InvalidReaction` |
| Только одна из пары полуреакций имеет стандартный потенциал | `InvalidReaction` |
| Приложенное напряжение < требуемого | `InvalidMixtureState` |
| Недостаточно реагентов для электролиза | `InvalidMixtureState` |
| Несовместимые среды двух полуреакций | `InvalidReaction` |

Константа: `FARADAY_CONSTANT_COULOMBS_PER_MOL = 96 485.332…`
Маркер свободного электрона: `ELECTRON_EXTERNAL_ID = "redox:electron"`

## Связи

- [[core-reaction|Reaction / ReactionBuilder]] — RedoxAnnotation встраивается в Reaction
- [[core-registry|ChemistryRegistry]] — хранит `RedoxHalfReaction`; регистрирует и валидирует `RedoxPair`
- [[core-simulation|simulation]] — `thermodynamic_rate_factor` включается в расчёт скорости реакции
- [[core-mixture|Mixture]] — `activity_of`, `apply_reaction_phase_deltas_by_index`, `temperature_kelvin`
- [[core-error|ChemistryError]] — `MassNotConserved`, `ChargeNotConserved`, `InvalidReaction`, `InvalidMixtureState`
