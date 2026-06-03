# Substance — описание химического вещества

Исходный код: `core/substance.rs`

## Назначение

Определяет единицу химического справочника: всё, что система знает о веществе
до появления его в реакторе. Сюда входят термодинамические параметры,
агрегатное поведение, растворимость, молекулярная структура и редокс-роли.
Типы этого модуля — строительные блоки для [[core-registry|реестра]] и
[[core-mixture|смеси]].

## Ключевые типы

### `SubstanceId` / `SubstanceTagId`

Newtype-обёртки над `String`. `SubstanceId` уникально идентифицирует вещество
(например, `"destroy:water"`). Пустая строка отвергается при создании.
`SubstanceTagId` используется для логических меток (например,
`"destroy:hypothetical"` блокирует добавление в смесь).

### `Substance`

Центральный тип данных. Хранит:

- `molar_mass_grams` — молярная масса (г/моль);
- `liquid_density_grams_per_bucket` / `solid_density_grams_per_bucket` — плотности для расчёта объёма смеси;
- `melting_point_kelvin` / `boiling_point_kelvin` — границы фаз;
- `molar_heat_capacity_j_per_mol_kelvin` — теплоёмкость для `heat()`;
- `fusion_heat_j_per_mol` — теплота плавления (латентная для твёрдой→жидкой);
- `latent_heat_j_per_mol` — теплота испарения;
- `critical_temperature_kelvin` / `critical_pressure_pascal` — критическая точка (оба или ни одного);
- `acentric_factor` — требует критической точки;
- `vapor_pressure_model` — явная модель давления пара (иначе используется Клаузиус–Клапейрон по точке кипения);
- `phase_properties: SubstancePhaseBehavior` — растворимость и фазовые предпочтения;
- `representation: SubstanceRepresentation` — тип материала;
- `molecular_structure` / `functional_groups` — граф молекулы и найденные функциональные группы;
- `redox_roles` / `explicit_oxidation_states` — окислительно-восстановительные роли.

### `SubstanceRepresentation`

Перечисление типов материала. Каждый вариант накладывает свои инварианты:

| Вариант | Ключевые ограничения |
|---|---|
| `Molecular` | нейтральные молекулы, молграф допустим |
| `Ion` | заряд ≠ 0, `preferred_liquid_phase` = Aqueous |
| `IonicSolid` | нейтральный, `can_precipitate = true`, не растворитель |
| `Metal` | нейтральный, не растворитель |
| `Oxide` | нейтральный, не растворитель |
| `Hydrate` | нейтральный, `water_count > 0` |
| `SurfaceMaterial` | нейтральный, не растворитель |
| `UnspecifiedMaterial` | требует непустого `reason` |

### `VaporPressureModel`

Два варианта модели давления пара:

- `ClausiusClapeyron` — экспоненциальная форма через `ΔH`, опорную температуру и давление;
- `Log10PressurePascalAntoine` — уравнение Антуана: `log₁₀(P) = A − B/(T+C)`.
  Диапазон `[min_T, max_T]` опционален, но обязан быть упорядоченным;
  знаменатель `T+C` должен оставаться положительным на всём диапазоне.

Если `vapor_pressure_model` не задана, `vapor_pressure_pascal()` автоматически
применяет Клаузиус–Клапейрон, используя `boiling_point_kelvin` как опорную
точку при `P = 101 325 Па`. Ионы и вещества выше критической температуры
возвращают `None`.

### `SubstancePhaseBehavior`

Управляет распределением по фазам в [[core-mixture|смеси]]:

- `preferred_liquid_phase: LiquidPhasePreference` — в какую жидкую фазу идёт вещество по умолчанию;
- `aqueous_solubility_mol_per_bucket` / `organic_solubility_mol_per_bucket` — лимиты растворимости (`None` = безграничная);
- `can_precipitate` — может ли избыток образовывать твёрдую фазу;
- `can_form_liquid_phase` — может ли вещество само формировать жидкость (обязательно для растворителей);
- `solvent_role: SolventRole` — `NotSolvent` / `KnownSolvent` / `ConservativePredictedSolvent`.

Конструкторы-помощники: `aqueous_unlimited()`, `aqueous_solvent()`,
`organic_unlimited(aq_sol)`, `organic_solute(aq_sol)`.

## Публичные входы

### `Substance::new(...)`

Минимальный конструктор. Автоматически выбирает `phase_properties` и
`representation` по `id` и `charge`:

- `id == "destroy:water"` → `aqueous_solvent()`;
- `charge == 0` → `organic_unlimited(0.05)`, `Molecular`;
- `charge ≠ 0` → ион с ограниченной водной растворимостью `10 моль/ведро`.

### Builder-методы `with_*()`

Каждый метод возвращает `Self`, позволяя цепочку вызовов. Важные:
`with_phase_properties`, `with_representation`, `with_critical_point`,
`with_molecular_structure` (автоматически определяет `functional_groups`),
`with_melting_point_kelvin`, `with_fusion_heat_j_per_mol`.

### `aggregate_state_at(T)`

Возвращает `Solid` / `Liquid` / `Gas` по сравнению с `melting_point_kelvin` и
`boiling_point_kelvin`. Требует T ≥ 0.

### `vapor_pressure_pascal(T)`

Рассчитывает давление пара по явной модели или по Клаузиус–Клапейрону.
Возвращает `None` для ионов и сверхкритических состояний.

### `validate()`

Полный контроль консистентности объекта. Проверяет физические числа, пары
(critical T, critical P), что молграф согласован по заряду и молярной массе
с допуском `1e-6 г/моль`, что степени окисления корректно назначены.

## Инварианты и ошибки

Все ошибки — `ChemistryError::InvalidSubstance { substance_id, reason }`.

Ключевые правила:

- `melting_point ≤ boiling_point`;
- `critical_temperature > boiling_point` (если задана);
- критические T и P должны задаваться парой;
- `acentric_factor` требует критической точки;
- заряженные вещества: только водная фаза, не растворители;
- расплавы (`MoltenMetal`, `MoltenSlag`): не являются растворителями в обычном смысле;
- `IonicSolid` не может иметь молграф и обязан иметь `can_precipitate = true`.

## Связи

- [[core-mixture|Mixture]] — потребляет `Substance` через реестр для фазового распределения и нагрева
- [[core-registry|Registry]] — хранит и индексирует все `Substance`
- [[core-thermodynamics|Thermodynamics]] — использует `vapor_pressure_pascal` для VLE
- [[core-redox|Redox]] — `redox_roles`, `explicit_oxidation_states`
- [[molecule-graph|MolecularStructure]] — хранится внутри `Substance`
- [[molecule-functional-group|FunctionalGroup]] — вычисляется при `with_molecular_structure`
