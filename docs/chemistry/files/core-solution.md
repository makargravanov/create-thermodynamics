# Равновесия раствора (solution)

Исходный код: `core/solution.rs`

## Назначение

Решает химические равновесия в растворе: общие равновесия (комплексы, растворимость),
осаждение из пересыщенных растворов и кислотно-основные равновесия. Вычисляет pH,
активность ионов и коэффициенты активности через модель Дэвиса.

## Ключевые типы

### Спецификации (описание, не состояние)

| Тип | Роль |
|-----|------|
| `EquilibriumSpec` | одно равновесие: реагенты, продукты, K, опционально ΔH и T_ref |
| `EquilibriumTerm` | одна сторона равновесия: вещество, коэффициент, фаза |
| `PrecipitationSpec` | равновесие осаждения: твёрдое тело, ионы, K_sp, опц. ΔH и T_ref |
| `AcidBaseSpec` | кислотно-основная пара: кислота, сопряжённое основание, pKa, протон |

### Индексированные (рантайм, с SubstanceIndex)

| Тип | Роль |
|-----|------|
| `IndexedEquilibrium` | `EquilibriumSpec` + термы с `SubstanceIndex` |
| `IndexedPrecipitation` | `PrecipitationSpec` + индекс твёрдого тела |

### Состояние и наблюдаемые

| Тип | Роль |
|-----|------|
| `SolutionState` | pH (однозначный водный) + состояние по фазам |
| `LiquidPhaseSolutionState` | ионная сила, активность протона, pH, коэффициенты активности |
| `AqueousEquilibriumSystem` | снимок водной фазы для кислотно-основного решателя |
| `ActivityModel` | только `Davies`; `coefficient(charge, I, T)` |

### Внутренние структуры кислотно-основного решателя

`AcidBaseNumericalSystem`, `AcidBaseNumericalPair`, `AcidBaseSolution`, `AcidBaseTarget`

## Публичные входы

- `EquilibriumSpec::new(id, reactants, products, K)` + `.with_enthalpy_change_kj_per_mol / with_reference_temperature_kelvin`
- `PrecipitationSpec::new(id, solid, ions, K_sp)` + аналогичные модификаторы
- `AcidBaseSpec::new(id, acid, conjugate_base, pKa)` + `.with_proton`
- `solution_state(registry, mixture)` — pH и коэффициенты активности для всех жидких фаз
- `aqueous_equilibrium_systems(registry, mixture)` — полный снимок для отладки и JNI
- `activity_of(registry, mixture, substance_id, phase)` — активность с поправкой Davies; вода = 1 если присутствует
- `equilibrate_solution_equilibria(registry, mixture)` — главный решатель, мутирует `Mixture`, возвращает максимальное Δ

## Поток данных / Алгоритм

### Температурная поправка констант (уравнение Вант-Гоффа)

```
K(T) = K_ref · exp(−ΔH·1000/R · (1/T − 1/T_ref))
```
Если `ΔH = 0`, K(T) = K_ref. Применяется к `EquilibriumSpec` и `PrecipitationSpec`.

### Модель активности Дэвиса

```
log₁₀ γ = −A · (T_ref/T)^0.5 · z² · (√I/(1+√I) − 0.3·I)
```
A = 0.509 при 298 K. Нейтральные вещества и вода имеют γ = 1. Значение обрезается в [0, 1].

### Главный решатель (`equilibrate_solution_equilibria`)

```
1. equilibrate_aqueous_acid_base_systems  (до 64 итераций активности)
2. До 256 проходов:
     для каждого IndexedEquilibrium (не кислотно-основные, не автоионизация воды):
         apply_equilibrium  → bisection(80 шагов) по степени превращения
     для каждого IndexedPrecipitation:
         apply_precipitation → bisection(80 шагов)
     если max_delta ≤ 1e-10 → сходимость
3. Если не сошлось → EquilibriumInvariantViolation
```

### Кислотно-основной решатель (специализированный)

Используется вместо общего итератора для систем с `AcidBaseSpec`, поскольку стандартная
итерация по равновесиям не обеспечивает нужной точности pH.

```
Для каждой итерации активности (до 64):
  1. Собрать AqueousEquilibriumSystem
  2. build_acid_base_numerical_system:
       - вычислить target_charge = Σ c_i · z_i (все виды)
       - fixed_charge = заряд видов, не участвующих в кислотно-основных равновесиях
  3. solve_log_h_for_charge_balance:
       bisection(96 шагов, log₁₀[H⁺] в [−16, 2])
       + Newton(12 шагов) с численным градиентом
  4. acid_base_targets_at_log_h → концентрации [H⁺], [OH⁻] и каждой кислотно-основной пары
  5. apply_aqueous_acid_base_solution → мутирует Mixture
  если delta ≤ 1e-8 → выход
```

Балансное уравнение: `charge_residual(log_h) = fixed_charge + z_H·[H] + z_OH·[OH] + Σ(z_a·[A] + z_b·[B]) − target_charge`

### Осаждение (`apply_precipitation`)

```
Q = Π a_i^ν_i  (ионное произведение)
если Q > Ksp → осаждение (ионы → твёрдое)
если Q < Ksp и solid > 0 → растворение (твёрдое → ионы)
```
Бисекция по степени превращения; проверяется `ln(Q/K) ≤ 1e-8`.

## Инварианты и ошибки

| Условие | Ошибка |
|---------|--------|
| `equilibrium_constant ≤ 0` | `InvalidReaction` при регистрации |
| Несохранение массы в `EquilibriumSpec` | `MassNotConserved` |
| Несохранение заряда в `AcidBaseSpec` | `ChargeNotConserved` |
| Несохранение массы в `PrecipitationSpec` | `MassNotConserved` |
| Несколько водных жидких фаз при запросе pH | `InvalidMixtureState` |
| Решатель не сошёлся за 256 проходов | `EquilibriumInvariantViolation` |
| К-б решатель не сошёлся за 64 итерации активности | `EquilibriumInvariantViolation` |
| Несколько водных фаз в к-б решателе | `EquilibriumInvariantViolation` |
| Коэффициент активности вне [0, 1] | `InvalidMixtureState` |

Числовые константы: `ACID_BASE_BISECTION_STEPS = 96`, `ACID_BASE_NEWTON_STEPS = 12`,
`EQUILIBRIUM_MAX_PASSES = 256`, `WATER_ION_PRODUCT = 1e-14`, `DAVIES_A_AT_298_K = 0.509`.

## Связи

- [[core-complex|complex]] — `ComplexSpec::to_equilibrium()` регистрируется как `EquilibriumSpec`
- [[core-registry|ChemistryRegistry]] — хранит `IndexedEquilibrium`, `IndexedPrecipitation`, `AcidBaseSpec`
- [[core-mixture|Mixture]] — `apply_reaction_phase_deltas_by_index`, `liquid_phase_ionic_strengths`, `apply_aqueous_targets_by_index`
- [[core-simulation|simulation]] — вызывает `equilibrate_solution_equilibria` после каждого шага реакции
- [[core-error|ChemistryError]]
