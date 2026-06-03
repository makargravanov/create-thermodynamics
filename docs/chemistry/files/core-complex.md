# Координационные комплексы (complex)

Исходный код: `core/complex.rs`

## Назначение

Описывает координационные комплексы металлов: центральный ион, лиганды с их дентатностью,
геометрию, константу устойчивости и лабильность обмена лигандов. Автоматически строит
`EquilibriumSpec` и `Substance` для регистрации в реестре.

## Ключевые типы

| Тип | Роль |
|-----|------|
| `ComplexLigand` | один тип лиганда: `substance_id`, `count`, `denticity`; `occupied_sites = count × denticity` |
| `ComplexGeometry` | `Linear(2) / Tetrahedral(4) / SquarePlanar(4) / Octahedral(6) / Unknown` |
| `LigandExchangeLability` | `Labile / Intermediate / Inert` — скорость обмена лигандов |
| `ComplexSpec` | полное описание комплекса: ион, лиганды, геометрия, заряд, `formation_constant`, фаза, цвет, теги |

## Публичные входы

- `ComplexSpec::new(id, central_ion, ligands, charge, formation_constant)` — конструктор; координационное число вычисляется автоматически
- `ComplexSpec::with_coordination_number / with_geometry / with_ligand_exchange_lability / with_phase / with_color_argb / with_tags` — строитель
- `ComplexSpec::validate_shape()` — полная структурная валидация
- `ComplexLigand::new(substance_id, count).with_denticity(d)` — задание лиганда
- `ComplexLigand::occupied_sites()` — `count × denticity` (насыщающее умножение)

## Поток данных / Алгоритм

### Регистрация в реестре

`ComplexSpec` не регистрируется напрямую как `Substance` или `Reaction` — реестр вызывает
два внутренних метода:

```
ComplexSpec::to_equilibrium() → EquilibriumSpec
  реагенты: central_ion(1) + каждый ligand(count) в phase
  продукты: complex_id(1) в phase
  K = formation_constant

ComplexSpec::to_substance(molar_mass, charge) → Substance
  проверяет charge == self.charge
  preferred_liquid_phase = Aqueous | Organic (по self.phase)
  aqueous_solubility = 10 моль/ведро (если Aqueous), иначе 0
  organic_solubility  = 0 (если Aqueous), иначе 10
  can_precipitate = true
```

Это значит, что равновесие образования комплекса решается тем же итерационным
решателем `equilibrate_solution_equilibria` из [[core-solution|solution]], что и обычные
химические равновесия.

### Валидация формы (`validate_shape`)

1. Непустые `id`, `central_ion`, `ligands`
2. `coordination_number > 0`
3. Если геометрия не `Unknown`: `coordination_number` обязан совпадать с ожидаемым
4. `formation_constant > 0` и конечна
5. Фаза не `Gas` или `Solid`
6. Каждый лиганд: `count > 0`, `denticity > 0`, непустой `substance_id`
7. `Σ occupied_sites == coordination_number`

## Инварианты и ошибки

| Условие | Ошибка |
|---------|--------|
| `coordination_number` не совпадает с геометрией | `InvalidReaction` |
| `Σ occupied_sites ≠ coordination_number` | `InvalidReaction` |
| `formation_constant ≤ 0` или не конечна | `InvalidReaction` |
| Фаза `Gas` или `Solid` | `InvalidReaction` |
| Заряд вещества ≠ `self.charge` при `to_substance` | `ChargeNotConserved` |

## Связи

- [[core-solution|solution]] — `ComplexSpec::to_equilibrium()` используется при построении `IndexedEquilibrium`; `AqueousComplexForm` описывает состояние комплекса в растворе
- [[core-substance|Substance]] — `to_substance()` генерирует вещество для реестра
- [[core-registry|ChemistryRegistry]] — хранит `ComplexSpec`, вызывает `to_equilibrium` и `to_substance`
- [[core-error|ChemistryError]]
