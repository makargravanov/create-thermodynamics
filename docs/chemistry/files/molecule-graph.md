# Граф молекулы

Исходный код: `molecule/graph.rs`

## Назначение

Центральный модуль молекулярного ядра. Определяет структуры данных для представления молекулы как помеченного мультиграфа, парсеры нескольких текстовых форматов (linear, topology, java-legacy), редактор молекул и логику вычисления явных водородов.

> **Ключевой инвариант:** водороды в этой модели ВСЕГДА явные атомы в графе. `hydrogen_count` просто считает H-атомы в списке соседей — скрытых H нет. Метод `add_missing_hydrogens` добавляет H как реальные вершины при финализации `StructureBuilder`.

## Ключевые типы

### MolecularAtom

| Поле | Тип | Смысл |
|------|-----|-------|
| `element` | `String` | Символ элемента (`"C"`, `"H"`, `"O"`, `"R"`, …) |
| `charge` | `f64` | Частичный или целочисленный заряд |
| `r_group_number` | `u8` | Номер R-группы (0 = не R-группа) |

### MolecularBond

| Поле | Тип | Смысл |
|------|-----|-------|
| `from` | `usize` | Индекс первого атома |
| `to` | `usize` | Индекс второго атома |
| `order` | `f64` | Порядок связи: 1.0, 2.0, 3.0, 1.5 (ароматическая) |

Граф **не ориентированный**: связь `from-to` эквивалентна `to-from`.

### Stereochemistry

```
Tetrahedral(TetrahedralStereo)   — хиральный центр, 4 заместителя, cw/ccw
DoubleBond(DoubleBondStereo)     — E/Z или cis/trans по двойной связи
Mixture { atoms, kind }          — смесь стереоизомеров
```

### MolecularStructure

Неизменяемое (clone-on-write) представление молекулы:

```rust
pub struct MolecularStructure {
    pub source_code: String,       // исходный FROWNS-код или "generated"
    pub atoms: Vec<MolecularAtom>,
    pub bonds: Vec<MolecularBond>,
    pub stereochemistry: Vec<Stereochemistry>,
}
```

Важные методы:
- `validate_structure()` — структурная целостность (индексы, связность, элементы), без валентности
- `validate()` — полная проверка: структура + валентность + стереохимия
- `summary()` — молярная масса и суммарный заряд
- `neighbors(i)` → `Vec<(usize, f64)>` — соседи атома i
- `explicit_hydrogen_count(i)` — число явных H-атомов в соседях (единственный способ считать H)
- `is_connected()` — BFS-проверка связности

### MolecularEditor

Изменяемая обёртка над `MolecularStructure`. Содержит флаг `modified`; при `finish()` если флаг поднят, `source_code` заменяется на `"generated"`.

| Метод | Действие |
|-------|----------|
| `add_atom(parent, el, charge, order)` | Добавляет атом + связь к parent |
| `add_group(parent, group, root, order)` | Вставляет подструктуру, смещает индексы |
| `add_bond(a, b, order)` | Добавляет связь (нельзя дублировать) |
| `remove_bond(a, b)` | Удаляет связь, валидирует результат |
| `remove_atom(i)` | Удаляет атом + все его связи, смещает индексы |
| `remove_atoms(slice)` | Пакетное удаление, возвращает маппинг старых→новых индексов |
| `replace_atom(i, el, charge)` | Меняет элемент/заряд без изменения топологии |
| `set_bond_order(a, b, order)` | Меняет порядок связи |
| `insert_bridging_atom(a, b, el, charge)` | Разрывает связь a-b, вставляет мостиковый атом |
| `finish()` | validate → aromatize → validate, возвращает `MolecularStructure` |
| `join_structures(s1, a1, s2, a2, order)` | Статический: соединяет две структуры |
| `split_at_bond(s, a, b)` | Статический: делит молекулу по мостовой связи |

Любая операция, затрагивающая стереохимический центр, очищает его запись (`clear_stereo_at_atom` / `clear_stereo_at_bond`).

## Поток данных / Алгоритм

### Парсинг linear FROWNS

`parse_linear_structure` → `StructureBuilder` — рекурсивный однопроходный парсер:
- символы элементов — `[A-Z][a-z]*`
- `=` / `#` / `~` / `-` устанавливают `pending_bond` для следующего атома
- `(` / `)` — стек для ветвлений
- `^` — заряд при атоме
- `R0..R9` — R-группы
- После построения вызывается `add_missing_hydrogens()`

### Добавление водородов

`hydrogens_to_add(element, bonds, charge)` = `next_lowest_valency - |charge| - bonds`, округлённое вниз. `next_lowest_valency` возвращает наименьшую целевую валентность ≥ текущему числу связей. Каждый добавленный H — отдельный `MolecularAtom { element: "H", charge: 0.0 }`.

### Парсинг topology FROWNS

`parse_topology_structure` строит шаблонный граф (`benzene`, `cubane`, `cyclohexene`, `diborane`, `octasulfur`, `tetraborate`, `anthracene`, `isohydrobenzofuran`) и прикрепляет линейные группы к точкам присоединения.

## Инварианты и ошибки

- Граф всегда связен (`is_connected()` = true); нарушение → `InvalidSubstance`.
- Индексы связей всегда в диапазоне `[0, atoms.len())`, `from ≠ to`.
- Заряд конечный (`is_finite()`).
- Порядок связи > 0 и конечный.
- Валентность: `sum(bond_orders) - |charge| ≤ max_valency + ε` (ароматическая 1.5 считается как 1.0).
- R-атом с ненулевым `r_group_number` должен иметь `element == "R"`.
- Стереохимия: заместители тетраэдра — 4 различных одиночных соседа центра; двойная связь — порядок 2.0 между `first` и `second`.
- После удаления/замены атомов, затрагивающих стереоцентр, запись стереохимии удаляется автоматически.

## Связи

- [[molecule-frowns|FROWNS]] — использует `parse_legacy_structure`, `aromatize`, `MolecularStructure`
- [[molecule-canonical|Канонизация]] — использует `MolecularStructure`
- [[molecule-aromatic-perception|Ароматическое восприятие]] — экспортирует `aromatize`
- [[molecule-functional-group|Функциональные группы]] — использует `neighbors`, `hydrogen_count`, `carbon_degree`
- [[molecule-reactive-site|Реактивные центры]] — использует `MolecularStructure`
