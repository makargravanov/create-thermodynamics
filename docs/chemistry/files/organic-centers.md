# Реактивные центры органики

Исходный код: `organic/centers.rs`

## Назначение

Типизированные обёртки над `SiteParticipant`, дающие именованный доступ к конкретным атомам реактивного сайта. Каждый центр — это `clone`-абельная структура с лайфтаймом `'a`, привязанным к исходному `Substance`/`MolecularStructure`. Весь модуль — `impl SiteParticipant<'a>`: набор методов-конструкторов, которые проверяют тип сайта и извлекают атомы по элементу и порядку связи.

## Ключевые типы

### Центры реактивных сайтов

| Тип | Поля (атомы) | Дополнительно |
|-----|-------------|---------------|
| `HalideSite` | `carbon`, `halogen` | `degree` (1/2/3) |
| `AlcoholSite` | `carbon`, `oxygen`, `hydrogen` | `degree` |
| `AlkoxideSite` | `oxygen` | — |
| `CarbonylSite` | `carbon`, `oxygen` | `is_ketone: bool` |
| `AlphaCarbonCenter` | `carbonyl_carbon`, `carbonyl_oxygen`, `alpha_carbon`, `alpha_hydrogens` | `carbonyl_kind`, `acidity`, `steric_class`, `conjugation` |
| `CarboxylicAcidSite` | `carbon`, `hydroxyl_oxygen`, `hydroxyl_hydrogen` | — |
| `EsterSite` | `carbon`, `carbonyl_oxygen`, `alkoxy_oxygen` | — |
| `AcylChlorideSite` | `carbon`, `chlorine` | — |
| `AmideSite` | `carbon`, `nitrogen`, `nitrogen_hydrogens` | — |
| `AmineSite` | `nitrogen`, `hydrogens` | — |
| `PhosphineSite` | `phosphorus` | требует 3 заместителя |
| `PhosphoniumSaltSite` | `phosphorus`, `alpha_carbon`, `alpha_hydrogens` | требует 3 заместителя на P |
| `PhosphorusYlideSite` | `phosphorus`, `alpha_carbon` | `stability: YlideStability` |
| `PhosphonateCarbanionSite` | `phosphorus`, `alpha_carbon` | требует C с зарядом < −0.1 |
| `SulfoneCarbanionSite` | `sulfur`, `alpha_carbon` | требует ≥ 2 сульфоновых кислорода |
| `NitrileSite` | `carbon`, `nitrogen` | связь C≡N (порядок 3.0) |
| `NitroSite` | `nitrogen`, `oxygens: [usize; 2]` | ровно 2 кислорода |
| `UnsaturatedBondSite` | `high_degree_carbon`, `low_degree_carbon` | `is_alkyne: bool` |
| `BoraneSite` | `carbon`, `boron` | — |
| `BorateEsterSite` | `oxygen`, `boron` | — |
| `IsocyanateSite` | `nitrogen`, `functional_carbon`, `oxygen` | N=C=O |
| `ArylHalideSite` | `carbon`, `halogen` | — |

### Центры защитных групп

| Тип | Поля |
|-----|------|
| `SilylEtherCenter` | `oxygen`, `silicon` |
| `AcetalCenter` | `acetal_carbon`, `oxygen_a`, `oxygen_b` |
| `BocCarbamateCenter` | `nitrogen`, `carbonyl_carbon`, `carbonyl_oxygen`, `alkoxy_oxygen`, `tert_butyl_carbon` |
| `CbzCarbamateCenter` | `nitrogen`, `carbonyl_carbon`, `carbonyl_oxygen`, `alkoxy_oxygen` |

### Вспомогательные перечисления

| Тип | Варианты |
|-----|---------|
| `AlphaCarbonylKind` | `Aldehyde`, `Ketone`, `Ester` |
| `AlphaAcidityClass` | `Ordinary`, `Activated` (два карбонила на α-углерод) |
| `AlphaStericClass` | `Primary`, `Secondary`, `Tertiary` |
| `AlphaConjugation` | `None`, `Allylic`, `Benzylic` |
| `YlideStability` | `Unstabilized`, `SemiStabilized`, `Stabilized` |

## Поток данных / Алгоритм

Каждый метод `SiteParticipant::*_site(self)` следует шаблону:

1. `require_kind(expected)` — проверяет `self.site.kind`, иначе `InvalidReaction`.
2. Ищет атомы через `site_atom_by_element(element)` или `bonded_site_atom(parent, element, order)`.
3. Проверяет инварианты (количество водородов, заряды, количество заместителей).
4. Возвращает типизированную структуру-центр.

### Логика `AlphaCarbonCenter`

- `carbonyl_kind`: если C=O соединён ещё с O через одинарную связь → Ester; если есть H на карбонильном C → Aldehyde; иначе Ketone.
- `acidity`: Activated если α-C имеет второй карбонильный сосед (β-дикарбонил).
- `steric_class`: по числу углеродных соседей α-C (0–1 → Primary, 2 → Secondary, 3+ → Tertiary).
- `conjugation`: Benzylic если сосед α-C имеет ароматическую связь (1.5); Allylic если сосед имеет C=C.

### Логика `YlideStability` (фосфорные илиды)

- Stabilized: α-C имеет двойную связь C=C.
- SemiStabilized: α-C рядом с ароматикой (1.5) или аллильным C=C.
- Unstabilized: иначе.

### Логика `CarbonylSite.is_ketone`

`is_ketone = true` если `site.kind == Ketone` ИЛИ (`kind == Carbonyl` И ≥ 2 углеродных соседей карбонильного C).

## Инварианты и ошибки

- Все методы потребляют `self` (move), возвращая `ChemistryResult<XxxSite<'a>>`.
- Ошибки формируются через `site_error(reason)` → `ChemistryError::InvalidReaction` с ID `typed_site/{substance}/{atoms}`.
- `alpha_hydrogens.is_empty()` → ошибка: явные водороды обязательны.
- `PhosphineSite`: ровно 3 заместителя на P (нейтральный третичный фосфин).
- `SulfoneCarbanionSite`: минимум 2 сульфоновых кислорода (order ≥ 1.5).
- `AcetalCenter`: минимум 2 однократно связанных кислорода на ацетальном C.

## Связи

- [[organic-space-resolver]] — `SiteParticipant` (базовый тип всех центров)
- [[organic-generators]] — все функции `generate_*` принимают XxxSite<'a>
- [[organic-engine]] — диспетчеризует по `ReactiveSiteKind` и вызывает методы центров
- [[molecule-reactive-site]] — `ReactiveSite`, `ReactiveSiteKind`
