# Пространство генерации и резолвер

Исходный код: `organic/space.rs`, `organic/resolver.rs`

## Назначение

`space.rs` строит индекс «все реактивные сайты всех веществ, сгруппированные по типу» — `OrganicGenerationSpace`. Это позволяет движку за O(1) получить всех участников любого типа сайта.

`resolver.rs` реализует `DerivedSubstanceResolver` — кэш FROWNS-строка → `SubstanceId` с автосозданием новых веществ для продуктов, не существующих в исходном реестре.

## Ключевые типы

### `GenerationScope`

```rust
pub(crate) struct GenerationScope { substances: BTreeSet<SubstanceId> }
```

Фильтр: только вещества из этого множества участвуют в генерации сайтов. Создаётся через `GenerationScope::all(registry)` или `from_substances(set)`.

### `SiteParticipant<'a>`

```rust
pub(crate) struct SiteParticipant<'a> {
    pub substance: &'a Substance,
    pub structure: &'a MolecularStructure,
    pub site: ReactiveSite,
}
```

Базовая единица генерации. `Clone`-абелен. Метод `is_seed(seeds)` возвращает true если `seeds == None` или вещество входит в seeds.

### `OrganicGenerationSpace<'a>`

```rust
pub(crate) struct OrganicGenerationSpace<'a> {
    pub all_substances: Vec<&'a Substance>,
    participants_by_site: BTreeMap<ReactiveSiteKind, Vec<SiteParticipant<'a>>>,
}
```

| Метод | Возвращает |
|-------|-----------|
| `new(substances, scope)` | `ChemistryResult<Self>` |
| `from_substances_for_scope(substances, scope)` | `ChemistryResult<Self>` |
| `site_participants()` | `impl Iterator<Item = SiteParticipant<'a>>` (все сайты) |
| `sites_of(kind)` | `impl Iterator<Item = SiteParticipant<'a>>` (только данного типа) |

### `DerivedSubstanceResolver`

```rust
pub(crate) struct DerivedSubstanceResolver {
    canonical_to_id: BTreeMap<String, SubstanceId>,
    generated_id_to_canonical: BTreeMap<SubstanceId, String>,
    pub substances: Vec<Substance>,
}
```

Константы по умолчанию для новых веществ:

| Параметр | Значение |
|----------|---------|
| `density` | 1000.0 г/л |
| `heat_capacity` | 100.0 Дж/(моль·К) |
| `latent_heat` | 20 000.0 Дж/моль |
| `boiling_point` | 1000.0 К (нейтральные) / `f64::MAX` (заряженные) |
| Каталожный цвет | `0x20FF_FFFF` |

## Поток данных / Алгоритм

### Построение `OrganicGenerationSpace`

```
substances (итератор &'a Substance)
  for each substance:
    all_substances.push(substance)
    if !scope.contains(id): skip
    if structure is None: skip
    try_find_reactive_sites(structure)
      → Vec<ReactiveSite>
    for each site:
      participants_by_site[site.kind].push(SiteParticipant {...})
```

`try_find_reactive_sites` — внешняя функция из [[molecule-reactive-site]].

### Резолвер: `resolve(structure)`

```
write_frowns(structure) → canonical: String
if canonical_to_id contains canonical:
    return existing id
id = SubstanceId::new(canonical)   // FROWNS сам является ID
check id collision in generated_id_to_canonical
Substance::new(id, ...) с дефолтными физическими параметрами
canonical_to_id.insert(canonical, id)
generated_id_to_canonical.insert(id, canonical)
substances.push(substance)
return id
```

FROWNS-строка одновременно служит идентификатором и ключом канонической формы.

## Инварианты и ошибки

- Вещество без `molecular_structure` не участвует в генерации сайтов, но попадает в `all_substances`.
- `resolve` возвращает `Err(InvalidSubstance)` при коллизии ID (одинаковый ID, разная FROWNS) — защита от ошибок канонизатора.
- Заряженные вещества получают `boiling_point = f64::MAX` (бесконечная температура кипения = всегда ионы в растворе).
- `generated_id_to_canonical` хранит только **вновь созданные** вещества, `canonical_to_id` — все (включая из исходного реестра).

## Связи

- [[organic-engine]] — создаёт `OrganicGenerationSpace` и `DerivedSubstanceResolver`, передаёт их в генераторы
- [[organic-centers]] — `SiteParticipant` является базой для всех XxxSite<'a>
- [[molecule-reactive-site]] — `try_find_reactive_sites`, `ReactiveSite`, `ReactiveSiteKind`
- [[molecule-frowns]] — `write_frowns`, `parse_frowns`
- [[core-substance]] — `Substance`, `SubstanceId`
