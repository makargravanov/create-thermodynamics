# Ароматическая логика и каталог органики

Исходный код: `organic/aromatic.rs`, `organic/catalog.rs`

## Назначение

`aromatic.rs` реализует модель электронных и стерических эффектов заместителей ароматического кольца. Используется генераторами EAS (электрофильное ароматическое замещение), FC-реакций и SNAr.

`catalog.rs` — минималистичный контейнер результата генерации: новые вещества + новые реакции.

## Ключевые типы

### `SubstituentClass`

Классификация заместителей по электронному эффекту:

| Вариант | Примеры групп |
|---------|--------------|
| `StronglyActivating` | −OH, −NH2 |
| `ModeratelyActivating` | −OR, −NHCOR, −SR |
| `WeaklyActivating` | алкил, −C(без EWG) |
| `WeaklyDeactivating` | −F, −Cl, −Br, −I |
| `ModeratelyDeactivating` | −CHO, −COR, −COOH, −SO3H |
| `StronglyDeactivating` | −NO2, −NR3+, −CN, −CCl3 |

Метод `effect_for_distance(distance: usize) -> f64`:
- distance 1 (ипсо-положение): наибольший эффект
- distance 2 (орто/пара): промежуточный
- иные расстояния: суммарный фоновый эффект

### `AromaticSubstituent`

```rust
pub struct AromaticSubstituent {
    pub ring_carbon: usize,
    pub substituent_atom: usize,
    pub class: SubstituentClass,
}
```

### `AromaticPosition`

Описывает конкретную позицию кольца: `atom`, `has_h`, `substituents`, `steric_penalty`, `electronic_delta`.

### `AromaticRingDescriptor<'a>`

Главный объект анализа кольца. Строится двумя способами:

- `from_start_carbon(structure, start)` — BFS по связям порядка 1.5 от одного атома
- `from_ring_atoms(structure, atoms)` — по готовому списку атомов кольца

### `GeneratedOrganicCatalog`

```rust
pub(crate) struct GeneratedOrganicCatalog {
    pub substances: Vec<Substance>,
    pub reactions: Vec<Reaction>,
}
```

`Default` — пустой каталог.

## Поток данных / Алгоритм

### `AromaticRingDescriptor::substituents()`

Обходит атомы кольца, ищет нециклические, не-H соседей с нероматической связью и классифицирует их через `classify_substituent`.

### `classify_substituent`

Классификация первого атома заместителя:
- O: StronglyActivating (есть H) или ModeratelyActivating (OR/ester)
- N: StronglyDeactivating (заряд > 0) или ModeratelyActivating (amide) или StronglyActivating (amine)
- C: ModeratelyDeactivating (C=O), StronglyDeactivating (C≡N или CCl3), иначе WeaklyActivating
- Halogens: WeaklyDeactivating
- S: ModeratelyDeactivating (≥2 S=O), иначе WeaklyActivating

### `compute_eas_activation_delta(target_carbon)`

Суммирует `effect_for_distance` всех заместителей кольца до target_carbon + стерический штраф за объёмные заместители (Br, I, третичный C) в ортоположении (+3.0).

Отрицательный delta → сайт активирован (реакция легче). Положительный → деактивирован.

### `compute_snar_activation_delta(halogen_carbon)`

Только для SNAr: суммирует вклад EWG-заместителей (ModeratelyDeactivating / StronglyDeactivating) на расстоянии 1 или 3 от атома галогена:
- StronglyDeactivating (орто/пара) → −12.0
- ModeratelyDeactivating → −6.0

Возвращает отрицательное число (активация SNAr); если ≥ 0 — реакция не генерируется.

### `is_deactivated_for_fc()`

Возвращает true при любом WeaklyDeactivating или сильнее — FC-реакции блокируются.

### `is_aromatic_ring_preserved(descriptor, product, mapping)`

Проверяет, что все связи порядка 1.5 между атомами кольца сохранились в продукте после атомного маппинга. Используется как постусловие в EAS-генераторах.

## Инварианты и ошибки

- BFS в `from_start_carbon` посещает только атомы с ароматической связью (1.5). Если кольцо не существует или атомов < 5, `ring_atoms` будет маленьким, и SNAr-генераторы вернут `None`.
- `effect_for_distance` не обрезает сумму — итоговый `activation_delta` может стать отрицательным (очень активированный сайт), минимальная EA зажата в `max(..., 5.0)` в генераторах.

## Связи

- [[organic-generators]] — `generate_aromatic_*`, `generate_fc_*`, `generate_aryl_halide_*` используют `AromaticRingDescriptor`
- [[organic-engine]] — собирает `GeneratedOrganicCatalog` как итог генерации
- [[molecule-graph]] — `MolecularStructure`, `bond_order_matches`
- [[core-reaction]] — `Reaction`
- [[core-substance]] — `Substance`
