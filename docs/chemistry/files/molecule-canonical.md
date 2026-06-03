# Канонизация молекулярного графа

Исходный код: `molecule/canonical.rs`

## Назначение

Вычисляет канонический код молекулы — детерминированную строку-идентификатор, не зависящую от порядка атомов во входном представлении. Используется для дедупликации веществ в реестре и как выходной формат `write_frowns`.

## Ключевые типы

### CanonicalizationJob

Пошаговый объект, реализующий алгоритм. Поддерживает как синхронный (`run_to_completion`), так и инкрементный (`step`) режимы.

```rust
pub struct CanonicalizationJob {
    graph: CanonicalGraph,
    stack: Vec<Vec<Vec<usize>>>,  // стек разбиений для поиска
    best_body: Option<String>,
    processed_branches: usize,
    stage: CanonicalizationStage,
    completed: bool,
    result: Option<String>,
}
```

### CanonicalizationHandle

Обёртка для фоновой канонизации в отдельном потоке (`spawn_canonicalization`). Предоставляет `progress()`, `result()`, `join()` через `Arc<Mutex<…>>`.

### CanonicalGraph (internal)

Рабочее представление: вершины — `CanonicalVertex { token, neutral_hydrogens, original_index }`, рёбра — `CanonicalEdge { from, to, order: String }`.

**Нейтральные водороды:** ненаряженные H-атомы, которые являются соседями включённой вершины, вместо того чтобы быть отдельными вершинами графа канонизации, складываются в `neutral_hydrogens`. Заряженные H и H в стереохимии остаются явными вершинами. Это ключевое отличие: в `MolecularStructure` все H явные, в `CanonicalGraph` нейтральные H сворачиваются для компактного кода.

## Публичные входы

```rust
pub fn canonical_structure_code(structure: &MolecularStructure) -> ChemistryResult<String>
pub fn spawn_canonicalization(structure: MolecularStructure) -> CanonicalizationHandle
```

## Поток данных / Алгоритм

### Деревья (ациклические молекулы без стереохимии)

```
validate → CanonicalGraph::new → has_cycle? NO
→ tree_centers() (обрезка листьев)
→ rooted_tree_code(center) для каждого центра
→ min по строковому сравнению
→ "destroy:linear:<tree_code>"
```

Код дерева строится рекурсивно: `atom_token + sorted(branches)`, ветви в скобках `(…)`. Пример: `C(C)(O)` → ветви отсортированы лексикографически.

### Циклы и/или стереохимия — алгоритм разбиения (Nauty-подобный)

```
initial_partition()           — группировка по initial_label
  (токен | h_count | degree | sorted bond orders)
↓
stack: [initial_partition]
↓
loop:
  pop partition
  refine_partition()          — итерация до стабильности
    refinement_signature = initial_label + sorted neighbor classes
  if discrete → graph_body_for_order(order) → update best_body
  else → split smallest non-trivial cell → push children to stack
↓
"destroy:graph:<best_body>"
```

`graph_body_for_order` строит:
```
atoms=<tok0>.<tok1>...;bonds=<i>-<ord>-<j>,...[;stereo=...]
```
Связи сортируются лексикографически. Для стереохимии применяется `permutation_is_odd` — если перестановка заместителей нечётна, дескриптор инвертируется (cw↔ccw).

### Токены

| Сущность | Токен |
|----------|-------|
| Атом без заряда | `C`, `N`, `O`, `R2` |
| Атом с зарядом | `C^-1`, `O^-0.5` |
| Порядок связи (graph) | `1`, `2`, `3`, `1.5` |
| Порядок связи (linear) | `` (пусто), `=`, `#`, `~` |

## Инварианты и ошибки

- Нейтральные H не становятся вершинами `CanonicalGraph`, но учитываются в `neutral_hydrogens` — поэтому изотопные/заряженные H влияют на код.
- Молекулы с одинаковой связностью, но разными стереодескрипторами дают **разные** коды (тест `graph_code_preserves_tetrahedral_stereo`).
- Молекулы с разной связностью **никогда** не дают одинаковый код (тест `different_connectivity_does_not_collapse`).
- Цикл без допустимого порядка — `ChemistryError`.
- Фоновый поток при панике → `ChemistryError` из `join()`.

## Связи

- [[molecule-frowns|FROWNS]] — экспортирует `canonical_structure_code` как `write_frowns`
- [[molecule-graph|Граф молекулы]] — `MolecularStructure`, `Stereochemistry`
- [[data-catalog|Каталог данных]] — использует канонический код как ключ дедупликации
