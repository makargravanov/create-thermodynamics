# Планировщик синтеза

Исходный код: `chemistry/synthesis.rs`

## Назначение

`SynthesisPlanner` реализует поиск мультистадийных синтетических маршрутов:
заданы начальные вещества и целевая структура — планировщик находит последовательности
реакций, ведущих к цели, с ограничением по числу шагов, маршрутов и политикой безопасности.

## Ключевые типы

### SynthesisPlanner

| Поле | По умолчанию | Смысл |
|---|---|---|
| `max_steps` | 6 | Максимальная длина одного маршрута |
| `max_routes` | 8 | Максимальное число маршрутов в результате |
| `allowed_reaction_prefixes` | пусто (все) | Белый список префиксов ID реакций |
| `safety_policy` | пусто | Список запрещённых веществ/структур |

### SynthesisSafetyPolicy

Запрещает вещества по `SubstanceId` или по канонической FROWNS-структуре.
Проверяется как для целевого вещества (до поиска), так и для каждого продукта
промежуточного шага.

### SynthesisRoute

```rust
pub struct SynthesisRoute {
    pub target: SubstanceId,
    pub steps: Vec<SynthesisStep>,
    pub estimated_yield: f64,   // произведение product_fraction по шагам
    pub score: f64,             // len(steps) + (1 - yield); меньше = лучше
}
```

### SynthesisStep

Один шаг: реакция, список реагентов и продуктов, доля выхода.
Для реакций с `product_distribution` — `product_fraction` равен максимуму
по вариантам; для `channels` — `1 / count(products)`.

## Публичные входы

| Метод | Описание |
|---|---|
| `SynthesisPlanner::new()` | Конструктор с умолчаниями |
| `.with_max_steps(n)` | Ограничение длины маршрута |
| `.with_max_routes(n)` | Ограничение числа результатов |
| `.allow_reaction_prefix(prefix)` | Добавить префикс в белый список |
| `.with_safety_policy(policy)` | Установить политику безопасности |
| `.find_routes(registry, starting, target)` | Запустить поиск |
| `SynthesisSafetyPolicy::deny_substance(id, reason)` | Запретить вещество |
| `SynthesisSafetyPolicy::deny_structure(structure, reason)` | Запретить по структуре |

## Поток данных / Алгоритм

`find_routes` реализует обход в ширину (BFS) по пространству состояний:

```mermaid
flowchart TD
    A[начальные вещества → FrontierState] --> B{очередь пуста?}
    B -- нет --> C[state = pop_front]
    C --> D{steps.len >= max_steps?}
    D -- да --> B
    D -- нет --> E[generate_reactions_for_substances(state.known, 1)]
    E --> F[кандидаты: реакции с разрешёнными префиксами,<br>где все реагенты входят в state.known]
    F --> G{шаг добавляет новые продукты?}
    G -- нет --> B
    G -- да --> H{продукты разрешены safety_policy?}
    H -- нет --> B
    H -- да --> I{next_known содержит target?}
    I -- да --> J[добавить route, обрезать до max_routes, сортировать по score]
    I -- нет --> K{seen_known_sets содержит ключ?}
    K -- да --> B
    K -- нет --> L[push FrontierState]
    L --> B
    B -- пусто --> M[вернуть routes]
```

Состояние фронтира (`FrontierState`) — это набор известных веществ на текущем шаге.
Дедупликация по `Vec<SubstanceId>` (ключ = отсортированный список).

На каждом шаге вызывается `generate_reactions_for_substances(..., 1)` — т.е. ровно
одна итерация генерации для текущего набора веществ. Это позволяет планировщику
открывать новые реакции по мере продвижения без полного раскрытия пространства.

## Инварианты и ошибки

- `max_steps == 0` → `InvalidReaction`
- `max_routes == 0` → `InvalidReaction`
- Целевое вещество, запрещённое политикой → `InvalidReaction` с
  `reaction_id = "<synthesis-policy>"` до начала поиска
- Если начальные вещества уже содержат цель — возвращается единственный маршрут
  с пустым `steps` и `score = 0.0`
- `score = steps.len + (1 - estimated_yield)` — маршруты сортируются от лучшего к худшему

## Связи

- [[dynamic|DynamicChemistryRegistry]] — источник реакций и разрешение структур
- [[molecule-frowns|FROWNS]] — каноническая форма для проверки политики структур
- [[core-reaction|Reaction]] — анализ `products`, `product_distribution`, `channels`

