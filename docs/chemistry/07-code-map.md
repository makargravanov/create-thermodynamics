# Карта кода

## Главный модуль

`chemistry/mod.rs`

Собирает публичные модули и дает главные входы:

- `destroy_registry_builder`
- `destroy_registry_with_generated_reactions_builder`

## Данные Destroy

- `data/catalog.rs`: вещества Destroy.
- `data/destroy_reactions.generated.rs`: перенесенные явные реакции.
- `data/reactions.rs`: подключение сгенерированного файла.

## Основная модель

- `core/substance.rs`: вещество и теги.
- `core/reaction.rs`: реакция, условия, построитель.
- `core/registry.rs`: проверенный реестр, индексы реакции.
- `core/mixture.rs`: смесь, нагрев, кипение, смешивание.
- `core/simulation.rs`: расчет реакции за тик и до равновесия.

## Молекулярный слой

- `molecule/graph.rs`: граф, редактор графа, разбор старых структур.
- `molecule/frowns.rs`: разбор и запись FROWNS.
- `molecule/functional_group.rs`: поиск функциональных групп.

## Органическая динамика

- `organic/mod.rs`: генераторы органических реакций.
- `dynamic/mod.rs`: динамические вещества, динамические реакции, область генерации.

## Где начинать отладку

Если не идет реакция:

1. Проверить, есть ли вещества в `Mixture`.
2. Проверить кандидатов через `reaction_candidates_for_substances`.
3. Проверить `context_allows_reaction`.
4. Проверить скорость в `reaction_rate_mol_per_bucket_per_tick_with_context`.
5. Проверить ограничение реагентами в `limit_by_reactants_and_context`.

Если не создается динамическая реакция:

1. Проверить структуру вещества через `resolve_frowns`.
2. Проверить найденные функциональные группы в `molecule/functional_group.rs`.
3. Проверить область генерации в `generate_reactions_from_scope`.
4. Проверить конкретный генератор в `organic/mod.rs`.
