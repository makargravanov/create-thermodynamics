# Карта кода

## Главный модуль

`mod.rs`

Собирает публичные модули и дает главные входы:

- `destroy_registry_builder`
- `destroy_registry_with_generated_reactions_builder`

## Данные Destroy

- `catalog.rs`: вещества Destroy.
- `destroy_reactions.generated.rs`: перенесенные явные реакции.
- `reactions.rs`: подключение сгенерированного файла.

## Основная модель

- `substance.rs`: вещество и теги.
- `reaction.rs`: реакция, условия, построитель.
- `registry.rs`: проверенный реестр, индексы реакции.
- `mixture.rs`: смесь, нагрев, кипение, смешивание.
- `simulation.rs`: расчет реакции за тик и до равновесия.

## Молекулярный слой

- `molecule.rs`: граф, редактор графа, разбор старых структур.
- `frowns.rs`: разбор и запись FROWNS.
- `functional_group.rs`: поиск функциональных групп.

## Органическая динамика

- `organic.rs`: генераторы органических реакций.
- `dynamic.rs`: динамические вещества, динамические реакции, область генерации.

## Где начинать отладку

Если не идет реакция:

1. Проверить, есть ли вещества в `Mixture`.
2. Проверить кандидатов через `reaction_candidates_for_substances`.
3. Проверить `context_allows_reaction`.
4. Проверить скорость в `reaction_rate_mol_per_bucket_per_tick_with_context`.
5. Проверить ограничение реагентами в `limit_by_reactants_and_context`.

Если не создается динамическая реакция:

1. Проверить структуру вещества через `resolve_frowns`.
2. Проверить найденные функциональные группы в `functional_group.rs`.
3. Проверить область генерации в `generate_reactions_from_scope`.
4. Проверить конкретный генератор в `organic.rs`.
