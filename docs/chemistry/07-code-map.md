# Карта кода

Эта страница - указатель по нативной (Rust) части. Каждый файл химического ядра имеет подробную страницу в папке `files/`. Здесь дана структура и порядок чтения; детали - по ссылкам.

Корень крейта: `native/thermodynamics-jni/`. Пути в подстраницах даны относительно `src/chemistry/`, кроме `lib.rs`/`build.rs`.

## Главные входы

`chemistry/mod.rs` собирает публичные модули и даёт точки входа:

- `destroy_registry_builder`
- `destroy_registry_with_generated_reactions_builder`

JVM-слой вызывает ядро через JNI: см. [[jni-bindings|JNI-привязки]].

## Основная модель (core/)

- [[core-substance|Вещество]] — `core/substance.rs`: `Substance`, теги, фазовые свойства, давление пара.
- [[core-reaction|Реакция]] — `core/reaction.rs`: уравнение, условия, построитель, каналы, внешние участники.
- [[core-registry|Реестр]] — `core/registry.rs`: проверенный реестр, индексы веществ и реакций, кандидаты.
- [[core-mixture|Смесь]] — `core/mixture.rs`: фазы, концентрации, нагрев, кипение, растворимость, газы.
- [[core-simulation|Симуляция]] — `core/simulation.rs`: тик реакции, контекст, ограничение реагентами, равновесие.
- [[core-condition|Условия]] — `core/condition.rs`: pH, давление, температура как условия реакции.
- [[core-thermodynamics|Термодинамика]] — `core/thermodynamics.rs`: ΔG, константа равновесия, термодинамический фактор скорости.
- [[core-kinetics|Кинетика каналов]] — `core/kinetics.rs`: каналы реакции, энергетическая модель, свет, изомеры.
- [[core-catalysis|Катализ]] — `core/catalysis.rs`: поверхности, активные центры, адсорбция, отравление.
- [[core-complex|Комплексы]] — `core/complex.rs`: комплексообразование, лиганды, геометрия, устойчивость.
- [[core-solution|Раствор]] — `core/solution.rs`: равновесия раствора, осаждение, кислотно-основные пары.
- [[core-redox|Редокс]] — `core/redox.rs`: степени окисления, полуреакции, электролиз, редокс-среда.
- [[core-error|Ошибки]] — `error.rs`: общий тип `ChemistryError`.

## Молекулярный слой (molecule/)

Центральное правило слоя: водороды всегда явные атомы в графе.

- [[molecule-graph|Граф молекулы]] — `molecule/graph.rs`: `MolecularStructure`, `MolecularEditor`, проверки валентности и массы.
- [[molecule-frowns|FROWNS]] — `molecule/frowns.rs`: разбор и запись формата FROWNS, восстановление явных водородов.
- [[molecule-canonical|Канонизация]] — `molecule/canonical.rs`: канонический код структуры для дедупликации веществ.
- [[molecule-functional-group|Функциональные группы]] — `molecule/functional_group.rs`: поиск групп по фактическому графу.
- [[molecule-aromatic-perception|Ароматичность]] — `molecule/aromatic_perception.rs`: распознавание ароматических систем, правило Хюккеля.
- [[molecule-reactive-site|Реактивные центры]] — `molecule/reactive_site.rs`: отображение групп в реактивные центры.

## Органическая динамика

- [[dynamic|Динамический слой]] — `dynamic/mod.rs`: `DynamicChemistryRegistry`, разрешение веществ, генерация реакций, область генерации.
- [[organic-engine|Движок органики]] — `organic/engine.rs`, `organic/mod.rs`: два прохода генерации, диспетчеризация по типу центра.
- [[organic-centers|Центры органики]] — `organic/centers.rs`: типизированные реактивные центры, α-углерод.
- [[organic-space-resolver|Пространство и резолвер]] — `organic/space.rs`, `organic/resolver.rs`: область генерации, вывод производных веществ.
- [[organic-aromatic-catalog|Ароматика и каталог]] — `organic/aromatic.rs`, `organic/catalog.rs`: эффекты заместителей, EAS/SNAr.
- [[organic-generators|Генераторы реакций]] — `organic/generators/*`: обзор всех генераторов и создаваемых ими реакций.

## Селективность кинетики (kinetics/selectivity/)

- [[selectivity-types|Типы селективности]] — `kinetics/selectivity/types.rs`: дескрипторы центров, профили, контекст.
- [[selectivity-engine|Движок селективности]] — `kinetics/selectivity/engine.rs`: оценка профиля, стерика, электроника, редокс-среда.
- [[selectivity-reactions|Классы реакций]] — присоединение к карбонилу, элиминирование, этерификация, нуклеофильное замещение.

## Данные Destroy (data/)

- [[data-catalog|Каталог веществ]] — `data/catalog.rs`: вещества Destroy, теги, свойства, категории.
- [[data-reactions|Явные реакции]] — `data/reactions.rs`, `data/destroy_reactions.generated.rs`: перенесённые реакции, кодогенерация.

## Граница с JVM

- [[jni-bindings|JNI-привязки]] — `lib.rs`, `build.rs`: JNI-экспорты, ABI-версия, кодогенерация реакций при сборке.
- [[synthesis|Планировщик синтеза]] — `synthesis.rs`: поиск путей синтеза по пространству состояний.

## Где начинать отладку

Если не идёт реакция:

1. Проверить, есть ли вещества в [[core-mixture|Mixture]].
2. Проверить кандидатов через `reaction_candidates_for_substances` ([[core-registry|registry]]).
3. Проверить `context_allows_reaction` и условия ([[core-condition|condition]]).
4. Проверить скорость в `reaction_rate_mol_per_bucket_per_tick_with_context` ([[core-simulation|simulation]]).
5. Проверить ограничение реагентами в `limit_by_reactants_and_context`.
6. Если реакция поверхностная, проверить `ReactionContext::surfaces` ([[core-catalysis|catalysis]]).

Если не создаётся динамическая реакция:

1. Проверить структуру вещества через `resolve_frowns` ([[dynamic|dynamic]]).
2. Проверить найденные функциональные группы ([[molecule-functional-group|functional_group]]).
3. Проверить область генерации в `generate_reactions_from_scope`.
4. Проверить конкретный генератор ([[organic-generators|generators]]).

## Правило чтения кода

Сначала читай публичные входы и типы данных, потом внутренние проверки. В этой модели важнее поток данных, чем отдельные функции.
