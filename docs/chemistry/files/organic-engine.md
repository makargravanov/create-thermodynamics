# Движок генерации органических реакций

Исходный код: `organic/engine.rs`, `organic/mod.rs`

## Назначение

Главная точка входа органического модуля. Принимает реестр зарегистрированных веществ, обходит их реактивные сайты и порождает полный набор реакций органической химии, включая сгорание. Результат — [[organic-generators|GeneratedOrganicCatalog]] (новые вещества + реакции), который вливается в [[core-registry|ChemistryRegistry]].

## Ключевые типы

| Тип | Откуда | Роль |
|-----|--------|------|
| `GeneratedOrganicCatalog` | [[organic-aromatic-catalog]] | Контейнер результата |
| `OrganicGenerationSpace` | [[organic-space-resolver]] | Индекс сайтов по типу |
| `DerivedSubstanceResolver` | [[organic-space-resolver]] | Разрешение/создание промежуточных веществ |
| `SiteParticipant` | [[organic-space-resolver]] | Вещество + структура + один сайт |
| `SelectivityContext` | [[selectivity-types]] | Контекст для фильтрации |

## Публичные входы

```rust
pub fn destroy_registry_with_generated_reactions_builder()
    -> ChemistryResult<ChemistryRegistryBuilder>
```
Строит полный реестр: базовый реестр (`destroy_registry_builder`) плюс все сгенерированные органические реакции.

```rust
pub(crate) fn generate_organic_reactions(registry: &ChemistryRegistry)
    -> ChemistryResult<GeneratedOrganicCatalog>
```
Генерация по всему реестру (используется в тестах и внутри модуля).

```rust
pub(crate) fn generate_organic_reactions_for_seed_substances(
    space, seeds, canonical_to_id, context) -> ChemistryResult<GeneratedOrganicCatalog>
```
Генерация только для заданных seed-веществ (используется в динамическом движке).

## Поток данных / Алгоритм

```
registry
  └─ GenerationScope::all(registry)           // множество ID в scope
  └─ OrganicGenerationSpace::new(...)         // индекс SiteParticipant по ReactiveSiteKind
       └─ generate_organic_reactions_with_space
            └─ canonical_to_id_from_substances  // FROWNS → SubstanceId
            └─ generate_organic_reactions_for_seed_substances
                 ├─ Проход 1: generate_organic_reactions_for_seed_participants
                 │    └─ для каждого seed-сайта:
                 │         ├─ одиночные реакции (по типу сайта)
                 │         └─ парные реакции (seed × весь space)
                 ├─ generate_complete_combustion  // для каждого seed-вещества
                 └─ Проход 2: generate_site_reactions_for_seed_participants
                      └─ Enol, AromaticRing, Organomet., Epoxide, ArylHalide
```

### Одиночные реакции (по типу сайта)

| `ReactiveSiteKind` | Генераторы |
|--------------------|-----------|
| `Halide` | hydroxide-sub, ammonia-sub, cyanide-sub |
| `Alcohol` | oxidation, dehydration, thionyl-Cl-sub, silyl-protection |
| `Alkoxide` | protonation |
| `Nitrile` | hydrolysis, hydrogenation |
| `Nitro` | hydrogenation |
| `AcylChloride` | hydrolysis |
| `CarboxylicAcid` | acyl-chloride-formation |
| `Aldehyde/Ketone/Carbonyl` | oxidation, CN-addition, NaBH4-reduction, Wolff–Kishner |
| `Amide` | hydrolysis |
| `Ester` | hydrolysis, LAH-reduction |
| `PrimaryAmine` | phosgenation, Boc-protection, Cbz-protection |
| `NonTertiaryAmine` | cyanamide-addition, Boc/Cbz-protection |
| `PhosphoniumSalt` | ylide-formation |
| `PhosphorusYlide` | Wittig (× все карбонилы в space) |
| `PhosphonateCarbanion` | HWE (× все карбонилы) |
| `SulfoneCarbanion` | Julia (× все карбонилы) |
| `Isocyanate` | hydrolysis |
| `Borane` | oxidation |
| `BorateEster` | hydrolysis |
| `Alkene` | электрофильное присоединение (7 spec), Michael (× Enol) |
| `Alkyne` | электрофильное присоединение (7 spec) |
| `SilylEther` | deprotection |
| `Acetal/Ketal` | deprotection |
| `BocCarbamate` | deprotection |
| `CbzCarbamate` | deprotection |

### Парные реакции (seed × space)

Сайт seed тоже обходится в `generate_pair_reactions_for_seed`. Кросс-реакции:
CarboxylicAcid × Alcohol → esterification; Alcohol × AcylChloride → esterification,
Alcohol × Carbonyl → acetal; Halide × Enol → enolate-alkylation, Halide × Amine → sub,
Halide/AcylChloride × AromaticRing → FC-alkylation/acylation; Phosphine × Halide → phosphonium;
NonTertiaryAmine × Halide, Carbonyl → enamine; Ester × Enol → Claisen;
Carbonyl × Alcohol → acetal, × PrimaryAmine → imine; AromaticRing × Halide/AcylChloride → FC.

### Проход 2 (site-реакции)

Дополнительные реакции поверх уже найденных: органометаллическое присоединение к карбонилам
(Grignard/RLi/Gilman), альдольное присоединение, α-галогенирование, дегидратация альдола,
ароматическое нитрование/хлорирование/бромирование/сульфонирование, гидролиз эпоксидов,
нуклеофильная ароматическая замена арилгалогенидов.

## Инварианты и ошибки

- Реакции дедублируются по строковому ID через `BTreeSet<String>` (`push_unique_reaction`).
- Ошибки стерео-распределения типа "unknown stereo distribution" пропускаются (`continue`),
  не прерывая генерацию.
- `canonical_to_id_from_generated` обогащает карту FROWNS→ID после первого прохода,
  чтобы второй проход не дублировал уже сгенерированные вещества.
- Вещества только с молекулярной структурой попадают в `canonical_to_id`.

## Связи

- [[organic-space-resolver]] — `OrganicGenerationSpace`, `DerivedSubstanceResolver`, `SiteParticipant`
- [[organic-centers]] — типизированные сайты (HalideSite, AlcoholSite, ...)
- [[organic-generators]] — все функции `generate_*`
- [[organic-aromatic-catalog]] — `GeneratedOrganicCatalog`
- [[core-registry]] — `ChemistryRegistry`, `ChemistryRegistryBuilder`
- [[selectivity-engine]] — `SelectivityContext`, `SelectivityProfile`
