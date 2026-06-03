# Каталитические поверхности (catalysis)

Исходный код: `core/catalysis.rs`

## Назначение

Описывает гетерогенный катализ: каталитические поверхности с именованными сайтами,
их состояние (свободные / занятые / отравленные), а также шаги реакции на поверхности
(адсорбция, десорбция, отравление, восстановление).

## Ключевые типы

| Тип | Роль |
|-----|------|
| `CatalystSurfaceId` | строковый идентификатор поверхности, не может быть пустым |
| `SurfaceSiteId` | идентификатор типа активного сайта; дефолт `"default"` |
| `CatalystSurfaceSpec` | описание поверхности: масса, заряд, фазы, сайт; два конструктора: `chemical` и `unchecked` |
| `CatalystSurfaceState` | рантайм-состояние: общее число сайтов, занятые сайты по `SurfaceSiteId`, отравленные |
| `SurfaceRequirement` | требование реакции к поверхности: сайт, фазы, сайтов на оборот |
| `SurfaceStep` | шаг реакции на поверхности: `Adsorb / Desorb / Poison / Restore` |

## Публичные входы

- `CatalystSurfaceSpec::chemical(id, mass, charge)` — поверхность с проверяемой массой
- `CatalystSurfaceSpec::unchecked(id, reason)` — поверхность без проверки массы (напр. металл-носитель)
- `CatalystSurfaceSpec::validate()` — проверяет согласованность (масса+заряд, ≥1 фаза)
- `CatalystSurfaceState::new(total_sites)` — создаёт состояние; `total_sites ≥ 0`
- `CatalystSurfaceState::free_sites()` — `total − occupied − poisoned`, не ниже 0
- `CatalystSurfaceState::validate()` — инварианты состояния (численные и ограничения)
- `SurfaceRequirement::validate(reaction_id)` — `sites_per_reaction > 0`, ≥1 фаза
- `SurfaceStep::validate(reaction_id)` — ненулевые id и `sites > 0`

## Поток данных / Алгоритм

### Жизненный цикл сайтов

```mermaid
stateDiagram-v2
    [*] --> Free : new(total_sites)
    Free --> Occupied : SurfaceStep::Adsorb
    Occupied --> Free : SurfaceStep::Desorb
    Occupied --> Poisoned : SurfaceStep::Poison
    Poisoned --> Free : SurfaceStep::Restore
```

`CatalystSurfaceState` хранит `occupied_sites: BTreeMap<SurfaceSiteId, f64>`, что позволяет
нескольким типам реакций занимать разные типы сайтов на одной поверхности одновременно.

`free_sites = total − Σ occupied[site] − poisoned`

### Связь с кинетикой

`ChannelConditionEffect::Surface { surface_id, multiplier }` в [[core-kinetics|kinetics]] читает
`ReactionContext::free_sites(surface_id)` — если свободных сайтов ≤ `TRACE`, канал гасится.

`SurfaceRequirement` используется [[core-simulation|simulation]] для проверки доступности
поверхности перед шагом реакции.

### Два конструктора `CatalystSurfaceSpec`

- `chemical(id, mass, charge)` — масса и заряд проверяются при валидации реестра; масса
  используется при проверке сохранения в шагах адсорбции/десорбции.
- `unchecked(id, reason)` — сохраняет строку-причину; масса и заряд не задаются; используется
  для пористых носителей и твёрдых фаз, где атомарная масса активного центра не определена.

## Инварианты и ошибки

| Условие | Ошибка |
|---------|--------|
| Пустой `CatalystSurfaceId` или `SurfaceSiteId` | `InvalidReaction` |
| Пустой список `accessible_phases` | `InvalidReaction` |
| Масса отрицательна или не конечна (`chemical`) | `InvalidReaction` |
| Ни масса+заряд, ни `unchecked_mass_reason` не заданы | `InvalidReaction` |
| `total_sites < 0` или нечисловые | `InvalidMixtureState` |
| `occupied + poisoned > total + TRACE` | `InvalidMixtureState` |
| `sites_per_reaction ≤ 0` в `SurfaceRequirement` | `InvalidReaction` |

## Связи

- [[core-kinetics|kinetics]] — `CatalystSurfaceId` в `ChannelConditionEffect::Surface` и `IsomerEnergy`
- [[core-simulation|simulation]] — обрабатывает `SurfaceStep` и `SurfaceRequirement`
- [[core-registry|ChemistryRegistry]] — регистрирует `CatalystSurfaceSpec`
- [[core-error|ChemistryError]]
