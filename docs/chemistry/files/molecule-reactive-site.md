# Реактивные центры

Исходный код: `molecule/reactive_site.rs`

## Назначение

Высокоуровневый слой поверх функциональных групп. Обогащает каждую группу ролями (`ReactiveRole`), определяет первичный атом, якорные атомы и уходящую группу. Является входным интерфейсом для движков органической химии и селективности.

## Ключевые типы

### ReactiveSiteKind

Надмножество `FunctionalGroupType` + дополнительные виды:

`Aldehyde`, `Ketone` (из `Carbonyl`), `AromaticCarbon`, `AromaticRing`, `ArylHalide`, `Azide`, `Diazonium`, `Enol`, `Enolate`, `Epoxide`, `Ether`, `Imine`, `Organocopper`, `Organolithium`, `Organomagnesium`, `Phenol`, `Sulfide`, `SulfonylChloride`, `Thiol`

### ReactiveRole

| Роль | Смысл |
|------|-------|
| `Electrophile` | принимает электроны |
| `Nucleophile` | отдаёт электроны |
| `AcidicProton` | способен к депротонированию |
| `LeavingGroup` | уходящая группа |
| `UnsaturatedBond` | π-система |
| `AlphaCarbon` | α-углерод при карбониле |
| `AromaticDirector` | управляет региоселективностью |
| `Oxidizable` / `Reducible` | окис./вос. |

### ReactiveSite

```rust
pub struct ReactiveSite {
    pub kind: ReactiveSiteKind,
    pub atoms: Vec<usize>,              // отсортированы, дедуплицированы
    pub roles: Vec<ReactiveRole>,       // отсортированы
    pub primary_atom: Option<usize>,    // главный атом сайта
    pub anchor_atoms: Vec<usize>,       // атомы реакционной связи
    pub leaving_atom: Option<usize>,    // уходящий атом (F/Cl/Br/I)
    pub bond_order: Option<i32>,        // порядок × 1000
    pub substitution_degree: Option<usize>,  // число C-соседей primary_atom
}
```

`ReactiveSiteKey` — ord/hash-образный снимок для дедупликации.

## Публичные входы

```rust
pub fn find_reactive_sites(structure: &MolecularStructure) -> Vec<ReactiveSite>
pub fn try_find_reactive_sites(structure: &MolecularStructure) -> ChemistryResult<Vec<ReactiveSite>>
```

`find_reactive_sites` паникует при невалидном сайте; `try_find_reactive_sites` возвращает ошибку.

## Поток данных / Алгоритм

```
find_functional_groups()
  ↓ маппинг FunctionalGroupType → (ReactiveSiteKind, roles)
remove_conflicting_protected_sites()   ← защищённые OH/NH не появляются
add_aromatic_sites()                   ← AromaticRing, AromaticCarbon, ArylHalide
add_oxygen_sites()                     ← Epoxide, Ether, Phenol
add_sulfur_sites()                     ← SulfonylChloride, SulfoneCarbanion, Thiol, Sulfide
add_nitrogen_sites()                   ← Diazonium, Azide, Imine
add_organometallic_sites()             ← Organomagnesium, Organolithium, Organocopper
add_alpha_sites()                      ← Enol (α-C при альдегиде/кетоне/эфире)
enrich_sites()                         ← primary_atom, anchor_atoms, leaving_atom, bond_order
deduplicate_sites()                    ← BTreeSet по ключу, сортировка
validate_against(structure)            ← проверка индексов и консистентности
```

### Логика обогащения (enrich_sites)

- `primary_atom` — первый атом в `atoms` для большинства сайтов; `None` для `AromaticRing`
- `anchor_atoms` — пара атомов связи для Alkene/Alkyne/Carbonyl/Enol/металлоорганики; все три атома для Epoxide; иначе список из `primary_atom`
- `leaving_atom` — первый F/Cl/Br/I в atoms (только если роль `LeavingGroup`)
- `bond_order` — порядок связи между первыми двумя anchor_atoms × 1000 (целое)
- `substitution_degree` — `carbon_degree(primary_atom)` если C

### Удаление конфликтов защитных групп

- SilylEther → соответствующий Alcohol по кислороду удаляется
- BocCarbamate/CbzCarbamate → PrimaryAmine/NonTertiaryAmine по азоту удаляются; Ester по карбонильному C удаляется

### Ароматические сайты

`AromaticRing` создаётся если ≥5 C с двумя ароматическими связями. Каждый такой C получает отдельный `AromaticCarbon` или `ArylHalide` (если есть одиночная связь с галогеном).

## Инварианты и ошибки

- `atoms` в каждом сайте отсортированы и дедуплицированы в конструкторе `ReactiveSite::new`.
- `leaving_atom` обязан быть в `atoms` и в ролях должна быть `LeavingGroup`.
- `bond_order` требует ≥2 `anchor_atoms`.
- `validate_against` проверяет все эти условия; нарушение → `ChemistryError`.
- Сайты дедуплицируются по `ReactiveSiteKey` — одна молекула, несколько путей распознавания, одна запись.

## Связи

- [[molecule-functional-group|Функциональные группы]] — источник сайтов
- [[molecule-graph|Граф молекулы]] — `MolecularStructure`, `bond_order_matches`
- [[selectivity-engine|Движок селективности]] — потребляет `ReactiveSite` для матчинга реакций
- [[organic-engine|Органический движок]] — использует сайты для генерации продуктов
