# Реакции Destroy

Исходный код: `data/reactions.rs`, `data/destroy_reactions.generated.rs`

## Назначение

`data/reactions.rs` состоит из единственной строки:

```rust
include!("destroy_reactions.generated.rs");
```

Весь контент — в сгенерированном файле `destroy_reactions.generated.rs`.
Файл содержит функцию `destroy_reactions_registry_builder`, регистрирующую все
явные химические реакции Destroy в [[core-registry|реестре химии]].

## Структура данных / Ключевые типы

Используется `Reaction::builder(id)` — паттерн строитель из [[core-reaction]].

Доступные вызовы на строителе:

| Метод | Смысл |
|---|---|
| `.reactant(id, stoich, catalyst_order)` | вещество-реагент; `catalyst_order` = 0 → катализатор |
| `.product(id, stoich)` | продукт реакции |
| `.catalyst_order(id, order)` | жидкофазный катализатор |
| `.external_reactant(java_expr, amount)` | реагент-твердотело (предмет из инвентаря), строка = Java-выражение |
| `.external_catalyst(java_expr, amount)` | твёрдотельный катализатор |
| `.pre_exponential_factor(A)` | фактор Аррениуса A |
| `.activation_energy_kj_per_mol(Ea)` | энергия активации |
| `.enthalpy_change_kj_per_mol(ΔH)` | тепловой эффект |
| `.requires_uv()` | требует УФ-облучения |
| `.reaction_result(java_expr, amount)` | Java-выражение побочного результата (осадок, предмет) |
| `.reverse_reaction_id(id)` | связь с обратной реакцией |
| `.display_as_reversible()` | показать в JEI как обратимую |
| `.show_in_jei(bool)` | спрятать из JEI (напр. авто-сгенерированные обратные) |
| `.redox_annotation(…)` | проверяемая/непроверяемая редокс-аннотация |
| `.allow_mass_imbalance()` | разрешить несбалансированность по массе |

## Публичные входы

### `destroy_reactions_registry_builder(builder) -> ChemistryResult<ChemistryRegistryBuilder>`

Принимает уже заполненный веществами `ChemistryRegistryBuilder` (из
`destroy_substances_registry_builder` в [[data-catalog]]) и добавляет реакции.
Возвращает builder; вызывающий код затем вызывает `.build()`.

### Константы

```rust
pub const DESTROY_EXPLICIT_REACTION_COUNT: usize = 120;
pub const DESTROY_REGISTERED_REACTION_COUNT: usize = 157;
```

`DESTROY_EXPLICIT_REACTION_COUNT` — число уникальных реакций верхнего уровня
(без авто-генерированных кислотно-основных пар).

`DESTROY_REGISTERED_REACTION_COUNT` — итоговое число записей в реестре,
включая обратные реакции и диссоциации кислот/оснований, сгенерированных
автоматически из `AcidBaseSpec` / `EquilibriumSpec` при сборке.

## Категории реакций

### Промышленные неорганические процессы

| Реакция | Описание |
|---|---|
| `haber_process` | N₂ + 3H₂ → 2NH₃; кат. Fe |
| `contact_process` | 2SO₂ + O₂ → 2SO₃; кат. Pt |
| `ostwald_process` | NH₃ + O₂ → HNO₃ + H₂O; кат. Rh |
| `andrussow_process` | CH₄ + NH₃ + O₂ → HCN + H₂O; кат. Pt |
| `steam_reformation` | CH₄ + H₂O ⇌ CO + 3H₂; кат. Ni |
| `methanol_synthesis` | CO + 2H₂ → CH₃OH; кат. Cu/Zn |
| `oleum_formation` / `.reverse` | H₂SO₄ + SO₃ ⇌ олеум |
| `sulfur_trioxide_hydration` / `.reverse` | SO₃ + H₂O ⇌ H₂SO₄ (ΔH = −200 кДж/моль) |

### Сера и хлор

`sulfur_oxidation`, `thionyl_chloride_synthesis`,
`phosgene_formation` (CO + Cl₂ → COCl₂; ΔH = −107.6),
`chlorine_solvation` / `.reverse` (требует УФ),
`hydrogen_chloride_synthesis` (H₂ + Cl₂ → 2HCl; УФ).

### Хлорирование метана (УФ-цепочка, 4 реакции)

`methane_uv_chlorination` → `chloromethane_uv_chlorination` →
`dichloromethane_uv_chlorination` → `chloroform_uv_chlorination`.
Каждый шаг требует УФ, Ea растёт от 22.5 до 30 кДж/моль.

### Фреоны и галогенированные углеводороды

`carbon_tetrachloride_fluorination` (CCl₄ + HF → CCl₂F₂ + CCl₃F + HCl),
`chloroform_fluorination` (CHCl₃ + 2HF → CHClF₂ + 2HCl),
`chlorodifluoromethane_pyrolysis` (2CHClF₂ → 2HCl + C₂F₄).

### Растворение металлов и руд

Попарные реакции (пыль / дроблёная руда) для Fe, Cu, Ni, Zn, Pb, Cr:
`<metal>_dissolution` (форж-тег) + `<metal>_ore_dissolution` (предмет).
Общая схема: 2H⁺ + Me → H₂ + Me²⁺/³⁺.

Также `lime_slaking` (CaO + H₂O → Ca²⁺ + 2OH⁻),
`fluorite_dissolution`, `chromium_ore_dissolution`, `crocoite_dissolution`,
`borax_dissolution`, `gold_dissolution` (царская водка → [AuCl₄]⁻).

### Хром

`chromium_iii_oxidation` (Cr³⁺ → CrO₄²⁻; H₂O₂, щелочь; редокс-аннотация),
`chromate_conversion` (2CrO₄²⁻ + 2H⁺ ⇌ Cr₂O₇²⁻ + H₂O).

### Бор

`diborane_hydrolysis` (B₂H₆ + 6H₂O → 2B(OH)₃ + 6H₂; ΔH = −466),
`brown_schlesinger_process` (NaBH₄ синтез),
`borohydride_iodine_oxidation`,
`borax_dissolution` / `borax_precipitation`,
`boric_acid_neutralization`,
`tetraborate_equilibrium` / `.reverse`.

### Азот и нитро-химия

`nitronium_formation` (HNO₃ + H₂SO₄ → NO₂⁺ + H₂O + HSO₄⁻),
`toluene_nitration` → `dinitrotoluene_nitration` → TNT,
`phenol_nitration` → пикриновая кислота,
`glycerol_nitration` → нитроглицерин,
`cellulose_nitration` → нитроцеллюлоза.

`haber_process`, `peroxide_process` (гидразин),
`frank_caro_process` (Ca₂C₂ + N₂ → цианамид),
`andrussow_process`, `ostwald_process`.

### Ксилолы и алкилбензолы (зеолитный цикл)

`toluene_transalkylation` (8Tol → 4Bz + xyl·3 + EB),
`ethylbenzene_transalkylation` / `metaxylene_transalkylation` /
`orthoxylene_transalkylation` / `paraxylene_transalkylation` — все через DestroyItems.ZEOLITE.
`benzene_ethylation` (Bz + этен → EB; H⁺),
`ethylbenzene_dehydrogenation` (EB → стирол + H₂; Fe³⁺/H₂O).

### Ароматика и красители

`orthoxylene_oxidation` → фталевый ангидрид (кат. Hg),
`ethylanthraquinone_synthesis` (ZEOLITE),
`cumene_process` (бензол + пропен + O₂ → фенол + ацетон),
`kolbe_schmitt_reaction` (CO₂ + фенол → салициловая кислота; Na⁺/H⁺).

### Антрахиноновый цикл (H₂O₂)

`anthraquinone_reduction` (этилантрахинон + H₂ → гидрохинон; Pd) →
`anthraquinone_process` (гидрохинон + O₂ → антрахинон + H₂O₂).

### Полимеризация (AIBN-катализ)

Аддитивные полимеры через `aibn` как катализатор нулевого порядка:
ABS, PAN, PVC, PE, PP, PTFE, полистирол, полистирол-бутадиен (SBR),
полиизопрен, PMMA.

Конденсационные: `nylon_polymerisation` (адипиновая кислота + гексаметилендиамин),
`urethane_hdi_polymerization` / `urethane_tdi_polymerization` (полиуретан).

### Натрий

`sodium_dissolution`, `sodium_ingot_dissolution`, `sodium_oxide_dissolution`,
`sodium_hydride_formation` / `sodium_hydride_hydrolysis`,
`sodium_amalgamization` / `.reverse` (ртутный цикл).

### Кислотно-основные авто-реакции (Hidden, show_in_jei: false)

Автоматически генерируются из `AcidBaseSpec` для каждой кислоты 3 реакции:
`dissociation` (кислота + H₂O → H⁺ + основание),
`neutralization` (кислота + OH⁻ → основание + H₂O),
`association` (основание + H⁺ → кислота).

Параметры рассчитываются из pKa: общая Ea = 2.477 кДж/моль.
Эти 37 реакций (9 пар × 3 + авто-ионизация) составляют большую часть
разрыва 120 → 157.

### Прочие

`baby_blue_precipitation` (метилсалицилат + Na⁺),
`creatine_precipitation`, `nhn_synthesis` (никелевый гидразинат),
`cordite_precipitation` (ацетон + нитроглицерин + нитроцеллюлоза),
`tatp` (ацетон + H₂O₂ → ацетонпероксид),
`mercury_fulmination`, `cisplatin_synthesis`,
`touch_powder_synthesis`, `kelp_dissolution`,
`carbon_capture` (Ca²⁺ + CO₂ → CaCO₃↓),
`tetraethyllead_synthesis`,
`naughty_reaction` (show_in_jei: false — пасхальное яйцо).

## Инварианты и ошибки

- Все реакции используют `.allow_mass_imbalance()` — Rust-движок не проверяет
  стехиометрический баланс; ответственность на авторе реакции.
- Обратная реакция, скрытая от JEI (`show_in_jei(false)`), сохраняет
  `.reverse_reaction_id` для двунаправленной кинетики.
- `RedoxAnnotation::checked` требует, чтобы заявленный перенос электронов
  соответствовал полуреакциям в каталоге; `unchecked` отключает проверку.

## Формат генерации

Файл `destroy_reactions.generated.rs` генерируется внешним скриптом из
Java-источников мода. Структура неизменна:

```
use ...;
pub const DESTROY_EXPLICIT_REACTION_COUNT: usize = N;
pub const DESTROY_REGISTERED_REACTION_COUNT: usize = M;
pub fn destroy_reactions_registry_builder(...) -> ... {
    builder = builder.reaction(Reaction::builder("destroy:…") … .build());
    …
    Ok(builder)
}
```

`reactions.rs` включает файл через `include!`, делая всё содержимое
частью модуля `data`.

## Связи

- [[core-reaction]] — `Reaction`, строитель, кинетические параметры
- [[core-registry]] — `ChemistryRegistryBuilder::reaction()`
- [[core-redox]] — `RedoxAnnotation`, `RedoxEnvironment`
- [[data-catalog]] — `destroy_substances_registry_builder` — поставщик builder'а
- [[core-kinetics]] — обработка Ea, A, ΔH при симуляции
