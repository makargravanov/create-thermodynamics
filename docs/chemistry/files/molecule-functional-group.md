# Функциональные группы

Исходный код: `molecule/functional_group.rs`

## Назначение

Поиск функциональных групп в молекуле по реальному графу с явными водородами. Никаких скрытых H — все водороды уже являются атомами, поэтому проверки вида «атом O имеет H-сосед» работают непосредственно через `neighbors()`.

## Ключевые типы

### FunctionalGroupType

Перечисление распознаваемых групп:

| Группа | Структурный паттерн |
|--------|---------------------|
| `CarboxylicAcid` | C(=O)O-H |
| `Ester` | C(=O)-O-C |
| `AcidAnhydride` | C(=O)-O-C(=O) |
| `AcylChloride` | C(=O)-Cl |
| `UnsubstitutedAmide` | C(=O)-N(H)(H) |
| `Carbonyl` | C=O (альдегид или кетон по `is_ketone`) |
| `Alcohol` | C-O-H |
| `Alkoxide` | C-O⁻ |
| `Halide` | C-X (X = Cl, Br, I) |
| `PrimaryAmine` | C-N(H)(H) |
| `NonTertiaryAmine` | C-N-H |
| `Alkene` | C=C |
| `Alkyne` | C≡C |
| `Nitrile` | C≡N |
| `Nitro` | C-N(~O)(~O) |
| `Isocyanate` | C-N=C=O |
| `Borane` / `NonTertiaryBorane` | C-B / C-B-H |
| `BoricAcid` | B-O-H |
| `BorateEster` | C-O-B |
| `Phosphine` | P(C)(C)(C) нейтральный |
| `PhosphonateCarbanion` | P(=O)(C⁻) |
| `PhosphoniumSalt` | P⁺-C-H |
| `PhosphorusYlide` | P⁺-C⁻ |
| `SulfoneCarbanion` | S(=O)(=O)-C⁻ |
| `SilylEther` | C-O-Si |
| `Acetal` | C(OR)(OR)(H) |
| `Ketal` | C(OR)(OR)(C)(C) |
| `BocCarbamate` | N-C(=O)-O-C(CH₃)₃ |
| `CbzCarbamate` | N-C(=O)-O-CH₂-Ar |

### FunctionalGroup

```rust
pub struct FunctionalGroup {
    pub group_type: FunctionalGroupType,
    pub atoms: Vec<usize>,       // индексы атомов группы
    pub degree: Option<usize>,   // степень замещения (C-атом)
    pub is_ketone: Option<bool>, // для Carbonyl: true=кетон, false=альдегид
}
```

## Публичные входы

```rust
pub fn find_functional_groups(structure: &MolecularStructure) -> Vec<FunctionalGroup>
```

## Поток данных / Алгоритм

Функция выполняет три независимых прохода по атомам:

### Проход 1: атомы углерода

Для каждого C (не в `carbonyl_carbons_to_ignore`) анализируются соседи:

1. **Карбонильный путь** (`carbonyl_oxygens.len() == 1`):
   - + одиночный O → Ester или AcidAnhydride (второй карбонил добавляется в `ignore`)
   - + O-H → CarboxylicAcid
   - + N(H)(H) → UnsubstitutedAmide
   - + Cl → AcylChloride
   - + 2 соседних C → Carbonyl(ketone=true)
   - иначе → Carbonyl(ketone=false)

2. **Не-карбонильный путь**: галогены, спирты, алкоксиды, амины, нитрилы, бораны

3. **Алкены/алкины:** дедупликация через `alkenes_to_ignore` / `alkynes_to_ignore`; каждая связь C=C записывается дважды если обе степени равны (позволяет реакциям выбирать любой конец).

### Проход 2: атомы бора → BoricAcid

### Проход 3: атомы фосфора → Phosphine / PhosphonateCarbanion / PhosphoniumSalt / PhosphorusYlide

### Защитные группы (`add_protecting_groups`)

Отдельный подпроход:
- **SilylEther:** O связан с C и Si
- **BocCarbamate:** N-C(=O)-O-C-(CH₃)₃ — три метила проверяются по числу H-соседей (3 явных H у каждого)
- **CbzCarbamate:** N-C(=O)-O-CH₂-C(ар.), ароматический C определяется через порядок связи 1.5
- **Acetal/Ketal:** C с двумя O-соседями без H на кислороде

## Инварианты и ошибки

- Функция чисто вычислительная, ошибок не бросает.
- Водороды ОБЯЗАНЫ быть явными атомами: `hydrogen_count(oxygen) == 1` проверяет H-соседа напрямую.
- Нитро-группа определяется через ароматические связи N-O (order 1.5), а не обычные двойные.
- `degree` для галогенидов/спиртов/аминов = число C-соседей атома углерода.

## Связи

- [[molecule-graph|Граф молекулы]] — `MolecularStructure`, `bond_order_matches`, `hydrogen_count`
- [[molecule-reactive-site|Реактивные центры]] — потребляет результат `find_functional_groups`
- [[organic-centers|Органические центры]] — использует группы для определения реакционных центров
