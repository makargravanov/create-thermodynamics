# Каталог веществ Destroy

Исходный код: `data/catalog.rs`

## Назначение

Статический каталог всех веществ мода Destroy — единственная точка сборки
[[core-registry|реестра химии]]. Описывает 170 явных веществ (`DESTROY_SUBSTANCES`),
металлургические ионы / металлы / материалы, теги веществ, таблицы фаз
(растворимость газов, смешиваемость растворителей), кислотно-основные пары,
равновесия, редокс-полуреакции и координационные комплексы.

## Структура данных / Ключевые типы

### RawSubstance

Промежуточный строительный блок для большинства органических и неорганических веществ.

| Поле | Тип | Смысл |
|---|---|---|
| `id` | `&str` | ключ без префикса (итог: `destroy:<id>`) |
| `structure_code` | `Option<&str>` | нативный SMILES-подобный код (`destroy:linear:…`, `destroy:benzene:…` и т.д.) |
| `java_structure_code` | `Option<&str>` | легаси-код Java API (`LegacyMolecularStructure.atom(…)`) |
| `boiling_point_celsius` / `boiling_point_kelvin` | `Option<f64>` | предпочтителен °C; К — для криогенных газов |
| `density` | `Option<f64>` | г/ведро; дефолт 1000 |
| `molar_heat_capacity` / `specific_heat_capacity` | `Option<f64>` | одно из двух |
| `color_argb` | `u32` | 0 → дефолт `0x20FFFFFF` |
| `tags` | `&[&str]` | теги без префикса (добавляется `destroy:`) |

### RawMetallurgyIon / RawMetallurgyMetal / RawMetallurgyMaterial

Три отдельных типа для металлургической ветки.

`RawMetallurgyMetal` — чистый металл; хранит `solid_density`, `melting_point_kelvin`,
`fusion_heat`, `vaporization_heat`. Фаза жидкости — `LiquidPhasePreference::MoltenMetal`.

`RawMetallurgyMaterial` — оксид или ионный твёрдый (сульфид, карбонат, силикат).
Молярная масса рассчитывается через `material_formula_mass` по формульным единицам.
Фаза жидкости — `LiquidPhasePreference::MoltenSlag`.

`RawMetallurgyIon` — заряженный металлический ион в растворе (напр. `Fe²⁺`, `Al³⁺`);
нет твёрдой фазы, нет кипения (`boiling = f64::MAX`).

## Публичные входы

### `destroy_substances_registry_builder() -> ChemistryResult<ChemistryRegistryBuilder>`

Главная функция сборки: регистрирует теги → явные вещества →
металлургические вещества → таблицы фаз (катализаторы, комплексы, растворимость газов,
смешиваемость растворителей, равновесия, кислотно-основные пары, редокс-полуреакции).

Результат передаётся в `destroy_reactions_registry_builder` (см. [[data-reactions]]),
после чего вызывается `ChemistryRegistryBuilder::build()`.

### `summarize_legacy_structure(structure_code) -> ChemistryResult<MolecularSummary>`

Публичная утилита: парсит легаси-код структуры и возвращает молярную массу + заряд
(используется в JNI-слое для валидации Java-сторон).

## Категории данных

### Неорганические вещества (~30 ед.)

Простые газы: `oxygen`, `nitrogen`, `hydrogen`, `carbon_dioxide`, `carbon_monoxide`,
`sulfur_dioxide`, `nitrogen_dioxide`, `ammonia`, `chlorine` и др.

Кислоты и основания: `sulfuric_acid`, `nitric_acid`, `hydrochloric_acid`,
`hydrofluoric_acid`, `hydrogen_cyanide`, `boric_acid`, `hypochlorous_acid`, `oleum` и др.

Ионы: `proton`, `hydroxide`, `chloride`, `fluoride`, `sulfate`, `hydrogensulfate`,
`nitrate`, `oxide`, `sulfide`, `carbonate`, `silicate`, `borohydride`,
`tetrahydroxyborate`, `nitronium`, `cyanide`, `hypochlorite`, `iodide`, `chromate`,
`dichromate`, `chloroaurate`, `cyanamide_ion`, `ammonium`, `acetate` и др.

Переходные металлы (ионы в растворе): `iron_ii`, `iron_iii`, `copper_i`, `copper_ii`,
`nickel_ion`, `zinc_ion`, `lead_ii`, `chromium_iii`, `calcium_ion`, `potassium_ion`,
`sodium_ion`.

Прочие: `water`, `mercury`, `iodine`, `sodium_metal`, `argon`, `octasulfur`,
`hydrogen_peroxide`, `cisplatin`, `tetraethyllead`.

### Органика — алифатические (~45 ед.)

Одноатомные/лёгкие: `methane`, `ethene`, `acetylene`, `propene`, `butadiene`,
`isoprene`, `methanol`, `ethanol`, `isopropanol`, `acetone`, `glycerol`,
`acetic_acid`, `acetic_anhydride`, `acetamide` и др.

Галогенированные: `chloromethane`, `dichloromethane`, `chloroform`, `carbon_tetrachloride`,
`chloroethane`, `chloroethene`, `chlorodifluoromethane`, `dichlorodifluoromethane`,
`trichlorofluoromethane`, `tetrafluoroethene`, `iodomethane`, `bromine`,
`hydrobromic_acid`, `hydroiodic_acid`, `hydrogen_iodide`.

Нитро и взрывчатые: `nitroglycerine`, `hydrazine`, `dinitrotoluene`, `tnt`,
`picric_acid`, `aibn`.

Прочие функциональные группы: `acrylonitrile`, `adipic_acid`, `adiponitrile`,
`hexanediamine`, `hexane_diisocyanate`, `vinyl_acetate`, `methyl_methacrylate`,
`methyl_acetate`, `methylamine`, `trimethylamine`, `cyanamide`, `creatine`.

Бор-органика: `diborane`, `trimethyl_borate`, `trimethylphosphine`.

### Органика — ароматические (~20 ед.)

`benzene`, `toluene`, `ethylbenzene`, `styrene`, `orthoxylene`, `metaxylene`,
`paraxylene`, `phenol`, `aspirin`, `salicylic_acid`, `phenylacetic_acid`,
`phenylacetone`, `methyl_salicylate`, `cubane`, `cubanedicarboxylic_acid`,
`phthalic_anhydride`, `ethylanthraquinone`, `ethylanthrahydroquinone`,
`benzyl_chloride`.

### Generics — шаблонные вещества (~17 ед.)

Вещества с `id` вида `generic_*` (напр. `generic_alcohol`, `generic_alkene`,
`generic_amine`, `generic_nitrile`, `generic_nitro` и т.д.) — не участвуют в реакциях
напрямую, служат образцами для органического движка ([[organic-engine]]) и JEI-отображения.

### Защитные группы и реагенты (~6 ед.)

`trimethylsilyl_chloride`, `trimethylsilyl_fluoride`, `di_tert_butyl_dicarbonate`,
`benzyl_chloroformate`, `tert_butanol`, `acetyl_chloride`.

Катализаторы Фриделя–Крафтса: `ferric_chloride`, `ferric_bromide`, `aluminum_trichloride`.

### Металлургические ионы (5)

`aluminum_iii`, `magnesium_ion`, `silicon_iv`, `carbonate`, `silicate`.
Регистрируются отдельно через `METALLURGY_IONS`.

### Металлы — расплавы (8)

Fe, Cu, Zn, Ni, Pb, Al, Ca, Mg — через `METALLURGY_METALS`.
Каждый имеет точку плавления, теплоту плавления/испарения,
фазу `MoltenMetal`, представление `SubstanceRepresentation::Metal`.

### Металлургические материалы — оксиды и сульфиды/карбонаты/силикаты (24)

**Оксиды:** FeO, Fe₂O₃, Fe₃O₄ (магнетит), Cu₂O, CuO, ZnO, NiO, PbO, Al₂O₃, CaO,
MgO, SiO₂ (кремнезём).

**Сульфиды:** FeS, Cu₂S, CuS, ZnS, PbS, NiS.

**Карбонаты:** CaCO₃, MgCO₃, ZnCO₃, FeCO₃, CuCO₃, PbCO₃.

**Силикаты:** Ca₂SiO₄, Mg₂SiO₄, Fe₂SiO₄.

Все имеют фазу `MoltenSlag`, представление `Oxide` или `IonicSolid`.

## Таблицы фаз

### Катализаторы (`register_phase_tables`)

Металлические порошки как контактные катализаторы через `CatalystSurfaceSpec::chemical`:
Ni, Pd, Pt, Rh, Fe, Cu, Zn (дублированные записи forge- и Destroy-тегов).
Также `DestroyItems.ZEOLITE` и `DestroyItems.SILICA` через `unchecked`.

### Координационные комплексы (5)

| Комплекс | Центр | Лиганды | Геометрия | lg Kf |
|---|---|---|---|---|
| `copper_ii_tetraammine` | Cu²⁺ | 4 × NH₃ | квадратная плоская | 13 |
| `nickel_tetraammine` | Ni²⁺ | 4 × NH₃ | тетраэдр | 8 |
| `ferric_hexacyanide` | Fe³⁺ | 6 × CN⁻ | октаэдр | 31 |
| `cuprous_dicyanide` | Cu⁺ | 2 × CN⁻ | линейный | 24 |
| `zinc_tetraammine` | Zn²⁺ | 4 × NH₃ | тетраэдр | 9 |
| `ferric_tetrachloride` | Fe³⁺ | 4 × Cl⁻ | тетраэдр | 2 |

### Растворимость газов (Henry, 10 газов)

O₂, N₂, CO, CO₂, Cl₂, NH₃, HCl, SO₂, NO₂, H₂.
Каждый описан `GasSolubilityModel::Henry` с константой Генри, поправкой на
«высаливание» и коэффициентом массопереноса. Все помечены `estimated: true`.

### Смешиваемость растворителей

`water ↔ ethanol` — полная; `water ↔ acetone` — полная;
`water ↔ chloroform` — частичная (0.1 моль/ведро).

### Кислотно-основные пары (9)

`acetic_acid/acetate` (pKa 4.76), `ammonium/ammonia` (9.25),
`HCl/Cl⁻` (−6.3), `HF/F⁻` (3.17), `HCN/CN⁻` (9.21),
`HI/I⁻` (−10), `HSO₄⁻/SO₄²⁻` (1.99), `HClO/ClO⁻` (7.53),
`HNO₃/NO₃⁻` (−1.4), `H₂SO₄/HSO₄⁻` (−3.0).

Также особое равновесие: авто-ионизация воды (Kw = 10⁻¹⁴),
гидролиз борной кислоты (Ka = 10⁻⁹·²⁴).

### Редокс-полуреакции (9 пар)

Fe²⁺/Fe³⁺ (E°=0.771 В), Cu⁺/Cu²⁺ (0.153 В), I⁻/I₂ (0.535 В),
H₂/H⁺ (0.0 В), O₂/H₂O (1.229 В), H₂O₂/H₂O (1.776 В),
H₂O₂/O₂ (−0.682 В), Cl₂/Cl⁻ (1.358 В), HClO/Cl⁻ (1.49 В),
Cr₂O₇²⁻/Cr³⁺ (1.33 В).

## Вспомогательные функции

`estimate_boiling_point(M)` — линейная аппроксимация по молярной массе
(`Tbp = 2.04·M + 178.2 K`) для веществ без явной точки кипения.

`estimate_phase_properties` — эвристика фазовых свойств:
`water` → растворитель; ионы (charge ≠ 0) → водная фаза, может осаждаться;
тег `solvent` → органический растворитель без предела; кислоты/аммиак → водная фаза
с ограниченной органической растворимостью; остальные → органические растворители.

## Инварианты и ошибки

- Каждое `RawSubstance` обязано иметь хотя бы один из `structure_code` / `java_structure_code`.
- Молярная масса вычисляется парсером структуры; невалидная масса (≤ 0 или NaN) → `InvalidSubstance`.
- `DESTROY_SUBSTANCE_COUNT` (170) и итоговый `registry.substance_count()` (224) проверяются в тесте.
  Разница — вещества, порождённые неявно (металлы + ионы + материалы = 54).

## Связи

- [[core-registry]] — `ChemistryRegistryBuilder`, куда регистрируются все данные
- [[core-substance]] — тип `Substance`, `SubstancePhaseProperties`, `SubstanceRepresentation`
- [[core-complex]] — `ComplexSpec`, `ComplexLigand`, `ComplexGeometry`
- [[core-redox]] — `RedoxHalfReaction`, `RedoxEnvironment`
- [[core-solution]] — `AcidBaseSpec`, `EquilibriumSpec`
- [[core-catalysis]] — `CatalystSurfaceSpec`
- [[data-reactions]] — функция, принимающая готовый builder
- [[molecule-graph]] — парсинг `structure_code` / `java_structure_code`
