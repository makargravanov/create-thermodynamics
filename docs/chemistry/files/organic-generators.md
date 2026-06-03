# Генераторы органических реакций

Исходный код: `organic/generators/mod.rs` и все подмодули

## Назначение

Каждый подмодуль реализует одну или несколько функций `generate_*`, принимающих типизированный реактивный центр из [[organic-centers]] и `DerivedSubstanceResolver` из [[organic-space-resolver]], и возвращающих `ChemistryResult<Reaction>` или `ChemistryResult<Option<Reaction>>`. Реакция строится через `Reaction::builder(id)` паттерном.

Общие утилиты (`common.rs`) используются всеми генераторами.

## Общие утилиты (`common.rs`)

| Функция | Назначение |
|---------|-----------|
| `add_hydroxyl(editor, parent)` | Добавляет −OH к parent |
| `add_addition_group(editor, parent, group)` | Добавляет Atom/Hydroxyl/Borane |
| `bonded_hydrogens(structure, parent)` | Список явных H у атома |
| `first_bonded_hydrogen(structure, atom)` | Первый H или None |
| `halide_ion(structure, halogen, ...)` | `"destroy:chloride"` / fluoride / iodide |
| `carbonyl_atoms_from_site(structure, site, role)` | (C, O) из сайта с C=O |
| `organometallic_atoms(structure, site)` | (organo_C, metal, residue_atoms) |
| `atom_mass_sum / atom_charge_sum` | Масса/заряд набора атомов |
| `mapped_atom(mapping, old_idx, role)` | Атом после remove_atoms, или ошибка |
| `generated_site_reaction_id(prefix, participant)` | `"prefix/substance_id/atoms"` |
| `generated_pair_site_reaction_id(prefix, a, b)` | ID для парных реакций |
| `sanitize_id(value)` | Заменяет не-ASCII-алфавитные символы на `_` |

---

## Генератор: `substitution.rs`

**Входные сайты:** `HalideSite`, `AlcoholSite`, `AlkoxideSite`, `AmineSite`

| Функция | Реагент | Продукт | Тип |
|---------|---------|---------|-----|
| `generate_halide_hydroxide_substitution` | `destroy:hydroxide` | R−OH + галогенид | SN2 |
| `generate_halide_ammonia_substitution` | `destroy:ammonia` (2 экв.) | R−NH2 + галогенид + аммоний | SN2 |
| `generate_halide_cyanide_substitution` | `destroy:cyanide` | R−CN + галогенид | SN2 |
| `generate_halide_amine_substitution` | `HalideSite` + `AmineSite` | R−NR'2 + галогенид + протон | SN2 |
| `generate_thionyl_chloride_substitution` | `destroy:thionyl_chloride` | R−Cl + HCl + SO2 | — |
| `generate_alkoxide_protonation` | `destroy:proton` | R−OH | — |

Все SN2-реакции: EA = 25 кДж/моль, `SelectivityProfile(SN2, halide_desc)`.
Галоген-замена работает только для F, Cl, I (Br не поддержан на момент написания кода).

---

## Генератор: `alcohol.rs`

**Входные сайты:** `AlcoholSite`

| Функция | Условие | Реагент | Продукт |
|---------|---------|---------|---------|
| `generate_alcohol_oxidation` | degree < 3, есть α-H | `destroy:dichromate` + 8H⁺ (3:1:8) | альдегид/кетон + Cr³⁺ + H2O |
| `generate_alcohol_dehydration` | есть β-H | `destroy:oleum` (n:n) | алкен(ы) + H2SO4 |

Окисление: 3 моля спирта на 1 дихромат, EA = 25 кДж/моль.
Дегидратация: перебирает все β-C соседей → может выдать несколько реакций Zaitsev.

---

## Генератор: `acid_derivatives.rs`

**Входные сайты:** `CarboxylicAcidSite`, `AcylChlorideSite`, `AmideSite`, `EsterSite`

| Функция | Реагент | Продукты | EA кДж/моль |
|---------|---------|---------|------------|
| `generate_carboxylic_acid_esterification` | кислота + спирт | эфир + H2O | 25 |
| `generate_acyl_chloride_formation` | RCOOH + `destroy:phosgene` | RCOCl + HCl + CO2 | — |
| `generate_acyl_chloride_hydrolysis` | RCOCl + H2O | RCOOH + HCl | — |
| `generate_acyl_chloride_esterification` | RCOCl + ROH | эфир + HCl | — |
| `generate_ester_hydrolysis` | эфир + H2O | кислота + спирт | 42 |
| `generate_lah_ester_reduction` | эфир + LiAlH4 (внешн.) | 2 спирта | 18 |
| `generate_amide_hydrolysis` | амид + H2O | кислота + NH3 | — |

Этерификация карбоновой кислоты: катализатор `destroy:sulfuric_acid`, профиль `EsterProtection`.
Гидролиз эфира: кислотные условия (`AcidityCondition::Acidic`), min water 0.35, профиль `EsterHydrolysis`.
LAH-восстановление: dry, max water 0.02, внешний реагент 4.04 г/моль (4H эквиваленты).
Функция `ester_alkoxy_branch` — DFS по атомам эфира, отсекая по карбонильному C: выделяет алкокси-ветвь для расщепления на кислоту + спирт.

---

## Генератор: `carbonyl.rs`

**Входные сайты:** `CarbonylSite`, `AmineSite`, `AlcoholSite`, `SiteParticipant` (organometallic)

| Функция | Вход | Реагент | Продукт | EA |
|---------|------|---------|---------|---|
| `generate_acetal_formation` | Carbonyl + Alcohol | 2 ROH | ацеталь + H2O | 25 |
| `generate_imine_formation` | Carbonyl + PrimaryAmine | амин (H2N−) | имин + H2O | 25 |
| `generate_organometallic_carbonyl_addition` | Carbonyl + Organomet. | донор H (внешн.) | спирт | 15 |
| `generate_aldehyde_oxidation` | CarbonylSite (!is_ketone) | `destroy:dichromate` + H⁺ | карбоновая кислота | 25 |
| `generate_cyanide_nucleophilic_addition` | CarbonylSite | `destroy:hydrogen_cyanide` | циангидрин | 20 |
| `generate_borohydride_carbonyl_reduction` | CarbonylSite | `destroy:borohydride` + H2O | спирт | 16 |
| `generate_wolff_kishner_reduction` | CarbonylSite | `destroy:hydrazine` | алкан + N2 + H2O | 30 |

Ацетальное образование: катализатор H⁺, кислые условия, max water 0.35. При наличии стереоцентра — расширение в два канала (CW/CCW) через `expand_stereo_product_distribution`.
Иминообразование: требует 2 явных H на амине (первичный амин); связь C=N порядок 2.0.
Органометаллическое присоединение: условие inert + dry (max H2O/O2 0.02); остаток металла выходит как `chemical_external_product`.
Восстановление NaBH4: `SelectivityProfile(CarbonylReduction, strong)`, never_suppress.
Реакция Вольфа–Кишнера: катализатор OH⁻, EA 30, удаляет кислород (→ CH2).

---

## Генератор: `addition.rs`

**Входные сайты:** `UnsaturatedBondSite` (алкен или алкин)

### Электрофильное присоединение

`electrophilic_addition_specs(alkyne: bool)` возвращает 7 спецификаций:

| Реакция | Электрофил | Марковников: high-C получает | low-C получает | Стерео (алкин) |
|---------|-----------|------------------------------|----------------|----------------|
| Хлорирование | Cl2 | Cl | Cl | Anti |
| Хлоргидрин | HOCl | OH | Cl | Anti |
| Гидратация | H2O | OH | H | смесь |
| Гидроборирование | B2H6 (2 экв.) | H | BH2 | смесь |
| Гидрохлорирование | HCl | Cl | H | смесь |
| Гидрогенизация | H2 | H | H | смесь |
| Гидройодирование | HI | I | H | смесь |
| Иодирование | I2 | I | I | Anti |

EA: алкин 10 кДж/моль, алкен 25 кДж/моль (кроме гидратации: 20).
Гидрогенизация: внешний катализатор `forge:dusts/nickel`.

При отсутствии стерео-дистрибуции (смесь) `expand_stereo_product_distribution` разворачивает в два канала (E/Z или CW/CCW).

### Вспомогательные функции

| Функция | Назначение |
|---------|-----------|
| `expand_stereo_product_distribution` | Разворачивает `Stereochemistry::Mixture` в несколько `StereoProductVariant` |
| `geometric_z_steric_penalty_kj_per_mol` | Штраф EA для Z-изомера по объёму заместителей: 1.5..8.0 кДж/моль |
| `z_pre_exponential_multiplier` | Предэкспонент для Z: 0.55..0.95 |
| `substituent_steric_bulk` | Рекурсивная мера стерического объёма заместителя |
| `atomic_stereo_priority` | CIP-приоритет по атомному номеру |
| `apply_alkyne_stereo_rule` | Anti-правило для алкинов (→ Trans) |

---

## Генератор: `enolate.rs`

**Входные сайты:** `AlphaCarbonCenter`, `EsterSite`, `AmineSite`, `UnsaturatedBondSite`, `CarbonylSite`

| Функция | Вход | Условие | Продукт | EA |
|---------|------|---------|---------|---|
| `generate_aldol_addition` | Enol + Carbonyl | basic, ≤ 323 K | β-гидрокси-карбонил | 28 |
| `generate_alpha_halogenation` | AlphaCarbonCenter | acid | α-хлор/иодкарбонил | 30 |
| `generate_aldol_dehydration` | AlphaCarbonCenter (β-OH) | acid | α,β-ненасыщенный карбонил | 32 |
| `generate_enamine_formation` | Carbonyl + NonTertiaryAmine + Alpha | acid, max water 0.5 | enamine + H2O | 30 |
| `generate_enolate_alkylation` | AlphaCarbonCenter + HalideSite | basic | α-алкилированный карбонил | 26 |
| `generate_michael_addition` | AlphaCarbonCenter + Alkene | basic | Michael-аддукт | 28 |
| `generate_claisen_condensation` | AlphaCarbonCenter (Ester) + EsterSite | basic | β-кетоэфир + ROH | 31 |

`generate_enamine_formation` требует ровно 1 H на амине (вторичный амин).
`generate_michael_addition` работает только для алкенов с сопряжённой карбонильной группой (`michael_acceptor_atoms` детектирует α/β-положения).
`generate_claisen_condensation` доступна только если `center.carbonyl_kind == Ester`.
`alpha_selectivity_profile`: строит `SiteDescriptor` с учётом `steric_class`, `acidity`, `conjugation`, `carbonyl_kind`.

---

## Генератор: `heteroatom.rs`

**Входные сайты:** `NitrileSite`, `NitroSite`, `AmineSite`, `IsocyanateSite`, `SiteParticipant` (epoxide)

| Функция | Реагент | Продукт |
|---------|---------|---------|
| `generate_nitrile_hydrolysis` | R−C≡N + H2O + H⁺ | R−CONH2 (амид) |
| `generate_nitrile_hydrogenation` | R−C≡N + 2H2 | R−CH2−NH2 + Ni-кат. |
| `generate_nitro_hydrogenation` | R−NO2 + 3H2 | R−NH2 + 2H2O + Pd-кат. |
| `generate_amine_phosgenation` | R−NH2 + COCl2 | изоцианат + 2HCl |
| `generate_cyanamide_addition` | R−NH + NCN | R−N−C(=NH)−NH2 |
| `generate_isocyanate_hydrolysis` | R−N=C=O + H2O | R−NH2 + CO2 |
| `generate_epoxide_hydrolysis` | эпоксид + H2O | диол |

Гидрирование нитро: `forge:dusts/palladium`.
Гидрирование нитрила: `forge:dusts/nickel`.
Гидролиз эпоксида: кислые условия (Acidic), min water 0.1; раскрытие по наименее замещённому C (удаляется связь O−C[0]).
`deprotonated_alcohol_fragment` — вспомогательная: возвращает фрагмент спирта без протона (для использования в других генераторах).

---

## Генератор: `aromatic.rs`

**Входные сайты:** `SiteParticipant` (AromaticRing, ArylHalide, Halide, AcylChloride)

### EAS-каркас: `generate_eas_reaction`

Общий шаблон для всех EAS:
1. Строит `AromaticRingDescriptor` из `site.atoms`.
2. Если FC-реакция и кольцо деактивировано (`is_deactivated_for_fc`) → `Ok(None)`.
3. Перебирает атомы кольца с H: удаляет H, вызывает `editor_transform(editor, carbon)`, проверяет сохранность кольца.
4. `compute_eas_activation_delta(carbon)` → поправка к EA.
5. При 1 варианте → одноканальная реакция; при N → N каналов `position_0..N`.

| Функция | Реагент | EA базовая | Катализатор |
|---------|---------|-----------|------------|
| `generate_aromatic_nitration` | HNO3 | 30 | H2SO4 (acid cond.) |
| `generate_aromatic_chlorination` | Cl2 | 28 | FeCl3 |
| `generate_aromatic_bromination` | Br2 | 28 | FeBr3 |
| `generate_aromatic_sulfonation` | H2SO4 | 32 | — |
| `generate_fc_alkylation` | арен + R−X | 35 | AlCl3 |
| `generate_fc_acylation` | арен + RCOCl | 30 | AlCl3 |

FC-реакции (алкилирование и ацилирование): `join_structures` для объединения фрагментов.

### SNAr-реакции: `generate_aryl_halide_*`

- Строит кольцо через `from_start_carbon`.
- Если `ring_atoms.len() < 5` → `None` (не настоящее ароматическое кольцо).
- `compute_snar_activation_delta` — если ≥ 0 → `None` (не активировано EWG).

| Функция | Нуклеофил | Продукт |
|---------|-----------|---------|
| `generate_aryl_halide_hydroxide_substitution` | OH⁻ | ArOH + галогенид |
| `generate_aryl_halide_ammonia_substitution` | NH3 (2 экв.) | ArNH2 + галогенид + NH4⁺ |

---

## Генератор: `phosphorus.rs`

**Входные сайты:** `HalideSite`, `PhosphineSite`, `PhosphoniumSaltSite`, `PhosphorusYlideSite`, `PhosphonateCarbanionSite`, `SulfoneCarbanionSite`, `CarbonylSite`

| Функция | Реагенты | Продукты | EA |
|---------|---------|---------|---|
| `generate_phosphonium_salt_formation` | R−X + PPh3 | [Ph3P−R]⁺X⁻ | 22 |
| `generate_phosphonium_ylide_formation` | [Ph3P−CH−R]⁺ + EtO⁻ | Ph3P=C(R) | 28 |
| `generate_wittig_olefination` | ylide + carbonyl | алкен + Ph3P=O | 24 |
| `generate_horner_wadsworth_emmons_olefination` | phosphonate + carbonyl | алкен + фосфат | 22 |
| `generate_julia_olefination` | sulfone + carbonyl | алкен + сульфонат | 22 |

**Реакция Виттига**: при наличии стерео — два канала (E/Z). Селективность определяется стабильностью илида:
- Unstabilized: предпочитает Z (−2.0 к EA_Z, ×1.25 к предэкспоненту)
- Stabilized: предпочитает E (−1.5 к EA_E)

**HWE**: E-предпочтение (strong, −2.0 / +1.5). Побочный продукт — фосфат (P с дополнительным O⁻).
**Julia**: умеренное E-предпочтение (−1.0 / +0.75). Побочный продукт — сульфонат (S с дополнительным O⁻).

Phosphonium salt formation: только для degree < 3 на галоген-C.

---

## Генератор: `boron.rs`

**Входные сайты:** `BoraneSite`, `BorateEsterSite`

| Функция | Реагент | Продукт |
|---------|---------|---------|
| `generate_borane_oxidation` | R−C−BH2 + H2O2 + OH⁻ | R−C−O−BH2 (→ спирт) |
| `generate_borate_ester_hydrolysis` | R−O−B + H2O | R−OH + бороновая кислота |

Боран-окисление: вставляет мостиковый O между C и B через `insert_bridging_atom`.
Гидролиз борного эфира: `split_at_bond(O, B)` делит молекулу, к O добавляется H, к B − OH.

---

## Генератор: `protecting_groups.rs`

**Входные сайты:** `AlcoholSite`, `SilylEtherCenter`, `AcetalCenter`, `AmineSite`, `BocCarbamateCenter`, `CbzCarbamateCenter`

### Силильная защита спирта (TMS)

| Реакция | Реагент | Продукт | EA |
|---------|---------|---------|---|
| Защита (TMS) | R−OH + Me3SiCl | R−OSiMe3 + HCl | 15 |
| Снятие защиты | R−OSiMe3 + F⁻ + H⁺ | R−OH + Me3SiF | 20 |

Защита: dry (max water 0.1). К кислороду добавляется Si с тремя метильными группами.
Снятие: удаляет Si и все метильные C (включая их H), добавляет H к O.

### Ацеталь/кеталь (защита карбонила)

Образование: `generate_acetal_formation` в `carbonyl.rs` (см. выше).

| Реакция | Реагент | Продукт | EA |
|---------|---------|---------|---|
| Снятие (гидролиз) | R2C(OR')2 + H2O | R2C=O + 2 R'OH | 25 |

Снятие: кислые условия, min water 0.35. `branch_atoms` (DFS) выделяет каждую ветвь, проверяет что ветви не пересекаются (иначе ошибка — циклический ацеталь).

### Boc-защита аминов

| Реакция | Реагент | Продукт | EA |
|---------|---------|---------|---|
| Защита | R2NH + Boc2O | R2NCO2tBu + tBuOH + CO2 | 20 |
| Снятие | R2NCO2tBu + H2O + H⁺ | R2NH + tBuOH + CO2 | 18 |

Защита: basic, max water 0.2. Функция `add_boc_group` добавляет карбамат с трет-бутилом.
Снятие: кислые условия, min water 0.2. Удаляет карбонильный C, O, алкокси-O, трет-бутильный C с тремя CH3.

### Cbz-защита аминов

| Реакция | Реагент | Продукт | EA |
|---------|---------|---------|---|
| Защита | R2NH + CbzCl | R2NCO2CH2Ph + HCl | 20 |
| Снятие | R2NCO2CH2Ph + H2 | R2NH + толуол + CO2 | 12 |

Снятие: Pd-катализатор, гидрогенолиз. `add_cbz_group` добавляет карбамат с бензильной группой (6-членное кольцо из C с bond order 1.5).

---

## Генератор: `combustion.rs`

**Вход:** `&Substance` (без сайтов)

Генерирует реакцию полного сгорания для органических CHO-соединений.

### Условия генерации

- `fuel.charge == 0`
- Вещество не помечено тегом `destroy:hypothetical`
- Структура содержит только C, H, O

### Алгоритм стехиометрии

Для формулы CₙHₘOₖ:

```
oxygen_quarters = 4n + m − 2k
multiplier = 4 / gcd(oxygen_quarters, 4)
→ fuel×multiplier + O2×(oxygen_quarters×multiplier/4)
  → CO2×(n×multiplier) + H2O×(m×multiplier/2)
```

### Оценка энтальпии

```
ΔH = −393.5 × n_CO2 − 241.8 × n_H2O  (кДж/моль реакции)
```

### Кинетика

Pre-exponential: 2.5×10¹¹. EA: 115 кДж/моль. Фаза: Gas для всех участников.

---

## Общие инварианты генераторов

- Все функции `generate_*` явно передают водороды в структуре (нет неявных H).
- `resolver.resolve(structure)` — единственный способ создать/найти вещество по структуре.
- Стереохимическая смесь (`Stereochemistry::Mixture`) разворачивается через `expand_stereo_product_distribution` → несколько каналов `ReactionChannel` с разными EA и предэкспонентами.
- `General` стерео-смесь → ошибка `"stereo mixture has no quantitative distribution"`, которая в движке перехватывается и пропускается.

## Связи

- [[organic-centers]] — все XxxSite<'a> как входные параметры
- [[organic-space-resolver]] — `DerivedSubstanceResolver`
- [[organic-aromatic-catalog]] — `AromaticRingDescriptor`
- [[organic-engine]] — вызывает все функции `generate_*`
- [[selectivity-engine]] — `SiteDescriptorBuilder`, `SelectivityProfile`
- [[selectivity-types]] — `ReactionType`, `NucleophileStrength`
- [[core-reaction]] — `Reaction`, `ReactionChannel`, `StoichiometricTerm`
- [[molecule-graph]] — `MolecularEditor`, `MolecularStructure`



