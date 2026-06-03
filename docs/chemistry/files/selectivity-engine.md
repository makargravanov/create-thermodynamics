# Движок оценки селективности

Исходный код: `kinetics/selectivity/engine.rs`

## Назначение

`SelectivityEngine` — единственная точка входа, через которую генераторы органических реакций получают числовые поправки к кинетике. Движок делегирует оценку специализированным функциям модулей и нормирует результат в `SelectivityRuntimeEffect`.

`SiteDescriptorBuilder` — вспомогательный конструктор дескрипторов: как параметрический (`build`), так и factory-методы (десятки готовых шаблонов: `primary_halide`, `aldehyde`, `boc_carbamate` и т.д.), и внутренние методы `from_*_site`, строящие дескриптор напрямую из данных молекулярного графа.

## Поток данных / Алгоритм

### SelectivityEngine::evaluate_profile

```
SelectivityProfile → SelectivityRuntimeEffect
```

1. По полю `mechanism` вызывается соответствующая функция оценки:
   - `SN2` → `evaluate_sn2` (возвращает `SelectivityResult`, берётся `primary`)
   - `SN1` → `evaluate_sn1`
   - `E2` / `E1` → `evaluate_e2` / `evaluate_e1`
   - `FischerEsterification` → `evaluate_fischer_esterification` (требует `secondary_site`)
   - `CarbonylAddition` → `evaluate_carbonyl_addition` + `apply_inorganic_environment_to_carbonyl_score`
   - `CarbonylReduction` → то же + проверки кислотной/окислительной среды
   - `Wittig` / `HWE` / `Julia` → `evaluate_carbonyl_addition` с `nucleophile_strength = VeryStrong`, затем штраф за стерику реагента и резонансную стабилизацию
   - Группа защитных групп → `evaluate_protecting_group_profile`
   - Группа alpha-carbon → `evaluate_alpha_carbon_profile`

2. На основе `score.value` определяется `recommendation`:
   - Для `SN2` и `FischerEsterification` recommendation берётся напрямую из `SelectivityResult`.
   - Иначе: value < 0.01 → `None`; < 0.1 → `Suppressed`; иначе → `Preferred`.

3. `suppressed = SuppressWhenDisfavored AND (Suppressed OR None)`

4. Возвращается `SelectivityRuntimeEffect`:
   - `rate_multiplier = 0.0` если suppressed, иначе `score.value.max(0.0)`
   - `activation_delta_kj_per_mol` из `score.activation_delta`
   - `pre_exp_multiplier` из `score.pre_exp_multiplier`

### apply_inorganic_environment_to_carbonyl_score

Применяется только к `VeryStrong` нуклеофилам (Grignard, RLi).

| Условие | rate × | ΔEa (кДж/моль) |
|---|---|---|
| water_activity ≥ 0.35 (богатая вода) | 0.02 | +24 |
| water_activity > 0.02 (следы воды) | 0.2 | +10 |
| is_oxygen_rich OR is_oxidizing | 0.25 | +8 |
| свободные комплексируемые металлы | 1/(1+act×20) | +(1−penalty)×12 |

Штрафы независимы и мультипликативны. Причины дописываются к `score.reason`.

### redox_competition_penalty

Используется при `CarbonylReduction` в окисляющей среде:

```
oxidizing_fraction = oxidizing / (oxidizing + reducing)
penalty = clamp(1.0 − oxidizing_fraction × 0.9, 0.05, 1.0)
```

### evaluate_alpha_carbon_profile

Общая логика для реакций через α-углерод: базовый `value = steric_accessibility`, затем:
- ≥ 2 EWG: value × 2, ΔEa − 4 кДж/моль («активированный α-углерод»)
- 0 EWG: value × 0.6, ΔEa + 3 кДж/моль

Дополнительные поправки по типу механизма:
- `AlphaHalogenation`: кислота или основание × 1.4 / иначе × 0.4 (нужна кислота или основание)
- `AldolAddition`, `EnolateAlkylation`, `MichaelAddition`, `ClaisenCondensation`: основание × 1.8; иначе × 0.35; вода × 0.65
- `AldolDehydration`: кислота или высокая T × 1.5
- `EnamineFormation`: кислота × 1.2; основание × 0.4

### evaluate_protecting_group_profile

Условия для каждого типа защитной группы:

| Тип | Требуется для активации |
|---|---|
| `SilylEtherFormation` | сухая среда (water_poor) + основание |
| `SilylEtherCleavage` | фторид (×4, ΔEa −10); без него ×0.02 |
| `AcetalFormation` | кислота + сухая среда |
| `AcetalHydrolysis` | кислота + вода |
| `CarbamateFormation` | основание + сухая среда |
| `CarbamateCleavage` | (кислота+вода) ИЛИ (H₂+Pd+поверхность) |
| `EsterProtection` | кислота + сухая среда |
| `EsterHydrolysis` | кислота/основание + вода |

### SiteDescriptorBuilder

Factory-методы без аргументов для стандартных центров (`primary_alcohol`, `secondary_alcohol`, `tertiary_alcohol`, `benzylic_alcohol`, `aldehyde`, `ketone`, `aromatic_aldehyde`, `primary_halide` … `benzylic_halide`, `carboxylic_acid`, `silyl_ether`, `acetal`, `boc_carbamate`, `cbz_carbamate`).

Внутренние методы `from_*_site` строят `SiteDescriptor` из молекулярного графа через `descriptor_from_carbon`:
1. Определяет степень замещения: приоритет бензильного (≥2 ароматических связи у соседа) → аллильного → 3° → 2° → 1°.
2. BFS до глубины 3 для подсчёта EDG/EWG. Правила: O/N/S с одиночной связью → донор; O/N/F/Cl/Br/I в общем случае → акцептор; C с =O/=N/≡N → акцептор; C иначе → донор.
3. «Громоздкий заместитель» (bulky): Br/I или C-сосед с ≥2 нехудожественными заместителями.
4. β-водород: одиночная связь к C, у которого есть H.

## Инварианты и ошибки

- Если для `FischerEsterification`, `Wittig`/`HWE`/`Julia` не передан `secondary_site`, возвращается эффект с `rate_multiplier=0`, `activation_delta=50 кДж/моль`, `suppressed=true`.
- `SelectivitySuppressionPolicy::NeverSuppress` обходит проверку рекомендации — реакция получает ненулевой `rate_multiplier` даже при `Suppressed`/`None`.
- Движок не вызывает `from_mixture`; всё заполнение `SelectivityContext` — ответственность вызывающего кода.

## Связи

- [[selectivity-types|Типы селективности]] — все входные и выходные типы.
- [[selectivity-reactions|Реакции]] — `evaluate_sn2`, `evaluate_e2`, `evaluate_fischer_esterification`, `evaluate_carbonyl_addition`.
- [[organic-generators|Генераторы]] — потребляют `evaluate_profile` и convenience-методы движка.
- [[organic-centers|Центры]] — передают `*Site` объекты в `SiteDescriptorBuilder::from_*_site`.
