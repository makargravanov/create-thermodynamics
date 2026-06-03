# Классы реакций: модели селективности

Исходный код: `kinetics/selectivity/nucleophilic_substitution.rs`, `kinetics/selectivity/elimination.rs`, `kinetics/selectivity/esterification.rs`, `kinetics/selectivity/carbonyl_addition.rs`

## Назначение

Четыре модуля реализуют расчётные функции для конкретных классов реакций. Каждый возвращает либо `ReactivityScore` (одна реакция), либо `SelectivityResult` (реакция с конкуренцией). Движок [[selectivity-engine|SelectivityEngine]] вызывает эти функции напрямую.

---

## Нуклеофильное замещение (SN1/SN2)

Файл `nucleophilic_substitution.rs`.

### evaluate_sn2

Возвращает `SelectivityResult` с конкурентом E2 (если есть β-водород).

Факторы `sn2_score = base × steric³ × electronic × temp × solvent`:

| Степень | base |
|---|---|
| Primary | 1.0 |
| Benzylic | 0.8 |
| Allylic | 0.7 |
| Secondary | 0.25 |
| Tertiary | 0.005 |

- Стерика входит **в кубе** (`steric_accessibility³`) — главный ингибирующий фактор SN2.
- Высокая T: ×0.7 (очень высокая) / ×0.85 (высокая).
- `AproticPolar`: ×1.5; `Protic`: ×0.6; `Basic`: ×0.8.
- ΔEa += steric_score × 20 кДж/моль; pre_exp ≥ 0.3.

E2-конкурент внутри (`evaluate_e2_competition`): base 3° 1.0 → 1° 0.05; умножается на `temp_bonus` (до ×2) и `base_bonus` (до ×2.5); для 3°/2° ΔEa −8 кДж/моль (E2 выгоднее по барьеру).

### evaluate_sn1

Используется напрямую только через `evaluate_profile(SN1, …)`.

| Степень | base |
|---|---|
| Tertiary | 1.0 |
| Benzylic | 0.9 |
| Allylic | 0.8 |
| Secondary | 0.2 |
| Primary | 0.001 |

`Protic`: ×1.5; `Acidic`: ×1.3; `AproticPolar`: ×0.5; `Basic`: ×0.1. EDG стабилизируют карбокатион (+0.15/ед.), резонанс +0.3. Высокая T: ΔEa −5 кДж/моль.

---

## Элиминирование (E1/E2)

Файл `elimination.rs`.

### evaluate_e2

Обязательное условие: `has_beta_hydrogen`. Без него возвращается `value=0`.

`score = base × temp × base_factor × solvent × (1 + steric_bonus)`:

| Степень | base |
|---|---|
| Tertiary | 1.0 |
| Secondary | 0.6 |
| Allylic | 0.5 |
| Benzylic | 0.4 |
| Primary | 0.15 |

- Очень высокая T: ×2.5, ΔEa −12 кДж/моль; высокая T: ×1.8, ΔEa −6 (энтропийный вклад).
- Основание pH > 8: ×3.0; `SolventType::Basic`: ×2.5; pH > 10: ×2.0; нейтрально: ×0.3 (E2 требует основания).
- `AproticPolar`: ×1.4; `Protic`: ×0.8.
- Стерика даёт бонус (×(1 + score×0.5)): объёмистый субстрат вытесняет SN2 в пользу E2.

### evaluate_e1

Обязательное условие: `has_beta_hydrogen`.

Tertiary 0.7, Benzylic 0.6, Allylic 0.5, Secondary 0.15, Primary ≈0. `Protic`: ×1.4; `Basic`: ×0.05. EDG и резонанс стабилизируют карбокатион.

### zaitsev_hofmann_ratio

Возвращает отношение Цайцева к Гофману. Базовое значение 4.0 (обычный нагрев). Объёмное основание (`SolventType::Basic`): ×0.5 → сдвиг к Гофману. Стерика субстрата снижает коэффициент (×(1 − score×0.5)). Высокая T слегка увеличивает долю термодинамического продукта Цайцева (×1.2).

---

## Этерификация Фишера

Файл `esterification.rs`.

### evaluate_fischer_esterification

Возвращает `SelectivityResult` с конкурентом E2 для вторичных/третичных спиртов.

`score = acid_factor × alcohol_factor × steric × electronic × temp_factor`:

| Спирт | alcohol_factor |
|---|---|
| Primary | 1.0 |
| Benzylic | 0.95 |
| Allylic | 0.9 |
| Secondary | 0.35 |
| Tertiary | 0.02 |

- Кислота: Primary 1.0, Secondary 0.8, Tertiary 0.4 (пиваловая кислота — сильное затруднение).
- Стерика: `(1 − acid.steric×0.5) × (1 − alc.steric×0.7)`.
- Температура для 2°: высокая T × 0.5, очень высокая × 0.2 (конкуренция E2).
- Для третичных третичных спиртов: temp_factor = 0.1 независимо от T.

Конкурент E2: запускается для 2°/3° с `solvent_type=Acidic`-клоном контекста. Если E2 побеждает — `dominant_competitor` не пустой.

Дополнительные утилиты: `esterification_competition` (ранжирует несколько спиртов по относительным скоростям), `will_esterify` (быстрый предикат: Primary/Benzylic → true; Secondary при не-очень-высокой T → true; Tertiary → false).

---

## Нуклеофильное присоединение к карбонилу

Файл `carbonyl_addition.rs`.

### evaluate_carbonyl_addition

`score = carbonyl_factor × steric × electronic × aromatic_penalty × nucleophile × solvent`:

| Степень | carbonyl_factor | Интерпретация |
|---|---|---|
| Primary | 1.0 | альдегид R-CHO |
| Secondary | 0.3 | кетон R₂C=O |

Стерика альдегида 0.95 (H мал); кетона: `max(0.4 − score×0.2, 0.15)`.

| Нуклеофил | ×factor |
|---|---|
| VeryStrong (RMgX) | 2.0 |
| Strong (NaBH₄) | 1.5 |
| Moderate (OH⁻) | 1.0 |
| Weak (вода) | 0.5 |

- Ароматический карбонил (бензальдегид): ×0.6 (резонанс с кольцом снижает электрофильность).
- EDG: −0.15 за группу; EWG: +0.2 за группу.
- `AproticPolar`: ×1.3; `Protic`: ×0.8.
- ΔEa += steric × 5 кДж/моль (альдегид) или × 15 кДж/моль (кетон).
- Сильный нуклеофил снижает чувствительность к стерике: `activation_delta × 0.5` (VeryStrong) или ×0.7 (Strong).

Вспомогательные: `aldehyde_vs_ketone` (пара скоростей), `chemoselectivity_aldehyde_over_ketone` (отношение rate_ald / rate_ket).

## Связи

- [[selectivity-engine|Движок]] — импортирует и оркестрирует все четыре модуля.
- [[selectivity-types|Типы]] — `SiteDescriptor`, `SelectivityContext`, `ReactivityScore`, `SelectivityResult`.
- [[organic-generators|Генераторы]] — строят профили и запускают оценку через `SelectivityEngine`.
