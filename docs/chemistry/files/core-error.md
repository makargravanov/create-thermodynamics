# Ошибки химии

Исходный код: `error.rs`

## Назначение

Единый тип ошибок химического ядра. Модель следует правилу: невозможное или
повреждённое состояние должно явно возвращать ошибку, а не молча давать ноль
или пропуск. Все публичные входы ядра возвращают `ChemistryResult<T>`.

## Ключевые типы

```rust
pub type ChemistryResult<T> = Result<T, ChemistryError>;
```

`ChemistryError` — перечисление. Каждый вариант несёт достаточно данных, чтобы
понять, какая реакция, вещество, генератор или равновесие вызвало проблему.

| Вариант | Когда возникает |
|---|---|
| `DuplicateSubstance(id)` | вещество зарегистрировано дважды |
| `DuplicateReaction(id)` | реакция зарегистрирована дважды |
| `UnknownSubstance { reaction_id, substance_id }` | реакция ссылается на отсутствующее вещество |
| `UnknownReaction(id)` | запрос несуществующей реакции |
| `InvalidSubstance { substance_id, reason }` | вещество не прошло `validate` |
| `InvalidReaction { reaction_id, reason }` | реакция не прошла `validate_shape` |
| `ChargeNotConserved { reaction_id, reactants, products }` | сумма зарядов не сходится |
| `MassNotConserved { reaction_id, reactants, products }` | сумма масс не сходится |
| `ReversibleThermodynamicsMismatch { reaction_id, reverse_id, reason }` | пара обратных реакций несогласована |
| `GenerationInvariantViolation { generator, substance_id, reason }` | нарушен инвариант динамической генерации |
| `EquilibriumInvariantViolation { equilibrium_id, reason }` | нарушен инвариант равновесия |
| `InvalidMixtureState(reason)` | смесь в недопустимом состоянии |

## Публичные входы

`ChemistryError` реализует `Display` (человекочитаемое сообщение для каждого
варианта) и `std::error::Error`. Это позволяет пробрасывать ошибку через `?`
и логировать на JVM-стороне.

## Связи

- [[core-registry|Реестр]] — источник большинства ошибок сборки (дубли, ссылки, масса, заряд)
- [[core-substance|Вещество]] — `InvalidSubstance` из `Substance::validate`
- [[core-reaction|Реакция]] — `InvalidReaction` из `validate_shape`
- [[core-mixture|Смесь]] — `InvalidMixtureState` из `Mixture::validate`
- [[dynamic|Динамический слой]] — `GenerationInvariantViolation`
- [[core-solution|Раствор]] — `EquilibriumInvariantViolation`
