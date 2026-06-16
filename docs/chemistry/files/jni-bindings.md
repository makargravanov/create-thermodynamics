# JNI-привязки

Исходный код: `lib.rs`, `build.rs`

## Назначение

`lib.rs` — точка входа нативной библиотеки. Экспортирует функции через JNI,
позволяя JVM-коду (пакет `dev.makargravanov.create_thermodynamics.common.rust`)
вызывать Rust-ядро напрямую через нативный интерфейс Java.

`build.rs` — скрипт сборки Cargo. Парсит Java-исходники устаревшего Destroy-мода
(`DestroyReactions.java`, `DestroyMolecules.java`) и кодогенерирует Rust-файл
`destroy_reactions.rs` с регистрацией всех реакций.

## Ключевые типы

Нативная библиотека компилируется как `cdylib` (динамическая библиотека для JNI).
Модуль `chemistry` реэкспортируется как `pub mod chemistry`.

## Публичные входы

### nativeIdealGasPressure

```rust
Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeIdealGasPressure(
    moles: jdouble,
    temperature_kelvin: jdouble,
    volume_cubic_meters: jdouble,
) -> jdouble
```

Вычисляет давление по уравнению идеального газа: `P = nRT/V`.
Константа `R = 8.314_462_618_153_24` Дж/(моль·К).
При `volume == 0.0` возвращает `NaN`.

### nativeAbiVersion

```rust
Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeAbiVersion(
) -> jint
```

Возвращает `2`. Используется JVM-стороной для проверки совместимости
нативной библиотеки при загрузке.

### nativeReplaceMinecraftItemChemicalBindings

```rust
Java_dev_makargravanov_create_1thermodynamics_common_rust_ThermodynamicsNative_nativeReplaceMinecraftItemChemicalBindings(
    item_ids: Array<String>,
    substance_ids: Array<String>,
    mol_per_items: DoubleArray,
)
```

Атомарно заменяет весь набор соответствий Minecraft-предметов химическим веществам.
JVM-сторона передаёт три массива одинаковой длины:

- `item_ids` — идентификаторы предметов Minecraft;
- `substance_ids` — устойчивые `SubstanceId` из химического каталога;
- `mol_per_items` — количество вещества в одном предмете.

Rust сначала строит новый `MinecraftChemicalRegistry` во временной структуре и
проверяет все записи. Если есть неизвестное вещество, дубль предмета, `NaN`,
бесконечность или неположительное количество, текущий набор биндингов не меняется,
а JVM получает `IllegalArgumentException`.

### nativeClearMinecraftItemChemicalBindings

Очищает все предметные соответствия.

### nativeMinecraftItemChemicalBindingCount

Возвращает количество зарегистрированных предметных соответствий.

### nativeHasMinecraftItemChemicalBinding

Проверяет, есть ли соответствие для конкретного идентификатора предмета.

## Поток данных / Алгоритм

### build.rs: кодогенерация реакций

1. Читает `DestroyReactions.java` и `DestroyMolecules.java` из смежного модуля Destroy
2. `split_reaction_blocks` — разбивает Java-источник на блоки по паттерну
   `= builder()` … `.build()`
3. `parse_molecule_id_map` — строит `HashMap<ConstantName, SubstanceId>` из блоков молекул
4. `parse_acid_entries` — извлекает строки `builder().acid(acid, base, pKa)`
5. `emit_reaction_builder` — транслирует каждый Java-блок в цепочку вызовов
   `Reaction::builder(...).reactant(...).product(...).build()`
6. `emit_acid_reactions` — для каждой кислоты генерирует три реакции:
   диссоциация, нейтрализация и ассоциация; скорость `0.5 × 10^(-pKa)`
7. Результат записывается в `$OUT_DIR/destroy_reactions.rs`

Файлы `DestroyReactions.java` и `DestroyMolecules.java` объявлены как зависимости
через `cargo:rerun-if-changed`, поэтому кодогенерация повторяется только при их изменении.

## Инварианты и ошибки

- Дублирующиеся ID реакций при кодогенерации → `panic!` в `build.rs`
- Несбалансированные скобки в Java-источнике → `panic!` в `matching_paren`
- ABI-версия `2` жёстко зашита; изменение нативного API требует её инкремента
  и синхронизации с JVM-стороной
- Предметные биндинги заменяются атомарно: ошибка в одной записи не оставляет
  частично применённый набор

## Связи

- [[dynamic|DynamicChemistryRegistry]] — основной объект химического ядра
- [[data-reactions|destroy_reactions.rs]] — генерируемый файл (каталог реакций Destroy)
- [[core-registry|ChemistryRegistry]] / [[core-registry|ChemistryRegistryBuilder]] —
  используются в сгенерированном коде для регистрации реакций
- [[core-thermodynamics|термодинамика]] — `nativeIdealGasPressure` реализует
  простейший газовый закон, остальные расчёты — в отдельных модулях
