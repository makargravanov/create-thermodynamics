# Система стилей интерфейсов

## Кратко

Нужно расширить `ui-dsl` так, чтобы стиль был не набором ручных цветов и прямоугольников в каждом экране, а проверяемой частью описания интерфейса. Экран должен описывать смысловые элементы: окно, панель, вкладку, кнопку, карточку, список, сетку слотов. Тема и таблица стилей должны определять, как эти элементы выглядят, а компилятор интерфейса должен раскрывать их в текущую низкоуровневую программу, проверять и оптимизировать до запуска Minecraft.

Главная цель: разработчик или ИИ-агент пишет интерфейс через устойчивые компоненты и стили, а не подгоняет пиксели. Если стиль неполный, не подходит элементу, плохо читается, не поддерживается целью вывода или создаёт невозможную отрисовку, это должно быть видно в тестах, отчёте сборки и предварительном просмотре.

## Принципы

- Никаких строковых имён стилей вроде `"primaryButton"`. Стиль должен быть типизированным значением Kotlin.
- Никаких скрытых значений по умолчанию для важных визуальных свойств. Цвет, отступ, шрифт, рамка и политика текста берутся из явно выбранной темы.
- Экран не должен руками собирать распространённые элементы из `box + background + padding`, если для них есть смысловой компонент.
- Состояния элемента являются частью стиля. Если кнопка может быть выбранной, стиль кнопки обязан описывать выбранное состояние.
- Предварительный просмотр и Minecraft получают один и тот же результат раскрытия стилей.
- Оптимизатор должен видеть стили как источник сведений о статичности: фон окна, рамки, панели и декоративные слои можно запекать, если они доказуемо постоянны.
- Ошибки стиля ловятся до игры. Во время игры экран не должен падать из-за переполнения строки или отсутствующего цвета.
- Тема не должна знать химию, реакторы или Minecraft-сущности. Она описывает интерфейс.

## Текущая опора

Уже есть:

- `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/foundation/UiElement.kt` - дерево интерфейса.
- `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/foundation/modifier/Modifier.kt` - модификаторы.
- `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/layout/UiLayoutResolver.kt` - расчёт раскладки.
- `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/PrimitiveScreenProgram.kt` - низкоуровневая программа.
- `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/PrimitiveScreenAnalysis.kt` - проверка программы.
- `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/PrimitiveScreenOptimizer.kt` - оптимизация, включая запекание текстур.
- `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/MinecraftScreenSourceGenerator.kt` - генерация Kotlin-кода под Minecraft.
- `modules/ui/src/uiSpec/kotlin/dev/makargravanov/create_thermodynamics/ui/reactor/ReactorControllerUi.kt` - текущий реакторный интерфейс, который надо перевести с ручных цветов и размеров на стиль.

## Целевая архитектура

Поток сборки:

```text
DSL-дерево
  -> StyleResolver
  -> UiLayoutResolver
  -> PrimitiveScreenProgram
  -> PrimitiveProgramAnalyzer
  -> PrimitiveScreenOptimizer
  -> MinecraftScreenSourceGenerator / предварительный просмотр
```

`StyleResolver` должен быть отдельным проходом между смысловым описанием и раскладкой. Он раскрывает стилевые элементы в уже существующие элементы `UiElement`, поэтому не ломает текущую низкоуровневую программу и генератор.

## Новые основные типы

Создать пакет:

```text
external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/
```

В нём:

```kotlin
data class UiTheme(
    val tokens: UiTokens,
    val styles: UiStyleSheet,
)
```

```kotlin
data class UiTokens(
    val colors: UiColorTokens,
    val spacing: UiSpacingTokens,
    val typography: UiTypographyTokens,
    val borders: UiBorderTokens,
    val textures: UiTextureTokens,
)
```

```kotlin
data class UiStyleSheet(
    val window: WindowStyle,
    val panel: PanelStyle,
    val button: ButtonStyle,
    val tab: TabStyle,
    val metricCard: MetricCardStyle,
    val listRow: ListRowStyle,
    val slotGrid: SlotGridStyle,
    val tooltip: TooltipStyle,
)
```

Стиль поверхности:

```kotlin
data class SurfaceStyle(
    val fill: StyleColor,
    val border: BorderStyle? = null,
    val texture: TextureStyle? = null,
    val padding: Insets = Insets.Zero,
    val bakeHint: BakeHint = BakeHint.Neutral,
)
```

Стиль текста:

```kotlin
data class TextStyle(
    val color: StyleColor,
    val alignment: TextAlignment = TextAlignment.Start,
    val overflow: TextOverflowPolicy,
    val lineHeight: Int,
)
```

Стиль управляющего элемента:

```kotlin
data class ControlStyle<S : Any>(
    val states: Map<S, SurfaceStyle>,
    val label: TextStyle,
)
```

Состояния лучше делать отдельными типами, а не одним общим перечислением:

```kotlin
enum class ButtonState {
    Normal,
    Hovered,
    Pressed,
    Selected,
    Disabled,
}

enum class TabState {
    Normal,
    Selected,
    Hovered,
    Disabled,
}
```

Цвета и текстуры:

```kotlin
sealed interface StyleColor {
    data class Constant(val color: Color) : StyleColor
    data class Dynamic(val expression: Value<Color>) : StyleColor
}
```

```kotlin
sealed interface TextureStyle {
    data class Resource(
        val namespace: String,
        val path: String,
    ) : TextureStyle
}
```

Запекание:

```kotlin
enum class BakeHint {
    Neutral,
    PreferBakedTexture,
    KeepPrimitiveCommands,
}
```

## Смысловые элементы

Создать пакет:

```text
external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/styled/
```

Минимальные элементы:

```kotlin
sealed interface StyledUiElement<out Action>
```

```kotlin
data class StyledWindow<Action>(
    val title: Value<String>,
    val children: List<StyledUiElement<Action>>,
) : StyledUiElement<Action>
```

```kotlin
data class StyledPanel<Action>(
    val style: PanelStyle,
    val children: List<StyledUiElement<Action>>,
) : StyledUiElement<Action>
```

```kotlin
data class StyledButton<Action>(
    val style: ButtonStyle,
    val state: Value<ButtonState>,
    val action: Value<Action?>,
    val label: Value<String>,
) : StyledUiElement<Action>
```

```kotlin
data class StyledTab<Action>(
    val state: Value<TabState>,
    val action: Value<Action?>,
    val label: Value<String>,
) : StyledUiElement<Action>
```

```kotlin
data class StyledMetricCard(
    val title: Value<String>,
    val lines: List<Value<String>>,
) : StyledUiElement<Nothing>
```

Публичный DSL должен выглядеть примерно так:

```kotlin
screen(
    size = ReactorControllerUiSize,
    theme = CreateLikeTheme,
) {
    window(title = "Reactor Controller") {
        tabs {
            tab("Overview", selected = stateValue(...), action = stateValue(...))
            tab("Zones", selected = stateValue(...), action = stateValue(...))
            tab("Mixture", selected = stateValue(...), action = stateValue(...))
        }

        metricGrid(columns = 3) {
            metric("State", stateValue(...))
            metric("Native", stateValue(...))
            metric("Zones", stateValue(...))
        }
    }
}
```

В пользовательском экране допускаются размеры сеток, число колонок и выбор компонента. Цвета, рамки, внутренние отступы и высоты типовых элементов должны жить в стиле.

## Проверка стилей

Создать:

```text
external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/StyleAnalysis.kt
```

Диагностики:

```kotlin
sealed interface StyleDiagnostic {
    val path: String

    data class MissingControlState(
        override val path: String,
        val style: String,
        val state: String,
    ) : StyleDiagnostic

    data class UnsupportedTextureForTarget(
        override val path: String,
        val target: String,
        val texture: String,
    ) : StyleDiagnostic

    data class LowTextContrast(
        override val path: String,
        val foreground: Color,
        val background: Color,
        val ratio: Double,
    ) : StyleDiagnostic

    data class UnsafeDynamicTextPolicy(
        override val path: String,
        val policy: TextOverflowPolicy,
    ) : StyleDiagnostic
}
```

Правила:

- `ControlStyle<ButtonState>` обязан содержать `Normal`, `Hovered`, `Pressed`, `Selected`, `Disabled`.
- `ControlStyle<TabState>` обязан содержать `Normal`, `Selected`, `Hovered`, `Disabled`.
- Динамический текст не может использовать `TextOverflowPolicy.FailInValidation`.
- Если цвет текста и фон постоянны, надо проверять контраст. Для Minecraft-пиксельного интерфейса порог можно сделать умеренным, но не отключать проверку молча.
- Если стиль использует текстуру, цель вывода должна её поддерживать.

## Раскрытие стилей

Создать:

```text
external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/StyleResolver.kt
```

Контракт:

```kotlin
data class StyleResolutionResult<Action>(
    val element: UiElement<Action>,
    val diagnostics: List<StyleDiagnostic>,
)
```

`StyleResolver` должен:

- принимать `StyledUiElement`;
- принимать `UiTheme`;
- возвращать обычный `UiElement`;
- сохранять пути элементов для последующей диагностики;
- не заниматься расчётом координат;
- не генерировать Minecraft-код;
- не читать состояние игры.

Пример раскрытия кнопки:

```kotlin
StyledButton(
    style = theme.styles.button.primary,
    state = state,
    action = action,
    label = label,
)
```

раскрывается в:

```kotlin
button(
    modifier =
        Modifier
            .fillMaxSize()
            .background(resolvedSurface.fill),
    action = action,
) {
    text(
        modifier = Modifier.fillMaxSize().textAlign(resolvedLabel.alignment),
        text = label,
        color = resolvedLabel.color,
    )
}
```

Если цвет зависит от состояния, `StyleResolver` должен сохранить это как выражение значения, а не вычислять во время сборки.

## Связь со значениями состояния

Для цветов, поверхностей и состояний нужны выражения, совместимые с текущим `Value`.

Минимальный путь:

- оставить `Value<T>` как внешний механизм;
- добавить стилевые функции, которые принимают `Value<ButtonState>` и возвращают `Value<Color>`;
- в низкоуровневой программе это должно опускаться в `PrimitiveValueExpression.StateField`, если значение берётся из состояния сгенерированного экрана.

Если текущий `PrimitiveValueExpression` не умеет выбирать цвет по состоянию, добавить тип:

```kotlin
data class Match(
    val subject: PrimitiveValueExpression,
    val cases: Map<Any, PrimitiveValueExpression>,
    val default: PrimitiveValueExpression,
) : PrimitiveValueExpression
```

Генератор Kotlin должен превращать это в `when`.

## Темы

Создать тему библиотеки:

```text
external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/CreateLikeTheme.kt
```

Она не должна ссылаться на мод `create-thermodynamics`. Это общая тема для Create-подобных интерфейсов.

Создать тему мода:

```text
modules/ui/src/uiSpec/kotlin/dev/makargravanov/create_thermodynamics/ui/style/ThermodynamicsUiTheme.kt
```

Она может переиспользовать `CreateLikeTheme`, но задавать цвета и компоненты конкретного мода.

## Отчёт

Расширить текущие отчёты генерации:

```text
modules/ui/build/generated-ui/reports/
```

Добавить:

```text
reactor-controller-style-report.txt
reactor-controller-style-report.json
```

Содержимое текстового отчёта:

```text
Theme:
  name: ThermodynamicsUiTheme

Styles:
  window: used 1
  tab: used 3
  metricCard: used 6
  panel: used 4

Diagnostics:
  errors: 0
  warnings: 0

Optimization hints:
  PreferBakedTexture: 3 surfaces
  KeepPrimitiveCommands: 0 surfaces
```

Отчёт не должен быть декоративным. Если есть диагностическая ошибка, задача генерации должна падать.

## Задачи реализации

### Задача 1. Типы темы и стилей

Файлы:

- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/UiTheme.kt`
- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/UiTokens.kt`
- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/UiStyles.kt`
- Создать `external/ui-dsl/src/test/kotlin/ru/lazyhat/kraftui/style/UiStyleModelTest.kt`

Проверки:

- стиль кнопки без состояния `Disabled` не проходит проверку;
- стиль вкладки без состояния `Selected` не проходит проверку;
- `SurfaceStyle` с отрицательным отступом невозможно создать;
- `TextStyle` с `lineHeight <= 0` невозможно создать.

Команда:

```text
./gradlew :ui-dsl:test --tests "*UiStyleModelTest"
```

### Задача 2. Смысловые элементы

Файлы:

- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/styled/StyledUiElement.kt`
- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/styled/StyledUiScope.kt`
- Создать `external/ui-dsl/src/test/kotlin/ru/lazyhat/kraftui/styled/StyledUiScopeTest.kt`

Проверки:

- `window { tabs { ... } }` создаёт устойчивое дерево без ручных цветов;
- `metricCard` хранит заголовок и строки как значения, а не как уже отрисованные элементы;
- действие кнопки остаётся типизированным значением, а не строкой.

Команда:

```text
./gradlew :ui-dsl:test --tests "*StyledUiScopeTest"
```

### Задача 3. Проверка стилей

Файлы:

- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/StyleAnalysis.kt`
- Создать `external/ui-dsl/src/test/kotlin/ru/lazyhat/kraftui/style/StyleAnalysisTest.kt`

Проверки:

- отсутствующее состояние управляющего элемента возвращает `MissingControlState`;
- динамический текст с `FailInValidation` возвращает `UnsafeDynamicTextPolicy`;
- постоянный текст с низким контрастом возвращает `LowTextContrast`;
- текстура, неподдержанная целью вывода, возвращает `UnsupportedTextureForTarget`.

Команда:

```text
./gradlew :ui-dsl:test --tests "*StyleAnalysisTest"
```

### Задача 4. Раскрытие стилей в обычные элементы

Файлы:

- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/StyleResolver.kt`
- Создать `external/ui-dsl/src/test/kotlin/ru/lazyhat/kraftui/style/StyleResolverTest.kt`

Проверки:

- `StyledMetricCard` раскрывается в панель, заголовок и строки;
- `StyledTab` раскрывается в кнопку с выбранным состоянием;
- пути элементов сохраняются так, чтобы диагностика указывала на исходный смысловой элемент;
- результат можно передать в `UiLayoutResolver`.

Команда:

```text
./gradlew :ui-dsl:test --tests "*StyleResolverTest"
```

### Задача 5. Выражения выбора по состоянию

Файлы:

- Изменить `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/PrimitiveScreenProgram.kt`
- Изменить `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/MinecraftScreenSourceGenerator.kt`
- Изменить `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/PrimitiveScreenAnalysis.kt`
- Добавить тесты в `external/ui-dsl/src/test/kotlin/ru/lazyhat/kraftui/program/PrimitiveValueExpressionTest.kt`

Проверки:

- выражение `Match` с постоянным значением сворачивается оптимизатором;
- выражение `Match` по полю состояния генерирует `when` в Kotlin-коде;
- неподдержанный тип ключа в `Match` даёт диагностическую ошибку генерации;
- цвет вкладки меняется через состояние без ручного `If` вокруг фона.

Команды:

```text
./gradlew :ui-dsl:test --tests "*PrimitiveValueExpressionTest"
./gradlew :ui-dsl:test --tests "*MinecraftScreenSourceGeneratorTest"
```

### Задача 6. Общая Create-подобная тема

Файлы:

- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/CreateLikeTheme.kt`
- Создать `external/ui-dsl/src/test/kotlin/ru/lazyhat/kraftui/style/CreateLikeThemeTest.kt`

Проверки:

- тема проходит `StyleAnalysis` без ошибок;
- все стили компонентов заполнены;
- все динамические тексты в компонентах используют безопасную политику переполнения;
- тема не импортирует классы мода `create-thermodynamics`.

Команда:

```text
./gradlew :ui-dsl:test --tests "*CreateLikeThemeTest"
```

### Задача 7. Тема мода и перенос реакторного экрана

Файлы:

- Создать `modules/ui/src/uiSpec/kotlin/dev/makargravanov/create_thermodynamics/ui/style/ThermodynamicsUiTheme.kt`
- Изменить `modules/ui/src/uiSpec/kotlin/dev/makargravanov/create_thermodynamics/ui/reactor/ReactorControllerUi.kt`
- Изменить `modules/ui/src/uiSpec/kotlin/dev/makargravanov/create_thermodynamics/ui/reactor/GenerateReactorControllerMinecraftUi.kt`
- Изменить тесты в `modules/ui/src/test/kotlin/dev/makargravanov/create_thermodynamics/ui/reactor/`

Цель переноса:

- убрать локальный объект `Colors` из `ReactorControllerUi.kt`;
- заменить ручные `box(...background(...))` для типовых элементов на `window`, `tabButton`, `metricCard`, `panel`;
- оставить в экране только структуру: вкладки, сетки, карточки, данные состояния и действия.

Проверки:

- реакторный экран генерируется;
- предварительный просмотр совпадает с ожидаемым снимком;
- выбранная вкладка визуально отличается от невыбранной без ручного `If` в экране;
- отчёт стилей не содержит ошибок.

Команды:

```text
./gradlew :modules:ui:test
./gradlew :modules:ui:generateReactorControllerMinecraftUi
./gradlew :modules:ui-preview:renderUiPreviews
```

### Задача 8. Отчёт стилей и связь с оптимизатором

Файлы:

- Создать `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/style/StyleReport.kt`
- Изменить `external/ui-dsl/src/main/kotlin/ru/lazyhat/kraftui/program/PrimitiveScreenOptimizer.kt`
- Изменить `modules/ui/src/uiSpec/kotlin/dev/makargravanov/create_thermodynamics/ui/reactor/GenerateReactorControllerMinecraftUi.kt`
- Добавить тесты `external/ui-dsl/src/test/kotlin/ru/lazyhat/kraftui/style/StyleReportTest.kt`

Проверки:

- отчёт перечисляет использованные стили и количество использований;
- отчёт показывает стилевые поверхности с `BakeHint.PreferBakedTexture`;
- если поверхность не запеклась, отчёт объясняет причину: динамический цвет, условная видимость, малое число команд или неподдержанная операция;
- задача генерации пишет отчёт в `modules/ui/build/generated-ui/reports/`.

Команды:

```text
./gradlew :ui-dsl:test --tests "*StyleReportTest"
./gradlew :modules:ui:generateReactorControllerMinecraftUi
```

### Задача 9. Предварительный просмотр стилевых сценариев

Файлы:

- Изменить `modules/ui-preview/src/main/kotlin/dev/makargravanov/create_thermodynamics/ui/preview/UiPreviewMain.kt`
- Изменить `modules/ui-preview/src/test/kotlin/dev/makargravanov/create_thermodynamics/ui/preview/UiPreviewRendererTest.kt`

Сценарии:

- обзор реактора с нормальными строками;
- обзор реактора с длинными названиями веществ;
- экран зон с несколькими зонами;
- экран смеси с пустой смесью;
- экран смеси с длинными идентификаторами веществ.

Проверки:

- каждый сценарий создаёт картинку;
- отчёт не содержит ошибок переполнения;
- динамические строки используют многоточие или перенос согласно стилю;
- предварительный просмотр и Minecraft-генератор используют один и тот же `UiTheme`.

Команда:

```text
./gradlew :modules:ui-preview:renderUiPreviews :modules:ui-preview:test
```

### Задача 10. Архитектурная защита от возврата к ручным стилям

Файлы:

- Добавить тест в `modules/ui/src/test/kotlin/dev/makargravanov/create_thermodynamics/ui/ReactorUiArchitectureTest.kt`

Проверки:

- `ReactorControllerUi.kt` не содержит локальный объект `Colors`;
- `ReactorControllerUi.kt` не вызывает `background(Color(...))` напрямую для типовых элементов;
- `ReactorControllerUi.kt` не содержит строковых имён стилей;
- действия остаются типизированными, а не строковыми.

Команда:

```text
./gradlew :modules:ui:test --tests "*ReactorUiArchitectureTest"
```

## Финальная проверка

После всех задач:

```text
./gradlew :ui-dsl:test
./gradlew :modules:ui:test
./gradlew :modules:ui:generateReactorControllerMinecraftUi
./gradlew :modules:ui-preview:renderUiPreviews
./gradlew :modules:v1_21_1:v1_21_1-neoforge:compileKotlin
./gradlew :modules:v1_21_1:v1_21_1-neoforge:checkThinLoaderBoundary
```

Ожидаемый результат:

- все проверки проходят;
- реакторный экран собирается из смысловых компонентов;
- стиль задаётся через тему;
- отчёт показывает использованные стили, диагностику и результат оптимизаций;
- предварительный просмотр продолжает отображать тот же экран, что и Minecraft;
- новая система не требует запускать игру, чтобы найти ошибки стиля.

## Что не делать

- Не добавлять CSS-подобные строки и селекторы.
- Не делать глобальную изменяемую тему.
- Не оставлять старый ручной путь реакторного экрана как запасной.
- Не скрывать отсутствующие стили значениями по умолчанию.
- Не смешивать настройки стиля с состоянием реактора.
- Не превращать стиль в обёртку над цветами. Стиль должен описывать компонент целиком.
- Не добавлять новую отдельную низкоуровневую программу, если текущую `PrimitiveScreenProgram` можно расширить.

## Следующие слои после этого плана

- Девять-срезов для рамок и панелей.
- Иконки и символьные кнопки.
- Таблицы с закреплёнными заголовками.
- Поля ввода и числовые настройки.
- Подсказки и подробные всплывающие карточки.
- Поддержка нескольких тем в одном моде.
- Проверка доступности: контраст, читаемость, минимальный размер областей нажатия.
