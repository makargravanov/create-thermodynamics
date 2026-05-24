# Каталог и реестр

## Каталог

Каталог - это исходные данные Destroy: вещества, теги, явные реакции и свойства.

Код:

- `catalog.rs`
- `destroy_reactions.generated.rs`
- `reactions.rs`

Каталог сам по себе не рассчитывает смесь. Он только наполняет построитель реестра.

## Реестр

Реестр - это проверенная, готовая к расчету форма данных.

Код: `registry.rs`

При сборке реестр:

- проверяет вещества;
- проверяет реакции;
- проверяет, что все ссылки на вещества существуют;
- проверяет сохранение массы и заряда;
- проверяет обратные пары реакций;
- строит индекс реакций по веществам.

## Поток сборки

```mermaid
sequenceDiagram
    participant M as mod.rs
    participant C as catalog.rs
    participant R as reactions.rs
    participant B as ChemistryRegistryBuilder
    participant V as ChemistryRegistry

    M->>C: destroy_substances_registry_builder()
    C-->>M: builder с веществами и тегами
    M->>R: destroy_reactions_registry_builder(builder)
    R-->>M: builder с реакциями
    M->>B: build()
    B->>V: validate_substance_tags()
    B->>V: validate_reactions()
    B->>V: build_reaction_index()
    V-->>M: готовый реестр
```

## Каталог против реестра

Каталог отвечает на вопрос: “что известно из Destroy?”

Реестр отвечает на вопрос: “какие проверенные данные можно использовать в расчете?”

Это разделение важно: ошибка в данных должна проявиться при сборке реестра, а не во время тика смеси.
