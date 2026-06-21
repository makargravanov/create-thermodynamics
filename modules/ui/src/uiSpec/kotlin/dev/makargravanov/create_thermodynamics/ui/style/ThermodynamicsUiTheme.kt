package dev.makargravanov.create_thermodynamics.ui.style

import ru.lazyhat.kraftui.foundation.Color
import ru.lazyhat.kraftui.foundation.modifier.TextOverflowPolicy
import ru.lazyhat.kraftui.style.CreateLikeTheme
import ru.lazyhat.kraftui.style.StyleColor
import ru.lazyhat.kraftui.style.TextStyle
import ru.lazyhat.kraftui.style.UiTheme

object ThermodynamicsUiTheme {
    val theme: UiTheme = CreateLikeTheme.theme

    val goodText: TextStyle =
        theme.styles.panel.body.copy(
            color = StyleColor.Constant(Color.rgb(35, 134, 78)),
            overflow = TextOverflowPolicy.Ellipsize,
        )
}
