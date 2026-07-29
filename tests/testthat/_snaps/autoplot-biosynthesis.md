# figure dimensions are validated

    Code
      ggplot2::autoplot(graph, width = 0)
    Condition
      Error in `ggplot2::autoplot()`:
      ! `width` and `height` must be larger than 0.

---

    Code
      ggplot2::autoplot(graph, width = -1)
    Condition
      Error in `autoplot.glyenzy_biosynthesis_network()`:
      ! Assertion on 'width' failed: Element 1 is not >= 0.

---

    Code
      ggplot2::autoplot(graph, height = Inf)
    Condition
      Error in `autoplot.glyenzy_biosynthesis_network()`:
      ! Assertion on 'height' failed: Must be finite.

---

    Code
      ggplot2::autoplot(graph, units = "px")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "in", "cm", "mm"

---

    Code
      ggplot2::autoplot(graph, max_panel_width = 5)
    Condition
      Error in `.validate_biosynthesis_glycan_dots()`:
      ! The old panel-size argument `max_panel_width` has been removed.
      i Use `width` and `height` to set the complete figure size.

---

    Code
      ggplot2::autoplot(graph, max_panel_height = 5)
    Condition
      Error in `.validate_biosynthesis_glycan_dots()`:
      ! The old panel-size argument `max_panel_height` has been removed.
      i Use `width` and `height` to set the complete figure size.

# unrenderable enzyme labels are hidden with a warning

    Code
      plot <- ggplot2::autoplot(graph)
    Condition
      Warning:
      Enzyme labels cannot be rendered at the fitted plot size and have been hidden.
      i Increase `width` or `height` to display enzyme labels.
