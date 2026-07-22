# virtual tracing only accepts enzymes when annotating

    Code
      trace_biosynthesis_virtual("Gal(b1-3)GalNAc(a1-", enzymes = "C1GALT1")
    Condition
      Error in `.validate_virtual_enzymes()`:
      ! `enzymes` must be `NULL` unless `annotate_enzymes` is `TRUE`.
      i Set `annotate_enzymes = TRUE` to match known enzyme rules.

# virtual tracing reports unreachable paths and has no search controls

    Code
      path_biosynthesis_virtual("Man(a1-3)GlcNAc(b1-", target)
    Condition
      Error in `.perform_virtual_synthesis()`:
      ! No synthesis path found for 1 target(s) within 1 steps.

