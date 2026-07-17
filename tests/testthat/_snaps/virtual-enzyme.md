# virtual tracing validates incompatible options

    Code
      trace_biosynthesis("Gal(b1-3)GalNAc(a1-", enzymes = "C1GALT1", method = "virtual")
    Condition
      Error in `.validate_virtual_enzymes()`:
      ! `enzymes` must be `NULL` when `method = "virtual"`.
      i Virtual-enzyme tracing does not use known enzyme rules.

# virtual tracing reports unreachable paths

    Code
      trace_biosynthesis(target, method = "virtual", max_steps = 1)
    Condition
      Error in `.perform_virtual_synthesis()`:
      ! No synthesis path found for 1 target(s) within 1 steps.

---

    Code
      trace_biosynthesis(target, method = "virtual", filter = function(glycans)
        length(glycans) == 0L)
    Condition
      Error in `.perform_virtual_synthesis()`:
      ! No synthesis path found for 1 target(s) within 20 steps.

---

    Code
      path_biosynthesis("Man(a1-3)GlcNAc(b1-", target, method = "virtual")
    Condition
      Error in `.perform_virtual_synthesis()`:
      ! No synthesis path found for 1 target(s) within 10 steps.
