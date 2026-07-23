# fallback counts every sulfate against both step budgets

    Code
      call(2L, 1L)
    Condition
      Error in `.abort_virtual_fallback()`:
      ! No synthesis path found for 1 target(s) within 2 steps.
      i Virtual fallback was limited to 1 step(s).

---

    Code
      call(1L, 2L)
    Condition
      Error in `.abort_virtual_fallback()`:
      ! No synthesis path found for 1 target(s) within 1 steps.
      i Virtual fallback was limited to 2 step(s).

# fallback respects the virtual-step limit

    Code
      call(1L)
    Condition
      Error in `.abort_virtual_fallback()`:
      ! No synthesis path found for 1 target(s) within 2 steps.
      i Virtual fallback was limited to 1 step(s).

# filters apply to virtual products

    Code
      path_biosynthesis("GalNAc(a1-", "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-", enzymes = "ST3GAL1",
        max_steps = 2, filter = reject_core_1, max_virtual_steps = 1)
    Condition
      Error in `.abort_virtual_fallback()`:
      ! No synthesis path found for 1 target(s) within 2 steps.
      i Virtual fallback was limited to 1 step(s).

# filters apply to virtual sulfate products

    Code
      path_biosynthesis("Gal(b1-3)GalNAc(a1-", "Gal3S(b1-3)GalNAc(a1-", enzymes = "C1GALT1",
        max_steps = 1, filter = reject_sulfate, max_virtual_steps = 1)
    Condition
      Error in `.abort_virtual_fallback()`:
      ! No synthesis path found for 1 target(s) within 1 steps.
      i Virtual fallback was limited to 1 step(s).
