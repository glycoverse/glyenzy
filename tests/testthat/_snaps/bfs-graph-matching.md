# prepared graph matching preserves duplicate-reject errors

    Code
      trace_biosynthesis("GlcNAc(b1-3)GalNAc(a1-", enzymes = list(enzyme), max_steps = 1)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! `motifs` cannot have duplications.
      x Duplicate motifs: "Fuc(a1-2)GalNAc(a1-".
      i Consider using `unique()`.
