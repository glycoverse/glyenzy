# abstract enzymes remain absent from existing discovery interfaces

    Code
      enzyme("ManI")
    Condition
      Error in `enzyme()`:
      ! Unknown enzyme: "ManI".

---

    Code
      .enzymes_from_arg(abstract_names)
    Condition
      Error in `.check_known_enzyme_names()`:
      ! Unknown enzymes: "ManI", "GnTI", "ManII", "GnTII", "GnTIII", "a6FucT", "GnTIV", "GnTV", "b4GalT", "iGnT", "a3SiaT", and "GlcH".

---

    Code
      abstract_enzymes(NA)
    Condition
      Error in `abstract_enzymes()`:
      ! Assertion on 'return_str' failed: May not be NA.

# printing abstract enzymes includes localization and rejects

    Code
      print(abstract_enzymes()$ManII)
    Message
      -- Enzyme: ManII ---------------------------------------------------------------
      i Type: "GH" (Glycoside hydrolase)
      i Species: "unspecified"
      i Glycan type: "N"
      i Abstract enzyme; localization: "medial"
      -- Rules (1) --
      > Rule 1: core alignment
      Acceptor:
      "GlcNAc(b1-2)Man(a1-3)[Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
      Product: "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
      Rejects:
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-3/6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

---

    Code
      print(abstract_enzymes()$iGnT)
    Message
      -- Enzyme: iGnT ----------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "unspecified"
      i Glycan type: "N"
      i Abstract enzyme; localization: "trans"
      -- Rules (1) --
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(?1-"
      Product: "GlcNAc(b1-3)Gal(b1-4)GlcNAc(?1-"
      Rejects:
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)Man(b1-"
      "Gal(b1-4)GlcNAc(b1-4)Man(a1-3)Man(b1-"
      "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-"

# abstract tracing retains precursor enzymes and rejects impossible routes

    Code
      trace_biosynthesis(cases$glycan[cases$name == "Man5"], enzymes = enzymes[names(
        enzymes) != "GlcH"], max_virtual_steps = 0L)
    Condition
      Error in `engine$run()`:
      ! No synthesis path found for 1 target(s) within 11 steps.

---

    Code
      trace_biosynthesis(wrong_arm, enzymes = enzymes, max_virtual_steps = 0L)
    Condition
      Error in `engine$run()`:
      ! No synthesis path found for 1 target(s) within 15 steps.
