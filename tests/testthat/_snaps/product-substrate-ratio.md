# product_substrate_ratio validates its inputs

    Code
      product_substrate_ratio(SummarizedExperiment::SummarizedExperiment(assays = list(
        abundance = matrix(1))), "ST6GAL1")
    Condition
      Error in `.product_substrate_exp_type()`:
      ! `exp` must be a <GlycomicSE> or <GlycoproteomicSE> object.

---

    Code
      product_substrate_ratio(missing_structure, "ST6GAL1")
    Condition
      Error in `.check_product_substrate_exp()`:
      ! `exp` must contain glycan_structure in its row data.

---

    Code
      product_substrate_ratio(exp, "UNKNOWN")
    Condition
      Error in `.check_known_enzyme_names()`:
      ! Unknown enzymes: "UNKNOWN".

---

    Code
      product_substrate_ratio(exp, c("ST6GAL1", "MAN2A1"))
    Condition
      Error in `.product_substrate_enzymes()`:
      ! `product_substrate_ratio()` only supports glycosyltransferases and sulfotransferases.
      x Unsupported enzymes: "MAN2A1".

---

    Code
      product_substrate_ratio(exp, list(enzyme("ST6GAL1"), enzyme("MAN2A1")))
    Condition
      Error in `.product_substrate_enzymes()`:
      ! `product_substrate_ratio()` only supports glycosyltransferases and sulfotransferases.
      x Unsupported enzymes: "MAN2A1".

---

    Code
      product_substrate_ratio(exp, list(empty_enzyme))
    Condition
      Error in `.product_substrate_enzymes()`:
      ! Enzymes must contain at least one rule: "EMPTY".

---

    Code
      product_substrate_ratio(exp, "GALNT1")
    Condition
      Error in `.product_substrate_enzymes()`:
      ! Enzymes with non-glycan substrates are not supported.
      x Their substrates cannot be quantified: "GALNT1".

