# enzymes_from_rnaseq() validates TPM input

    Code
      enzymes_from_rnaseq(matrix(1, dimnames = list(NULL, "sample")))
    Condition
      Error in `enzymes_from_rnaseq()`:
      ! Assertion on 'rownames(mat)' failed: Must have names.

---

    Code
      enzymes_from_rnaseq(matrix(-1, dimnames = list("FUT8", "sample")))
    Condition
      Error in `enzymes_from_rnaseq()`:
      ! Assertion on 'as.vector(mat)' failed: Element 1 is not >= 0.

---

    Code
      enzymes_from_rnaseq(matrix(1, dimnames = list("FUT8", "sample")), threshold = -
      1)
    Condition
      Error in `enzymes_from_rnaseq()`:
      ! Assertion on 'threshold' failed: Element 1 is not >= 0.

