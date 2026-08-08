test_that("enzymes_from_rnaseq() filters built-in enzymes by mean TPM", {
  mat <- rbind(
    FUT8 = c(0.5, 0.5),
    ST3GAL3 = c(1, 1),
    MOGS = c(0, 2),
    NOT_IN_GLYENZY = c(10, 10)
  )

  result <- enzymes_from_rnaseq(mat)

  expect_type(result, "list")
  expect_setequal(names(result), c("ST3GAL3", "MOGS"))
  expect_equal(
    unname(purrr::map_chr(result, "name")),
    names(result)
  )
  expect_true(all(purrr::map_lgl(result, inherits, "glyenzy_enzyme")))
})

test_that("enzymes_from_rnaseq() respects a custom threshold", {
  mat <- matrix(
    c(0.5, 1),
    ncol = 1,
    dimnames = list(c("FUT8", "ST3GAL3"), "sample")
  )

  expect_setequal(
    names(enzymes_from_rnaseq(mat, threshold = 0.5)),
    c("FUT8", "ST3GAL3")
  )
  expect_length(enzymes_from_rnaseq(mat, threshold = 2), 0)
})

test_that("enzymes_from_rnaseq() validates TPM input", {
  expect_snapshot(
    error = TRUE,
    enzymes_from_rnaseq(matrix(1, dimnames = list(NULL, "sample")))
  )
  expect_snapshot(
    error = TRUE,
    enzymes_from_rnaseq(matrix(-1, dimnames = list("FUT8", "sample")))
  )
  expect_snapshot(
    error = TRUE,
    enzymes_from_rnaseq(
      matrix(1, dimnames = list("FUT8", "sample")),
      threshold = -1
    )
  )
})
