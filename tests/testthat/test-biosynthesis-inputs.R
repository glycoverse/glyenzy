test_that("trace_biosynthesis rejects incompatible glycans", {
  expect_snapshot(
    trace_biosynthesis("Hex(b1-3)GalNAc(a1-"),
    error = TRUE
  )
  expect_snapshot(
    trace_biosynthesis(c(
      "Gal(b1-3)GalNAc(a1-",
      "Hex(b1-3)HexNAc(a1-"
    )),
    error = TRUE
  )
  expect_snapshot(
    trace_biosynthesis("Gal(b1-?)GalNAc(a1-"),
    error = TRUE
  )
  expect_snapshot(
    trace_biosynthesis(NA_character_),
    error = TRUE
  )
})

test_that("trace_biosynthesis_virtual rejects mixed structure levels", {
  expect_snapshot(
    trace_biosynthesis_virtual(c(
      "Gal(b1-3)GalNAc(a1-",
      "Gal(??-?)GalNAc(??-"
    )),
    error = TRUE
  )
})

test_that("trace_biosynthesis accepts homogeneous generic intact glycans", {
  glycans <- c(
    "HexNAc(a1-",
    "Hex(b1-3)HexNAc(a1-"
  )
  path <- suppressWarnings(trace_biosynthesis(
    glycans,
    enzymes = "C1GALT1",
    max_steps = 1
  ))

  expect_equal(igraph::vcount(path), 2L)
  expect_identical(igraph::V(path)$target, rep(TRUE, 2L))
})

test_that("path_biosynthesis requires compatible endpoint structures", {
  expect_snapshot(
    path_biosynthesis("GalNAc(a1-", "HexNAc(a1-"),
    error = TRUE
  )
  expect_snapshot(
    path_biosynthesis("GalNAc(?1-", "Gal(b1-3)GalNAc(?1-"),
    error = TRUE
  )
})

test_that("path_biosynthesis_virtual requires compatible endpoints", {
  expect_snapshot(
    path_biosynthesis_virtual("GalNAc(a1-", "GalNAc(??-"),
    error = TRUE
  )
})
