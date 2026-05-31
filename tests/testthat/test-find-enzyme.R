# ===== Input type =====
test_that("find_enzyme works for `glyrepr_structure`", {
  glycan <- glyparse::parse_iupac_condensed(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_true("MGAT1" %in% find_enzyme(glycan))
})

test_that("find_enzyme works for characters", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_true("MGAT1" %in% find_enzyme(glycan))
})

test_that("find_enzyme works for vectorized inputs", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_false("MGAT2" %in% find_enzyme(glycans)[[1]])
  expect_true("MGAT2" %in% find_enzyme(glycans)[[2]])
})

test_that("find_enzyme rejects invalid inputs", {
  expect_error(find_enzyme(123), "`glycans` must be")
})

test_that("find_enzyme handles return_list = TRUE correctly", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_true("MGAT1" %in% find_enzyme(glycan, return_list = TRUE)[[1]])
})

test_that("find_enzyme raises error when return_list = FALSE and length(glycans) > 1", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_error(find_enzyme(glycans, return_list = FALSE), "must have length 1")
})

test_that("find_enzyme considers starter GTs", {
  glycan <- glyrepr::n_glycan_core()
  expect_true("DPAGT1" %in% find_enzyme(glycan))
})