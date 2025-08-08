# ===== Input type =====
test_that("get_involved_enzymes works for `glyrepr_structure`", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
  expect_true("MGAT1" %in% get_involved_enzymes(glycan))
})

test_that("get_involved_enzymes works for characters", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_true("MGAT1" %in% get_involved_enzymes(glycan))
})

test_that("get_involved_enzymes works for vectorized inputs", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_false("MGAT2" %in% get_involved_enzymes(glycans)[[1]])
  expect_true("MGAT2" %in% get_involved_enzymes(glycans)[[2]])
})

test_that("get_involved_enzymes rejects invalid inputs", {
  expect_error(get_involved_enzymes(123), "`glycans` must be")
})

test_that("get_involved_enzymes handles return_list = TRUE correctly", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_true("MGAT1" %in% get_involved_enzymes(glycan, return_list = TRUE)[[1]])
})

test_that("get_involved_enzymes raises error when return_list = FALSE and length(glycans) > 1", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_error(get_involved_enzymes(glycans, return_list = FALSE), "must have length 1")
})
