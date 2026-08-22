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

test_that("find_enzyme batched motif matching agrees with enzyme-wise matching", {
  glycans <- glyparse::auto_parse(c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-8)Neu5Ac(a2-3)Gal(b1-4)Glc(b1-",
    "Gal3S6S(b1-4)GlcNAc(b1-"
  ))
  enzymes <- c("B4GALT1", "ST8SIA5", "GAL3ST2", "MGAT1", "MAN1A1", "ALG10")

  found <- find_enzyme(glycans, return_list = TRUE)
  observed <- vapply(
    enzymes,
    \(enzyme) vapply(found, \(x) enzyme %in% x, logical(1)),
    logical(length(glycans))
  )
  expected <- vapply(
    enzymes,
    \(enzyme) have_enzyme(glycans, enzyme),
    logical(length(glycans))
  )

  expect_identical(observed, expected)
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

test_that("find_enzyme path method includes N-glycan precursor GTs", {
  path_res <- suppressMessages(find_enzyme(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    method = "path"
  ))

  expect_true(all(c("ALG2", "ALG6", "ALG10") %in% path_res))
})

test_that("find_enzyme can use trace-derived path enzymes", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]GalNAc(a1-"

  motif_res <- find_enzyme(glycan)
  path_res <- suppressMessages(find_enzyme(glycan, method = "path"))

  expect_true(all(c("FUT3", "FUT4", "FUT9") %in% motif_res))
  expect_false(any(c("FUT3", "FUT4", "FUT9", "GALNT1") %in% path_res))
  expect_true(all(c("FUT5", "FUT6", "FUT7") %in% path_res))
})
