# Input type
test_that("is_synthesized_by works for `glyrepr_structure` and `glyenzy_enzyme`", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
  enzyme <- enzyme("ST3GAL3")
  expect_true(is_synthesized_by(glycan, enzyme), TRUE)
})

test_that("is_synthesized_by works for characters", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  enzyme <- "ST3GAL3"
  expect_true(is_synthesized_by(glycan, enzyme), TRUE)
})

test_that("is_synthesized_by works vectorizedly", {
  glycans <- c(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    "Gal(b1-3)GlcNAc(b1-"
  )
  expect_equal(is_synthesized_by(glycans, "ST3GAL3"), c(TRUE, FALSE))
})