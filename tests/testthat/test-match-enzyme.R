test_that("match_enzyme works for glyrepr_structure and glyenzy_enzyme", {
  glycan <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  )
  enzyme <- enzyme("ST3GAL3")

  expect_equal(match_enzyme(glycan, enzyme), list(1L))
})

test_that("match_enzyme doesn't work for characters", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  expect_error(match_enzyme(glycan, "ST3GAL3"))
})

test_that("match_enzyme works vectorizedly and preserves names", {
  glycans <- glyrepr::as_glycan_structure(c(
    sialylated = "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    unsialylated = "Gal(b1-3)GlcNAc(b1-"
  ))

  expect_equal(
    match_enzyme(glycans, "ST3GAL3"),
    list(sialylated = 1L, unsialylated = integer())
  )
})

test_that("match_enzyme pairs product and acceptor matches in the same glycan", {
  glycan <- glyrepr::as_glycan_structure(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )

  expect_equal(match_enzyme(glycan, "B4GALT1"), list(c(1L, 4L)))
})

test_that("match_enzyme only accepts glycosyltransferases", {
  glycan <- glyrepr::as_glycan_structure(
    "Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_error(
    match_enzyme(glycan, "MOGS"),
    "only supports glycosyltransferases"
  )
})
