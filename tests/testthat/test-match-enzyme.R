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

test_that("match_enzyme uses product alignment derived from enzyme rules", {
  enz <- make_enzyme(
    name = "TEST_CORE_GT",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "GalNAc(a1-",
      acceptor_alignment = "core",
      rejects = NULL,
      product = "Gal(b1-3)GalNAc(a1-"
    ))
  )
  glycan <- glyrepr::as_glycan_structure(
    "Gal(b1-3)GalNAc(a1-3)GalNAc(a1-"
  )

  expect_equal(match_enzyme(glycan, enz), list(integer()))
})

test_that("match_enzyme filters acceptor matches rejected by enzyme rules", {
  enz <- make_enzyme(
    name = "TEST_REJECT_GT",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-4)GlcNAc(b1-",
      acceptor_alignment = "terminal",
      rejects = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
      product = "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-"
    ))
  )
  glycan <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )

  expect_equal(match_enzyme(glycan, enz), list(integer()))
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

test_that("match_enzyme works with starter GTs", {
  glycan <- glyrepr::n_glycan_core()
  expect_equal(match_enzyme(glycan, "DPAGT1"), list(5L))
})

test_that("match_enzyme works for ALG13/ALG14", {
  glycan <- glyrepr::n_glycan_core()
  expect_equal(match_enzyme(glycan, "ALG13"), list(4L))
  expect_equal(match_enzyme(glycan, "ALG14"), list(4L))
})

test_that("match_enzyme works for ALG1", {
  glycan <- glyrepr::n_glycan_core()
  expect_equal(match_enzyme(glycan, "ALG1"), list(3L))
})

test_that("match_enzyme works for ALG2", {
  glycan <- glyrepr::n_glycan_core()
  expect_equal(match_enzyme(glycan, "ALG2"), list(c(1L, 2L)))
})

test_that("match_enzyme works for ALG11", {
  glycan <- glyrepr::as_glycan_structure(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(match_enzyme(glycan, "ALG11"), list(c(1L, 2L)))
})

test_that("match_enzyme works for ALG3", {
  glycan <- glyrepr::as_glycan_structure(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(match_enzyme(glycan, "ALG3"), list(3L))
})

test_that("match_enzyme works for ALG12", {
  glycan <- glyrepr::as_glycan_structure(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(match_enzyme(glycan, "ALG12"), list(4L))
})

test_that("match_enzyme works for ALG9", {
  glycan <- glyrepr::as_glycan_structure(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(match_enzyme(glycan, "ALG9"), list(c(7L, 9L)))
})

test_that("match_enzyme works for ALG6", {
  glycan <- glyrepr::as_glycan_structure(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(match_enzyme(glycan, "ALG6"), list(3L))
})

test_that("match_enzyme works for ALG8", {
  glycan <- glyrepr::as_glycan_structure(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(match_enzyme(glycan, "ALG8"), list(2L))
})

test_that("match_enzyme works for ALG10", {
  glycan <- glyrepr::as_glycan_structure(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(match_enzyme(glycan, "ALG10"), list(1L))
})
