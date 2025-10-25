test_that("functions check glycan monosaccharides", {
  glycan <- "Hex(b1-3)HexNAc(a1-"
  expect_error(
    is_synthesized_by(glycan, "B3GNT6"),
    "All glycans must have concrete monosaccharides"
  )
})

test_that("functions check glycan linkages", {
  glycan <- "Gal(b1-?)GalNAc(a1-"
  expect_error(
    is_synthesized_by(glycan, "B3GNT6"),
    "All linkages must be intact"
  )

  glycan <- "Gal(b1-3)GalNAc(?1-"
  expect_error(
    is_synthesized_by(glycan, "B3GNT6"),
    "All linkages must be intact"
  )
})

test_that("functions check glycan substituents", {
  glycan <- "Gal3Me(b1-3)GalNAc(a1-"
  expect_error(
    is_synthesized_by(glycan, "B3GNT6"),
    "Glycans with substituents are not supported"
  )
})