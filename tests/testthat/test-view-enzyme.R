test_that("view_enzyme highlights residues matched by an enzyme", {
  glycan <- glyrepr::as_glycan_structure(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  )

  result <- suppressMessages(view_enzyme(glycan, "ST3GAL3"))
  expected <- glydraw::draw_cartoon(
    glycan,
    highlight = match_enzyme(glycan, "ST3GAL3")[[1]]
  )

  expect_s3_class(result, "glydraw_cartoon")
  expect_s3_class(result, "ggplot")
  expect_equal(
    ggplot2::ggplot_build(result)$data,
    ggplot2::ggplot_build(expected)$data
  )
})

test_that("view_enzyme returns an unhighlighted plot when no match is found", {
  glycan <- glyrepr::as_glycan_structure("Gal(b1-3)GlcNAc(b1-")

  expect_snapshot(result <- view_enzyme(glycan, "ST3GAL3"))

  expect_s3_class(result, "glydraw_cartoon")
  expect_s3_class(result, "ggplot")
})

test_that("view_enzyme accepts structure strings", {
  result <- suppressMessages(view_enzyme(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    "ST3GAL3"
  ))

  expect_s3_class(result, "glydraw_cartoon")
  expect_s3_class(result, "ggplot")
})

test_that("view_enzyme rejects multiple glycans", {
  glycans <- glyrepr::as_glycan_structure(c(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    "Gal(b1-3)GlcNAc(b1-"
  ))

  expect_error(
    view_enzyme(glycans, "ST3GAL3"),
    "Only one glycan can be visualized at a time"
  )
})
