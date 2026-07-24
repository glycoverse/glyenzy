test_that("functions check glycan monosaccharides", {
  glycan <- "Hex(b1-3)HexNAc(a1-"
  expect_error(
    have_enzyme(glycan, "B3GNT6"),
    "All glycans must have concrete monosaccharides"
  )
})

test_that("functions warn on non-intact glycan structures", {
  glycan <- "Gal(b1-?)GalNAc(a1-"
  expect_warning(
    expect_true(have_enzyme(glycan, "C1GALT1")),
    "non-intact glycan structures"
  )

  glycan <- "Gal(b1-3)GalNAc(?1-"
  expect_warning(
    expect_true(have_enzyme(glycan, "C1GALT1")),
    "non-intact glycan structures"
  )
})

test_that("functions allow sulfates and reject other substituents", {
  sulfate <- .process_glycans_arg("Gal6S(b1-3)GalNAc(a1-")
  expect_equal(as.character(sulfate), "Gal6S(b1-3)GalNAc(a1-")

  unsupported <- c(
    methyl = "Gal3Me(b1-3)GalNAc(a1-",
    phosphate = "Gal6P(b1-3)GalNAc(a1-",
    acetyl = "Gal3Ac(b1-3)GalNAc(a1-"
  )
  purrr::walk(unsupported, function(glycan) {
    expect_error(
      have_enzyme(glycan, "B3GNT6"),
      "Only sulfate substituents are supported"
    )
  })
})

test_that("substituent-subset matching is topology aware", {
  target <- glyparse::auto_parse(c(
    "Gal3S6S(b1-4)GlcNAc(b1-",
    "Gal3S(b1-4)GlcNAc(b1-",
    "Gal6S(b1-3)GlcNAc(b1-"
  ))

  expect_equal(
    .have_motif_substituent_subset(
      target,
      "Gal6S(b1-4)GlcNAc(b1-",
      alignment = "whole"
    ),
    c(TRUE, FALSE, FALSE)
  )
  expect_equal(
    .have_motif_substituent_subset(
      target,
      "Gal(b1-4)GlcNAc(b1-",
      alignment = "whole"
    ),
    c(TRUE, TRUE, FALSE)
  )
})

test_that("substituent token containment supports lenient positions", {
  expect_true(.substituent_tokens_contained("3S,6S", "6S", "strict"))
  expect_false(.substituent_tokens_contained("3S", "6S", "strict"))
  expect_true(.substituent_tokens_contained("?S", "6S", "lenient"))
  expect_false(.substituent_tokens_contained("?S", "6S", "strict"))
})

test_that("enzyme list processing accepts supported input forms", {
  expect_equal(
    names(.process_enzymes_arg(NULL)),
    NULL
  )

  enzyme_names <- c("B3GNT6", "C1GALT1")
  by_name <- .process_enzymes_arg(enzyme_names)
  expect_equal(
    purrr::map_chr(by_name, "name"),
    enzyme_names
  )

  by_object <- .process_enzymes_arg(purrr::map(enzyme_names, enzyme))
  expect_equal(
    purrr::map_chr(by_object, "name"),
    enzyme_names
  )
})

test_that("enzyme list processing rejects unknown enzyme names", {
  expect_error(
    .process_enzymes_arg(c("B3GNT6", "UNKNOWN_ENZYME")),
    "Unknown enzymes"
  )
})

test_that("enzyme list processing can prefilter against target glycans", {
  glycan <- glyparse::auto_parse("Gal(b1-3)GalNAc(a1-")

  filtered <- .process_enzymes_arg(
    c("C1GALT1", "ST6GAL1"),
    glycan = glycan,
    apply_prefilter = TRUE
  )
  unfiltered <- .process_enzymes_arg(
    c("C1GALT1", "ST6GAL1"),
    glycan = glycan,
    apply_prefilter = FALSE
  )

  expect_equal(purrr::map_chr(filtered, "name"), "C1GALT1")
  expect_equal(purrr::map_chr(unfiltered, "name"), c("C1GALT1", "ST6GAL1"))
})

test_that("enzyme prefiltering vectorizes targets with a scalar fallback", {
  calls <- integer()
  have_method <- function(glycans, enzyme) {
    calls <<- c(calls, length(glycans))
    rep(TRUE, length(glycans))
  }
  rlang::local_bindings(
    .have_enzyme_motif.test_prefilter_enzyme = have_method,
    .env = globalenv()
  )
  custom <- enzyme("C1GALT1")
  class(custom) <- c("test_prefilter_enzyme", class(custom))
  targets <- glyparse::auto_parse(c(
    "Gal(b1-3)GalNAc(a1-",
    "GlcNAc(b1-3)GalNAc(a1-"
  ))

  expect_equal(.can_enzymes_contribute(list(custom), targets), TRUE)
  expect_equal(calls, 2L)

  calls <- integer()
  fallback_method <- function(glycans, enzyme) {
    calls <<- c(calls, length(glycans))
    if (length(glycans) > 1L) {
      stop("scalar only")
    }
    as.character(glycans) == "GlcNAc(b1-3)GalNAc(a1-"
  }
  rlang::local_bindings(
    .have_enzyme_motif.test_prefilter_enzyme = fallback_method,
    .env = globalenv()
  )

  expect_equal(.can_enzymes_contribute(list(custom), targets), TRUE)
  expect_equal(calls, c(2L, 1L, 1L))
})

test_that("enzyme prefiltering rejects scalar-length vectorized results", {
  calls <- integer()
  scalar_method <- function(glycans, enzyme) {
    calls <<- c(calls, length(glycans))
    as.character(glycans[1]) == "GlcNAc(b1-3)GalNAc(a1-"
  }
  rlang::local_bindings(
    .have_enzyme_motif.test_scalar_prefilter_enzyme = scalar_method,
    .env = globalenv()
  )
  custom <- enzyme("C1GALT1")
  class(custom) <- c("test_scalar_prefilter_enzyme", class(custom))
  targets <- glyparse::auto_parse(c(
    "Gal(b1-3)GalNAc(a1-",
    "GlcNAc(b1-3)GalNAc(a1-"
  ))

  expect_equal(.can_enzymes_contribute(list(custom), targets), TRUE)
  expect_equal(calls, c(2L, 1L, 1L))
})

test_that("enzyme list processing errors when no enzyme can contribute", {
  glycan <- glyparse::auto_parse("Gal(b1-3)GalNAc(a1-")

  expect_error(
    .process_enzymes_arg("ST6GAL1", glycan = glycan, apply_prefilter = TRUE),
    "No enzymes are predicted to contribute"
  )
})
