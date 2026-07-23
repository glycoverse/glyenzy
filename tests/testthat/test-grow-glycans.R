test_that("grow_glycans_step works for `glyrepr_structure` and `glyenzy_enzyme`", {
  glycan <- glyparse::parse_iupac_condensed(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzyme <- list(glyenzy::enzyme("MGAT3")) # Pass as list of enzyme objects
  result <- suppressMessages(grow_glycans_step(glycan, enzyme))
  expect_s3_class(result, "glyrepr_structure")
})

test_that("grow_glycans_step works for characters", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzyme <- "MGAT3"
  result <- suppressMessages(grow_glycans_step(glycan, enzyme))
  expect_s3_class(result, "glyrepr_structure")
})

test_that("grow_glycans works with progress bar", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")

  # Test that function works with progress bar (suppress messages to keep test clean)
  result <- suppressMessages(grow_glycans(glycans, enzymes, n_steps = 3))
  expect_s3_class(result, "glyrepr_structure")
  expect_true(length(result) > 0)
})

test_that("grow_glycans returns unique glycans", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")

  result <- suppressMessages(grow_glycans(glycans, enzymes, n_steps = 5))

  # Check that all glycans are unique
  result_strings <- as.character(result)
  expect_equal(length(result_strings), length(unique(result_strings)))
})

test_that("grow_glycans handles early termination", {
  # Use a glycan that cannot be further modified by the given enzyme
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-" # O-glycan
  enzyme <- "MGAT1" # N-glycan specific enzyme

  result <- suppressMessages(grow_glycans(glycan, enzyme, n_steps = 10))

  # Should return only the original glycan since no modifications are possible
  expect_equal(length(result), 1)
  expect_equal(as.character(result), glycan)
})

test_that("grow_glycans handles starter GTs", {
  glycan <- glyrepr::n_glycan_core()

  step_result <- suppressMessages(grow_glycans_step(glycan, "DPAGT1"))
  full_result <- suppressMessages(grow_glycans(glycan, "DPAGT1", n_steps = 1))

  expect_equal(step_result, glyrepr::glycan_structure())
  expect_equal(as.character(full_result), as.character(glycan))
})

# ===== Filter Parameter Tests =====
test_that("grow_glycans works without filter (NULL)", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")

  result <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 2,
    filter = NULL
  ))
  expect_s3_class(result, "glyrepr_structure")
  expect_true(length(result) >= length(glycans)) # Should include original + new glycans
})

test_that("grow_glycans works with function filter", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")

  # Filter to keep only structures with an N-glycan core motif.
  filter_func <- function(x) glymotif::have_motif(x, "GlcNAc(b1-4)GlcNAc(b1-")
  result <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 2,
    filter = filter_func
  ))

  expect_s3_class(result, "glyrepr_structure")
  expect_true(all(filter_func(result)))
})

test_that("grow_glycans works with formula filter", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")

  # Filter to keep only structures with an N-glycan core motif using formula syntax.
  filter_fn <- function(x) glymotif::have_motif(x, "GlcNAc(b1-4)GlcNAc(b1-")
  result <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 2,
    filter = ~ filter_fn(.x)
  ))

  expect_s3_class(result, "glyrepr_structure")
  expect_true(all(filter_fn(result)))
})

test_that("grow_glycans works with have_enzyme filter", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3", "MGAT2")

  # Filter to keep only glycans that could be synthesized by MGAT2
  result <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 3,
    filter = ~ have_enzyme(.x, "MGAT2")
  ))

  expect_s3_class(result, "glyrepr_structure")
  expect_true(length(result) > 0)
  # All results should be synthesized by MGAT2
  expect_true(all(have_enzyme(result, "MGAT2")))
})

test_that("grow_glycans filter reduces result size", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")

  # Run without filter
  result_no_filter <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 3
  ))

  # Run with restrictive filter (keep only glycans with exactly 3 mannose residues)
  result_with_filter <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 3,
    filter = ~ glyrepr::count_mono(.x, "Man") == 3
  ))

  expect_s3_class(result_no_filter, "glyrepr_structure")
  expect_s3_class(result_with_filter, "glyrepr_structure")
  expect_true(length(result_with_filter) <= length(result_no_filter))

  # All filtered results should have exactly 3 mannose residues
  if (length(result_with_filter) > 0) {
    expect_true(all(glyrepr::count_mono(result_with_filter, "Man") == 3))
  }
})

test_that("grow_glycans filter works with empty results", {
  glycans <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzymes <- "B4GALT1"

  # Use a filter that excludes everything (no glycans have 100 mannose residues)
  result <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 2,
    filter = ~ glyrepr::count_mono(.x, "Man") == 100
  ))

  expect_s3_class(result, "glyrepr_structure")
  expect_equal(length(result), 0)
})

test_that("grow_glycans filter is applied after each step", {
  glycans <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzymes <- c("B4GALT1", "ST3GAL3")

  # Filter to keep only glycans without sialic acid
  # This should prevent accumulation of sialylated glycans
  result <- suppressMessages(grow_glycans(
    glycans,
    enzymes,
    n_steps = 3,
    filter = ~ glyrepr::count_mono(.x, "Neu5Ac") == 0
  ))

  expect_s3_class(result, "glyrepr_structure")
  expect_true(length(result) > 0)
  # No result should contain sialic acid
  expect_true(all(glyrepr::count_mono(result, "Neu5Ac") == 0))
})

test_that("grow_glycans supports sequential ST actions", {
  st6 <- make_enzyme(
    name = "TEST_ST6",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-",
      acceptor_alignment = "whole",
      rejects = NULL,
      product = "Gal6S(b1-"
    ))
  )
  st3 <- make_enzyme(
    name = "TEST_ST3",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "Gal6S(b1-",
      acceptor_alignment = "whole",
      rejects = NULL,
      product = "Gal3S6S(b1-"
    ))
  )

  one_step <- grow_glycans_step("Gal(b1-", list(st6, st3))
  result <- suppressMessages(grow_glycans(
    "Gal(b1-",
    list(st6, st3),
    n_steps = 2
  ))

  expect_equal(as.character(one_step), "Gal6S(b1-")
  expect_true("Gal3S6S(b1-" %in% as.character(result))
})
