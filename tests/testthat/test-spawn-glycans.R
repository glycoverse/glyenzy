test_that("spawn_glycans_step works for `glyrepr_structure` and `glyenzy_enzyme`", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
  enzyme <- list(glyenzy::enzyme("MGAT3"))  # Pass as list of enzyme objects
  result <- suppressMessages(spawn_glycans_step(glycan, enzyme))
  expect_s3_class(result, "glyrepr_structure")
})

test_that("spawn_glycans_step works for characters", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzyme <- "MGAT3"
  result <- suppressMessages(spawn_glycans_step(glycan, enzyme))
  expect_s3_class(result, "glyrepr_structure")
})

test_that("spawn_glycans works with progress bar", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")
  
  # Test that function works with progress bar (suppress messages to keep test clean)
  result <- suppressMessages(spawn_glycans(glycans, enzymes, n_steps = 3))
  expect_s3_class(result, "glyrepr_structure")
  expect_true(length(result) > 0)
})

test_that("spawn_glycans returns unique glycans", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzymes <- c("B4GALT1", "ST3GAL3")
  
  result <- suppressMessages(spawn_glycans(glycans, enzymes, n_steps = 5))
  
  # Check that all glycans are unique
  result_strings <- as.character(result)
  expect_equal(length(result_strings), length(unique(result_strings)))
})

test_that("spawn_glycans handles early termination", {
  # Use a glycan that cannot be further modified by the given enzyme
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"  # O-glycan
  enzyme <- "MGAT1"  # N-glycan specific enzyme
  
  result <- suppressMessages(spawn_glycans(glycan, enzyme, n_steps = 10))
  
  # Should return only the original glycan since no modifications are possible
  expect_equal(length(result), 1)
  expect_equal(as.character(result), glycan)
})
