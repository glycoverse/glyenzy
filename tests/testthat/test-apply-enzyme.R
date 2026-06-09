# ===== Input Type =====
test_that("apply_enzyme takes `glyrepr_structure` and `glyenzy_enzyme` as inputs", {
  glycan <- glyparse::parse_iupac_condensed(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzyme <- enzyme("MGAT3")
  expect_s3_class(apply_enzyme(glycan, enzyme), "glyrepr_structure")
})

test_that("apply_enzyme takes characters as inputs", {
  glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzyme <- "MGAT3"
  expect_s3_class(apply_enzyme(glycan, enzyme), "glyrepr_structure")
})

test_that("apply_enzyme rejects invalid inputs", {
  expect_error(apply_enzyme(123, "MGAT3"), "`glycans` must be")
  expect_error(
    apply_enzyme(
      "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      123
    ),
    "`enzyme` must be"
  )
})

test_that("apply_enzyme works vectorizedly", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_type(apply_enzyme(glycans, "MGAT3"), "list")
})

# ===== Normal Cases for GT Enzymes =====
test_that("apply_enzyme works for MGAT3", {
  glycans <- c(
    # GOOD: a perfect substrate
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # GOOD: also a substrate
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # BAD: does not have the acceptor
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # BAD: has the acceptor but also has a reject
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # BAD: already has the product
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  res <- apply_enzyme(glycans, "MGAT3")

  expect_equal(length(res), 5L)
  expect_equal(
    as.character(res[[1]]),
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[[2]]),
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(res[[3]], glyrepr::glycan_structure())
  expect_equal(res[[4]], glyrepr::glycan_structure())
  expect_equal(res[[5]], glyrepr::glycan_structure())
})

test_that("apply_enzyme works for B4GALT1", {
  glycans <- c(
    # zero antenna
    "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # one antenna
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # two antennas
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  res <- apply_enzyme(glycans, "B4GALT1")

  expect_equal(length(res), 3L)
  expect_equal(res[[1]], glyrepr::glycan_structure())
  expect_equal(
    as.character(res[[2]]),
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(
    sort(as.character(res[[3]])),
    sort(c(
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-6)[GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ))
  )
})

test_that("apply_enzyme works with special reject rules", {
  # The N-glycan has two galactosylated antennas, one of which has a a1-3 fucose
  glycan <- "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  # And we define an enzyme that only afucosylated Gal can be sialylated
  rule <- new_enzyme_rule(
    acceptor = glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-"),
    product = glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"),
    acceptor_alignment = "terminal",
    rejects = glyrepr::as_glycan_structure("Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  )
  enzyme <- new_enzyme("TEST_ENZYME", list(rule), "GT", "human")
  enzyme <- enhance_enzyme(enzyme)
  res <- apply_enzyme(glycan, enzyme)
  # In the result, only one glycan should be generated, with the afucosylated Gal being sialylated
  expect_equal(
    as.character(res),
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
})

# ===== Normal Cases for GH Enzymes =====
test_that("apply_enzyme works for MAN2A1", {
  glycans <- c(
    # GOOD: a perfect substrate
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # GOOD: also a substrate
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # BAD: does not have the acceptor
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # BAD: has the reject
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  res <- apply_enzyme(glycans, "MAN2A1")

  expect_equal(length(res), 4L)
  expect_equal(
    sort(as.character(res[[1]])),
    sort(c(
      "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ))
  )
  expect_equal(
    as.character(res[[2]]),
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(res[[3]], glyrepr::glycan_structure())
  expect_equal(res[[4]], glyrepr::glycan_structure())
})

# ===== Multiple Rules =====
test_that("apply_enzyme works for ST6GAL2", {
  glycan <- "GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  res <- apply_enzyme(glycan, "ST6GAL2")
  expect_equal(
    sort(as.character(res)),
    sort(c(
      "Neu5Ac(a2-6)GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)[GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ))
  )
})

# ===== Starter GTs =====
test_that("apply_enzyme is compatible with starter GTs", {
  glycan <- glyrepr::n_glycan_core()
  res <- apply_enzyme(glycan, "DPAGT1")
  # Applying a starter GT on any glycan should always result in nothing
  expect_equal(res, glyrepr::glycan_structure())
})

test_that("apply_enzyme keeps return shape for starter GTs", {
  glycans <- glyparse::parse_iupac_condensed(c(
    "GlcNAc(b1-",
    "GalNAc(a1-"
  ))

  res <- apply_enzyme(glycans, "DPAGT1", return_list = TRUE)

  expect_type(res, "list")
  expect_length(res, length(glycans))
  purrr::walk(res, ~ expect_equal(.x, glyrepr::glycan_structure()))
  expect_error(
    apply_enzyme(glycans, "DPAGT1", return_list = FALSE),
    "must have length 1"
  )
})

# ===== N-glycan precursor GTs =====
test_that("apply_enzyme is compatible with N-glycan precursor GTs", {
  glycan <- glyrepr::n_glycan_core()
  res <- apply_enzyme(glycan, "ALG2")
  # Applying a starter GT on any glycan should always result in nothing
  expect_equal(res, glyrepr::glycan_structure())
})

# ===== FUT 3/4/5/6/7/9 =====
test_that("fucosylation of A-antigen (Type 1)", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-3)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), character(0))
  expect_equal(as.character(res[["FUT5"]]), character(0))
  expect_equal(as.character(res[["FUT6"]]), character(0))
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

test_that("fucosylation of H (Type 1)", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Fuc(a1-2)Gal(b1-3)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Fuc(a1-2)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), character(0))
  expect_equal(as.character(res[["FUT5"]]), character(0))
  expect_equal(as.character(res[["FUT6"]]), character(0))
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

test_that("fucosylation of B-antigen (Type 1)", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Fuc(a1-2)[Gal(a1-3)]Gal(b1-3)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Fuc(a1-2)[Gal(a1-3)]Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), character(0))
  expect_equal(as.character(res[["FUT5"]]), character(0))
  expect_equal(as.character(res[["FUT6"]]), character(0))
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

test_that("fucosylation of Type 1 chain", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Gal(b1-3)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), character(0))
  expect_equal(as.character(res[["FUT5"]]), character(0))
  expect_equal(as.character(res[["FUT6"]]), character(0))
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

test_that("fucosylation of sialylated Type 1 chain", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), character(0))
  expect_equal(as.character(res[["FUT5"]]), character(0))
  expect_equal(as.character(res[["FUT6"]]), character(0))
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

test_that("fucosylation of Type 2 chain", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Gal(b1-4)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT5"]]), "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT6"]]), "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-")
})

test_that("fucosylation of H (Type 2)", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Fuc(a1-2)Gal(b1-4)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT5"]]), "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT6"]]), "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
})

test_that("fucosylation of A antigen (Type 2)", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), character(0))
  expect_equal(as.character(res[["FUT5"]]), character(0))
  expect_equal(as.character(res[["FUT6"]]), character(0))
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

test_that("fucosylation of poly-LacNAc", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT4"]]), "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT5"]]), "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT6"]]), "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

test_that("fucosylation of sialyl-LacNAc", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(as.character(res[["FUT3"]]), character(0))
  expect_equal(as.character(res[["FUT4"]]), character(0))
  expect_equal(as.character(res[["FUT5"]]), "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT6"]]), "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT7"]]), "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-")
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

# ===== Regression Tests =====
test_that("apply_enzyme regression: GH enzymes do not create invalid out-tree structures", {
  # Test that GH enzymes only remove terminal residues and do not break tree structure
  # This prevents the "Glycan structure must be an out tree" error

  # Use a high-mannose structure where some GH enzymes might try to remove non-terminal residues
  glycan <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  # Test various GH enzymes that trim mannose residues
  gh_enzymes <- c(
    "MAN1A1",
    "MAN1A2",
    "MAN1C1",
    "MAN2A1",
    "MAN2A2",
    "GANAB",
    "MAN1B1"
  )

  for (enzyme_name in gh_enzymes) {
    # Should not throw "out tree" validation errors
    expect_no_error({
      res <- suppressMessages(apply_enzyme(glycan, enzyme_name))
      # Result should be a valid glyrepr_structure (possibly empty if no valid products)
      expect_s3_class(res, "glyrepr_structure")
    })
  }
})
