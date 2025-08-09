# ===== Input Type =====
test_that("apply_enzyme takes `glyrepr_structure` and `glyenzy_enzyme` as inputs", {
  glycan <- glyparse::parse_iupac_condensed("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
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
  expect_error(apply_enzyme("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", 123), "`enzyme` must be")
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
  expect_equal(as.character(res[[1]]), "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
  expect_equal(as.character(res[[2]]), "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
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
  expect_equal(as.character(res[[2]]), "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
  expect_equal(
    sort(as.character(res[[3]])),
    sort(c(
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-6)[GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ))
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
  expect_equal(as.character(res[[2]]), "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
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