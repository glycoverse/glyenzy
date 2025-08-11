# Test enzyme function
test_that("enzyme() returns a glyenzy_enzyme object", {
  expect_s3_class(enzyme("FUT8"), "glyenzy_enzyme")
})

# Test new_enzyme_rule function
test_that("new_enzyme_rule creates valid enzyme rule objects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))
  rejects_alignment <- character(0)

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects, rejects_alignment)

  expect_s3_class(rule, "glyenzy_enzyme_rule")
  expect_equal(rule$acceptor, acceptor)
  expect_equal(rule$product, product)
  expect_equal(rule$acceptor_alignment, "terminal")
  expect_equal(rule$rejects, rejects)
  expect_equal(rule$rejects_alignment, rejects_alignment)
})

# Test de novo synthesis (DPAGT1-like)
test_that("new_enzyme_rule handles de novo synthesis", {
  acceptor <- glyparse::parse_iupac_condensed(character(0))
  product <- glyparse::parse_iupac_condensed("GlcNAc(b1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))
  rejects_alignment <- character(0)

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, NULL, rejects, rejects_alignment)

  expect_s3_class(rule, "glyenzy_enzyme_rule")
  expect_equal(length(rule$acceptor), 0)
  expect_equal(rule$product, product)
  expect_null(rule$acceptor_alignment)
})

# Test enzyme with rejects
test_that("new_enzyme_rule handles rejects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(c("Gal(b1-4)GalNAc(a1-", "Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"))
  rejects_alignment <- c("terminal", "terminal")

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects, rejects_alignment)

  expect_s3_class(rule, "glyenzy_enzyme_rule")
  expect_equal(rule$rejects, rejects)
  expect_equal(rule$rejects_alignment, rejects_alignment)
})

# Test new_enzyme function
test_that("new_enzyme creates valid enzyme objects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))
  rejects_alignment <- character(0)

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects, rejects_alignment)
  enzyme <- glyenzy:::new_enzyme("ST3GAL2", list(rule), "GT", "human")

  expect_s3_class(enzyme, "glyenzy_enzyme")
  expect_equal(enzyme$name, "ST3GAL2")
  expect_equal(enzyme$rules, list(rule))
  expect_equal(enzyme$type, "GT")
  expect_equal(enzyme$species, "human")
})

# Test print method for glyenzy_enzyme
test_that("print.glyenzy_enzyme works correctly", {
  enzyme_obj <- enzyme("ST3GAL3")
  expect_snapshot(print(enzyme_obj))
})

# Test DPAGT1 enzyme (de novo synthesis)
test_that("DPAGT1 enzyme works correctly", {
  dpagt1 <- enzyme("DPAGT1")
  expect_s3_class(dpagt1, "glyenzy_enzyme")
  expect_equal(dpagt1$name, "DPAGT1")
  expect_equal(dpagt1$type, "GT")
  expect_snapshot(print(dpagt1))
})

# Test product_idx field for GT enzymes
test_that("product_idx is correctly calculated for GT enzymes", {
  # Test a regular GT enzyme (ST3GAL3)
  st3gal3 <- enzyme("ST3GAL3")
  rule <- st3gal3$rules[[1]]

  # Check that product_idx is set
  expect_true("product_idx" %in% names(rule))
  expect_type(rule$product_idx, "integer")
  expect_true(rule$product_idx > 0)

  # For ST3GAL3, the product should have one more residue than acceptor
  acceptor_size <- igraph::vcount(glyrepr::get_structure_graphs(rule$acceptor))
  product_size <- igraph::vcount(glyrepr::get_structure_graphs(rule$product))
  expect_equal(product_size, acceptor_size + 1)
})

# Test product_idx field for de novo synthesis
test_that("product_idx is correctly set for de novo synthesis", {
  dpagt1 <- enzyme("DPAGT1")
  rule <- dpagt1$rules[[1]]

  # For de novo synthesis, product_idx should be 1 (only one residue)
  expect_equal(rule$product_idx, 1)
  expect_equal(rule$acceptor_idx, 0)
})

# Test product_idx field for GH enzymes
test_that("product_idx is NULL for GH enzymes", {
  # Create a mock GH enzyme rule for testing
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))
  rejects_alignment <- character(0)

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects, rejects_alignment)
  enhanced_rule <- glyenzy:::enhance_enzyme_rule(rule, "GH")

  # For GH enzymes, product_idx should be NULL
  expect_null(enhanced_rule$product_idx)
  expect_type(enhanced_rule$acceptor_idx, "integer")
})