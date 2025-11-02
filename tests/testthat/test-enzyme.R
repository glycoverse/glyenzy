# Test enzyme function
test_that("enzyme() returns a glyenzy_enzyme object", {
  expect_s3_class(enzyme("FUT8"), "glyenzy_enzyme")
})

test_that("enzyme() throws error for unknown enzyme", {
  expect_error(
    enzyme("UNKNOWN_ENZYME"),
    "Unknown enzyme"
  )
})

# Test new_enzyme_rule function
test_that("new_enzyme_rule creates valid enzyme rule objects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_s3_class(rule, "glyenzy_enzyme_rule")
  expect_equal(rule$acceptor, acceptor)
  expect_equal(rule$product, product)
  expect_equal(rule$acceptor_alignment, "terminal")
  expect_equal(rule$rejects, rejects)
})

# Test enzyme with rejects
test_that("new_enzyme_rule handles rejects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(c("Gal(b1-4)GalNAc(a1-", "Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_s3_class(rule, "glyenzy_enzyme_rule")
  expect_equal(rule$rejects, rejects)
})

# Test new_enzyme function
test_that("new_enzyme creates valid enzyme objects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)
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

# Test acceptor_idx field for GT enzymes
test_that("acceptor_idx is correctly calculated for GT enzymes", {
  # Test a regular GT enzyme (ST3GAL3)
  st3gal3 <- enzyme("ST3GAL3")
  rule <- st3gal3$rules[[1]]

  # Check that acceptor_idx is set
  expect_true("acceptor_idx" %in% names(rule))
  expect_type(rule$acceptor_idx, "integer")
  expect_true(rule$acceptor_idx > 0)

  # acceptor_idx should be within the range of acceptor nodes
  acceptor_size <- igraph::vcount(glyrepr::get_structure_graphs(rule$acceptor, return_list = FALSE))
  expect_true(rule$acceptor_idx <= acceptor_size)
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
  acceptor_size <- igraph::vcount(glyrepr::get_structure_graphs(rule$acceptor, return_list = FALSE))
  product_size <- igraph::vcount(glyrepr::get_structure_graphs(rule$product, return_list = FALSE))
  expect_equal(product_size, acceptor_size + 1)

  # product_idx should be within the range of product nodes
  expect_true(rule$product_idx <= product_size)
})

# Test acceptor_idx and product_idx fields for GH enzymes
test_that("acceptor_idx and product_idx are correctly set for GH enzymes", {
  # Create a mock GH enzyme rule for testing
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)
  enhanced_rule <- glyenzy:::enhance_enzyme_rule(rule, "GH")

  # For GH enzymes, acceptor_idx should point to the removed residue
  expect_type(enhanced_rule$acceptor_idx, "integer")
  expect_true(enhanced_rule$acceptor_idx > 0)

  # acceptor_idx should be within the range of acceptor nodes
  acceptor_size <- igraph::vcount(glyrepr::get_structure_graphs(acceptor, return_list = FALSE))
  expect_true(enhanced_rule$acceptor_idx <= acceptor_size)

  # For GH enzymes, product_idx should be NULL (no new residue added)
  expect_null(enhanced_rule$product_idx)
})

# Test all_enzymes function
test_that("all_enzymes() returns enzyme list by default", {
  enzymes <- all_enzymes()
  expect_type(enzymes, "list")
  expect_true(all(purrr::map_lgl(enzymes, ~ inherits(.x, "glyenzy_enzyme"))))
  expect_true(length(enzymes) > 0)
})

test_that("all_enzymes(return_str = TRUE) returns character vector", {
  enzyme_names <- all_enzymes(return_str = TRUE)
  expect_type(enzyme_names, "character")
  expect_true(length(enzyme_names) > 0)
  expect_true("ST3GAL3" %in% enzyme_names)
})

test_that("all_enzymes(return_str = FALSE) returns enzyme list", {
  enzymes <- all_enzymes(return_str = FALSE)
  expect_type(enzymes, "list")
  expect_true(all(purrr::map_lgl(enzymes, ~ inherits(.x, "glyenzy_enzyme"))))
})

# Test validate_enzyme function
test_that("validate_enzyme works with valid enzyme", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)
  enzyme_obj <- glyenzy:::new_enzyme("TEST_ENZYME", list(rule), "GT", "human")

  expect_invisible(glyenzy:::validate_enzyme(enzyme_obj))
})

test_that("validate_enzyme fails with invalid class", {
  expect_error(
    glyenzy:::validate_enzyme(list(name = "test")),
    "Assertion on 'x' failed"
  )
})

# Test validate_enzyme_rule function
test_that("validate_enzyme_rule works with valid GT rule", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_invisible(glyenzy:::validate_enzyme_rule(rule, "GT"))
})

test_that("validate_enzyme_rule works with valid GH rule", {
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_invisible(glyenzy:::validate_enzyme_rule(rule, "GH"))
})

test_that("validate_enzyme_rule fails with multiple acceptors", {
  acceptor <- glyparse::parse_iupac_condensed(c("Gal(b1-3)GalNAc(a1-", "Gal(b1-4)GalNAc(a1-"))
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_error(
    glyenzy:::validate_enzyme_rule(rule, "GT"),
    "acceptor.*must be a single structure"
  )
})

test_that("validate_enzyme_rule fails with multiple products", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed(c("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-", "Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"))
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_error(
    glyenzy:::validate_enzyme_rule(rule, "GT"),
    "product.*must be a single structure"
  )
})

test_that("validate_enzyme_rule fails when acceptor is not substructure of product for GT", {
  acceptor <- glyparse::parse_iupac_condensed("Man(a1-3)GalNAc(a1-")  # Different structure
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_error(
    glyenzy:::validate_enzyme_rule(rule, "GT"),
    "acceptor.*must be a substructure"
  )
})

test_that("validate_enzyme_rule fails when product is not substructure of acceptor for GH", {
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Man(a1-3)GalNAc(a1-")  # Different structure
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)

  expect_error(
    glyenzy:::validate_enzyme_rule(rule, "GH"),
    "product.*must be a substructure"
  )
})

# Test enhance_enzyme function
test_that("enhance_enzyme enhances all rules", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)
  enzyme_obj <- glyenzy:::new_enzyme("TEST_ENZYME", list(rule), "GT", "human")

  enhanced_enzyme <- glyenzy:::enhance_enzyme(enzyme_obj)

  res_rule <- enhanced_enzyme$rules[[1]]
  expect_s3_class(enhanced_enzyme, "glyenzy_enzyme")
  expect_equal(res_rule$acceptor_idx, 1)
  expect_equal(res_rule$product_idx, 1)
  expect_equal(res_rule$new_residue, "Neu5Ac")
  expect_equal(res_rule$new_linkage, "a2-3")
})

# Test enhance_enzyme_rule function
test_that("enhance_enzyme_rule works for GT enzymes", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)
  enhanced_rule <- glyenzy:::enhance_enzyme_rule(rule, "GT")

  expect_true("acceptor_idx" %in% names(enhanced_rule))
  expect_true("product_idx" %in% names(enhanced_rule))
  expect_true("new_residue" %in% names(enhanced_rule))
  expect_true("new_linkage" %in% names(enhanced_rule))
  expect_type(enhanced_rule$acceptor_idx, "integer")
  expect_type(enhanced_rule$product_idx, "integer")
  expect_type(enhanced_rule$new_residue, "character")
  expect_type(enhanced_rule$new_linkage, "character")
})

test_that("enhance_enzyme_rule works for GH enzymes", {
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(character(0))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)
  enhanced_rule <- glyenzy:::enhance_enzyme_rule(rule, "GH")

  expect_true("acceptor_idx" %in% names(enhanced_rule))
  expect_type(enhanced_rule$acceptor_idx, "integer")
  # For GH enzymes, these fields should be NULL
  expect_null(enhanced_rule$product_idx)  # GH enzymes don't add residues
  expect_null(enhanced_rule$new_residue)
  expect_null(enhanced_rule$new_linkage)
})

# Test print method edge cases
test_that("print.glyenzy_enzyme works with enzyme with no rules", {
  enzyme_obj <- glyenzy:::new_enzyme("EMPTY_ENZYME", list(), "GT", "human")
  expect_snapshot(print(enzyme_obj))
})

test_that("print.glyenzy_enzyme works with enzyme with rejects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  rejects <- glyparse::parse_iupac_condensed(c("Gal(b1-4)GalNAc(a1-", "Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"))

  rule <- glyenzy:::new_enzyme_rule(acceptor, product, "terminal", rejects)
  enzyme_obj <- glyenzy:::new_enzyme("TEST_ENZYME", list(rule), "GT", "human")

  expect_snapshot(print(enzyme_obj))
})