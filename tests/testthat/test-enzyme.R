# Test new_enzyme_rule function
test_that("new_enzyme_rule creates valid enzyme rule objects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")

  expect_s3_class(rule, "glyenzy_enzyme_rule")
  expect_equal(rule$acceptor, acceptor)
  expect_equal(rule$product, product)
  expect_equal(rule$acceptor_alignment, "terminal")
  expect_equal(rule$type, "GT")
})

test_that("new_enzyme_rule validates input parameters", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")

  # Invalid acceptor class
  expect_error(
    new_enzyme_rule("invalid", product, "terminal", "GT"),
    "Assertion on 'acceptor' failed"
  )

  # Invalid product class
  expect_error(
    new_enzyme_rule(acceptor, "invalid", "terminal", "GT"),
    "Assertion on 'product' failed"
  )

  # Invalid alignment
  expect_error(
    new_enzyme_rule(acceptor, product, "invalid", "GT"),
    "Assertion on 'acceptor_alignment' failed"
  )

  # Invalid type
  expect_error(
    new_enzyme_rule(acceptor, product, "terminal", "invalid"),
    "Assertion on 'type' failed"
  )
})

test_that("validate_enzyme_rule fails for multiple structures", {
  acceptor <- glyparse::parse_iupac_condensed(c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)GalNAc(a1-"))
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  expect_error(
    validate_enzyme_rule(rule),
    "The `acceptor` must be a single structure"
  )

  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed(c("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-", "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"))

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  expect_error(
    validate_enzyme_rule(rule),
    "The `product` must be a single structure"
  )
})

test_that("validate_enzyme_rule fails with invalid acceptor-product pairs for GTs", {
  # Case 1: the acceptor is not a substructure of the product
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  expect_error(
    validate_enzyme_rule(rule),
    "`acceptor` must be a substructure of `product`"
  )

  # Case 2: the acceptor is the same as the product
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  expect_error(
    validate_enzyme_rule(rule),
    "`product` must have exactly one more residue than `acceptor`"
  )

  # Case 3: the acceptor is two residues shorter than the product
  acceptor <- glyparse::parse_iupac_condensed("GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  expect_error(
    validate_enzyme_rule(rule),
    "`product` must have exactly one more residue than `acceptor`"
  )
})

test_that("validate_enzyme_rule fails with invalid acceptor-product pairs for GDs", {
  # Case 1: the product is not a substructure of the acceptor
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Gal(b1-4)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GD")
  expect_error(
    validate_enzyme_rule(rule),
    "`product` must be a substructure of `acceptor`"
  )

  # Case 2: the product is the same as the acceptor
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GD")
  expect_error(
    validate_enzyme_rule(rule),
    "`acceptor` must have exactly one more residue than `product`"
  )

  # Case 3: the product is two residues shorter than the acceptor
  acceptor <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GD")
  expect_error(
    validate_enzyme_rule(rule),
    "`acceptor` must have exactly one more residue than `product`"
  )
})

# Test new_motif_set function
test_that("new_motif_set creates valid motif set objects", {
  motifs <- glyparse::parse_iupac_condensed(c("Gal(b1-3)GalNAc(a1-", "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"))
  alignments <- c("terminal", "substructure")

  motif_set <- new_motif_set(motifs, alignments)

  expect_s3_class(motif_set, "glyenzy_motif_set")
  expect_equal(motif_set$motifs, motifs)
  expect_equal(motif_set$alignments, alignments)
})

test_that("new_motif_set validates input parameters", {
  motifs <- glyparse::parse_iupac_condensed(c("Gal(b1-3)GalNAc(a1-", "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"))

  # Invalid motifs class
  expect_error(
    new_motif_set("invalid", c("terminal", "substructure")),
    "Assertion on 'motifs' failed"
  )

  # Invalid alignment
  expect_error(
    new_motif_set(motifs, c("invalid", "substructure")),
    "Assertion on 'alignments' failed"
  )
})

test_that("validate_motif_set fails for mismatched lengths", {
  motifs <- glyparse::parse_iupac_condensed(c("Gal(b1-3)GalNAc(a1-", "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"))
  alignments <- c("terminal")  # Only one alignment for two motifs

  motif_set <- new_motif_set(motifs, alignments)
  expect_error(
    validate_motif_set(motif_set),
    "The length of `motifs` must be equal to the length of `alignments`"
  )
})

# Test new_enzyme function
test_that("new_enzyme creates valid enzyme objects", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  motifs <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  rejects <- new_motif_set(motifs, "terminal")
  markers <- new_motif_set(motifs, "terminal")

  enzyme <- new_enzyme("ST3GAL2", list(rule), rejects, markers, "GT", "human")

  expect_s3_class(enzyme, "glyenzy_enzyme")
  expect_equal(enzyme$name, "ST3GAL2")
  expect_equal(enzyme$rules, list(rule))
  expect_equal(enzyme$rejects, rejects)
  expect_equal(enzyme$markers, markers)
  expect_equal(enzyme$type, "GT")
  expect_equal(enzyme$species, "human")
})

test_that("new_enzyme validates input parameters", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  motifs <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  rejects <- new_motif_set(motifs, "terminal")
  markers <- new_motif_set(motifs, "terminal")

  # Invalid name
  expect_error(
    new_enzyme(123, list(rule), rejects, markers, "GT", "human"),
    "Assertion on 'name' failed"
  )

  # Invalid rules
  expect_error(
    new_enzyme("ST3GAL2", "invalid", rejects, markers, "GT", "human"),
    "Assertion on 'rules' failed"
  )

  # Invalid rejects
  expect_error(
    new_enzyme("ST3GAL2", list(rule), "invalid", markers, "GT", "human"),
    "Assertion on 'rejects' failed"
  )

  # Invalid markers
  expect_error(
    new_enzyme("ST3GAL2", list(rule), rejects, "invalid", "GT", "human"),
    "Assertion on 'markers' failed"
  )

  # Invalid type
  expect_error(
    new_enzyme("ST3GAL2", list(rule), rejects, markers, "invalid", "human"),
    "Assertion on 'type' failed"
  )

  # Invalid species
  expect_error(
    new_enzyme("ST3GAL2", list(rule), rejects, markers, "GT", 123),
    "Assertion on 'species' failed"
  )
})

test_that("validate_enzyme fails for mismatched rule types", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  motifs <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  gt_rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  gd_rule <- new_enzyme_rule(product, acceptor, "terminal", "GD")
  rejects <- new_motif_set(motifs, "terminal")
  markers <- new_motif_set(motifs, "terminal")

  # Enzyme type is GT but contains GD rule
  enzyme <- new_enzyme("TEST", list(gt_rule, gd_rule), rejects, markers, "GT", "human")
  expect_error(
    validate_enzyme(enzyme),
    "All rules must have the same type as the enzyme"
  )
})

# Test print method for glyenzy_enzyme
test_that("print.glyenzy_enzyme works correctly", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  motifs <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  rejects <- new_motif_set(motifs, "terminal")
  markers <- new_motif_set(motifs, "terminal")
  enzyme <- new_enzyme("ST3GAL2", list(rule), rejects, markers, "GT", "human")

  expect_snapshot(print(enzyme))
})

test_that("print.glyenzy_enzyme handles empty motif sets", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  empty_motifs <- glyparse::parse_iupac_condensed(character(0))

  rule <- new_enzyme_rule(acceptor, product, "terminal", "GT")
  empty_rejects <- new_motif_set(empty_motifs, character(0))
  empty_markers <- new_motif_set(empty_motifs, character(0))
  enzyme <- new_enzyme("TEST", list(rule), empty_rejects, empty_markers, "GT", "human")

  expect_snapshot(print(enzyme))
})

test_that("print.glyenzy_enzyme handles multiple rules", {
  acceptor1 <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product1 <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  acceptor2 <- glyparse::parse_iupac_condensed("Gal(b1-4)GlcNAc(b1-")
  product2 <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-")

  rule1 <- new_enzyme_rule(acceptor1, product1, "terminal", "GT")
  rule2 <- new_enzyme_rule(acceptor2, product2, "terminal", "GT")
  empty_motifs <- glyparse::parse_iupac_condensed(character(0))
  empty_set <- new_motif_set(empty_motifs, character(0))

  enzyme <- new_enzyme("ST3GAL1", list(rule1, rule2), empty_set, empty_set, "GT", "human")

  expect_snapshot(print(enzyme))
})