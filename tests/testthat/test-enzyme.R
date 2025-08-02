test_that("create_enzyme works", {
  acceptor <- glyparse::parse_iupac_condensed("Gal(b1-3)GalNAc(a1-")
  product <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  enzyme <- create_enzyme("ST3GAL2", acceptor, product, "GT", "human")
  expect_snapshot(print(enzyme))
})

test_that("create_enzyme works with structure strings", {
  acceptor <- "Gal(b1-3)GalNAc(a1-"
  product <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  enzyme <- create_enzyme("ST3GAL2", acceptor, product, "GT", "human")
  expect_snapshot(print(enzyme))
})

test_that("create_enzyme fails for invalid structure strings", {
  acceptor <- "invalid"
  product <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  expect_error(create_enzyme("ST3GAL2", acceptor, product, "GT", "human"))
})

test_that("create_enzyme fails for multiple structures", {
  # Case 1: `acceptor` has multiple structures
  acceptor <- c("Gal(b1-3)GalNAc(a1-", "Gal(b1-3)GalNAc(a1-")
  product <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  expect_error(
    create_enzyme("ST3GAL2", acceptor, product, "GT", "human"),
    "`acceptor` must be a single structure"
  )

  # Case 2: `product` has multiple structures
  acceptor <- "Gal(b1-3)GalNAc(a1-"
  product <- c("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-", "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-")
  expect_error(
    create_enzyme("ST3GAL2", acceptor, product, "GT", "human"),
    "`product` must be a single structure"
  )
})

test_that("create_enzyme fails with invalid acceptor-product pairs for GTs", {
  # Case 1: the acceptor is not a substructure of the product
  acceptor <- "Gal(b1-3)GalNAc(a1-"
  product <- "Neu5Ac(a2-3)Gal(b1-4)GalNAc(a1-"
  expect_error(
    create_enzyme("ST3GAL2", acceptor, product, "GT", "human"),
    "`acceptor` must be a substructure of `product`"
  )

  # Case 2: the acceptor is the same of the product
  acceptor <- "Gal(b1-3)GalNAc(a1-"
  product <- "Gal(b1-3)GalNAc(a1-"
  expect_error(
    create_enzyme("ST3GAL2", acceptor, product, "GT", "human"),
    "`product` must have exactly one more residue than `acceptor`"
  )

  # Case 3: the acceptor is two residues shorter than the product
  acceptor <- "GalNAc(a1-"
  product <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  expect_error(
    create_enzyme("ST3GAL2", acceptor, product, "GT", "human"),
    "`product` must have exactly one more residue than `acceptor`"
  )

  # Case 4: the extra residue not at the terminal
  acceptor <- "Gal(b1-3)GalNAc(a1-"
  product <- "Gal(b1-3)GalNAc(a1-4)Neu5Ac(a2-"
  expect_error(
    create_enzyme("ST3GAL2", acceptor, product, "GT", "human"),
    "The extra residue in `product` must be at the terminal"
  )
})

test_that("create_enzyme fails with invalid acceptor-product pairs for GDs", {
  # Case 1: the product is not a substructure of the acceptor
  acceptor <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  product <- "Gal(b1-4)GalNAc(a1-"
  expect_error(
    create_enzyme("NEU1", acceptor, product, "GD", "human"),
    "`product` must be a substructure of `acceptor`"
  )

  # Case 2: the product is the same of the acceptor
  acceptor <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  product <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  expect_error(
    create_enzyme("NEU1", acceptor, product, "GD", "human"),
    "`acceptor` must have exactly one more residue than `product`"
  )

  # Case 3: the product is two residues longer than the acceptor
  acceptor <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  product <- "GalNAc(a1-"
  expect_error(
    create_enzyme("NEU1", acceptor, product, "GD", "human"),
    "`acceptor` must have exactly one more residue than `product`"
  )

  # Case 4: the extra residue not at the terminal
  acceptor <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  product <- "Neu5Ac(a2-3)Gal(b1-"
  expect_error(
    create_enzyme("NEU1", acceptor, product, "GD", "human"),
    "The extra residue in `acceptor` must be at the terminal"
  )
})