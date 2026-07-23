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

test_that("apply_enzyme uses lenient matching for non-intact glycans", {
  glycan <- "GalNAc(?1-"

  expect_warning(
    expect_equal(
      as.character(apply_enzyme(glycan, "C1GALT1")),
      "Gal(b1-3)GalNAc(?1-"
    ),
    "non-intact glycan structures"
  )
})

test_that("apply_enzyme can return reduced-level structures", {
  topological_res <- suppressWarnings(apply_enzyme(
    "GalNAc(??-",
    "C1GALT1",
    structure_level = "topological"
  ))
  basic_res <- suppressWarnings(apply_enzyme(
    "HexNAc(??-",
    "C1GALT1",
    structure_level = "basic"
  ))

  expect_equal(glyrepr::get_structure_level(topological_res), "topological")
  expect_equal(as.character(topological_res), "Gal(??-?)GalNAc(??-")
  expect_equal(glyrepr::get_structure_level(basic_res), "basic")
  expect_equal(as.character(basic_res), "Hex(??-?)HexNAc(??-")
})

test_that("apply_enzyme rejects structure_level lower than input structures", {
  expect_error(
    apply_enzyme("GalNAc(a1-", "C1GALT1", structure_level = "topological"),
    "lower than"
  )
  expect_error(
    apply_enzyme("GalNAc(a1-", "C1GALT1", structure_level = "basic"),
    "lower than"
  )
  expect_error(
    suppressWarnings(
      apply_enzyme("GalNAc(?1-", "C1GALT1", structure_level = "topological")
    ),
    "lower than"
  )
  expect_error(
    suppressWarnings(
      apply_enzyme("GalNAc(??-", "C1GALT1", structure_level = "basic")
    ),
    "lower than"
  )
})

test_that("apply_enzyme validates structure_level", {
  expect_error(
    apply_enzyme("GalNAc(a1-", "C1GALT1", structure_level = "partial"),
    "structure_level"
  )
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

test_that("trusted graph construction deduplicates symmetric basic products", {
  glycan <- glyrepr::reduce_structure_level(
    glyrepr::as_glycan_structure(
      "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ),
    "basic"
  )

  products <- suppressWarnings(apply_enzyme(
    glycan,
    "B4GALT1",
    structure_level = "basic"
  ))

  expect_length(products, 1L)
  expect_silent(
    glyrepr::validate_glycan_graph(glyrepr::get_structure_graphs(products))
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

test_that("apply_enzyme does not reuse an occupied acceptor carbon", {
  enzyme <- make_enzyme(
    name = "TEST_C3",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "GalNAc(a1-",
      acceptor_alignment = "core",
      rejects = NULL,
      product = "GlcNAc(b1-3)GalNAc(a1-"
    ))
  )

  result <- apply_enzyme("Gal(b1-3)GalNAc(a1-", enzyme)

  expect_equal(result, glyrepr::glycan_structure())
})

test_that("apply_enzyme adds one sulfate while preserving existing sulfates", {
  enz <- make_enzyme(
    name = "TEST_ST",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "Gal6S(b1-4)GlcNAc(b1-",
      acceptor_alignment = "terminal",
      rejects = NULL,
      product = "Gal3S6S(b1-4)GlcNAc(b1-"
    ))
  )

  result <- apply_enzyme("Gal6S(b1-4)GlcNAc(b1-", enz)
  expect_equal(as.character(result), "Gal3S6S(b1-4)GlcNAc(b1-")
  expect_equal(apply_enzyme(result, enz), glyrepr::glycan_structure())
})

test_that("apply_enzyme returns one unique ST product per modifiable site", {
  rule <- list(
    acceptor = "Gal(b1-4)GlcNAc(b1-",
    acceptor_alignment = "terminal",
    rejects = NULL,
    product = "Gal6S(b1-4)GlcNAc(b1-"
  )
  enz <- make_enzyme(
    name = "TEST_ST",
    type = "ST",
    species = "human",
    rules = list(rule, rule)
  )
  glycan <- "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-"

  result <- apply_enzyme(glycan, enz)

  expect_length(result, 2L)
  expect_true(all(glyrepr::count_mono(result, "Gal") == 2L))
  expect_true(all(.has_substituents(result)))
})

test_that("ST requirements are OR-valued with sulfate-subset matching", {
  enz <- make_enzyme(
    name = "TEST_REQUIRED_ST",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-",
      acceptor_alignment = "substructure",
      rejects = NULL,
      requires = list(
        list(motif = "GalNAc(a1-", alignment = "core"),
        list(
          motif = "GlcNAc(b1-4)GlcNAc(b1-",
          alignment = "core"
        )
      ),
      product = "Gal6S(b1-"
    ))
  )
  glycans <- c(
    "Gal(b1-3)GalNAc6S(a1-",
    "Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-3)GlcNAc(b1-4)GlcNAc(b1-"
  )

  result <- apply_enzyme(glycans, enz)

  expect_equal(
    as.character(result[[1]]),
    "Gal6S(b1-3)GalNAc6S(a1-"
  )
  expect_equal(result[[2]], glyrepr::glycan_structure())
  expect_equal(
    as.character(result[[3]]),
    "Gal6S(b1-3)GlcNAc(b1-4)GlcNAc(b1-"
  )
})

test_that("ST requirements use subset matching while acceptors and rejects stay exact", {
  enz <- make_enzyme(
    name = "TEST_EXACT_ST",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "Gal6S(b1-",
      acceptor_alignment = "substructure",
      rejects = "GlcNAc(b1-4)Gal6S(b1-",
      requires = list(list(
        motif = "GlcNAc(b1-4)Gal6S(b1-",
        alignment = "substructure"
      )),
      product = "Gal3S6S(b1-"
    ))
  )

  expect_equal(
    as.character(apply_enzyme("GlcNAc6S(b1-4)Gal6S(b1-", enz)),
    "GlcNAc6S(b1-4)Gal3S6S(b1-"
  )
  expect_equal(
    apply_enzyme("GlcNAc6S(b1-4)Gal4S6S(b1-", enz),
    glyrepr::glycan_structure()
  )
})

test_that("ST action does not sulfate an occupied carbon", {
  enz <- make_enzyme(
    name = "TEST_ST",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-",
      acceptor_alignment = "substructure",
      rejects = NULL,
      product = "Gal6S(b1-"
    ))
  )

  result <- apply_enzyme("Neu5Ac(a2-6)Gal(b1-", enz)

  expect_equal(result, glyrepr::glycan_structure())
})

test_that("existing GT rules can act on sulfated substrates", {
  result <- apply_enzyme("GlcNAc6S(b1-", "B4GALT4")

  expect_equal(as.character(result), "Gal(b1-4)GlcNAc6S(b1-")
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
  expect_equal(
    as.character(res[["FUT3"]]),
    "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-"
  )
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
  expect_equal(
    as.character(res[["FUT3"]]),
    "Fuc(a1-2)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-"
  )
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
  expect_equal(
    as.character(res[["FUT3"]]),
    "Fuc(a1-2)[Gal(a1-3)]Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-"
  )
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
  expect_equal(
    as.character(res[["FUT3"]]),
    "Neu5Ac(a2-3)Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-"
  )
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
  expect_equal(
    as.character(res[["FUT3"]]),
    "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT4"]]),
    "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT5"]]),
    "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT6"]]),
    "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(as.character(res[["FUT7"]]), character(0))
  expect_equal(
    as.character(res[["FUT9"]]),
    "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
})

test_that("fucosylation of A antigen (Type 2)", {
  enzymes <- c("FUT3", "FUT4", "FUT5", "FUT6", "FUT7", "FUT9")
  glycan <- "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)GlcNAc(b1-"
  res <- purrr::map(enzymes, ~ apply_enzyme(glycan, .x))
  names(res) <- enzymes
  expect_equal(
    as.character(res[["FUT3"]]),
    "Fuc(a1-2)[GalNAc(a1-3)]Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
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
  expect_equal(
    as.character(res[["FUT3"]]),
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT4"]]),
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT5"]]),
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT6"]]),
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
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
  expect_equal(
    as.character(res[["FUT5"]]),
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT6"]]),
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(
    as.character(res[["FUT7"]]),
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
  )
  expect_equal(as.character(res[["FUT9"]]), character(0))
})

# ===== N-Glycan Trimming =====
test_that("Trimming on Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1B1_result <- apply_enzyme(glycan, "MAN1B1")
  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  MAN1B1_expected <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1B1_result), MAN1B1_expected)

  MAN1A1_expected <- c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(sort(as.character(MAN1A1_result)), sort(MAN1A1_expected))
  expect_equal(sort(as.character(MAN1A2_result)), sort(MAN1A1_expected))
  expect_equal(sort(as.character(MAN1C1_result)), sort(MAN1A1_expected))
})

test_that("Trimming on Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  MAN1A1_expected <- "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1A1_result), MAN1A1_expected)
  expect_equal(as.character(MAN1A2_result), MAN1A1_expected)

  MAN1C1_expected <- c(
    "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(sort(as.character(MAN1C1_result)), sort(MAN1C1_expected))
})

test_that("Trimming on Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1A1_result), expected)
  expect_equal(as.character(MAN1A2_result), expected)
  expect_equal(as.character(MAN1C1_result), expected)
})

test_that("Trimming on Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(sort(as.character(MAN1A1_result)), sort(expected))
  expect_equal(sort(as.character(MAN1A2_result)), sort(expected))
  expect_equal(sort(as.character(MAN1C1_result)), sort(expected))
})

test_that("Trimming on Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expect_equal(MAN1A1_result, glyrepr::glycan_structure())
  expect_equal(MAN1A2_result, glyrepr::glycan_structure())

  MAN1C1_expected <- "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1C1_result), MAN1C1_expected)
})

test_that("Trimming on Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1A1_result), expected)
  expect_equal(as.character(MAN1A2_result), expected)
  expect_equal(as.character(MAN1C1_result), expected)
})

test_that("Trimming on Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1A1_result), expected)
  expect_equal(as.character(MAN1A2_result), expected)
  expect_equal(as.character(MAN1C1_result), expected)
})

test_that("Trimming on Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)[Man(a1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- c(
    "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(sort(as.character(MAN1A1_result)), sort(expected))
  expect_equal(sort(as.character(MAN1A2_result)), sort(expected))
  expect_equal(sort(as.character(MAN1C1_result)), sort(expected))
})

test_that("Trimming on Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1A1_result), expected)
  expect_equal(as.character(MAN1A2_result), expected)
  expect_equal(as.character(MAN1C1_result), expected)
})

test_that("Trimming on Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1A1_result), expected)
  expect_equal(as.character(MAN1A2_result), expected)
  expect_equal(as.character(MAN1C1_result), expected)
})

test_that("Trimming on Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", {
  glycan <- "Man(a1-2)Man(a1-6)[Man(a1-3)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  MAN1A1_result <- apply_enzyme(glycan, "MAN1A1")
  MAN1A2_result <- apply_enzyme(glycan, "MAN1A2")
  MAN1C1_result <- apply_enzyme(glycan, "MAN1C1")

  expected <- "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_equal(as.character(MAN1A1_result), expected)
  expect_equal(as.character(MAN1A2_result), expected)
  expect_equal(as.character(MAN1C1_result), expected)
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
