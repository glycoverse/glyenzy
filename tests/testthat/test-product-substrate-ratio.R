test_that("product_substrate_ratio calculates glycomics rule sums", {
  structures <- glyparse::auto_parse(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
    "Fuc(a1-2)Gal(b1-4)GlcNAc(b1-"
  ))
  abundance <- matrix(
    c(2, 4, 6, 8, 4, 2),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("G1", "G2", "G3"), c("S1", "S2"))
  )
  col_data <- S4Vectors::DataFrame(
    group = c("control", "case"),
    row.names = colnames(abundance)
  )
  exp <- glyexp::real_experiment2[seq_len(3), seq_len(2)]
  rownames(exp) <- rownames(abundance)
  colnames(exp) <- colnames(abundance)
  SummarizedExperiment::assay(exp) <- abundance
  SummarizedExperiment::rowData(exp)$glycan_structure <- structures
  SummarizedExperiment::colData(exp) <- col_data
  S4Vectors::metadata(exp) <- list(glycan_type = "N", source = "test")
  enzyme <- make_enzyme(
    name = "TEST_GT",
    rules = list(
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = NULL,
        product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      ),
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = NULL,
        product = "Fuc(a1-2)Gal(b1-4)GlcNAc(b1-"
      ),
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = NULL,
        product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      )
    ),
    type = "GT",
    species = "human"
  )

  result <- product_substrate_ratio(exp, list(enzyme, enzyme))

  expect_s4_class(result, "SummarizedExperiment")
  expect_identical(methods::is(result, "GlycomicSE"), FALSE)
  expect_identical(
    SummarizedExperiment::assayNames(result),
    "product_substrate_ratio"
  )
  expect_equal(
    SummarizedExperiment::assay(result),
    matrix(
      rep(c(8 / 3, 1.5), each = 2),
      nrow = 2,
      dimnames = list(c("TEST_GT", "TEST_GT.1"), c("S1", "S2"))
    )
  )
  expect_identical(
    SummarizedExperiment::rowData(result)$enzyme,
    c("TEST_GT", "TEST_GT")
  )
  expect_identical(SummarizedExperiment::colData(result), col_data)
  expect_identical(S4Vectors::metadata(result), S4Vectors::metadata(exp))
})

test_that("product_substrate_ratio uses lenient matching for non-intact data", {
  exp <- glyexp::real_experiment2
  structure_levels <- unique(
    glyrepr::get_structure_level(
      SummarizedExperiment::rowData(exp)$glycan_structure
    )
  )
  expect_true(all(structure_levels %in% c("partial", "topological")))

  by_name <- product_substrate_ratio(exp, c("MGAT3", "FUT8"))
  by_object <- product_substrate_ratio(
    exp,
    list(enzyme("MGAT3"), enzyme("FUT8"))
  )

  expect_identical(by_name, by_object)
  expect_identical(
    rowSums(is.na(SummarizedExperiment::assay(by_name))),
    c(MGAT3 = 0, FUT8 = 0)
  )
  expect_gt(min(SummarizedExperiment::assay(by_name)), 0)
})

test_that("product_substrate_ratio is site specific and keeps site order", {
  structures <- glyparse::auto_parse(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
  ))
  abundance <- matrix(
    c(2, 4, 6, 12, 4, 2, 8, 6),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(paste0("G", 1:4), c("S1", "S2"))
  )
  exp <- glyexp::real_experiment[seq_len(4), seq_len(2)]
  rownames(exp) <- rownames(abundance)
  colnames(exp) <- colnames(abundance)
  SummarizedExperiment::assay(exp) <- abundance
  SummarizedExperiment::rowData(exp)$protein <- c("P2", "P2", "P1", "P1")
  SummarizedExperiment::rowData(exp)$protein_site <- c(
    20L,
    20L,
    NA_integer_,
    NA_integer_
  )
  SummarizedExperiment::rowData(exp)$glycan_structure <- structures
  SummarizedExperiment::colData(exp) <- S4Vectors::DataFrame(
    batch = c("A", "B"),
    row.names = colnames(abundance)
  )
  S4Vectors::metadata(exp) <- list(glycan_type = "N", source = "test")
  enzyme <- make_enzyme(
    name = "TEST_GT",
    rules = list(list(
      acceptor = "Gal(b1-4)GlcNAc(b1-",
      acceptor_alignment = "terminal",
      rejects = NULL,
      product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
    )),
    type = "GT",
    species = "human"
  )

  result <- product_substrate_ratio(exp, list(enzyme))

  expect_equal(
    SummarizedExperiment::assay(result),
    matrix(
      c(3, 3, 2, 3),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("P2-20-TEST_GT", "P1-NA-TEST_GT"), c("S1", "S2"))
    )
  )
  result_rows <- SummarizedExperiment::rowData(result)
  expect_identical(
    colnames(result_rows),
    c("enzyme", "protein", "protein_site")
  )
  expect_identical(result_rows$protein, c("P2", "P1"))
  expect_identical(result_rows$protein_site, c(20L, NA_integer_))
  expect_identical(S4Vectors::metadata(result), S4Vectors::metadata(exp))
})

test_that("product_substrate_ratio ignores rule context and handles zero", {
  structures <- glyparse::auto_parse(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
  ))
  abundance <- matrix(
    c(2, 6),
    ncol = 1,
    dimnames = list(c("G1", "G2"), "S1")
  )
  exp <- glyexp::real_experiment2[seq_len(2), 1, drop = FALSE]
  rownames(exp) <- rownames(abundance)
  colnames(exp) <- colnames(abundance)
  SummarizedExperiment::assay(exp) <- abundance
  SummarizedExperiment::rowData(exp)$glycan_structure <- structures
  S4Vectors::metadata(exp) <- list(glycan_type = "N")
  plain <- make_enzyme(
    name = "PLAIN",
    rules = list(list(
      acceptor = "Gal(b1-4)GlcNAc(b1-",
      acceptor_alignment = "terminal",
      rejects = NULL,
      product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
    )),
    type = "GT",
    species = "human"
  )
  contextual <- make_enzyme(
    name = "CONTEXTUAL",
    rules = list(list(
      acceptor = "Gal(b1-4)GlcNAc(b1-",
      acceptor_alignment = "terminal",
      rejects = "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-",
      requires = list(list(motif = "GalNAc(a1-", alignment = "core")),
      product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
    )),
    type = "GT",
    species = "human"
  )

  result <- product_substrate_ratio(exp, list(plain, contextual))
  zero_result <- product_substrate_ratio(exp[2, ], list(plain))

  expect_equal(
    unname(SummarizedExperiment::assay(result)[, 1]),
    c(3, 3)
  )
  expect_identical(
    unname(SummarizedExperiment::assay(zero_result)[1, 1]),
    NA_real_
  )
})

test_that("product_substrate_ratio uses substituent-subset motif matching", {
  structures <- glyparse::auto_parse(c(
    "Gal3S6S(b1-",
    "Gal3S(b1-"
  ))
  exp <- glyexp::real_experiment2[seq_len(2), 1, drop = FALSE]
  rownames(exp) <- c("G1", "G2")
  colnames(exp) <- "S1"
  SummarizedExperiment::assay(exp) <- matrix(
    c(6, 2),
    ncol = 1,
    dimnames = list(c("G1", "G2"), "S1")
  )
  SummarizedExperiment::rowData(exp)$glycan_structure <- structures
  sulfate_transferase <- make_enzyme(
    name = "TEST_ST6",
    rules = list(list(
      acceptor = "Gal(b1-",
      acceptor_alignment = "whole",
      rejects = NULL,
      product = "Gal6S(b1-"
    )),
    type = "ST",
    species = "human"
  )

  result <- product_substrate_ratio(exp, list(sulfate_transferase))

  expect_equal(
    unname(SummarizedExperiment::assay(result)[1, 1]),
    6 / 8
  )
})

test_that("product_substrate_ratio validates its inputs", {
  structures <- glyparse::auto_parse("Gal(b1-4)GlcNAc(b1-")
  exp <- glyexp::real_experiment2[1, 1, drop = FALSE]
  rownames(exp) <- "G1"
  colnames(exp) <- "S1"
  SummarizedExperiment::assay(exp) <- matrix(
    1,
    nrow = 1,
    dimnames = list("G1", "S1")
  )
  SummarizedExperiment::rowData(exp)$glycan_structure <- structures
  S4Vectors::metadata(exp) <- list(glycan_type = "N")
  missing_structure <- exp
  SummarizedExperiment::rowData(missing_structure)$glycan_structure <- NULL
  empty_enzyme <- new_enzyme("EMPTY", list(), "GT", "human")

  expect_snapshot(error = TRUE, {
    product_substrate_ratio(
      SummarizedExperiment::SummarizedExperiment(
        assays = list(abundance = matrix(1))
      ),
      "ST6GAL1"
    )
  })
  expect_snapshot(error = TRUE, {
    product_substrate_ratio(missing_structure, "ST6GAL1")
  })
  expect_snapshot(error = TRUE, {
    product_substrate_ratio(exp, "UNKNOWN")
  })
  expect_snapshot(error = TRUE, {
    product_substrate_ratio(exp, c("ST6GAL1", "MAN2A1"))
  })
  expect_snapshot(error = TRUE, {
    product_substrate_ratio(
      exp,
      list(enzyme("ST6GAL1"), enzyme("MAN2A1"))
    )
  })
  expect_snapshot(error = TRUE, {
    product_substrate_ratio(exp, list(empty_enzyme))
  })
  expect_snapshot(error = TRUE, {
    product_substrate_ratio(exp, "GALNT1")
  })
})
