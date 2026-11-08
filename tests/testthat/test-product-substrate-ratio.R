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
  S4Vectors::metadata(exp) <- list(
    glycan_type = "N",
    source = "test",
    scale = "raw"
  )
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
  expect_identical(
    S4Vectors::metadata(result),
    S4Vectors::metadata(exp)
  )
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
  S4Vectors::metadata(exp) <- list(
    glycan_type = "N",
    source = "test",
    scale = "raw"
  )
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
  expect_identical(
    S4Vectors::metadata(result),
    S4Vectors::metadata(exp)
  )
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
  S4Vectors::metadata(exp) <- list(glycan_type = "N", scale = "raw")
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
  S4Vectors::metadata(exp) <- list(glycan_type = "N", scale = "raw")
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

test_that("batched ratio counts preserve graph matching semantics", {
  glycans <- glyparse::auto_parse(c(
    "Gal(b1-4)GlcNAc(b1-",
    "Gal3S6S(b1-4)GlcNAc(b1-",
    "Gal?S(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-3)[Gal(b1-6)]GalNAc(a1-",
    "Gal(b1-3)[Gal(b1-6)]GalNAc-ol(a1-"
  ))
  motifs <- glyparse::auto_parse(c(
    "Gal(b1-",
    "Gal3S(b1-",
    "Gal?S(b1-",
    "Gal3S6S(b1-",
    "Gal(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-6)Gal(b1-",
    "Gal(?1-?)[Gal(?1-?)]GalNAc(a1-",
    "Hex(b1-"
  ))
  names(motifs) <- paste0("M", seq_along(motifs))
  cases <- list(
    glycans,
    glyrepr::remove_linkages(glycans),
    glyrepr::convert_to_generic(glycans),
    glyparse::auto_parse(c("Gal(?1-4)GlcNAc(b1-", "Gal(b1-4)GlcNAc(b1-"))
  )
  for (x in cases) {
    x <- x[c(seq_along(x), 1L, 1L)]
    names(x) <- paste0("G", seq_along(x))
    for (alignment in c("substructure", "terminal", "core", "whole")) {
      native <- .product_substrate_native_counts(x, motifs, alignment)
      expect_type(native, "integer")
      reference <- .count_motifs_substituent_subset(x, motifs, alignment)
      expect_identical(unname(native), unname(reference))
      expect_identical(
        .product_substrate_counts(x, motifs, alignment),
        reference
      )
    }
  }
})

test_that("ratio counting falls back for floating and additional metadata", {
  floating <- glyparse::auto_parse("{Fuc(a1-2)}Gal(b1-4)GlcNAc(b1-")
  extra <- glyparse::auto_parse("Gal(b1-4)GlcNAc(b1-")
  graphs <- glyrepr::get_structure_graphs(extra, return_list = TRUE)
  graphs[[1]] <- igraph::set_graph_attr(graphs[[1]], "source", "test")
  extra <- glyrepr::new_glycan_structure(
    as.character(extra),
    stats::setNames(graphs, as.character(extra))
  )
  motifs <- glyparse::auto_parse(c("Gal(b1-", "Fuc(a1-"))
  for (x in list(floating, extra)) {
    x <- x[c(1, 1)]
    expect_null(.product_substrate_native_counts(x, motifs, "terminal"))
    expect_identical(
      .product_substrate_counts(x, motifs, "terminal"),
      .count_motifs_substituent_subset(x, motifs, "terminal")
    )
  }
})

test_that("batched ratios preserve real-data site and enzyme ordering", {
  source <- glyexp::real_experiment
  indices <- c(1:12, 1:4)
  abundance <- SummarizedExperiment::assay(source)[indices, 1:3]
  rownames(abundance) <- paste0("G", seq_along(indices))
  row_data <- SummarizedExperiment::rowData(source)[indices, , drop = FALSE]
  rownames(row_data) <- rownames(abundance)
  exp <- glyexp::GlycoproteomicSE(
    abundance,
    rowData = row_data,
    colData = SummarizedExperiment::colData(source)[1:3, , drop = FALSE],
    metadata = S4Vectors::metadata(source)
  )
  expect_true(methods::validObject(exp))
  enzymes <- c("ST6GAL1", "FUT8", "ST6GAL1")
  actual <- product_substrate_ratio(exp, enzymes)
  local_mocked_bindings(
    .product_substrate_counts = .count_motifs_substituent_subset
  )
  expect_identical(actual, product_substrate_ratio(exp, enzymes))
})

ratio_scale_experiment <- function(type) {
  exp <- if (type == "glycomics") {
    glyexp::real_experiment2[seq_len(4), seq_len(6)]
  } else {
    glyexp::real_experiment[seq_len(4), seq_len(6)]
  }
  rownames(exp) <- paste0("G", seq_len(4))
  colnames(exp) <- paste0("S", seq_len(6))
  SummarizedExperiment::assay(exp) <- matrix(
    c(
      2,
      8,
      4,
      12,
      0,
      2,
      0,
      4,
      NA,
      2,
      4,
      NA,
      2,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      NA,
      0.25,
      0,
      0.5,
      0.125
    ),
    nrow = 4,
    dimnames = dimnames(exp)
  )
  SummarizedExperiment::rowData(exp)$glycan_structure <- glyparse::auto_parse(
    rep(c("Gal(b1-4)GlcNAc(b1-", "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"), 2)
  )
  if (type == "glycoproteomics") {
    SummarizedExperiment::rowData(exp)$protein <- c("P2", "P2", "P1", "P1")
    SummarizedExperiment::rowData(exp)$protein_site <- c(20L, 20L, 10L, 10L)
  }
  S4Vectors::metadata(exp) <- list(
    glycan_type = "N",
    source = "scale regression",
    custom = list(batch = 1L),
    scale = "raw"
  )
  exp
}

test_that("ratios use linear abundances across reversible scales", {
  transforms <- list(raw = identity, log = log, log2 = log2, log10 = log10)
  for (type in c("glycomics", "glycoproteomics")) {
    raw <- ratio_scale_experiment(type)
    expected <- product_substrate_ratio(raw, "ST6GAL1")
    expected_values <- if (type == "glycomics") {
      matrix(c(10 / 3, NA, 0.5, 0, NA, 1 / 6), nrow = 1)
    } else {
      rbind(c(4, NA, NA, 0, NA, 0), c(3, NA, 0, NA, NA, 0.25))
    }
    expect_equal(unname(SummarizedExperiment::assay(expected)), expected_values)

    for (scale in names(transforms)) {
      input <- raw
      S4Vectors::metadata(input)$scale <- scale
      SummarizedExperiment::assay(input) <- transforms[[scale]](
        SummarizedExperiment::assay(raw)
      )
      expect_true(methods::validObject(input))
      original_assay <- SummarizedExperiment::assay(input)
      original_metadata <- S4Vectors::metadata(input)

      result <- product_substrate_ratio(input, "ST6GAL1")

      expect_equal(result, expected)
      expect_identical(S4Vectors::metadata(result), S4Vectors::metadata(raw))
      expect_identical(
        SummarizedExperiment::colData(result),
        SummarizedExperiment::colData(input)
      )
      expect_identical(SummarizedExperiment::assay(input), original_assay)
      expect_identical(S4Vectors::metadata(input), original_metadata)
    }
  }
})

test_that("legacy ratios default only absent scale metadata to raw", {
  for (type in c("glycomics", "glycoproteomics")) {
    raw <- ratio_scale_experiment(type)
    legacy <- raw
    S4Vectors::metadata(legacy)$scale <- NULL
    original_assay <- SummarizedExperiment::assay(legacy)
    original_metadata <- S4Vectors::metadata(legacy)

    expect_equal(
      product_substrate_ratio(legacy, "ST6GAL1"),
      product_substrate_ratio(raw, "ST6GAL1")
    )
    expect_identical(SummarizedExperiment::assay(legacy), original_assay)
    expect_identical(S4Vectors::metadata(legacy), original_metadata)

    for (scale in list(
      NULL,
      NA_character_,
      "",
      "unknown",
      "ln",
      c("raw", "log2")
    )) {
      invalid <- raw
      S4Vectors::metadata(invalid)["scale"] <- list(scale)
      expect_error(product_substrate_ratio(invalid, "ST6GAL1"), "scale")
    }
    logx <- raw
    S4Vectors::metadata(logx)$scale <- "logx"
    expect_error(
      product_substrate_ratio(logx, "ST6GAL1"),
      "recoverable abundances.*logx"
    )
  }
})

test_that("ratios reject invalid or unrepresentable linear abundances", {
  for (type in c("glycomics", "glycoproteomics")) {
    raw <- ratio_scale_experiment(type)
    for (value in c(-1, Inf, -Inf)) {
      invalid <- raw
      # Deliberately model a corrupted historical raw object. Modern assay
      # replacement rejects negative raw values before glyenzy can see them.
      S4Vectors::metadata(invalid)$scale <- "log2"
      SummarizedExperiment::assay(invalid)[1, 1] <- value
      S4Vectors::metadata(invalid)$scale <- "raw"
      expect_error(
        product_substrate_ratio(invalid, "ST6GAL1"),
        "finite, non-negative abundances"
      )
    }
    for (scale in c("log", "log2", "log10")) {
      invalid <- raw
      S4Vectors::metadata(invalid)$scale <- scale
      for (value in c(Inf, 1e4)) {
        SummarizedExperiment::assay(invalid)[1, 1] <- value
        expect_error(
          product_substrate_ratio(invalid, "ST6GAL1"),
          "finite, non-negative abundances"
        )
      }
      SummarizedExperiment::assay(invalid)[1, 1] <- -1e4
      expect_error(product_substrate_ratio(invalid, "ST6GAL1"), "underflow")
    }
  }
})
