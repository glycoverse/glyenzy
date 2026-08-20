test_that("built-in enzymes reproduce curated human GGDB reactions", {
  cases <- utils::read.csv(
    test_path("fixtures", "ggdb-curated-reactions.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(cases), 119L)
  expect_equal(sum(startsWith(cases$gene_symbol, "ALG")), 0L)

  actual <- vapply(seq_len(nrow(cases)), function(i) {
    products <- suppressWarnings(apply_enzyme(
      cases$acceptor_structure[[i]],
      cases$gene_symbol[[i]]
    ))
    paste(sort(unique(as.character(products))), collapse = ";")
  }, character(1))
  expected <- cases$expected_product
  case_names <- paste(cases$source_row, cases$gene_symbol, sep = ":")
  names(actual) <- case_names
  names(expected) <- case_names

  expect_equal(actual, expected)
})

test_that("ambiguous lipid/free acceptors are enabled explicitly", {
  cases <- utils::read.csv(
    test_path("fixtures", "ggdb-curated-reactions.csv"),
    stringsAsFactors = FALSE
  )
  genes <- unique(cases$gene_symbol[cases$requires_free])

  metadata <- lapply(genes, function(gene) enzyme(gene)$glycan_type)
  names(metadata) <- genes
  typed <- metadata[!vapply(metadata, is.null, logical(1))]

  expect_true(length(typed) > 0L)
  expect_true(all(vapply(typed, function(types) "free" %in% types, logical(1))))
})
