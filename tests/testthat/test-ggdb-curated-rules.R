test_that("built-in enzymes reproduce curated human GGDB reactions", {
  cases <- utils::read.csv(
    test_path("fixtures", "ggdb-curated-reactions.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(cases), 119L)
  expect_equal(sum(startsWith(cases$gene_symbol, "ALG")), 0L)

  actual <- lapply(seq_len(nrow(cases)), function(i) {
    suppressWarnings(apply_enzyme(
      cases$acceptor_structure[[i]],
      cases$gene_symbol[[i]]
    ))
  })
  matched <- vapply(seq_len(nrow(cases)), function(i) {
    products <- as.character(actual[[i]])
    expected <- cases$expected_product[[i]]
    if (identical(expected, "")) {
      length(products) == 0L
    } else {
      expected %in% products
    }
  }, logical(1))
  case_names <- paste(cases$source_row, cases$gene_symbol, sep = ":")
  failures <- vapply(which(!matched), function(i) {
    sprintf(
      "%s expected <%s>, got <%s>",
      case_names[[i]],
      cases$expected_product[[i]],
      paste(as.character(actual[[i]]), collapse = ";")
    )
  }, character(1))

  expect_length(failures, 0L)
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

  expect_gt(length(typed), 0L)
  missing_free <- names(typed)[!vapply(
    typed,
    function(types) "free" %in% types,
    logical(1)
  )]
  expect_equal(missing_free, character())
})
