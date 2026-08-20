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
  matched <- vapply(
    seq_len(nrow(cases)),
    function(i) {
      products <- as.character(actual[[i]])
      expected <- cases$expected_product[[i]]
      if (identical(expected, "")) {
        length(products) == 0L
      } else {
        expected %in% products
      }
    },
    logical(1)
  )
  case_names <- paste(cases$source_row, cases$gene_symbol, sep = ":")
  failures <- vapply(
    which(!matched),
    function(i) {
      sprintf(
        "%s expected <%s>, got <%s>",
        case_names[[i]],
        cases$expected_product[[i]],
        paste(as.character(actual[[i]]), collapse = ";")
      )
    },
    character(1)
  )

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
  missing_free <- names(typed)[
    !vapply(
      typed,
      function(types) "free" %in% types,
      logical(1)
    )
  ]
  expect_equal(missing_free, character())
})

test_that("generalized rules retain curated context boundaries", {
  cases <- data.frame(
    gene = c(
      "B3GALT5",
      "B3GNT2",
      "B3GNT3",
      "MGAT4A",
      "ST8SIA6",
      "ST8SIA6",
      "ST8SIA6"
    ),
    acceptor = c(
      "Gal(b1-4)GlcNAc(b1-3)GalNAc(a1-",
      "Gal(b1-3)GlcNAc(b1-3)[Neu5Ac(a2-6)]GalNAc(a1-",
      "Gal(b1-4)GlcNAc(b1-",
      paste0(
        "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)]",
        "[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
      ),
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-",
      "Neu5Ac(a2-6)Gal(b1-4)Glc(b1-",
      paste0(
        "Gal(b1-3)GalNAc(b1-4)[Neu5Ac(a2-3)]",
        "Gal(b1-4)Glc(b1-"
      )
    )
  )

  products <- Map(
    function(acceptor, gene) {
      suppressWarnings(apply_enzyme(acceptor, gene))
    },
    cases$acceptor,
    cases$gene
  )
  failures <- sprintf(
    "%s on <%s> produced <%s>",
    cases$gene,
    cases$acceptor,
    vapply(products, paste, collapse = ";", FUN.VALUE = character(1))
  )

  expect_length(failures[lengths(products) > 0L], 0L)
})
