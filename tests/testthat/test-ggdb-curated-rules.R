test_that("built-in enzymes reproduce curated human GGDB reactions", {
  cases <- utils::read.csv(
    test_path("fixtures", "ggdb-curated-reactions.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(cases), 125L)
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

test_that("B3GNT2 excludes curated N-glycan acceptors", {
  expect_identical(enzyme("B3GNT2")$glycan_type, c("O", "free"))

  acceptors <- c(
    paste0(
      "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
      "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
      "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ),
    paste0(
      "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)",
      "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
      "Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
    ),
    paste0(
      "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)",
      "[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]",
      "Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
    )
  )
  products <- lapply(
    acceptors,
    function(acceptor) suppressWarnings(apply_enzyme(acceptor, "B3GNT2"))
  )

  expect_identical(lengths(products), rep(0L, length(acceptors)))
})

test_that("MGAT branching enzymes remain N-glycan-specific", {
  genes <- c("MGAT1", "MGAT2", "MGAT4A", "MGAT4B")
  types <- lapply(genes, function(gene) enzyme(gene)$glycan_type)

  expect_true(all(vapply(types, identical, logical(1), "N")))
  expect_length(
    suppressWarnings(apply_enzyme("Man(a1-3)Man(a1-", "MGAT1")),
    0L
  )
})

test_that("MGAT2 requires the complete N-glycan core context", {
  truncated <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-"
  expect_length(suppressWarnings(apply_enzyme(truncated, "MGAT2")), 0L)

  acceptor <- paste0(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expected <- paste0(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )

  expect_identical(
    as.character(suppressWarnings(apply_enzyme(acceptor, "MGAT2"))),
    expected
  )
})

test_that("B3GALT2 supports curated low-activity acceptors", {
  expect_length(
    suppressWarnings(apply_enzyme("Gal(b1-4)Glc(b1-", "B3GALT2")),
    0L
  )
  expect_identical(
    as.character(suppressWarnings(apply_enzyme(
      "GlcNAc(b1-3)Gal(b1-4)Glc(b1-",
      "B3GALT2"
    ))),
    "Gal(b1-3)GlcNAc(b1-3)Gal(b1-4)Glc(b1-"
  )

  acceptors <- c("Gal(b1-4)GlcNAc(b1-", "Gal(b1-", "Glc(b1-")
  expected <- c(
    "Gal(b1-3)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-3)Gal(b1-",
    "Gal(b1-3)Glc(b1-"
  )
  products <- Map(
    function(acceptor) {
      as.character(suppressWarnings(apply_enzyme(acceptor, "B3GALT2")))
    },
    acceptors
  )

  expect_identical(unlist(products, use.names = FALSE), expected)
})

test_that("GAL3ST2 excludes Lewis a and Lewis x acceptors", {
  acceptors <- c(
    "Gal(b1-3)[Fuc(a1-4)]GlcNAc(b1-3)Gal(b1-4)Glc(b1-",
    "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-",
    "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-3)Gal(b1-4)Glc(b1-"
  )
  products <- lapply(
    acceptors,
    function(acceptor) suppressWarnings(apply_enzyme(acceptor, "GAL3ST2"))
  )

  expect_identical(lengths(products), rep(0L, length(acceptors)))
  expect_identical(
    as.character(suppressWarnings(apply_enzyme(
      "Gal(b1-4)GlcNAc(b1-3)GalNAc(a1-",
      "GAL3ST2"
    ))),
    "Gal3S(b1-4)GlcNAc(b1-3)GalNAc(a1-"
  )
})

test_that("B3GNT4 excludes long type-2 poly-LacNAc acceptors", {
  acceptors <- c(
    "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-",
    paste0(
      "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)",
      "Gal(b1-4)GlcNAc(b1-"
    ),
    paste0(
      "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)",
      "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)",
      "Gal(b1-4)GlcNAc(b1-"
    )
  )
  products <- lapply(
    acceptors,
    function(acceptor) suppressWarnings(apply_enzyme(acceptor, "B3GNT4"))
  )

  expect_identical(lengths(products), rep(0L, length(acceptors)))
  expect_identical(
    as.character(suppressWarnings(apply_enzyme(
      "Gal(b1-4)GlcNAc(b1-",
      "B3GNT4"
    ))),
    "GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-"
  )
})

test_that("B4GALNT2 requires alpha2-3-sialylated Gal", {
  expect_length(
    suppressWarnings(apply_enzyme("Gal(b1-3)GlcNAc(b1-", "B4GALNT2")),
    0L
  )
  expect_identical(
    as.character(suppressWarnings(apply_enzyme(
      "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
      "B4GALNT2"
    ))),
    "Neu5Ac(a2-3)[GalNAc(b1-4)]Gal(b1-3)GlcNAc(b1-"
  )
})

test_that("GCNT1 distinguishes core-1 from core-3", {
  expect_length(
    suppressWarnings(apply_enzyme("GlcNAc(b1-3)GalNAc(a1-", "GCNT1")),
    0L
  )
  expect_identical(
    as.character(suppressWarnings(apply_enzyme(
      "Gal(b1-3)GalNAc(a1-",
      "GCNT1"
    ))),
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  )
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
