test_that("built-in glycolipid enzymes have stable types and rule counts", {
  added_symbols <- c("UGCG", "UGT8", "B4GALT5", "B4GALT6", "ST3GAL5")
  changed_symbols <- c(
    "B3GALT4",
    "B4GALNT1",
    "ST6GALNAC3",
    "ST8SIA1",
    "ST8SIA5"
  )
  enzymes <- db_enzymes()[c(added_symbols, changed_symbols)]

  expect_identical(names(enzymes), c(added_symbols, changed_symbols))
  expect_identical(
    unname(vapply(enzymes, \(x) x$type, character(1))),
    rep("GT", length(enzymes))
  )
  expect_identical(
    unname(vapply(enzymes, \(x) x$species, character(1))),
    rep("human", length(enzymes))
  )
  expect_identical(
    unname(vapply(enzymes, \(x) length(x$rules), integer(1))),
    c(1L, 1L, 1L, 1L, 2L, 4L, 4L, 2L, 3L, 2L)
  )
})

test_that("UGCG and UGT8 are glycosphingolipid starter GTs", {
  expected_products <- c(UGCG = "Glc(b1-", UGT8 = "Gal(b1-")

  for (symbol in names(expected_products)) {
    enzyme_obj <- enzyme(symbol)
    rule <- enzyme_obj$rules[[1]]

    expect_s3_class(enzyme_obj, "glyenzy_starter_gt_enzyme")
    expect_length(rule$acceptor, 0L)
    expect_identical(as.character(rule$product), expected_products[[symbol]])
  }
})

test_that("each context-free glycolipid rule creates its exact product", {
  symbols <- c("B4GALT5", "B4GALT6", "ST3GAL5", "B3GALT4", "B4GALNT1")
  enzymes <- db_enzymes()[symbols]

  for (symbol in names(enzymes)) {
    for (i in seq_along(enzymes[[symbol]]$rules)) {
      rule <- enzymes[[symbol]]$rules[[i]]
      single_rule_enzyme <- enzymes[[symbol]]
      single_rule_enzyme$rules <- list(rule)

      expect_identical(
        as.character(apply_enzyme(rule$acceptor, single_rule_enzyme)),
        as.character(rule$product),
        info = paste(symbol, "rule", i)
      )
    }
  }
})

test_that("B4GALT5 and B4GALT6 synthesize the LacCer headgroup", {
  acceptor <- "Glc(b1-"
  product <- "Gal(b1-4)Glc(b1-"

  for (symbol in c("B4GALT5", "B4GALT6")) {
    expect_identical(as.character(apply_enzyme(acceptor, symbol)), product)
  }
})

test_that("ST3GAL5 supports the GM3 and GM4 headgroup routes", {
  acceptors <- c("Gal(b1-4)Glc(b1-", "Gal(b1-")
  products <- c(
    "Neu5Ac(a2-3)Gal(b1-4)Glc(b1-",
    "Neu5Ac(a2-3)Gal(b1-"
  )

  actual <- vapply(
    acceptors,
    \(acceptor) as.character(apply_enzyme(acceptor, "ST3GAL5")),
    character(1)
  )

  expect_identical(unname(actual), products)
})

test_that("B4GALNT1 and B3GALT4 cover a-, b-, and c-series gangliosides", {
  b4galnt1_acceptors <- c(
    "Gal(b1-4)Glc(b1-",
    "Neu5Ac(a2-3)Gal(b1-4)Glc(b1-",
    "Neu5Ac(a2-8)Neu5Ac(a2-3)Gal(b1-4)Glc(b1-",
    "Neu5Ac(a2-8)Neu5Ac(a2-8)Neu5Ac(a2-3)Gal(b1-4)Glc(b1-"
  )
  b3galt4_products <- c(
    "Gal(b1-3)GalNAc(b1-4)Gal(b1-4)Glc(b1-",
    "Gal(b1-3)GalNAc(b1-4)[Neu5Ac(a2-3)]Gal(b1-4)Glc(b1-",
    "Gal(b1-3)GalNAc(b1-4)[Neu5Ac(a2-8)Neu5Ac(a2-3)]Gal(b1-4)Glc(b1-",
    paste0(
      "Gal(b1-3)GalNAc(b1-4)",
      "[Neu5Ac(a2-8)Neu5Ac(a2-8)Neu5Ac(a2-3)]Gal(b1-4)Glc(b1-"
    )
  )

  b4galnt1_products <- vapply(
    b4galnt1_acceptors,
    \(acceptor) as.character(apply_enzyme(acceptor, "B4GALNT1")),
    character(1)
  )
  actual <- vapply(
    b4galnt1_products,
    \(acceptor) as.character(apply_enzyme(acceptor, "B3GALT4")),
    character(1)
  )

  expect_identical(
    unname(actual),
    as.character(glyparse::auto_parse(b3galt4_products))
  )
})

test_that("ST6GALNAC3 uses the ganglioside GalNAc linkage", {
  acceptor <- paste0(
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(b1-4)",
    "Gal(b1-4)Glc(b1-"
  )
  product <- paste0(
    "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-6)]GalNAc(b1-4)",
    "Gal(b1-4)Glc(b1-"
  )

  expect_identical(as.character(apply_enzyme(acceptor, "ST6GALNAC3")), product)
})

test_that("ST8SIA1 and ST8SIA5 recognize the ganglioside motif", {
  gm1b <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(b1-4)Gal(b1-4)Glc(b1-"
  gd1c <- paste0(
    "Neu5Ac(a2-8)Neu5Ac(a2-3)Gal(b1-3)GalNAc(b1-4)",
    "Gal(b1-4)Glc(b1-"
  )

  for (symbol in c("ST8SIA1", "ST8SIA5")) {
    expect_identical(as.character(apply_enzyme(gm1b, symbol)), gd1c)
  }
})

test_that("ST3GAL2 existing rule covers branched gangliosides", {
  gm1 <- "Gal(b1-3)GalNAc(b1-4)[Neu5Ac(a2-3)]Gal(b1-4)Glc(b1-"
  gd1a <- paste0(
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(b1-4)",
    "[Neu5Ac(a2-3)]Gal(b1-4)Glc(b1-"
  )

  expect_identical(as.character(apply_enzyme(gm1, "ST3GAL2")), gd1a)
})

test_that("trace_biosynthesis follows GlcCer and GalCer headgroup routes", {
  gm3 <- "Neu5Ac(a2-3)Gal(b1-4)Glc(b1-"
  gm3_path <- trace_biosynthesis(
    gm3,
    enzymes = c("B4GALT5", "ST3GAL5")
  )
  gm4 <- "Neu5Ac(a2-3)Gal(b1-"
  gm4_path <- trace_biosynthesis(gm4, enzymes = "ST3GAL5")

  expect_identical(igraph::E(gm3_path)$enzyme, c("B4GALT5", "ST3GAL5"))
  expect_identical(igraph::E(gm4_path)$enzyme, "ST3GAL5")
  expect_identical(
    igraph::V(gm3_path)$name[igraph::V(gm3_path)$target],
    gm3
  )
  expect_identical(
    igraph::V(gm4_path)$name[igraph::V(gm4_path)$target],
    gm4
  )
})
