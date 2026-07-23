test_that("built-in sulfotransferases have stable order and rule counts", {
  st_symbols <- c(
    "CHST1",
    "CHST2",
    "CHST3",
    "CHST4",
    "CHST5",
    "CHST6",
    "CHST8",
    "CHST9",
    "CHST10",
    "GAL3ST2",
    "GAL3ST3",
    "GAL3ST4"
  )
  enzymes <- db_enzymes()
  expected_counts <- c(4L, 4L, 2L, 2L, 3L, 3L, 1L, 1L, 1L, 6L, 4L, 2L)

  expect_identical(tail(names(enzymes), length(st_symbols)), st_symbols)
  expect_identical(
    unname(vapply(enzymes[st_symbols], \(x) length(x$rules), integer(1))),
    expected_counts
  )
  expect_identical(
    unname(vapply(enzymes[st_symbols], \(x) x$type, character(1))),
    rep("ST", length(st_symbols))
  )
  expect_identical(sum(expected_counts), 33L)
})

test_that("each context-free sulfotransferase rule creates its exact product", {
  st_symbols <- c(
    "CHST1",
    "CHST2",
    "CHST3",
    "CHST4",
    "CHST5",
    "CHST6",
    "CHST8",
    "CHST9",
    "CHST10",
    "GAL3ST2",
    "GAL3ST3",
    "GAL3ST4"
  )
  enzymes <- db_enzymes()[st_symbols]

  for (symbol in names(enzymes)) {
    for (i in seq_along(enzymes[[symbol]]$rules)) {
      rule <- enzymes[[symbol]]$rules[[i]]
      if (length(rule$requires) > 0L) {
        next
      }
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

test_that("CHST8 and CHST9 require an O-GalNAc or N-glycan core", {
  o_glycan <- "GalNAc(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  o_product <- "GalNAc4S(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  n_glycan <- paste0(
    "GalNAc(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  n_product <- paste0(
    "GalNAc4S(b1-4)GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )

  for (symbol in c("CHST8", "CHST9")) {
    expect_identical(as.character(apply_enzyme(o_glycan, symbol)), o_product)
    expect_identical(as.character(apply_enzyme(n_glycan, symbol)), n_product)
    expect_length(apply_enzyme("GalNAc(b1-4)GlcNAc(b1-", symbol), 0L)
  }
})

test_that("CHST3 excludes its glycosaminoglycan reaction", {
  gag <- "GlcA(b1-3)GalNAc(b1-4)GlcA(b1-"

  expect_length(apply_enzyme(gag, "CHST3"), 0L)
})

test_that("find_enzyme identifies built-in sulfotransferases", {
  found <- find_enzyme("Gal3S(b1-3)GalNAc(a1-")

  expect_identical("GAL3ST4" %in% found, TRUE)
})

test_that("B3GNT7 retains its sulfate-dependent rule", {
  acceptor <- paste0(
    "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)",
    "GlcNAc6S(b1-"
  )
  product <- paste0(
    "GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)",
    "Gal(b1-4)GlcNAc6S(b1-"
  )

  expect_identical(as.character(apply_enzyme(acceptor, "B3GNT7")), product)
})
