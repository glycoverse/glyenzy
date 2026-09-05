.normalization_fixture <- function(states, core_fuc = FALSE, bisect = FALSE) {
  antenna <- c(
    "GlcNAc(??-?)",
    "Gal(??-?)GlcNAc(??-?)",
    "Neu5Ac(??-?)Gal(??-?)GlcNAc(??-?)"
  )[states + 1L]
  left <- paste0(antenna[1], "[", antenna[2], "]Man(??-?)")
  right <- if (length(states) == 3L) {
    paste0(antenna[3], "Man(??-?)")
  } else {
    paste0(antenna[3], "[", antenna[4], "]Man(??-?)")
  }
  glyparse::auto_parse(paste0(
    left,
    "[",
    right,
    "]",
    if (bisect) "[GlcNAc(??-?)]",
    "Man(??-?)GlcNAc(??-?)",
    if (core_fuc) "[Fuc(??-?)]",
    "GlcNAc(??-"
  ))
}

test_that("all simple antenna distributions follow the documented convention", {
  representatives <- list(
    `3/1/0` = c(0, 0, 1),
    `3/1/1` = c(0, 0, 2),
    `3/2/0` = c(1, 0, 1),
    `3/2/1` = c(1, 0, 2),
    `3/2/2` = c(2, 0, 2),
    `3/3/1` = c(1, 1, 2),
    `3/3/2` = c(2, 1, 2),
    `4/2/0` = c(1, 0, 1, 0),
    `4/2/1` = c(2, 0, 1, 0),
    `4/2/2` = c(2, 0, 2, 0),
    `4/3/1` = c(1, 1, 2, 0),
    `4/3/2` = c(2, 1, 2, 0),
    `4/4/2` = c(2, 1, 2, 1)
  )
  for (n in 3:4) {
    states <- expand.grid(rep(list(0:2), n))
    for (i in seq_len(nrow(states))) {
      state <- as.integer(states[i, ])
      x <- .normalization_fixture(state)
      key <- paste(n, sum(state > 0), sum(state == 2), sep = "/")
      target <- representatives[[key]]
      expected <- if (is.null(target)) x else .normalization_fixture(target)
      result <- normalize_n_glycan(x)
      expect_identical(
        as.character(result),
        as.character(expected),
        info = paste(state, collapse = ",")
      )
      expect_identical(
        as.character(normalize_n_glycan(result)),
        as.character(result)
      )
      counts <- function(z) {
        sort(igraph::vertex_attr(glyrepr::get_structure_graphs(z), "mono"))
      }
      expect_identical(counts(result), counts(x))
      expect_identical(glyrepr::get_structure_level(result), "topological")
      graph <- glyrepr::get_structure_graphs(x)
      permuted <- igraph::permute(graph, rev(seq_len(igraph::vcount(graph))))
      expect_identical(
        as.character(glyrepr::as_glycan_structure(.normalize_n_glycan_graph(
          permuted
        ))),
        as.character(result)
      )
    }
  }
})

test_that("normalization retains core fucose and bisecting GlcNAc", {
  for (fuc in c(FALSE, TRUE)) {
    for (bisect in c(FALSE, TRUE)) {
      for (state in list(c(2, 1, 1), c(2, 2, 1, 1))) {
        expected <- if (length(state) == 3L) c(1, 1, 2) else c(2, 1, 2, 1)
        x <- .normalization_fixture(state, fuc, bisect)
        expect_identical(
          as.character(normalize_n_glycan(x)),
          as.character(.normalization_fixture(expected, fuc, bisect))
        )
      }
    }
  }
})

test_that("unsupported topological structures are unchanged", {
  base <- as.character(.normalization_fixture(c(2, 1, 0)))
  unusual <- c(
    sub("Neu5Ac", "Neu5Gc", base, fixed = TRUE),
    sub("Gal(??-?)", "Gal3S(??-?)", base, fixed = TRUE),
    sub("Gal(??-?)", "Gal(??-?)[Fuc(??-?)]", base, fixed = TRUE),
    sub("Gal(??-?)", "Gal(??-?)GlcNAc(??-?)Gal(??-?)", base, fixed = TRUE),
    paste0("{Fuc(??-?)}", base),
    paste0("{?S}", base),
    "Gal(??-?)GalNAc(??-",
    "Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-",
    "GlcNAc(??-?)Man(??-?)[GlcNAc(??-?)Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-"
  )
  for (value in unusual) {
    x <- glyparse::auto_parse(value)
    expect_identical(normalize_n_glycan(x), x)
  }
})

test_that("normalization preserves vector shape, names and missing values", {
  a <- .normalization_fixture(c(2, 0, 1))
  b <- .normalization_fixture(c(2, 2, 0, 0))
  x <- c(as.character(a), NA_character_, as.character(b), as.character(a))
  attr(x, "names") <- c("first", "", NA_character_, "first")
  result <- normalize_n_glycan(x)
  expect_s3_class(result, "glyrepr_structure")
  expect_length(result, 4L)
  expect_identical(names(result), names(x))
  expect_identical(is.na(result), is.na(glyparse::auto_parse(x)))
  expect_identical(
    unname(as.character(result)),
    c(
      as.character(normalize_n_glycan(a)),
      NA_character_,
      as.character(normalize_n_glycan(b)),
      as.character(normalize_n_glycan(a))
    )
  )
  expect_identical(
    as.character(normalize_n_glycan(result)),
    as.character(result)
  )
  expect_s3_class(normalize_n_glycan(character()), "glyrepr_structure")
  expect_length(normalize_n_glycan(character()), 0L)
  expect_length(normalize_n_glycan(glyrepr::glycan_structure()), 0L)
  expect_identical(
    as.character(normalize_n_glycan(c(NA_character_, NA_character_))),
    c(NA_character_, NA_character_)
  )
})

test_that("normalization rejects invalid input types", {
  expect_snapshot(error = TRUE, normalize_n_glycan(1))
})

test_that("normalization rejects generic monosaccharides", {
  expect_snapshot(error = TRUE, normalize_n_glycan("Hex(??-?)HexNAc(??-"))
})

test_that("normalization rejects intact structures", {
  expect_snapshot(error = TRUE, normalize_n_glycan("Gal(b1-4)GlcNAc(b1-"))
})

test_that("normalization rejects partial structures", {
  expect_snapshot(error = TRUE, normalize_n_glycan("Gal(??-?)GlcNAc(b1-"))
})
