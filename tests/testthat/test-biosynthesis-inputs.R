test_that("trace_biosynthesis rejects incompatible glycans", {
  expect_snapshot(
    trace_biosynthesis("Hex(b1-3)GalNAc(a1-"),
    error = TRUE
  )
  expect_snapshot(
    trace_biosynthesis(c(
      "Gal(b1-3)GalNAc(a1-",
      "Hex(b1-3)HexNAc(a1-"
    )),
    error = TRUE
  )
  expect_snapshot(
    trace_biosynthesis("Gal(b1-?)GalNAc(a1-"),
    error = TRUE
  )
  expect_snapshot(
    trace_biosynthesis(NA_character_),
    error = TRUE
  )
})

test_that("trace_biosynthesis_virtual rejects mixed structure levels", {
  expect_snapshot(
    trace_biosynthesis_virtual(c(
      "Gal(b1-3)GalNAc(a1-",
      "Gal(??-?)GalNAc(??-"
    )),
    error = TRUE
  )
})

test_that("trace_biosynthesis accepts homogeneous generic intact glycans", {
  glycans <- c(
    "HexNAc(a1-",
    "Hex(b1-3)HexNAc(a1-"
  )
  path <- suppressWarnings(trace_biosynthesis(
    glycans,
    enzymes = "C1GALT1",
    max_steps = 1
  ))

  expect_equal(igraph::vcount(path), 2L)
  expect_identical(igraph::V(path)$target, rep(TRUE, 2L))
})

test_that("path_biosynthesis requires compatible endpoint structures", {
  expect_snapshot(
    path_biosynthesis("GalNAc(a1-", "HexNAc(a1-"),
    error = TRUE
  )
  expect_snapshot(
    path_biosynthesis("GalNAc(?1-", "Gal(b1-3)GalNAc(?1-"),
    error = TRUE
  )
})

test_that("path_biosynthesis_virtual requires compatible endpoints", {
  expect_snapshot(
    path_biosynthesis_virtual("GalNAc(a1-", "GalNAc(??-"),
    error = TRUE
  )
})

test_that("direct condensed parsing preserves names and auto-parser results", {
  strings <- c(core1 = "Gal(b1-3)GalNAc(a1-", core3 = "GlcNAc(b1-3)GalNAc(a1-")
  expected <- .process_biosynthesis_glycans_arg(glyparse::auto_parse(strings))
  actual <- .process_biosynthesis_glycans_arg(strings)
  expect_identical(as.character(actual), as.character(expected))
  expect_identical(
    lapply(
      glyrepr::get_structure_graphs(actual, return_list = TRUE),
      .bfs_native_graph
    ),
    lapply(
      glyrepr::get_structure_graphs(expected, return_list = TRUE),
      .bfs_native_graph
    )
  )
  # Unsupported syntax and malformed condensed input keep auto_parse's errors.
  for (input in c("invalid", "NotAResidue(a1-")) {
    old <- tryCatch(glyparse::auto_parse(input), error = conditionMessage)
    actual <- tryCatch(
      .process_biosynthesis_glycans_arg(input),
      error = conditionMessage
    )
    expect_identical(actual, old)
  }
})

test_that("structure level fast path agrees with the public classifier", {
  glycans <- glyrepr::as_glycan_structure(c(
    intact = "Gal(b1-3)GalNAc(a1-",
    single = "GalNAc(a1-",
    ambiguous = "Gal(b1-3/6)GalNAc(a1-",
    unknown = "Gal(b1-?)GalNAc(?1-",
    topology = "Gal(??-?)GalNAc(??-",
    missing = NA_character_
  ))
  for (x in c(
    list(glycans, glycans[0]),
    lapply(seq_along(glycans), function(i) glycans[i])
  )) {
    expect_identical(
      .glycan_structure_levels(x),
      suppressWarnings(glyrepr::get_structure_level(x))
    )
  }
  graph <- glyrepr::get_structure_graphs(glycans[1])
  graph <- igraph::set_graph_attr(
    graph,
    "floating_parts",
    list(
      list(mono = "Fuc", sub = "", linkage = "a1-?", parents = 1:2)
    )
  )
  floating <- glyrepr::new_glycan_structure("floating", list(floating = graph))
  expect_identical(
    .glycan_structure_levels(floating),
    glyrepr::get_structure_level(floating)
  )
})
