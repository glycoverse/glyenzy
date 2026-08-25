# Tests for `trace_biosynthesis()`.
# The concrete results of this function is difficult to test,
# so we only validate some properties of the resulting igraph.

test_that("N-glycan starting structures depend on the tracing function", {
  target <- glyparse::auto_parse(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )

  expect_equal(
    as.character(.decide_starting_glycan(target)),
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(
    as.character(.decide_virtual_starting_glycan(target)),
    as.character(glyrepr::n_glycan_core())
  )
})

test_that("non-N-glycan starting structures use the reducing-end core", {
  glycans <- glyparse::auto_parse(c(
    "Gal(b1-3)GalNAc(a1-",
    "GlcNAc(b1-3)GalNAc(a1-",
    "Glc(a1-2)Gal(?1-",
    "Man(a1-3)GlcNAc(b1-"
  ))

  starts <- purrr::map_chr(glycans, ~ as.character(.decide_starting_glycan(.x)))

  expect_identical(
    starts,
    c("GalNAc(a1-", "GalNAc(a1-", "Gal(b1-", "GlcNAc(b1-")
  )
})

test_that("trace_biosynthesis infers max steps from target sizes", {
  n_glycans <- glyparse::auto_parse(c(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  o_glycans <- glyparse::auto_parse(c(
    "GalNAc(a1-",
    "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  ))
  sulfated_glycan <- glyparse::auto_parse("Gal3S(b1-3)GalNAc(a1-")

  expect_null(formals(trace_biosynthesis)$max_steps)
  expect_identical(.infer_trace_max_steps(n_glycans), 10L)
  expect_identical(.infer_trace_max_steps(o_glycans), 3L)
  expect_identical(.infer_trace_max_steps(sulfated_glycan), 2L)
})

test_that("trace_biosynthesis includes substituent transfers in max steps", {
  glycan <- "Gal3S(b1-3)GalNAc(a1-"
  path <- trace_biosynthesis(
    glycan,
    enzymes = c("C1GALT1", "GAL3ST4")
  )

  expect_equal(igraph::E(path)$enzyme, c("C1GALT1", "GAL3ST4"))
  expect_equal(igraph::V(path)$name[igraph::V(path)$target], glycan)
})

test_that("trace_biosynthesis returns starting glycans without enzyme steps", {
  path <- trace_biosynthesis("GalNAc(a1-")

  vertices <- igraph::as_data_frame(path, what = "vertices")
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_s3_class(path, "glyenzy_biosynthesis_network")
  expect_equal(vertices$name, "GalNAc(a1-")
  expect_equal(vertices$target, TRUE)
  expect_equal(nrow(edges), 0L)
})

test_that("trace_biosynthesis works for a high-mannose N-glycan", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  path <- trace_biosynthesis(glycan)

  expect_s3_class(path, "glyenzy_biosynthesis_network")
  expect_s3_class(path, "igraph")

  # The path starts with the N-glycan precursor
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(
    root_node$name,
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )

  # The path ends with the input glycan
  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  expect_equal(length(end_node), 1L)
  expect_equal(end_node$name, glycan)
})

test_that("trace_biosynthesis works for a complex N-glycan", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  path <- trace_biosynthesis(glycan)

  # The path starts with the N-glycan precursor
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(
    root_node$name,
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )

  # The path ends with the input glycan
  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  expect_equal(length(end_node), 1L)
  expect_equal(end_node$name, glycan)

  # The path includes GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
  H3N3 <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_true(H3N3 %in% igraph::V(path)$name)
})

test_that("trace_biosynthesis works for an O-GalNAc glycan", {
  glycan <- "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  path <- trace_biosynthesis(glycan)

  # The path starts with GalNAc(a1-
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(root_node$name, "GalNAc(a1-")

  # The path ends with the input glycan
  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  expect_equal(length(end_node), 1L)
  expect_equal(end_node$name, glycan)

  # More than one possible routes are found
  expect_gt(igraph::ecount(path), igraph::vcount(path) - 1)
})

test_that("trace_biosynthesis reaches topological target glycans", {
  glycan <- glyrepr::remove_linkages(glyrepr::o_glycan_core_2())
  expect_warning(
    path <- trace_biosynthesis(
      glycan,
      enzymes = c("C1GALT1", "GCNT1"),
      max_steps = 2
    ),
    "non-intact glycan structures"
  )

  vertices <- igraph::as_data_frame(path, what = "vertices")
  expect_true(as.character(glycan) %in% vertices$name)
  expect_true(all(
    glyrepr::get_structure_level(glyparse::auto_parse(vertices$name)) ==
      "topological"
  ))
})

test_that("trace_biosynthesis reaches generic topological target glycans", {
  glycan <- glyrepr::remove_linkages(
    glyrepr::convert_to_generic(glyrepr::o_glycan_core_2())
  )
  path <- suppressWarnings(trace_biosynthesis(
    glycan,
    enzymes = c("C1GALT1", "GCNT1"),
    max_steps = 2
  ))

  vertices <- igraph::as_data_frame(path, what = "vertices")
  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  end_glycans <- glyparse::auto_parse(end_node$name)
  expect_identical(end_node$target, TRUE)
  expect_true(any(glymotif::have_motif(
    end_glycans,
    glycan,
    alignment = "whole",
    mode = "lenient"
  )))
  expect_true(all(
    glyrepr::get_structure_level(glyparse::auto_parse(vertices$name)) ==
      "topological"
  ))
})

test_that("trace_biosynthesis supports custom enzyme objects", {
  enz <- make_enzyme(
    name = "TEST_ST3GAL",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-3)GalNAc(a1-",
      acceptor_alignment = "core",
      rejects = NULL,
      product = "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
    ))
  )

  path <- trace_biosynthesis(
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-",
    enzymes = list(enzyme("C1GALT1"), enz),
    max_steps = 3
  )

  expect_true("TEST_ST3GAL" %in% igraph::E(path)$enzyme)
})

test_that("trace_biosynthesis supports custom ST objects", {
  enz <- make_enzyme(
    name = "TEST_ST",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "GalNAc(a1-",
      acceptor_alignment = "core",
      rejects = NULL,
      product = "GalNAc6S(a1-"
    ))
  )
  target <- "GalNAc6S(a1-"

  path <- trace_biosynthesis(
    target,
    enzymes = list(enz),
    max_steps = 1
  )

  expect_true("TEST_ST" %in% igraph::E(path)$enzyme)
  expect_true(target %in% igraph::V(path)$name)
})

test_that("trace_biosynthesis works with O-Man glycans", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-"
  path <- trace_biosynthesis(glycan)

  # The path starts with Man(a1-
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(root_node$name, "Man(a1-")

  # The path ends with the input glycan
  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  expect_equal(length(end_node), 1L)
  expect_equal(end_node$name, glycan)
})

test_that("trace_biosynthesis works with complete graph (all paths)", {
  glycan <- "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  path <- trace_biosynthesis(glycan)
  # Should return complete graph with multiple paths
  expect_gte(igraph::ecount(path), igraph::vcount(path) - 1)
})

test_that("trace_biosynthesis works with multiple target glycans", {
  glycans <- c(
    "Gal(b1-3)GalNAc(a1-", # core 1 O-glycan
    "GlcNAc(b1-3)GalNAc(a1-" # core 3 O-glycan
  )
  path <- trace_biosynthesis(
    glycans,
    max_steps = 3,
    enzymes = c("B3GNT6", "C1GALT1")
  )

  # The path should start with GalNAc(a1-
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(root_node$name, "GalNAc(a1-")

  # Both target glycans should be present as vertices
  all_vertices <- igraph::V(path)$name
  expect_true(all(glycans %in% all_vertices))
  expect_identical(
    igraph::V(path)$target,
    all_vertices %in% glycans
  )

  # Should be a valid connected graph
  expect_true(igraph::is_connected(path, mode = "weak"))
})

test_that("trace_biosynthesis works with overlapping synthesis paths", {
  # Test case where some targets are intermediates of others
  glycans <- c(
    "Gal(b1-3)GalNAc(a1-", # core 1 O-glycan, intermediate
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-" # core 2 O-glycan, extension of the first
  )
  path <- trace_biosynthesis(
    glycans,
    max_steps = 10,
    enzymes = c("GCNT1", "C1GALT1")
  )

  # Both targets should be reachable
  all_vertices <- igraph::V(path)$name
  expect_true(all(glycans %in% all_vertices))

  # The intermediate should be on the path to the final target
  expect_true("Gal(b1-3)GalNAc(a1-" %in% all_vertices)
  expect_true("Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-" %in% all_vertices)
})

test_that("trace_biosynthesis works for multiple targets (complex situation)", {
  # core 1, 2, 3, 4
  glycans <- c(
    "Gal(b1-3)GalNAc(a1-", # core 1 O-glycan
    "GlcNAc(b1-3)GalNAc(a1-", # core 3 O-glycan
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-", # core 2 O-glycan, extension of the first
    "GlcNAc(b1-3)[GlcNAc(b1-6)]GalNAc(a1-" # core 4 O-glycan, extension of the first
  )
  path <- trace_biosynthesis(glycans, max_steps = 3) # no enzymes provided

  # The path should start with GalNAc(a1-
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(root_node$name, "GalNAc(a1-")

  # Both target glycans should be present as vertices
  all_vertices <- igraph::V(path)$name
  expect_true(all(glycans %in% all_vertices))

  # Should be a valid connected graph
  expect_true(igraph::is_connected(path, mode = "weak"))
})
