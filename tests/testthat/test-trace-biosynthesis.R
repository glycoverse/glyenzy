# Tests for `trace_biosynthesis()`.
# The concrete results of this function is difficult to test,
# so we only validate some properties of the resulting igraph.

test_that("trace_biosynthesis works for a high-mannose N-glycan", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
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

test_that("trace_biosynthesis reaches basic target glycans", {
  glycan <- glyrepr::reduce_structure_level(glyrepr::o_glycan_core_2(), "basic")
  path <- suppressWarnings(trace_biosynthesis(
    glycan,
    enzymes = c("C1GALT1", "GCNT1"),
    max_steps = 2
  ))

  vertices <- igraph::as_data_frame(path, what = "vertices")
  expect_true(as.character(glycan) %in% vertices$name)
  expect_true(all(
    glyrepr::get_structure_level(glyparse::auto_parse(vertices$name)) == "basic"
  ))
})

test_that("trace_biosynthesis matches partial targets with whole alignment", {
  glycan <- glyparse::auto_parse("Gal(b1-3)[GlcNAc(b1-?)]GalNAc(a1-")
  path <- suppressWarnings(trace_biosynthesis(
    glycan,
    enzymes = c("C1GALT1", "GCNT1"),
    max_steps = 2
  ))

  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  end_glycans <- glyparse::auto_parse(end_node$name)
  expect_true(any(glymotif::have_motif(
    end_glycans,
    glycan,
    alignment = "whole",
    mode = "lenient"
  )))
})

test_that("trace_biosynthesis records repeated initial endpoints", {
  glycans <- glyparse::auto_parse(c(
    "GalNAc(?1-",
    "GalNAc(a?-",
    "Gal(b1-3)GalNAc(?1-"
  ))

  path <- suppressWarnings(trace_biosynthesis(
    glycans,
    enzymes = "C1GALT1",
    max_steps = 1
  ))

  vertices <- igraph::as_data_frame(path, what = "vertices")
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_false(anyNA(vertices$name))
  expect_true("GalNAc(a1-" %in% vertices$name)
  expect_true("Gal(b1-3)GalNAc(a1-" %in% vertices$name)
  expect_equal(edges$enzyme, "C1GALT1")
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
