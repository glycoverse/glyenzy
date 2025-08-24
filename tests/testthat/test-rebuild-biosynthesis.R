# Tests for `rebuild_biosynthesis()`.
# The concrete results of this function is difficult to test,
# so we only validate some properties of the resulting igraph.

test_that("rebuild_biosynthesis works for a high-mannose N-glycan", {
  glycan <- "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  path <- rebuild_biosynthesis(glycan)

  # The path starts with the N-glycan precursor
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(root_node$name, "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")

  # The path ends with the input glycan
  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  expect_equal(length(end_node), 1L)
  expect_equal(end_node$name, glycan)
})

test_that("rebuild_biosynthesis works for a complex N-glycan", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  path <- rebuild_biosynthesis(glycan)

  # The path starts with the N-glycan precursor
  root_node <- igraph::V(path)[igraph::degree(path, mode = "in") == 0]
  expect_equal(length(root_node), 1L)
  expect_equal(root_node$name, "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")

  # The path ends with the input glycan
  end_node <- igraph::V(path)[igraph::degree(path, mode = "out") == 0]
  expect_equal(length(end_node), 1L)
  expect_equal(end_node$name, glycan)

  # The path includes GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
  H3N3 <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  expect_true(H3N3 %in% igraph::V(path)$name)
})

test_that("rebuild_biosynthesis works for an O-GalNAc glycan", {
  glycan <- "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  path <- rebuild_biosynthesis(glycan)

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

test_that("rebuild_biosynthesis works with complete graph (all paths)", {
  glycan <- "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  path <- rebuild_biosynthesis(glycan)
  # Should return complete graph with multiple paths
  expect_gte(igraph::ecount(path), igraph::vcount(path) - 1)
})

test_that("rebuild_biosynthesis works with multiple target glycans", {
  glycans <- c(
    "Gal(b1-3)GalNAc(a1-",  # core 1 O-glycan
    "GlcNAc(b1-3)GalNAc(a1-"  # core 3 O-glycan
  )
  path <- rebuild_biosynthesis(glycans, max_steps = 3, enzymes = c("B3GNT6", "C1GALT1"))

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

test_that("rebuild_biosynthesis works with overlapping synthesis paths", {
  # Test case where some targets are intermediates of others
  glycans <- c(
    "Gal(b1-3)GalNAc(a1-",  # core 1 O-glycan, intermediate
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"  # core 2 O-glycan, extension of the first
  )
  path <- rebuild_biosynthesis(glycans, max_steps = 10, enzymes = c("GCNT1", "C1GALT1"))
  
  # Both targets should be reachable
  all_vertices <- igraph::V(path)$name
  expect_true(all(glycans %in% all_vertices))
  
  # The intermediate should be on the path to the final target
  expect_true("Gal(b1-3)GalNAc(a1-" %in% all_vertices)
  expect_true("Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-" %in% all_vertices)
})

test_that("rebuild_biosynthesis works for multiple targets (complex situation)", {
  # core 1, 2, 3, 4
  glycans <- c(
    "Gal(b1-3)GalNAc(a1-",  # core 1 O-glycan
    "GlcNAc(b1-3)GalNAc(a1-",  # core 3 O-glycan
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",  # core 2 O-glycan, extension of the first
    "GlcNAc(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"  # core 4 O-glycan, extension of the first
  )
  path <- rebuild_biosynthesis(glycans, max_steps = 3)  # no enzymes provided

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