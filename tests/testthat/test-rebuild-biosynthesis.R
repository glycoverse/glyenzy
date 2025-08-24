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