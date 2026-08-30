test_that("all biosynthesis functions return one shared network schema", {
  from <- "GalNAc(a1-"
  to <- "Gal(b1-3)GalNAc(a1-"
  networks <- list(
    trace_biosynthesis(to, enzymes = "C1GALT1"),
    path_biosynthesis(from, to, enzymes = "C1GALT1", max_steps = 1),
    trace_biosynthesis_virtual(to),
    path_biosynthesis_virtual(from, to)
  )

  expect_equal(
    lapply(networks, igraph::vertex_attr_names),
    rep(list(c("name", "target")), 4L)
  )
  expect_equal(
    lapply(networks, igraph::edge_attr_names),
    rep(list(c("enzyme", "is_virtual", "step", "enzymes")), 4L)
  )
  expect_equal(
    vapply(networks, igraph::any_multiple, logical(1)),
    rep(FALSE, 4L)
  )
})

test_that("known isoenzymes share one substrate-product edge", {
  from <- "GlcNAc(b1-"
  to <- "Gal(b1-4)GlcNAc(b1-"
  candidates <- c("B4GALT1", "B4GALT2", "B4GALT3")
  network <- path_biosynthesis(
    from,
    to,
    enzymes = candidates,
    max_steps = 1
  )
  edges <- igraph::as_data_frame(network, what = "edges")

  expect_equal(nrow(edges), 1L)
  expect_equal(edges$enzyme, paste(candidates, collapse = " / "))
  expect_equal(edges$enzymes, list(candidates))
  expect_identical(edges$is_virtual, FALSE)
  expect_identical(edges$step, 1L)
  expect_equal(igraph::V(network)$name[igraph::V(network)$target], to)
})

test_that("zero-edge results retain the shared network schema", {
  glycan <- "GalNAc(a1-"
  networks <- list(
    trace_biosynthesis(glycan),
    path_biosynthesis(glycan, glycan),
    trace_biosynthesis_virtual(glycan),
    path_biosynthesis_virtual(glycan, glycan)
  )

  expect_equal(
    lapply(networks, igraph::vertex_attr_names),
    rep(list(c("name", "target")), 4L)
  )
  expect_equal(
    lapply(networks, igraph::edge_attr_names),
    rep(list(c("enzyme", "is_virtual", "step", "enzymes")), 4L)
  )
  expect_equal(vapply(networks, igraph::ecount, numeric(1)), rep(0, 4L))
})

test_that("known evidence supersedes a duplicate virtual edge", {
  from <- "GalNAc(a1-"
  to <- "Gal(b1-3)GalNAc(a1-"
  graph <- igraph::graph_from_data_frame(
    data.frame(
      from = c(from, from),
      to = c(to, to),
      enzyme = c("C1GALT1", "b3GalT"),
      is_virtual = c(FALSE, TRUE),
      step = c(1L, 1L)
    ),
    directed = TRUE
  )

  network <- .finalize_biosynthesis_network(graph, to)
  edges <- igraph::as_data_frame(network, what = "edges")

  expect_equal(nrow(edges), 1L)
  expect_equal(edges$enzyme, "C1GALT1")
  expect_equal(edges$enzymes, list("C1GALT1"))
  expect_identical(edges$is_virtual, FALSE)
})
