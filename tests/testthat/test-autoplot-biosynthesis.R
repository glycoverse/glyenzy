test_that("biosynthesis autoplot draws typed networks as glycan trees", {
  skip_if_not_installed("ggraph")
  root <- "GalNAc(a1-"
  core_1 <- "Gal(b1-3)GalNAc(a1-"
  branch <- "GlcNAc(b1-6)GalNAc(a1-"
  core_2 <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  target <- "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  graph <- igraph::graph_from_data_frame(
    data.frame(
      from = c(root, root, root, core_1, branch, core_2, core_2),
      to = c(core_1, core_1, branch, core_2, core_2, target, target),
      enzyme = c(
        "C1GALT1",
        "ALT-C1",
        "GCNT1",
        "GCNT1",
        "C1GALT1",
        "B4GALT1",
        "B4GALT2"
      )
    ),
    directed = TRUE
  )
  graph <- .new_biosynthesis_network(graph)

  plot <- ggplot2::autoplot(graph)
  collapsed <- .collapse_biosynthesis_reactions(graph)

  expect_s3_class(plot, "ggraph")
  expect_s3_class(plot$layers[[3]]$geom, "GeomGlycan")
  expect_s3_class(plot$layers[[3]]$stat, "StatFilter")
  expect_equal(igraph::ecount(collapsed), 5L)
  expect_true(
    "C1GALT1 / ALT-C1" %in%
      igraph::edge_attr(collapsed, "enzyme")
  )
  expect_equal(
    unname(plot$data$y[plot$data$name == root]),
    max(plot$data$y)
  )
  expect_equal(
    unname(plot$data$y[plot$data$name == target]),
    min(plot$data$y)
  )
  expect_named(
    attr(plot, "glyenzy_panel_size_in"),
    c("width", "height")
  )
})

test_that("biosynthesis tree geometry keeps nodes and edges apart", {
  skip_if_not_installed("ggraph")
  root <- "GalNAc(a1-"
  core_1 <- "Gal(b1-3)GalNAc(a1-"
  branch <- "GlcNAc(b1-6)GalNAc(a1-"
  core_2 <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  target <- "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  graph <- igraph::graph_from_data_frame(
    data.frame(
      from = c(root, root, core_1, branch, core_2),
      to = c(core_1, branch, core_2, core_2, target),
      enzyme = c(
        "C1GALT1",
        "GCNT1",
        "GCNT1",
        "C1GALT1",
        "B4GALT1"
      )
    ),
    directed = TRUE
  )
  graph <- .new_biosynthesis_network(graph)
  layout <- ggplot2::autoplot(graph)$data

  node_pairs <- combn(seq_len(nrow(layout)), 2)
  node_overlaps <- apply(node_pairs, 2, function(pair) {
    horizontal <- abs(diff(layout$x[pair])) < sum(layout$.node_width[pair]) / 2
    vertical <- abs(diff(layout$y[pair])) < sum(layout$.node_height[pair]) / 2
    horizontal && vertical
  })
  expect_equal(node_overlaps, rep(FALSE, ncol(node_pairs)))

  cap <- unclass(layout$.node_cap)
  expect_equal(
    cap$width,
    layout$.node_width + 2 * .biosynthesis_edge_padding
  )
  expect_equal(
    cap$height,
    layout$.node_height + 2 * .biosynthesis_edge_padding
  )

  edges <- igraph::as_data_frame(graph, what = "edges")
  edge_hits <- logical()
  for (edge_index in seq_len(nrow(edges))) {
    from <- match(edges$from[[edge_index]], layout$name)
    to <- match(edges$to[[edge_index]], layout$name)
    other_nodes <- setdiff(seq_len(nrow(layout)), c(from, to))
    along <- seq(0, 1, length.out = 1001)
    edge_x <- layout$x[from] +
      along * (layout$x[to] - layout$x[from])
    edge_y <- layout$y[from] +
      along * (layout$y[to] - layout$y[from])
    for (node in other_nodes) {
      edge_hits <- c(
        edge_hits,
        any(
          abs(edge_x - layout$x[node]) < layout$.node_width[node] / 2 &
            abs(edge_y - layout$y[node]) < layout$.node_height[node] / 2
        )
      )
    }
  }
  expect_equal(edge_hits, rep(FALSE, length(edge_hits)))
})

test_that("biosynthesis tree geometry keeps enzyme labels apart", {
  skip_if_not_installed("ggraph")
  graph <- igraph::graph_from_data_frame(
    data.frame(
      from = c(
        "GalNAc(a1-",
        "GalNAc(a1-",
        "Gal(b1-3)GalNAc(a1-",
        "GlcNAc(b1-6)GalNAc(a1-"
      ),
      to = c(
        "Gal(b1-3)GalNAc(a1-",
        "GlcNAc(b1-6)GalNAc(a1-",
        "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",
        "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
      ),
      enzyme = c("C1GALT1", "GCNT1", "GCNT1", "C1GALT1")
    ),
    directed = TRUE
  )
  graph <- .new_biosynthesis_network(graph)
  plot <- ggplot2::autoplot(graph)
  labels <- .biosynthesis_edge_labels(graph, plot$data)

  label_pairs <- combn(seq_len(nrow(labels)), 2)
  label_overlaps <- apply(label_pairs, 2, function(pair) {
    horizontal <- abs(diff(labels$x[pair])) < sum(labels$.label_width[pair]) / 2
    vertical <- abs(diff(labels$y[pair])) < sum(labels$.label_height[pair]) / 2
    horizontal && vertical
  })
  expect_equal(label_overlaps, rep(FALSE, ncol(label_pairs)))

  label_node_overlaps <- unlist(lapply(
    seq_len(nrow(labels)),
    function(label) {
      abs(labels$x[label] - plot$data$x) <
        (labels$.label_width[label] + plot$data$.node_width) / 2 &
        abs(labels$y[label] - plot$data$y) <
          (labels$.label_height[label] + plot$data$.node_height) / 2
    }
  ))
  expect_equal(
    unname(label_node_overlaps),
    rep(FALSE, length(label_node_overlaps))
  )
})

test_that("all biosynthesis functions return networks that can be plotted", {
  skip_if_not_installed("ggraph")
  concrete_target <- "Gal(b1-3)GalNAc(a1-"
  networks <- list(
    suppressMessages(trace_biosynthesis(
      concrete_target,
      enzymes = "C1GALT1"
    )),
    suppressMessages(path_biosynthesis(
      "GalNAc(a1-",
      concrete_target,
      enzymes = "C1GALT1",
      max_steps = 1
    )),
    trace_biosynthesis_virtual(concrete_target),
    path_biosynthesis_virtual("GalNAc(a1-", concrete_target)
  )

  plots <- lapply(networks, ggplot2::autoplot)
  builds <- lapply(plots, ggplot2::ggplot_build)

  expect_length(plots, 4L)
  expect_equal(
    vapply(plots, inherits, logical(1), "ggraph"),
    rep(TRUE, 4L)
  )
  expect_equal(
    vapply(builds, inherits, logical(1), "ggplot_built"),
    rep(TRUE, 4L)
  )
})

test_that("biosynthesis autoplot handles single-node networks", {
  skip_if_not_installed("ggraph")
  graph <- igraph::make_empty_graph(n = 1, directed = TRUE) |>
    igraph::set_vertex_attr("name", value = "GalNAc(a1-") |>
    .new_biosynthesis_network(virtual = TRUE)

  plot <- ggplot2::autoplot(graph, show_enzyme = FALSE)

  expect_s3_class(plot, "ggraph")
  expect_length(plot$layers, 1L)
  expect_s3_class(plot$layers[[1]]$geom, "GeomGlycan")
  expect_equal(nrow(plot$data), 1L)
})
