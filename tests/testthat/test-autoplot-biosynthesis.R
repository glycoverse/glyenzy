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

test_that("long labels do not overlap on converging O-glycan routes", {
  skip_if_not_installed("ggraph")
  target <- paste0(
    "Neu5Ac(a2-3)Gal(b1-3)",
    "[Gal(b1-4)GlcNAc(b1-6)]GalNAc(a1-"
  )
  graphs <- list(
    concrete = trace_biosynthesis(target),
    virtual = trace_biosynthesis_virtual(target)
  )
  excess_clearances <- lapply(graphs, function(graph) {
    collapsed <- .collapse_biosynthesis_reactions(graph)
    plot <- ggplot2::autoplot(graph)
    scale <- attr(plot, "glyenzy_layout_scale")
    labels <- .biosynthesis_edge_labels(
      collapsed,
      plot$data,
      scale = scale
    )
    same_rank <- combn(seq_len(nrow(labels)), 2)
    same_rank <- same_rank[,
      labels$.rank[same_rank[1, ]] == labels$.rank[same_rank[2, ]],
      drop = FALSE
    ]
    apply(same_rank, 2, function(pair) {
      abs(diff(labels$x[pair])) -
        sum(labels$.label_width[pair]) / 2 -
        .biosynthesis_label_gap * scale
    })
  })

  expect_gte(min(unlist(excess_clearances)), -1e-8)
})

test_that("large N-glycan networks fit the default panel limits", {
  skip_if_not_installed("ggraph")
  target <- paste0(
    "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  graph <- trace_biosynthesis(target)
  plot <- ggplot2::autoplot(graph)
  natural <- ggplot2::autoplot(
    graph,
    max_panel_width = Inf,
    max_panel_height = Inf
  )
  output <- tempfile(fileext = ".png")
  on.exit(unlink(output), add = TRUE)
  ggplot2::ggsave(
    output,
    plot,
    width = 7,
    height = 7,
    units = "in",
    dpi = 72
  )

  panel_size <- attr(plot, "glyenzy_panel_size_in")
  natural_size <- attr(natural, "glyenzy_panel_size_in")
  expect_lte(unname(panel_size[["width"]]), 6)
  expect_lte(unname(panel_size[["height"]]), 6)
  expect_gt(unname(natural_size[["height"]]), 6)
  expect_lt(attr(plot, "glyenzy_layout_scale"), 1)
  expect_gt(file.info(output)$size, 0)
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

test_that("biosynthesis edges default to dark grey with typed line styles", {
  skip_if_not_installed("ggraph")
  graph <- path_biosynthesis(
    "GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-",
    enzymes = "ST3GAL1",
    max_steps = 2,
    max_virtual_steps = 1
  )
  built <- ggplot2::ggplot_build(ggplot2::autoplot(graph))

  expect_identical(
    unique(built$data[[1]]$edge_colour),
    .biosynthesis_default_edge_color
  )
  expect_setequal(
    unique(built$data[[1]]$edge_linetype),
    c("solid", "22")
  )
  expect_identical(
    unique(built$data[[2]]$colour),
    .biosynthesis_default_edge_color
  )
})

test_that("colored biosynthesis edges match added residue colors", {
  skip_if_not_installed("ggraph")
  mixed <- path_biosynthesis(
    "GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-",
    enzymes = "ST3GAL1",
    max_steps = 2,
    max_virtual_steps = 1
  )
  mixed_graph <- mixed |>
    .collapse_biosynthesis_reactions() |>
    .color_biosynthesis_reactions()
  mixed_edges <- igraph::as_data_frame(mixed_graph, what = "edges")

  expect_equal(mixed_edges$.added_residue, c("Gal", "Neu5Ac"))
  expect_equal(mixed_edges$.edge_colour, c("#FFD400", "#A54399"))

  built <- ggplot2::autoplot(mixed, color_edge = TRUE) |>
    ggplot2::ggplot_build()
  label_colors <- setNames(
    built$data[[2]]$colour,
    built$data[[2]]$label
  )
  expect_setequal(
    unique(built$data[[1]]$edge_colour),
    c("#FFD400", "#A54399")
  )
  expect_setequal(
    unique(built$data[[1]]$edge_linetype),
    c("solid", "22")
  )
  expect_equal(
    unname(label_colors[c("b3GalT", "ST3GAL1")]),
    c("#FFD400", "#A54399")
  )

  glcnac <- path_biosynthesis_virtual(
    "GalNAc(a1-",
    "GlcNAc(b1-6)GalNAc(a1-"
  ) |>
    .collapse_biosynthesis_reactions() |>
    .color_biosynthesis_reactions()
  expect_identical(
    igraph::edge_attr(glcnac, ".added_residue"),
    "GlcNAc"
  )
  expect_identical(
    igraph::edge_attr(glcnac, ".edge_colour"),
    "#0072BC"
  )

  sulfate <- path_biosynthesis_virtual(
    "Gal(b1-3)GalNAc(a1-",
    "Gal6S(b1-3)GalNAc(a1-"
  ) |>
    .collapse_biosynthesis_reactions() |>
    .color_biosynthesis_reactions()
  expect_identical(
    igraph::edge_attr(sulfate, ".added_residue"),
    "S"
  )
  expect_identical(
    igraph::edge_attr(sulfate, ".edge_colour"),
    "black"
  )

  removal <- igraph::graph_from_data_frame(
    data.frame(
      from = "Gal(b1-3)GalNAc(a1-",
      to = "GalNAc(a1-",
      enzyme = "GH"
    ),
    directed = TRUE
  ) |>
    .new_biosynthesis_network() |>
    .collapse_biosynthesis_reactions() |>
    .color_biosynthesis_reactions()
  expect_true(is.na(igraph::edge_attr(removal, ".added_residue")))
  expect_identical(
    igraph::edge_attr(removal, ".edge_colour"),
    .biosynthesis_default_edge_color
  )
})

test_that("biosynthesis residue colors stay aligned with glyrepr", {
  glyrepr_color <- get(
    "get_mono_color",
    envir = asNamespace("glyrepr")
  )
  monosaccharides <- glyrepr::available_monosaccharides()

  expect_identical(
    vapply(
      monosaccharides,
      .biosynthesis_residue_color,
      character(1)
    ),
    vapply(monosaccharides, glyrepr_color, character(1))
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
