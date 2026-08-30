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
  expect_s3_class(plot$layers[[4]]$geom, "GeomGlycan")
  expect_s3_class(plot$layers[[4]]$stat, "StatFilter")
  expect_equal(igraph::ecount(collapsed), 5L)
  expect_true(
    "C1GALT1, ALT-C1" %in%
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

test_that("biosynthesis autoplot supports all glydraw orientations", {
  skip_if_not_installed("ggraph")
  graph <- trace_biosynthesis("Gal(b1-3)GalNAc(a1-")
  orientations <- c("left", "right", "up", "down")

  plots <- lapply(orientations, function(orient) {
    plot <- ggplot2::autoplot(
      graph,
      orient = orient,
      show_enzyme = FALSE
    )
    ggplot2::ggplot_build(plot)
    plot
  })

  expect_equal(
    vapply(
      plots,
      function(plot) plot$layers[[3]]$geom_params$orient,
      character(1)
    ),
    orientations
  )
})

test_that("biosynthesis autoplot highlights targets in multi-target networks", {
  skip_if_not_installed("ggraph")
  targets <- c(
    "GlcNAc(b1-4)Gal(b1-",
    "Fuc(a1-2)Gal(b1-"
  )
  graph <- trace_biosynthesis_virtual(targets)
  highlighted <- ggplot2::autoplot(graph)
  plain <- ggplot2::autoplot(graph, highlight_target = FALSE)
  highlighted_build <- ggplot2::ggplot_build(highlighted)
  plain_build <- ggplot2::ggplot_build(plain)

  expect_length(highlighted$layers, 5L)
  expect_identical(
    highlighted$layers[[4]]$geom_params$highlight,
    integer()
  )
  expect_s3_class(
    highlighted$layers[[5]]$geom,
    "GeomCustomAnn"
  )
  expect_length(
    highlighted$layers[[5]]$geom_params$grob$children,
    2L
  )
  expect_equal(
    highlighted_build$data[1:3],
    plain_build$data[1:3]
  )

  gtable <- ggplot2::ggplotGrob(highlighted)
  panel <- gtable$grobs[[which(gtable$layout$name == "panel")]]
  target_layer <- panel$children[[
    grep(
      "^biosynthesis_target_glycans",
      grid::childNames(panel),
      value = TRUE
    )
  ]]
  alphas <- lapply(
    list(
      panel$children[["geom_glycan"]],
      target_layer
    ),
    function(layer) {
      unname(unlist(lapply(
        layer$children,
        \(grob) unique(grob$polygon_coor$alpha)
      )))
    }
  )
  expect_equal(alphas, list(0.3, c(1, 1)))
})

test_that("target highlighting defaults off for single-target networks", {
  skip_if_not_installed("ggraph")
  graph <- trace_biosynthesis_virtual("GlcNAc(b1-4)Gal(b1-")

  default <- ggplot2::autoplot(graph)
  highlighted <- ggplot2::autoplot(graph, highlight_target = TRUE)

  expect_length(default$layers, 4L)
  expect_null(default$layers[[4]]$geom_params$highlight)
  expect_length(highlighted$layers, 5L)
  expect_length(
    highlighted$layers[[5]]$geom_params$grob$children,
    1L
  )
})

test_that("target highlighting supports path networks", {
  skip_if_not_installed("ggraph")
  target <- "Gal(b1-3)GalNAc(a1-"
  networks <- list(
    concrete = suppressMessages(path_biosynthesis(
      "GalNAc(a1-",
      target,
      enzymes = "C1GALT1",
      max_steps = 1
    )),
    virtual = path_biosynthesis_virtual("GalNAc(a1-", target)
  )

  expect_equal(
    vapply(
      networks,
      \(network) sum(igraph::vertex_attr(network, "target")),
      integer(1)
    ),
    c(concrete = 1, virtual = 1)
  )
  plots <- lapply(
    networks,
    ggplot2::autoplot,
    highlight_target = TRUE
  )
  expect_equal(
    vapply(plots, \(plot) length(plot$layers), integer(1)),
    c(concrete = 5L, virtual = 5L)
  )
})

test_that("figure dimensions support equivalent physical units", {
  skip_if_not_installed("ggraph")
  graph <- trace_biosynthesis("Gal(b1-3)GalNAc(a1-")
  sizes <- list(
    inches = c(width = 8, height = 5),
    centimetres = c(width = 20.32, height = 12.7),
    millimetres = c(width = 203.2, height = 127)
  )
  plots <- list(
    ggplot2::autoplot(
      graph,
      width = sizes$inches[["width"]],
      height = sizes$inches[["height"]],
      units = "in"
    ),
    ggplot2::autoplot(
      graph,
      width = sizes$centimetres[["width"]],
      height = sizes$centimetres[["height"]],
      units = "cm"
    ),
    ggplot2::autoplot(
      graph,
      width = sizes$millimetres[["width"]],
      height = sizes$millimetres[["height"]],
      units = "mm"
    )
  )

  for (plot in plots[-1]) {
    expect_equal(plot$data$x, plots[[1]]$data$x)
    expect_equal(plot$data$y, plots[[1]]$data$y)
    expect_equal(plot$data$.node_width, plots[[1]]$data$.node_width)
    expect_equal(plot$data$.node_height, plots[[1]]$data$.node_height)
    expect_equal(
      attr(plot, "glyenzy_panel_size_in"),
      attr(plots[[1]], "glyenzy_panel_size_in")
    )
    expect_equal(
      attr(plot, "glyenzy_layout_scale"),
      attr(plots[[1]], "glyenzy_layout_scale")
    )
    expect_equal(plot$theme$plot.margin, plots[[1]]$theme$plot.margin)
  }
})

test_that("figure dimensions are validated", {
  skip_if_not_installed("ggraph")
  graph <- trace_biosynthesis("Gal(b1-3)GalNAc(a1-")

  expect_snapshot(
    ggplot2::autoplot(graph, width = 0),
    error = TRUE
  )
  expect_snapshot(
    ggplot2::autoplot(graph, width = -1),
    error = TRUE
  )
  expect_snapshot(
    ggplot2::autoplot(graph, height = Inf),
    error = TRUE
  )
  expect_snapshot(
    ggplot2::autoplot(graph, units = "px"),
    error = TRUE
  )
  expect_snapshot(
    ggplot2::autoplot(graph, max_panel_width = 5),
    error = TRUE
  )
  expect_snapshot(
    ggplot2::autoplot(graph, max_panel_height = 5),
    error = TRUE
  )
})

test_that("graphics plot dispatches to the biosynthesis network method", {
  skip_if_not_installed("ggraph")
  graph <- trace_biosynthesis("Gal(b1-3)GalNAc(a1-")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  plot <- graphics::plot(graph)

  expect_s3_class(plot, "ggraph")
})

test_that("biosynthesis enzyme labels support full and condensed styles", {
  enzymes <- c(
    "B4GALT1",
    "B4GALT2",
    "B4GALT3",
    "B3GALT3",
    "B3GALT4"
  )

  expect_identical(
    .biosynthesis_enzyme_label(enzymes, style = "full"),
    paste(enzymes, collapse = " / ")
  )
  expect_identical(
    .biosynthesis_enzyme_label(enzymes, style = "condensed"),
    "B4GALT1/2/3, B3GALT3/4"
  )
  expect_identical(
    .biosynthesis_enzyme_label(enzymes),
    "B4GALT1/2/3, B3GALT3/4"
  )
  expect_identical(
    .biosynthesis_enzyme_label(
      c("b3GalT", "B4GALT1", "B4GALT1", "B4GALT2", "3SulfoT"),
      style = "condensed"
    ),
    "b3GalT, B4GALT1/2, 3SulfoT"
  )
})

test_that("canonical candidate lists do not repeat display labels", {
  candidates <- c("B4GALT1", "B4GALT2", "B4GALT3")
  graph <- path_biosynthesis(
    "GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-",
    enzymes = candidates,
    max_steps = 1
  )

  full <- .collapse_biosynthesis_reactions(
    graph,
    enzyme_label_style = "full"
  )
  condensed <- .collapse_biosynthesis_reactions(graph)

  expect_identical(
    igraph::E(full)$enzyme,
    paste(candidates, collapse = " / ")
  )
  expect_identical(igraph::E(condensed)$enzyme, "B4GALT1/2/3")
})

test_that("biosynthesis autoplot uses condensed multi-enzyme labels", {
  skip_if_not_installed("ggraph")
  enzymes <- c(
    "B4GALT1",
    "B4GALT2",
    "B4GALT3",
    "B3GALT3",
    "B3GALT4"
  )
  graph <- igraph::graph_from_data_frame(
    data.frame(
      from = rep("GalNAc(a1-", length(enzymes)),
      to = rep("Gal(b1-3)GalNAc(a1-", length(enzymes)),
      enzyme = enzymes
    ),
    directed = TRUE
  ) |>
    .new_biosynthesis_network()

  full <- ggplot2::autoplot(
    graph,
    enzyme_label_style = "full",
    width = 50,
    height = 50
  )
  condensed <- ggplot2::autoplot(
    graph,
    width = 50,
    height = 50
  )
  full_label <- ggplot2::ggplot_build(full)$data[[3]]$label
  condensed_label <- ggplot2::ggplot_build(condensed)$data[[3]]$label

  expect_identical(full_label, paste(enzymes, collapse = " / "))
  expect_identical(
    condensed_label,
    "B4GALT1/2/3, B3GALT3/4"
  )
  expect_lt(
    attr(condensed, "glyenzy_panel_size_in")[["width"]],
    attr(full, "glyenzy_panel_size_in")[["width"]]
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

test_that("large N-glycan networks fit the default figure", {
  skip_if_not_installed("ggraph")
  target <- paste0(
    "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  graph <- trace_biosynthesis(target)
  plot <- ggplot2::autoplot(graph)
  natural <- ggplot2::autoplot(
    graph,
    width = 50,
    height = 50
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

test_that("unrenderable enzyme labels are hidden with a warning", {
  skip_if_not_installed("ggraph")
  target <- paste0(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  graph <- trace_biosynthesis(target)

  expect_snapshot(
    plot <- ggplot2::autoplot(graph)
  )

  expect_length(plot$layers, 3L)
  expect_equal(
    unname(vapply(
      plot$layers,
      \(layer) inherits(layer$geom, "GeomLabel"),
      logical(1)
    )),
    rep(FALSE, 3L)
  )

  larger <- expect_no_warning(ggplot2::autoplot(
    graph,
    width = 12,
    height = 12
  ))
  expect_s3_class(larger$layers[[3]]$geom, "GeomLabel")
})

test_that("converging N-glycan networks use a balanced layered layout", {
  skip_if_not_installed("ggraph")
  target <- paste0(
    "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
    "[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
  )
  graph <- trace_biosynthesis(target)
  plot <- ggplot2::autoplot(
    graph,
    width = 50,
    height = 50
  )
  dimensions <- .biosynthesis_node_dimensions(
    .collapse_biosynthesis_reactions(graph),
    size = 0.4,
    show_linkage = FALSE,
    orient = "left",
    dots = list()
  )
  rank_midpoints <- vapply(
    split(plot$data$x, plot$data$y),
    function(x) mean(range(x)),
    numeric(1)
  )
  root_x <- plot$data$x[which.max(plot$data$y)]
  x_spacing <- max(dimensions$width) + 0.25

  expect_lte(max(abs(rank_midpoints - root_x)), 4 * x_spacing + 1e-8)
})

test_that("figure size includes and centers the complete panel", {
  skip_if_not_installed("ggraph")
  target <- "Gal(b1-3)GalNAc(a1-"
  plot <- ggplot2::autoplot(
    trace_biosynthesis(target),
    show_enzyme = FALSE,
    width = 8,
    height = 5
  )

  grDevices::pdf(NULL, width = 8, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  gtable <- ggplot2::ggplotGrob(plot)
  complete_size <- c(
    width = grid::convertWidth(
      sum(gtable$widths),
      "in",
      valueOnly = TRUE
    ),
    height = grid::convertHeight(
      sum(gtable$heights),
      "in",
      valueOnly = TRUE
    )
  )

  plot_margin <- plot$theme$plot.margin
  expect_equal(complete_size, c(width = 8, height = 5))
  expect_equal(attr(plot, "glyenzy_layout_scale"), 1)
  expect_equal(plot_margin[[1]], plot_margin[[3]])
  expect_equal(plot_margin[[2]], plot_margin[[4]])
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

test_that("biosynthesis edges use typed line styles without a legend", {
  skip_if_not_installed("ggraph")
  graph <- path_biosynthesis(
    "GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-",
    enzymes = "ST3GAL1",
    max_steps = 2,
    max_virtual_steps = 1
  )
  plot <- ggplot2::autoplot(graph)
  built <- ggplot2::ggplot_build(plot)

  expect_identical(
    unique(built$data[[1]]$edge_colour),
    .biosynthesis_default_edge_color
  )
  expect_setequal(
    unique(built$data[[1]]$edge_linetype),
    c("solid", "22")
  )
  expect_identical(
    unique(built$data[[2]]$edge_linetype),
    "solid"
  )
  expect_identical(
    unique(built$data[[2]]$edge_width),
    0
  )
  expect_s3_class(plot$layers[[2]]$geom_params$arrow, "arrow")
  expect_identical(plot$layers[[1]]$show.legend, FALSE)
  expect_identical(
    plot$scales$get_scales("edge_linetype")$guide,
    "none"
  )
  expect_length(built$plot$guides$guides, 0L)
  expect_identical(
    unique(built$data[[3]]$colour),
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
    built$data[[3]]$colour,
    built$data[[3]]$label
  )
  expect_setequal(
    unique(built$data[[1]]$edge_colour),
    c("#FFD400", "#A54399")
  )
  expect_setequal(
    unique(built$data[[1]]$edge_linetype),
    c("solid", "22")
  )
  expect_setequal(
    unique(built$data[[2]]$edge_colour),
    c("#FFD400", "#A54399")
  )
  expect_identical(
    unique(built$data[[2]]$edge_linetype),
    "solid"
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
