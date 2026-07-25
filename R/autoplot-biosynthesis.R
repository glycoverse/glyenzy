.biosynthesis_edge_padding <- 0.04
.biosynthesis_plot_padding <- 0.15
.biosynthesis_label_size <- 3
.biosynthesis_label_padding <- 0.06
.biosynthesis_label_gap <- 0.08
.biosynthesis_points_per_mm <- 72.27 / 25.4

#' Plot a glycan biosynthesis network
#'
#' `autoplot.glyenzy_biosynthesis_network()` draws the typed networks returned
#' by [trace_biosynthesis()], [trace_biosynthesis_virtual()],
#' [path_biosynthesis()], and [path_biosynthesis_virtual()]. A spanning tree
#' determines the node positions without discarding converging biosynthesis
#' routes. The original network edges are retained and parallel enzyme names
#' are combined on one reaction edge.
#'
#' Glycan dimensions are measured before the tree is laid out. Nodes at the
#' same rank are separated by their rendered bounds, while rectangular edge
#' caps stop arrows outside the source and target glycans. The plot uses a
#' fixed-size panel so these clearances remain physical rather than changing
#' with the coordinate range.
#'
#' @param object A `glyenzy_biosynthesis_network` object returned by a
#'   biosynthesis function.
#' @param show_enzyme Logical. Whether to label reaction edges with enzyme or
#'   virtual-enzyme names.
#' @param size Positive numeric whole-cartoon scale multiplier passed to
#'   [glydraw::geom_node_glycan()]. Defaults to `0.22`.
#' @param node_gap Non-negative physical clearance, in inches, between glycan
#'   nodes at the same tree rank.
#' @param level_gap Non-negative physical clearance, in inches, between
#'   adjacent tree ranks. Increase this when enzyme labels are unusually long.
#' @param max_panel_width,max_panel_height Positive maximum panel dimensions,
#'   in inches. Networks larger than either limit are scaled proportionally so
#'   the complete plot fits a standard graphics device. Use `Inf` to preserve
#'   the network's natural size along one dimension.
#' @param show_linkage Logical. Whether glycan linkage annotations are shown.
#'   Defaults to `FALSE` for a compact network.
#' @param orient Glycan drawing orientation passed to
#'   [glydraw::geom_node_glycan()].
#' @param ... Additional glycan appearance arguments accepted by
#'   [glydraw::glycanGrob()], such as `node_size`, `colors`, or `style`.
#'
#' @returns A ggraph/ggplot object with a fixed-size, collision-aware tree
#'   panel.
#'
#' @examples
#' \dontrun{
#' network <- trace_biosynthesis_virtual(
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
#' )
#' ggplot2::autoplot(network)
#' }
#' @importFrom rlang .data
#' @exportS3Method ggplot2::autoplot
autoplot.glyenzy_biosynthesis_network <- function(
  object,
  show_enzyme = TRUE,
  size = 0.22,
  node_gap = 0.25,
  level_gap = 0.6,
  max_panel_width = 6,
  max_panel_height = 6,
  show_linkage = FALSE,
  orient = c("H", "V"),
  ...
) {
  rlang::check_installed(
    "ggraph",
    reason = "to plot glycan biosynthesis networks"
  )
  checkmate::assert_flag(show_enzyme)
  checkmate::assert_number(size, lower = 0, finite = TRUE)
  checkmate::assert_number(node_gap, lower = 0, finite = TRUE)
  checkmate::assert_number(level_gap, lower = 0, finite = TRUE)
  checkmate::assert_number(max_panel_width, lower = 0)
  checkmate::assert_number(max_panel_height, lower = 0)
  checkmate::assert_flag(show_linkage)
  if (size == 0) {
    cli::cli_abort("{.arg size} must be larger than 0.")
  }
  if (max_panel_width == 0 || max_panel_height == 0) {
    cli::cli_abort(
      "{.arg max_panel_width} and {.arg max_panel_height} must be larger than 0."
    )
  }
  orient <- match.arg(orient)
  dots <- rlang::list2(...)
  .validate_biosynthesis_glycan_dots(dots)
  .validate_biosynthesis_network(object)

  graph <- .collapse_biosynthesis_reactions(object)
  dimensions <- .biosynthesis_node_dimensions(
    graph,
    size = size,
    show_linkage = show_linkage,
    orient = orient,
    dots = dots
  )
  geometry <- .biosynthesis_plot_geometry(
    graph,
    dimensions,
    node_gap = node_gap,
    level_gap = level_gap,
    show_enzyme = show_enzyme
  )
  fit_scale <- .biosynthesis_fit_scale(
    geometry$bounds,
    max_width = max_panel_width,
    max_height = max_panel_height
  )
  if (fit_scale < 1) {
    geometry <- .biosynthesis_plot_geometry(
      graph,
      dimensions,
      node_gap = node_gap,
      level_gap = level_gap,
      show_enzyme = show_enzyme,
      scale = fit_scale
    )
  }
  layout <- geometry$layout
  labels <- geometry$labels
  bounds <- geometry$bounds
  plot <- ggraph::ggraph(layout)

  if (igraph::ecount(graph) > 0L) {
    plot <- plot +
      .biosynthesis_edge_layer(fit_scale) +
      ggraph::scale_edge_colour_manual(
        values = c(
          Enzyme = "#536779",
          `Virtual enzyme` = "#C46B2D"
        )
      ) +
      ggraph::scale_edge_linetype_manual(
        values = c(
          Enzyme = "solid",
          `Virtual enzyme` = "22"
        )
      )
  }

  if (!is.null(labels) && nrow(labels) > 0L) {
    plot <- plot +
      ggplot2::geom_label(
        data = labels,
        mapping = ggplot2::aes(
          x = .data$x,
          y = .data$y,
          label = .data$enzyme,
          colour = .data$.reaction_type
        ),
        fill = "#FFFFFFF2",
        linewidth = 0.25 * fit_scale,
        size = .biosynthesis_label_size * fit_scale,
        label.padding = grid::unit(0.12, "lines"),
        label.r = grid::unit(0.08, "lines"),
        lineheight = 0.9,
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = c(
          Enzyme = "#536779",
          `Virtual enzyme` = "#C46B2D"
        ),
        guide = "none"
      )
  }

  plot <- plot +
    rlang::exec(
      glydraw::geom_node_glycan,
      mapping = ggplot2::aes(structure = .data$name),
      size = size * fit_scale,
      show_linkage = show_linkage,
      orient = orient,
      !!!dots
    ) +
    ggplot2::coord_fixed(
      xlim = bounds$x,
      ylim = bounds$y,
      expand = FALSE,
      clip = "off"
    ) +
    ggraph::theme_graph(base_family = "") +
    ggplot2::theme(
      panel.widths = grid::unit(diff(bounds$x), "in"),
      panel.heights = grid::unit(diff(bounds$y), "in"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )

  attr(plot, "glyenzy_panel_size_in") <- c(
    width = diff(bounds$x),
    height = diff(bounds$y)
  )
  attr(plot, "glyenzy_layout_scale") <- fit_scale
  plot
}

.validate_biosynthesis_glycan_dots <- function(dots) {
  if (length(dots) == 0L) {
    return(invisible(dots))
  }
  if (is.null(names(dots)) || any(names(dots) == "")) {
    cli::cli_abort(
      "Glycan appearance arguments in {.arg ...} must be named."
    )
  }

  allowed <- setdiff(
    names(formals(glydraw::glycanGrob)),
    c("structure", "...")
  )
  unknown <- setdiff(names(dots), allowed)
  if (length(unknown) > 0L) {
    cli::cli_abort(c(
      "Unsupported glycan appearance argument{?s}: {.arg {unknown}}.",
      "i" = paste(
        "Use the dedicated {.arg size}, {.arg show_linkage}, and",
        "{.arg orient} arguments for network layout."
      )
    ))
  }
  invisible(dots)
}

.validate_biosynthesis_network <- function(graph) {
  checkmate::assert_class(graph, "igraph")
  if (!igraph::is_directed(graph)) {
    cli::cli_abort("A biosynthesis network must be directed.")
  }
  if (!igraph::is_dag(graph)) {
    cli::cli_abort("A biosynthesis network must be acyclic.")
  }
  if (!igraph::is_connected(graph, mode = "weak")) {
    cli::cli_abort("A biosynthesis network must be weakly connected.")
  }
  if (is.null(igraph::vertex_attr(graph, "name"))) {
    cli::cli_abort("Biosynthesis network nodes must have a {.field name}.")
  }
  if (
    igraph::ecount(graph) > 0L &&
      is.null(igraph::edge_attr(graph, "enzyme"))
  ) {
    cli::cli_abort(
      "Biosynthesis network edges must have an {.field enzyme} attribute."
    )
  }

  roots <- which(igraph::degree(graph, mode = "in") == 0L)
  if (length(roots) != 1L) {
    cli::cli_abort("A biosynthesis network must have exactly one root node.")
  }
  invisible(graph)
}

.collapse_biosynthesis_reactions <- function(graph) {
  vertices <- igraph::as_data_frame(graph, what = "vertices")
  edges <- igraph::as_data_frame(graph, what = "edges")
  if (nrow(edges) == 0L) {
    empty_edges <- tibble::tibble(
      from = character(),
      to = character(),
      enzyme = character(),
      .reaction_type = character()
    )
    return(igraph::graph_from_data_frame(
      empty_edges,
      directed = TRUE,
      vertices = vertices
    ))
  }

  virtual <- if (
    inherits(
      graph,
      "glyenzy_virtual_biosynthesis_network"
    )
  ) {
    rep(TRUE, nrow(edges))
  } else {
    igraph::edge_attr(graph, "is_virtual")
  }
  if (is.null(virtual)) {
    virtual <- rep(FALSE, nrow(edges))
  }
  virtual[is.na(virtual)] <- FALSE

  keys <- paste(edges$from, edges$to, sep = "\r")
  groups <- split(
    seq_along(keys),
    factor(keys, levels = unique(keys))
  )
  collapsed <- purrr::map_dfr(groups, function(indices) {
    enzymes <- unique(as.character(edges$enzyme[indices]))
    reaction_type <- if (any(virtual[indices])) {
      "Virtual enzyme"
    } else {
      "Enzyme"
    }
    tibble::tibble(
      from = edges$from[indices[[1]]],
      to = edges$to[indices[[1]]],
      enzyme = paste(enzymes, collapse = " / "),
      .reaction_type = reaction_type
    )
  })
  collapsed$.reaction_type <- factor(
    collapsed$.reaction_type,
    levels = c("Enzyme", "Virtual enzyme")
  )

  igraph::graph_from_data_frame(
    collapsed,
    directed = TRUE,
    vertices = vertices
  )
}

.biosynthesis_node_dimensions <- function(
  graph,
  size,
  show_linkage,
  orient,
  dots
) {
  opened_device <- grDevices::dev.cur() == 1L
  if (opened_device) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  structures <- igraph::V(graph)$name
  unique_structures <- unique(structures)
  unique_dimensions <- lapply(unique_structures, function(structure) {
    grob <- rlang::exec(
      glydraw::glycanGrob,
      structure = structure,
      show_linkage = show_linkage,
      orient = orient,
      !!!dots
    )
    grob$glydraw_border_px <- 0
    grob$glydraw_background <- FALSE
    content <- grid::makeContent(grob)
    child <- content$children[[1]]
    c(
      width = grid::convertWidth(
        grid::grobWidth(child),
        "in",
        valueOnly = TRUE
      ),
      height = grid::convertHeight(
        grid::grobHeight(child),
        "in",
        valueOnly = TRUE
      )
    ) *
      size
  })
  unique_dimensions <- do.call(rbind, unique_dimensions)
  dimensions <- unique_dimensions[
    match(structures, unique_structures),
    ,
    drop = FALSE
  ]

  tibble::tibble(
    name = structures,
    width = unname(dimensions[, "width"]),
    height = unname(dimensions[, "height"])
  )
}

.biosynthesis_tree_layout <- function(
  graph,
  dimensions,
  node_gap,
  level_gap,
  edge_padding
) {
  root <- which(igraph::degree(graph, mode = "in") == 0L)
  search <- igraph::bfs(
    graph,
    root = root,
    mode = "out",
    unreachable = FALSE,
    parent = TRUE,
    dist = TRUE
  )
  if (any(!is.finite(search$dist))) {
    cli::cli_abort("Every biosynthesis node must be reachable from the root.")
  }

  edges <- igraph::as_edgelist(graph, names = FALSE)
  if (
    nrow(edges) > 0L &&
      any(search$dist[edges[, 2]] - search$dist[edges[, 1]] != 1)
  ) {
    cli::cli_abort(
      "Biosynthesis reactions must connect adjacent tree ranks."
    )
  }

  tree <- igraph::make_empty_graph(
    n = igraph::vcount(graph),
    directed = TRUE
  )
  igraph::V(tree)$name <- igraph::V(graph)$name
  parents <- as.integer(search$parent)
  children <- which(!is.na(parents))
  if (length(children) > 0L) {
    tree <- igraph::add_edges(
      tree,
      as.vector(rbind(parents[children], children))
    )
  }
  tree_layout <- ggraph::create_layout(
    tree,
    layout = "tree",
    root = root
  )
  tree_order <- match(
    igraph::V(graph)$name,
    tree_layout$name
  )
  raw_x <- tree_layout$x[tree_order]

  same_rank_spacing <- unlist(lapply(
    split(raw_x, search$dist),
    function(x) diff(sort(unique(x)))
  ))
  same_rank_spacing <- same_rank_spacing[same_rank_spacing > 0]
  min_tree_spacing <- if (length(same_rank_spacing) == 0L) {
    1
  } else {
    min(same_rank_spacing)
  }
  x_spacing <- max(dimensions$width) + node_gap
  y_spacing <- max(dimensions$height) + level_gap
  x <- raw_x / min_tree_spacing * x_spacing
  x <- x - mean(range(x))
  y <- -search$dist * y_spacing
  y <- y - mean(range(y))

  layout <- ggraph::create_layout(
    graph,
    layout = "manual",
    x = x,
    y = y
  )
  dimension_order <- match(layout$name, dimensions$name)
  layout$.node_width <- dimensions$width[dimension_order]
  layout$.node_height <- dimensions$height[dimension_order]
  layout$.node_cap <- ggraph::rectangle(
    width = layout$.node_width + 2 * edge_padding,
    height = layout$.node_height + 2 * edge_padding,
    width_unit = "in",
    height_unit = "in"
  )
  layout
}

.biosynthesis_plot_geometry <- function(
  graph,
  dimensions,
  node_gap,
  level_gap,
  show_enzyme,
  scale = 1
) {
  scaled_dimensions <- dimensions
  scaled_dimensions$width <- scaled_dimensions$width * scale
  scaled_dimensions$height <- scaled_dimensions$height * scale
  layout <- .biosynthesis_tree_layout(
    graph,
    scaled_dimensions,
    node_gap = node_gap * scale,
    level_gap = level_gap * scale,
    edge_padding = .biosynthesis_edge_padding * scale
  )
  labels <- if (show_enzyme) {
    .biosynthesis_edge_labels(graph, layout, scale = scale)
  } else {
    NULL
  }
  bounds <- .biosynthesis_plot_bounds(
    layout,
    labels,
    padding = .biosynthesis_plot_padding * scale
  )

  list(layout = layout, labels = labels, bounds = bounds)
}

.biosynthesis_fit_scale <- function(bounds, max_width, max_height) {
  min(
    1,
    max_width / diff(bounds$x),
    max_height / diff(bounds$y)
  )
}

.biosynthesis_edge_labels <- function(graph, layout, scale = 1) {
  edges <- igraph::as_data_frame(graph, what = "edges")
  from <- match(edges$from, layout$name)
  to <- match(edges$to, layout$name)
  label_dimensions <- .biosynthesis_label_dimensions(
    edges$enzyme,
    scale = scale
  )
  labels <- tibble::tibble(
    x = (layout$x[from] + layout$x[to]) / 2,
    y = (layout$y[from] + layout$y[to]) / 2,
    enzyme = edges$enzyme,
    .reaction_type = edges$.reaction_type,
    .label_width = label_dimensions[, "width"],
    .label_height = label_dimensions[, "height"],
    .rank = pmin(layout$y[from], layout$y[to])
  )

  rank_groups <- split(
    seq_len(nrow(labels)),
    factor(labels$.rank, levels = unique(labels$.rank))
  )
  for (indices in rank_groups) {
    order_in_rank <- indices[order(labels$x[indices])]
    packed <- labels$x[order_in_rank]
    if (length(packed) > 1L) {
      for (i in 2:length(packed)) {
        previous <- order_in_rank[[i - 1L]]
        current <- order_in_rank[[i]]
        separation <- (labels$.label_width[previous] +
          labels$.label_width[current]) /
          2 +
          .biosynthesis_label_gap * scale
        packed[[i]] <- max(
          packed[[i]],
          packed[[i - 1L]] + separation
        )
      }
      packed <- packed - mean(packed - labels$x[order_in_rank])
    }
    labels$x[order_in_rank] <- packed
  }
  labels
}

.biosynthesis_label_dimensions <- function(labels, scale) {
  if (length(labels) == 0L) {
    return(matrix(
      numeric(),
      nrow = 0L,
      ncol = 2L,
      dimnames = list(NULL, c("width", "height"))
    ))
  }

  opened_device <- grDevices::dev.cur() == 1L
  if (opened_device) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  dimensions <- lapply(labels, function(label) {
    text <- grid::textGrob(
      label,
      gp = grid::gpar(
        fontsize = .biosynthesis_label_size *
          .biosynthesis_points_per_mm *
          scale,
        lineheight = 0.9
      )
    )
    c(
      width = grid::convertWidth(
        grid::grobWidth(text),
        "in",
        valueOnly = TRUE
      ),
      height = grid::convertHeight(
        grid::grobHeight(text),
        "in",
        valueOnly = TRUE
      )
    ) +
      .biosynthesis_label_padding * scale
  })
  do.call(rbind, dimensions)
}

.biosynthesis_plot_bounds <- function(
  layout,
  labels = NULL,
  padding = .biosynthesis_plot_padding
) {
  x_values <- c(
    layout$x - layout$.node_width / 2,
    layout$x + layout$.node_width / 2
  )
  y_values <- c(
    layout$y - layout$.node_height / 2,
    layout$y + layout$.node_height / 2
  )
  if (!is.null(labels) && nrow(labels) > 0L) {
    x_values <- c(
      x_values,
      labels$x - labels$.label_width / 2,
      labels$x + labels$.label_width / 2
    )
    y_values <- c(
      y_values,
      labels$y - labels$.label_height / 2,
      labels$y + labels$.label_height / 2
    )
  }

  list(
    x = c(
      min(x_values) - padding,
      max(x_values) + padding
    ),
    y = c(
      min(y_values) - padding,
      max(y_values) + padding
    )
  )
}

.biosynthesis_edge_layer <- function(scale = 1) {
  ggraph::geom_edge_link(
    mapping = ggplot2::aes(
      start_cap = .data$node1..node_cap,
      end_cap = .data$node2..node_cap,
      edge_colour = .data$.reaction_type,
      edge_linetype = .data$.reaction_type
    ),
    data = ggraph::get_edges("short", "all"),
    arrow = grid::arrow(
      angle = 25,
      length = grid::unit(2.4 * scale, "mm"),
      type = "closed"
    ),
    linewidth = 0.65 * scale,
    alpha = 0.9,
    lineend = "round",
    show.legend = TRUE
  )
}
