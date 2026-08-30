.biosynthesis_edge_padding <- 0.04
.biosynthesis_plot_padding <- 0.15
.biosynthesis_label_size <- 3
.biosynthesis_label_padding <- 0.06
.biosynthesis_label_gap <- 0.08
.biosynthesis_points_per_mm <- 72.27 / 25.4
.biosynthesis_min_label_size_pt <- 4
.biosynthesis_default_edge_color <- "#4D4D4D"
.biosynthesis_residues_by_color <- strsplit(
  c(
    "#0072BC" = "Glc GlcNAc GlcN GlcA Qui QuiNAc Oli Bac Api",
    "#00A651" = "Man ManNAc ManN ManA Rha RhaNAc Tyv Ara Kdn Pse LDmanHep Fru",
    "#FFD400" = "Gal GalNAc GalN GalA Lyx Leg Kdo Tag",
    "#F47920" = "Gul GulNAc GulN GulA 6dGul Abe Xyl Dha Sor",
    "#F69EA1" = "Alt AltNAc AltN AltA 6dAlt 6dAltNAc Par Rib Aci DDmanHep Psi",
    "#A54399" = "All AllNAc AllN AllA Dig Neu5Ac MurNAc",
    "#8FCCE9" = "Tal TalNAc TalN TalA 6dTal 6dTalNAc Col Neu5Gc 4eLeg MurNGc",
    "#A17A4D" = "Ido IdoNAc IdoN IdoA Neu Mur",
    "#ED1C24" = "Fuc FucNAc Sia"
  ),
  " ",
  fixed = TRUE
)

#' Plot a glycan biosynthesis network
#'
#' `autoplot.glyenzy_biosynthesis_network()` draws the typed networks returned
#' by [trace_biosynthesis()], [trace_biosynthesis_virtual()],
#' [path_biosynthesis()], and [path_biosynthesis_virtual()]. A layered DAG
#' layout determines the node positions while accounting for converging
#' biosynthesis routes. Candidate enzyme names stored on each reaction edge
#' are combined into one label.
#'
#' Glycan dimensions are measured before the network is laid out. Nodes at the
#' same rank are separated by their rendered bounds, while rectangular edge
#' caps stop arrows outside the source and target glycans. The plot uses a
#' fixed-size panel inside a requested figure canvas so these clearances remain
#' physical rather than changing with the coordinate range.
#'
#' Concrete enzyme reactions use solid edge shafts, while virtual enzyme
#' reactions use dashed shafts; arrowheads remain solid in both cases. This
#' convention is documented rather than repeated in an in-plot legend.
#'
#' @param object A `glyenzy_biosynthesis_network` object returned by a
#'   biosynthesis function.
#' @param show_enzyme Logical. Whether to label reaction edges with enzyme or
#'   virtual-enzyme names. Labels that would be too small to render after panel
#'   fitting are hidden with a warning.
#' @param size Positive numeric whole-cartoon scale multiplier passed to
#'   [glydraw::geom_node_glycan()]. Defaults to `0.4`.
#' @param node_gap Non-negative physical clearance, in inches, between glycan
#'   nodes at the same rank.
#' @param level_gap Non-negative physical clearance, in inches, between
#'   adjacent ranks. Increase this when enzyme labels are unusually long.
#' @param width,height Positive, finite figure dimensions. Networks larger than
#'   the available figure are scaled proportionally, while smaller networks
#'   retain their natural size and are centered in the figure.
#' @param units Units for `width` and `height`. One of `"in"` (the default),
#'   `"cm"`, or `"mm"`.
#' @param show_linkage Logical. Whether glycan linkage annotations are shown.
#'   Defaults to `FALSE` for a compact network.
#' @param orient Glycan drawing orientation passed to
#'   [glydraw::geom_node_glycan()]. One of `"left"`, `"right"`, `"up"`, or
#'   `"down"`; the default is `"left"`.
#' @param color_edge Logical. Whether reaction edges and enzyme labels use the
#'   SNFG color of the residue added by each reaction. Defaults to `FALSE`,
#'   which draws them in dark grey. Concrete and virtual reactions use solid
#'   and dashed lines, respectively, in both modes. Reactions without one
#'   unambiguous added residue remain dark grey.
#' @param enzyme_label_style How candidate enzymes on one reaction edge are
#'   labeled. `"condensed"` (the default) groups names with the same prefix and
#'   terminal number, for example, `"B4GALT1/2/3, B3GALT3/4"`. `"full"` keeps
#'   complete names separated by `" / "`.
#' @param highlight_target Logical or `NULL`. When `TRUE`, target glycans are
#'   drawn normally and all other glycans are semi-transparent. The default,
#'   `NULL`, highlights targets in multi-target networks but not in
#'   single-target networks. Reaction edges and enzyme labels are unchanged.
#'   This argument cannot be `TRUE` for networks returned by
#'   [path_biosynthesis()] or [path_biosynthesis_virtual()].
#' @param ... Additional glycan appearance arguments accepted by
#'   [glydraw::glycanGrob()], such as `node_size`, `colors`, or `style`.
#'
#' @returns A ggraph/ggplot object with a fixed-size, collision-aware layered
#'   panel centered in the requested figure dimensions.
#'
#' @section Figure size:
#' `width`, `height`, and `units` describe the complete base figure returned by
#' this method. Use the same dimensions with [ggplot2::ggsave()] or the
#' corresponding knitr/Quarto figure options. Titles, legends, or margins added
#' after this method returns may require additional output space.
#'
#' @examples
#' \dontrun{
#' network <- trace_biosynthesis_virtual(
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
#' )
#' network_plot <- ggplot2::autoplot(network, width = 8, height = 5)
#' ggplot2::ggsave(
#'   "biosynthesis-network.png",
#'   network_plot,
#'   width = 8,
#'   height = 5,
#'   units = "in"
#' )
#' }
#' @importFrom rlang .data
#' @importFrom ggplot2 autoplot
#' @export
autoplot.glyenzy_biosynthesis_network <- function(
  object,
  show_enzyme = TRUE,
  size = 0.4,
  node_gap = 0.25,
  level_gap = 0.6,
  width = 6,
  height = 6,
  units = c("in", "cm", "mm"),
  show_linkage = FALSE,
  orient = c("left", "right", "up", "down"),
  color_edge = FALSE,
  enzyme_label_style = c("condensed", "full"),
  highlight_target = NULL,
  ...
) {
  rlang::check_installed(
    "ggraph",
    reason = "to plot glycan biosynthesis networks"
  )
  checkmate::assert_flag(show_enzyme)
  checkmate::assert_flag(color_edge)
  checkmate::assert_flag(highlight_target, null.ok = TRUE)
  checkmate::assert_number(size, lower = 0, finite = TRUE)
  checkmate::assert_number(node_gap, lower = 0, finite = TRUE)
  checkmate::assert_number(level_gap, lower = 0, finite = TRUE)
  checkmate::assert_number(width, lower = 0, finite = TRUE)
  checkmate::assert_number(height, lower = 0, finite = TRUE)
  checkmate::assert_flag(show_linkage)
  if (size == 0) {
    cli::cli_abort("{.arg size} must be larger than 0.")
  }
  if (width == 0 || height == 0) {
    cli::cli_abort("{.arg width} and {.arg height} must be larger than 0.")
  }
  units <- match.arg(units)
  orient <- match.arg(orient)
  enzyme_label_style <- match.arg(enzyme_label_style)
  dots <- rlang::list2(...)
  .validate_biosynthesis_glycan_dots(dots)
  .validate_biosynthesis_network(object)
  highlight_target <- .resolve_biosynthesis_highlight_target(
    object,
    highlight_target
  )
  figure_size <- .biosynthesis_figure_size_in(width, height, units)

  graph <- .collapse_biosynthesis_reactions(
    object,
    enzyme_label_style = enzyme_label_style
  )
  if (color_edge) {
    graph <- .color_biosynthesis_reactions(graph)
  }
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
    max_width = figure_size[["width"]],
    max_height = figure_size[["height"]]
  )
  label_size_pt <- .biosynthesis_label_size *
    .biosynthesis_points_per_mm *
    fit_scale
  if (
    show_enzyme &&
      !is.null(geometry$labels) &&
      nrow(geometry$labels) > 0L &&
      label_size_pt < .biosynthesis_min_label_size_pt
  ) {
    cli::cli_warn(c(
      "Enzyme labels cannot be rendered at the fitted plot size and have been hidden.",
      "i" = paste(
        "Increase {.arg width} or {.arg height} to display enzyme",
        "labels."
      )
    ))
    show_enzyme <- FALSE
    geometry <- .biosynthesis_plot_geometry(
      graph,
      dimensions,
      node_gap = node_gap,
      level_gap = level_gap,
      show_enzyme = show_enzyme
    )
    fit_scale <- .biosynthesis_fit_scale(
      geometry$bounds,
      max_width = figure_size[["width"]],
      max_height = figure_size[["height"]]
    )
  }
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
  if (highlight_target) {
    layout$.structure <- glyparse::auto_parse(layout$name)
  }
  labels <- geometry$labels
  bounds <- geometry$bounds
  panel_size <- c(
    width = diff(bounds$x),
    height = diff(bounds$y)
  )
  figure_padding <- pmax(figure_size - panel_size, 0)
  plot <- ggraph::ggraph(layout)

  if (igraph::ecount(graph) > 0L) {
    plot <- plot +
      .biosynthesis_edge_layer(fit_scale, color_edge = color_edge) +
      .biosynthesis_arrow_layer(fit_scale, color_edge = color_edge) +
      ggraph::scale_edge_linetype_manual(
        values = c(
          Enzyme = "solid",
          `Virtual enzyme` = "22"
        ),
        guide = "none"
      )
    if (color_edge) {
      plot <- plot + ggraph::scale_edge_colour_identity()
    }
  }

  if (!is.null(labels) && nrow(labels) > 0L) {
    plot <- plot +
      .biosynthesis_label_layer(
        labels,
        scale = fit_scale,
        color_edge = color_edge
      )
    if (color_edge) {
      plot <- plot + ggplot2::scale_colour_identity()
    }
  }

  glycan_args <- c(
    list(
      size = size * fit_scale,
      show_linkage = show_linkage,
      orient = orient
    ),
    dots
  )
  if (highlight_target) {
    dimmed_args <- glycan_args
    dimmed_args$highlight <- integer()
    target_grob <- .biosynthesis_target_grob(
      layout,
      size = size * fit_scale,
      show_linkage = show_linkage,
      orient = orient,
      dots = dots,
      bounds = bounds
    )
    plot <- plot +
      rlang::exec(
        glydraw::geom_node_glycan,
        mapping = ggplot2::aes(
          structure = .data$.structure,
          filter = !.data$target
        ),
        !!!dimmed_args
      ) +
      ggplot2::annotation_custom(
        target_grob,
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf
      )
  } else {
    plot <- plot +
      rlang::exec(
        glydraw::geom_node_glycan,
        mapping = ggplot2::aes(structure = .data$name),
        !!!glycan_args
      )
  }

  plot <- plot +
    ggplot2::coord_fixed(
      xlim = bounds$x,
      ylim = bounds$y,
      expand = FALSE,
      clip = "off"
    ) +
    ggraph::theme_graph(base_family = "") +
    ggplot2::theme(
      panel.widths = grid::unit(panel_size[["width"]], "in"),
      panel.heights = grid::unit(panel_size[["height"]], "in"),
      plot.margin = ggplot2::margin(
        t = figure_padding[["height"]] / 2,
        r = figure_padding[["width"]] / 2,
        b = figure_padding[["height"]] / 2,
        l = figure_padding[["width"]] / 2,
        unit = "in"
      )
    )

  attr(plot, "glyenzy_panel_size_in") <- panel_size
  attr(plot, "glyenzy_layout_scale") <- fit_scale
  plot
}

.biosynthesis_target_grob <- function(
  layout,
  size,
  show_linkage,
  orient,
  dots,
  bounds
) {
  target_rows <- which(layout$target)
  dots$highlight <- NULL
  grobs <- lapply(seq_along(target_rows), function(index) {
    row <- target_rows[[index]]
    grob <- rlang::exec(
      glydraw::glycanGrob,
      structure = layout$.structure[[row]],
      show_linkage = show_linkage,
      orient = orient,
      !!!dots
    )
    grob$name <- paste0("target_glycan.", index)
    grob$glydraw_scale <- size
    grob$glydraw_hjust <- 0.5
    grob$glydraw_vjust <- 0.5
    grob$glydraw_angle <- 0
    grob$glydraw_border_px <- 0
    grob$glydraw_background <- FALSE
    grob$vp <- grid::viewport(
      x = grid::unit(
        (layout$x[[row]] - bounds$x[[1]]) / diff(bounds$x),
        "npc"
      ),
      y = grid::unit(
        (layout$y[[row]] - bounds$y[[1]]) / diff(bounds$y),
        "npc"
      ),
      just = c("center", "center")
    )
    grob
  })
  grid::grobTree(
    children = rlang::exec(grid::gList, !!!grobs),
    name = "biosynthesis_target_glycans"
  )
}

.resolve_biosynthesis_highlight_target <- function(
  graph,
  highlight_target
) {
  targets <- igraph::vertex_attr(graph, "target")
  if (is.null(targets)) {
    if (identical(highlight_target, TRUE)) {
      cli::cli_abort(
        "The network does not have a {.field target} vertex attribute.",
        call = rlang::caller_env()
      )
    }
    return(FALSE)
  }
  if (
    !is.logical(targets) ||
      length(targets) != igraph::vcount(graph) ||
      anyNA(targets)
  ) {
    cli::cli_abort(
      "The network {.field target} vertex attribute must contain one non-missing logical value per glycan.",
      call = rlang::caller_env()
    )
  }
  if (is.null(highlight_target)) {
    return(sum(targets) > 1L)
  }
  highlight_target
}

#' @importFrom graphics plot
#' @export
plot.glyenzy_biosynthesis_network <- function(x, ...) {
  print(autoplot(x, ...))
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
  removed_size_args <- intersect(
    names(dots),
    c("max_panel_width", "max_panel_height")
  )
  if (length(removed_size_args) > 0L) {
    cli::cli_abort(c(
      "The old panel-size argument{?s} {.arg {removed_size_args}} {?has/have} been removed.",
      "i" = paste(
        "Use {.arg width} and {.arg height} to set the complete",
        "figure size."
      )
    ))
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

.biosynthesis_figure_size_in <- function(width, height, units) {
  unit_scale <- switch(
    units,
    "in" = 1,
    "cm" = 1 / 2.54,
    "mm" = 1 / 25.4
  )
  c(width = width, height = height) * unit_scale
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

.collapse_biosynthesis_reactions <- function(
  graph,
  enzyme_label_style = c("condensed", "full")
) {
  enzyme_label_style <- match.arg(enzyme_label_style)
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

  virtual <- if ("is_virtual" %in% names(edges)) {
    igraph::edge_attr(graph, "is_virtual")
  } else if (
    inherits(
      graph,
      "glyenzy_virtual_biosynthesis_network"
    )
  ) {
    rep(TRUE, nrow(edges))
  } else {
    rep(FALSE, nrow(edges))
  }
  virtual[is.na(virtual)] <- FALSE

  candidates <- if ("enzymes" %in% names(edges)) {
    lapply(edges$enzymes, .valid_biosynthesis_enzyme_names)
  } else {
    rep(list(character()), nrow(edges))
  }

  keys <- paste(edges$from, edges$to, sep = "\r")
  groups <- split(
    seq_along(keys),
    factor(keys, levels = unique(keys))
  )
  collapsed <- purrr::map_dfr(groups, function(indices) {
    reaction_is_virtual <- all(virtual[indices])
    enzymes <- if (reaction_is_virtual) {
      unique(as.character(edges$enzyme[indices]))
    } else {
      unique(c(
        unlist(candidates[indices], use.names = FALSE),
        as.character(edges$enzyme[indices[!virtual[indices]]])
      ))
    }
    enzymes <- .valid_biosynthesis_enzyme_names(enzymes)
    reaction_type <- if (reaction_is_virtual) {
      "Virtual enzyme"
    } else {
      "Enzyme"
    }
    tibble::tibble(
      from = edges$from[indices[[1]]],
      to = edges$to[indices[[1]]],
      enzyme = .biosynthesis_enzyme_label(
        enzymes,
        style = enzyme_label_style
      ),
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

.biosynthesis_enzyme_label <- function(
  enzymes,
  style = c("condensed", "full")
) {
  style <- match.arg(style)
  enzymes <- unique(as.character(enzymes))
  if (style == "full") {
    return(paste(enzymes, collapse = " / "))
  }

  has_terminal_number <- grepl("[0-9]+$", enzymes)
  prefixes <- sub("[0-9]+$", "", enzymes)
  suffixes <- sub("^.*?([0-9]+)$", "\\1", enzymes, perl = TRUE)
  group_keys <- ifelse(
    has_terminal_number,
    paste0("numbered\r", prefixes),
    paste0("literal\r", seq_along(enzymes))
  )
  groups <- split(
    seq_along(enzymes),
    factor(group_keys, levels = unique(group_keys))
  )
  labels <- vapply(
    groups,
    function(indices) {
      first <- indices[[1]]
      if (has_terminal_number[[first]] && length(indices) > 1L) {
        paste0(
          prefixes[[first]],
          paste(suffixes[indices], collapse = "/")
        )
      } else {
        enzymes[[first]]
      }
    },
    character(1)
  )

  paste(labels, collapse = ", ")
}

.color_biosynthesis_reactions <- function(graph) {
  edges <- igraph::as_data_frame(graph, what = "edges")
  if (nrow(edges) == 0L) {
    return(graph)
  }

  structures <- igraph::V(graph)$name
  parsed <- glyparse::auto_parse(structures)
  nodes <- glyrepr::structure_nodes(parsed)
  residue_counts <- lapply(
    seq_along(structures),
    .biosynthesis_residue_counts,
    nodes = nodes
  )
  from <- match(edges$from, structures)
  to <- match(edges$to, structures)
  added_residues <- purrr::map2_chr(
    residue_counts[from],
    residue_counts[to],
    .biosynthesis_added_residue
  )
  colors <- purrr::map_chr(
    added_residues,
    .biosynthesis_residue_color
  )

  graph <- igraph::set_edge_attr(
    graph,
    ".added_residue",
    value = added_residues
  )
  igraph::set_edge_attr(
    graph,
    ".edge_colour",
    value = colors
  )
}

.biosynthesis_residue_counts <- function(glycan_id, nodes) {
  glycan_nodes <- nodes[nodes$glycan_id == glycan_id, , drop = FALSE]
  substituents <- unlist(
    strsplit(
      glycan_nodes$sub[nzchar(glycan_nodes$sub)],
      ",",
      fixed = TRUE
    ),
    use.names = FALSE
  )
  substituents <- sub("^[0-9?]+", "", substituents)

  list(
    monosaccharides = table(glycan_nodes$mono),
    substituents = table(substituents)
  )
}

.biosynthesis_added_residue <- function(from, to) {
  mono_delta <- .biosynthesis_count_delta(
    from$monosaccharides,
    to$monosaccharides
  )
  substituent_delta <- .biosynthesis_count_delta(
    from$substituents,
    to$substituents
  )
  if (any(c(mono_delta, substituent_delta) < 0L)) {
    return(NA_character_)
  }

  added_monos <- names(mono_delta)[mono_delta > 0L]
  if (
    length(added_monos) == 1L &&
      unname(mono_delta[added_monos]) == 1L
  ) {
    return(added_monos)
  }

  added_substituents <- names(substituent_delta)[substituent_delta > 0L]
  if (
    length(added_monos) == 0L &&
      length(added_substituents) == 1L &&
      unname(substituent_delta[added_substituents]) == 1L
  ) {
    return(added_substituents)
  }

  NA_character_
}

.biosynthesis_count_delta <- function(from, to) {
  residues <- union(names(from), names(to))
  if (length(residues) == 0L) {
    return(integer())
  }

  from_counts <- integer(length(residues))
  names(from_counts) <- residues
  to_counts <- from_counts
  from_counts[names(from)] <- as.integer(from)
  to_counts[names(to)] <- as.integer(to)
  to_counts - from_counts
}

.biosynthesis_residue_color <- function(residue) {
  if (is.na(residue)) {
    return(.biosynthesis_default_edge_color)
  }

  # Keep the palette aligned with every residue form supported by glyrepr,
  # including generic and furanose residues added after this package's palette.
  glyrepr_color <- get("get_mono_color", envir = asNamespace("glyrepr"))
  color <- glyrepr_color(residue)
  if (length(color) == 1L && !identical(color, "black")) {
    return(color)
  }

  # Keep this palette aligned with glyrepr's internal SNFG color helper.
  matches <- vapply(
    .biosynthesis_residues_by_color,
    function(residues) residue %in% residues,
    logical(1)
  )
  colors <- names(.biosynthesis_residues_by_color)[matches]
  if (length(colors) == 1L) colors else "black"
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

.biosynthesis_layered_layout <- function(
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
      "Biosynthesis reactions must connect adjacent ranks."
    )
  }

  raw_x <- igraph::layout_with_sugiyama(
    graph,
    layers = as.integer(search$dist) + 1L
  )$layout[, 1]

  same_rank_spacing <- unlist(lapply(
    split(raw_x, search$dist),
    function(x) diff(sort(unique(x)))
  ))
  same_rank_spacing <- same_rank_spacing[same_rank_spacing > 0]
  min_layout_spacing <- if (length(same_rank_spacing) == 0L) {
    1
  } else {
    min(same_rank_spacing)
  }
  x_spacing <- max(dimensions$width) + node_gap
  y_spacing <- max(dimensions$height) + level_gap
  x <- raw_x / min_layout_spacing * x_spacing
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
  layout <- .biosynthesis_layered_layout(
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
    .edge_colour = edges$.edge_colour,
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

.biosynthesis_edge_layer <- function(scale = 1, color_edge = FALSE) {
  mapping <- if (color_edge) {
    ggplot2::aes(
      start_cap = .data$node1..node_cap,
      end_cap = .data$node2..node_cap,
      edge_colour = .data$.edge_colour,
      edge_linetype = .data$.reaction_type
    )
  } else {
    ggplot2::aes(
      start_cap = .data$node1..node_cap,
      end_cap = .data$node2..node_cap,
      edge_linetype = .data$.reaction_type
    )
  }
  arguments <- list(
    mapping = mapping,
    data = ggraph::get_edges("short", "all"),
    linewidth = 0.65 * scale,
    alpha = 0.9,
    lineend = "round",
    show.legend = FALSE
  )
  if (!color_edge) {
    arguments$edge_colour <- .biosynthesis_default_edge_color
  }
  rlang::exec(ggraph::geom_edge_link, !!!arguments)
}

.biosynthesis_arrow_layer <- function(scale = 1, color_edge = FALSE) {
  mapping <- if (color_edge) {
    ggplot2::aes(
      start_cap = .data$node1..node_cap,
      end_cap = .data$node2..node_cap,
      edge_colour = .data$.edge_colour
    )
  } else {
    ggplot2::aes(
      start_cap = .data$node1..node_cap,
      end_cap = .data$node2..node_cap
    )
  }
  arguments <- list(
    mapping = mapping,
    data = ggraph::get_edges("short", "all"),
    arrow = grid::arrow(
      angle = 25,
      length = grid::unit(2.4 * scale, "mm"),
      type = "closed"
    ),
    linewidth = 0,
    edge_linetype = "solid",
    alpha = 0.9,
    lineend = "round",
    show.legend = FALSE
  )
  if (!color_edge) {
    arguments$edge_colour <- .biosynthesis_default_edge_color
  }
  rlang::exec(ggraph::geom_edge_link, !!!arguments)
}

.biosynthesis_label_layer <- function(labels, scale = 1, color_edge = FALSE) {
  mapping <- if (color_edge) {
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      label = .data$enzyme,
      colour = .data$.edge_colour
    )
  } else {
    ggplot2::aes(
      x = .data$x,
      y = .data$y,
      label = .data$enzyme
    )
  }
  arguments <- list(
    data = labels,
    mapping = mapping,
    fill = "#FFFFFFF2",
    linewidth = 0.25 * scale,
    size = .biosynthesis_label_size * scale,
    label.padding = grid::unit(0.12, "lines"),
    label.r = grid::unit(0.08, "lines"),
    lineheight = 0.9,
    show.legend = FALSE,
    inherit.aes = FALSE
  )
  if (!color_edge) {
    arguments$colour <- .biosynthesis_default_edge_color
  }
  rlang::exec(ggplot2::geom_label, !!!arguments)
}
