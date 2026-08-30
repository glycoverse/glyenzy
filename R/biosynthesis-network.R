.new_biosynthesis_network <- function(graph, virtual = FALSE) {
  checkmate::assert_class(graph, "igraph")
  checkmate::assert_flag(virtual)

  network_classes <- c(
    if (virtual) "glyenzy_virtual_biosynthesis_network",
    "glyenzy_biosynthesis_network"
  )
  class(graph) <- unique(c(network_classes, class(graph)))
  graph
}

.finalize_biosynthesis_network <- function(
  graph,
  targets,
  match = c("key", "whole"),
  virtual = FALSE
) {
  match <- match.arg(match)
  graph |>
    .set_biosynthesis_targets(targets, match = match) |>
    .collapse_biosynthesis_edges(virtual = virtual) |>
    .new_biosynthesis_network(virtual = virtual)
}

.collapse_biosynthesis_edges <- function(graph, virtual = FALSE) {
  checkmate::assert_class(graph, "igraph")
  checkmate::assert_flag(virtual)

  vertices <- igraph::as_data_frame(graph, what = "vertices")
  edges <- igraph::as_data_frame(graph, what = "edges")
  if (nrow(edges) == 0L) {
    graph <- igraph::graph_from_data_frame(
      tibble::tibble(
        from = character(),
        to = character(),
        enzyme = character(),
        is_virtual = logical(),
        step = integer()
      ),
      directed = TRUE,
      vertices = vertices
    )
    return(igraph::set_edge_attr(
      graph,
      "enzymes",
      value = list()
    ))
  }

  edge_is_virtual <- if (virtual) {
    rep(TRUE, nrow(edges))
  } else if ("is_virtual" %in% names(edges)) {
    as.logical(edges$is_virtual)
  } else {
    rep(FALSE, nrow(edges))
  }
  edge_is_virtual[is.na(edge_is_virtual)] <- FALSE

  edge_steps <- if ("step" %in% names(edges)) {
    suppressWarnings(as.integer(edges$step))
  } else {
    rep(NA_integer_, nrow(edges))
  }
  edge_labels <- if ("enzyme" %in% names(edges)) {
    as.character(edges$enzyme)
  } else {
    rep(NA_character_, nrow(edges))
  }
  edge_enzymes <- if ("enzymes" %in% names(edges)) {
    lapply(edges$enzymes, .valid_biosynthesis_enzyme_names)
  } else {
    rep(list(character()), nrow(edges))
  }
  concrete_enzymes <- if ("concrete_enzymes" %in% names(edges)) {
    lapply(edges$concrete_enzymes, .valid_biosynthesis_enzyme_names)
  } else {
    rep(list(character()), nrow(edges))
  }

  keys <- paste(edges$from, edges$to, sep = "\r")
  groups <- split(
    seq_len(nrow(edges)),
    factor(keys, levels = unique(keys))
  )
  collapsed <- lapply(groups, function(indices) {
    keep <- if (any(!edge_is_virtual[indices])) {
      indices[!edge_is_virtual[indices]]
    } else {
      indices
    }
    is_virtual <- all(edge_is_virtual[keep])
    enzymes <- unique(c(
      unlist(edge_enzymes[keep], use.names = FALSE),
      unlist(concrete_enzymes[keep], use.names = FALSE),
      if (!is_virtual) edge_labels[keep]
    ))
    enzymes <- .valid_biosynthesis_enzyme_names(enzymes)
    labels <- if (is_virtual) {
      .valid_biosynthesis_enzyme_names(edge_labels[keep])
    } else {
      enzymes
    }
    steps <- edge_steps[keep]
    steps <- steps[!is.na(steps)]

    list(
      from = edges$from[keep[[1]]],
      to = edges$to[keep[[1]]],
      enzyme = paste(labels, collapse = " / "),
      enzymes = enzymes,
      is_virtual = is_virtual,
      step = if (length(steps) == 0L) NA_integer_ else min(steps)
    )
  })

  edge_data <- tibble::tibble(
    from = vapply(collapsed, `[[`, character(1), "from"),
    to = vapply(collapsed, `[[`, character(1), "to"),
    enzyme = vapply(collapsed, `[[`, character(1), "enzyme"),
    is_virtual = vapply(collapsed, `[[`, logical(1), "is_virtual"),
    step = vapply(collapsed, `[[`, integer(1), "step")
  )
  graph <- igraph::graph_from_data_frame(
    edge_data,
    directed = TRUE,
    vertices = vertices
  )
  igraph::set_edge_attr(
    graph,
    "enzymes",
    value = lapply(collapsed, `[[`, "enzymes")
  )
}

.valid_biosynthesis_enzyme_names <- function(x) {
  x <- as.character(x)
  unique(x[!is.na(x) & nzchar(x)])
}

.set_biosynthesis_targets <- function(
  graph,
  targets,
  match = c("key", "whole")
) {
  match <- match.arg(match)
  vertex_names <- igraph::vertex_attr(graph, "name")

  if (identical(match, "key")) {
    is_target <- vertex_names %in% as.character(targets)
  } else {
    vertex_graphs <- glyrepr::get_structure_graphs(
      glyparse::auto_parse(vertex_names),
      return_list = TRUE
    )
    target_graphs <- glyrepr::get_structure_graphs(
      targets,
      return_list = TRUE
    )
    is_target <- purrr::map_lgl(
      vertex_graphs,
      function(vertex_graph) {
        purrr::some(
          target_graphs,
          function(target_graph) {
            glymotif::.g_have_motif(
              vertex_graph,
              target_graph,
              alignment = "whole",
              mode = "lenient"
            )
          }
        )
      }
    )
  }

  igraph::set_vertex_attr(graph, "target", value = is_target)
}
