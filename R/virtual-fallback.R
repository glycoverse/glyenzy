# Target-directed virtual fallback for enzymatic searches -----------------

.perform_virtual_fallback_synthesis <- function(
  from_g,
  to_gs,
  enzymes,
  max_steps,
  filter,
  max_virtual_steps,
  structure_level,
  target_match,
  strict_result
) {
  from_key <- as.character(from_g)[[1]]
  strict_edges <- .offset_concrete_edges(strict_result$all_edges, 0L)
  found_keys <- setdiff(unique(strict_result$found_keys), from_key)
  strict_paths <- lapply(
    found_keys,
    function(found_key) {
      .build_virtual_fallback_result_graph(
        strict_edges,
        from_key,
        found_key,
        max_steps
      )
    }
  )

  to_keys <- as.character(to_gs)
  missing_idx <- which(to_keys %in% strict_result$missing_target_keys)
  fallback_paths <- lapply(
    missing_idx,
    function(i) {
      .search_target_with_virtual_fallback(
        from_g = from_g,
        to_g = to_gs[i],
        enzymes = enzymes,
        max_steps = max_steps,
        filter = filter,
        max_virtual_steps = max_virtual_steps,
        structure_level = structure_level,
        target_match = target_match,
        strict_result = strict_result
      )
    }
  )

  .combine_virtual_fallback_paths(
    c(strict_paths, fallback_paths),
    from_key
  )
}

.search_target_with_virtual_fallback <- function(
  from_g,
  to_g,
  enzymes,
  max_steps,
  filter,
  max_virtual_steps,
  structure_level,
  target_match,
  strict_result = NULL
) {
  from_key <- as.character(from_g)[[1]]
  to_key <- as.character(to_g)[[1]]
  if (is.null(strict_result)) {
    strict_result <- bfs_synthesis_search(
      from_g = from_g,
      to_gs = to_g,
      enzymes = enzymes,
      max_steps = max_steps,
      filter = filter,
      from_key = from_key,
      to_keys = to_key,
      structure_level = structure_level,
      target_match = target_match,
      allow_partial = TRUE
    )
  }

  if (length(strict_result$missing_target_keys) == 0L) {
    path <- build_synthesis_result_graph(
      strict_result,
      from_key,
      unique(strict_result$found_keys),
      max_steps
    )
    return(.mark_concrete_path(path))
  }

  candidates <- .prepare_virtual_fallback_candidates(to_g, enzymes)
  if (nrow(candidates) == 0L) {
    .abort_virtual_fallback(max_steps, max_virtual_steps)
  }

  anchors <- list(list(
    glycan = from_g,
    key = from_key,
    steps = 0L,
    result = strict_result
  ))
  all_edges <- list()

  for (virtual_steps_used in 0:max_virtual_steps) {
    segments <- vector("list", length(anchors))
    found_keys <- character()

    for (i in seq_along(anchors)) {
      anchor <- anchors[[i]]
      result <- anchor$result
      if (is.null(result)) {
        result <- bfs_synthesis_search(
          from_g = anchor$glycan,
          to_gs = to_g,
          enzymes = enzymes,
          max_steps = max_steps - anchor$steps,
          filter = filter,
          from_key = anchor$key,
          to_keys = to_key,
          structure_level = structure_level,
          target_match = target_match,
          allow_partial = TRUE
        )
      }

      concrete_edges <- .offset_concrete_edges(
        result$all_edges,
        anchor$steps
      )
      if (length(concrete_edges) > 0L) {
        all_edges <- c(all_edges, concrete_edges)
      }
      if (length(result$missing_target_keys) == 0L) {
        found_keys <- c(found_keys, result$found_keys)
      }

      segments[[i]] <- list(anchor = anchor, result = result)
    }

    if (length(found_keys) > 0L) {
      return(.build_virtual_fallback_result_graph(
        all_edges,
        from_key,
        unique(found_keys),
        max_steps
      ))
    }
    if (virtual_steps_used == max_virtual_steps) {
      break
    }

    next_anchors <- list()
    for (segment in segments) {
      reached_steps <- .concrete_segment_reached_steps(
        segment$anchor,
        segment$result
      )
      candidate_idx <- which(candidates$from %in% names(reached_steps))

      for (j in candidate_idx) {
        candidate <- candidates[j, , drop = FALSE]
        source_steps <- unname(reached_steps[[candidate$from]])
        product_steps <- source_steps + 1L
        if (product_steps > max_steps) {
          next
        }

        product <- glyparse::auto_parse(candidate$to)
        if (!.virtual_product_passes_filter(product, filter)) {
          next
        }

        all_edges[[length(all_edges) + 1L]] <- list(
          from = candidate$from,
          to = candidate$to,
          enzyme = candidate$enzyme,
          step = product_steps,
          is_virtual = TRUE
        )
        next_anchors[[length(next_anchors) + 1L]] <- list(
          glycan = product,
          key = candidate$to,
          steps = product_steps,
          result = NULL
        )
      }
    }

    anchors <- .deduplicate_virtual_anchors(next_anchors)
    if (length(anchors) == 0L) {
      break
    }
  }

  .abort_virtual_fallback(max_steps, max_virtual_steps)
}

.prepare_virtual_fallback_candidates <- function(to_g, enzymes) {
  from_g <- .decide_virtual_starting_glycan(to_g)
  result <- .virtual_trim_target(from_g, to_g)
  if (!result$found || length(result$all_edges) == 0L) {
    return(tibble::tibble(
      from = character(),
      to = character(),
      enzyme = character()
    ))
  }

  search_result <- list(
    found_keys = result$endpoint_key,
    all_edges = result$all_edges,
    parent = rlang::env(),
    parent_enzyme = rlang::env(),
    parent_step = rlang::env()
  )
  path <- build_synthesis_result_graph(
    search_result,
    as.character(from_g)[[1]],
    result$endpoint_key,
    .virtual_glycan_size(to_g) - .virtual_glycan_size(from_g)
  )
  path <- .amplify_virtual_edges(path, enzymes)
  edges <- igraph::as_data_frame(path, what = "edges")
  unsupported <- lengths(edges$concrete_enzymes) == 0L

  tibble::as_tibble(edges[unsupported, c("from", "to", "enzyme")])
}

.offset_concrete_edges <- function(edges, step_offset) {
  lapply(
    edges,
    function(edge) {
      list(
        from = edge$from,
        to = edge$to,
        enzyme = edge$enzyme,
        step = step_offset + edge$step,
        is_virtual = FALSE
      )
    }
  )
}

.concrete_segment_reached_steps <- function(anchor, result) {
  if (length(result$all_edges) == 0L) {
    return(stats::setNames(anchor$steps, anchor$key))
  }

  edges <- do.call(
    rbind,
    lapply(
      result$all_edges,
      \(edge) data.frame(from = edge$from, to = edge$to)
    )
  )
  vertices <- data.frame(name = unique(c(anchor$key, edges$from, edges$to)))
  graph <- igraph::graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = vertices
  )
  distances <- igraph::distances(
    graph,
    v = anchor$key,
    mode = "out",
    weights = NA
  )[1, ]
  distances <- distances[is.finite(distances)]
  stats::setNames(
    anchor$steps + as.integer(distances),
    names(distances)
  )
}

.virtual_product_passes_filter <- function(product, filter) {
  if (is.null(filter)) {
    return(TRUE)
  }

  keep <- filter(product)
  checkmate::assert_logical(
    keep,
    len = length(product),
    any.missing = FALSE
  )
  keep[[1]]
}

.deduplicate_virtual_anchors <- function(anchors) {
  if (length(anchors) < 2L) {
    return(anchors)
  }

  keys <- vapply(anchors, `[[`, character(1), "key")
  steps <- vapply(anchors, `[[`, integer(1), "steps")
  anchors <- anchors[order(steps)]
  keys <- keys[order(steps)]
  anchors[!duplicated(keys)]
}

.build_virtual_fallback_result_graph <- function(
  all_edges,
  from_key,
  endpoint_keys,
  max_steps
) {
  edges <- do.call(
    rbind,
    lapply(
      all_edges,
      function(edge) {
        tibble::tibble(
          from = edge$from,
          to = edge$to,
          enzyme = edge$enzyme,
          step = edge$step,
          is_virtual = edge$is_virtual
        )
      }
    )
  )
  edges <- edges[order(edges$step), ]
  edge_keys <- paste(
    edges$from,
    edges$to,
    edges$enzyme,
    edges$is_virtual,
    sep = "\r"
  )
  edges <- edges[!duplicated(edge_keys), ]

  graph <- igraph::graph_from_data_frame(
    edges[c("from", "to", "enzyme", "is_virtual")],
    directed = TRUE
  )
  endpoint_keys <- intersect(endpoint_keys, igraph::V(graph)$name)
  if (length(endpoint_keys) == 0L) {
    cli::cli_abort(
      "No synthesis path found within {.val {max_steps}} steps."
    )
  }

  virtual_penalty <- max_steps + 1L
  weights <- 1 + virtual_penalty * as.integer(igraph::E(graph)$is_virtual)
  from_id <- which(igraph::V(graph)$name == from_key)
  endpoint_ids <- match(endpoint_keys, igraph::V(graph)$name)
  from_distances <- igraph::distances(
    graph,
    v = from_id,
    mode = "out",
    weights = weights
  )[1, ]
  best_distance <- min(from_distances[endpoint_ids])
  best_endpoint_ids <- endpoint_ids[
    from_distances[endpoint_ids] == best_distance
  ]
  to_distances <- igraph::distances(
    graph,
    v = best_endpoint_ids,
    to = igraph::V(graph),
    mode = "in",
    weights = weights
  )
  if (is.null(dim(to_distances))) {
    to_distances <- matrix(to_distances, nrow = 1L)
  }
  to_distances <- apply(to_distances, 2L, min)

  edge_ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  keep <- from_distances[edge_ends[, 1L]] +
    weights +
    to_distances[edge_ends[, 2L]] ==
    best_distance
  path <- igraph::subgraph_from_edges(
    graph,
    igraph::E(graph)[keep],
    delete.vertices = TRUE
  )

  path_ends <- igraph::ends(path, igraph::E(path), names = TRUE)
  path_to_distances <- from_distances[
    match(path_ends[, 2L], igraph::V(graph)$name)
  ]
  virtual_counts <- path_to_distances %/% virtual_penalty
  igraph::set_edge_attr(
    path,
    "step",
    value = as.integer(
      path_to_distances - virtual_counts * virtual_penalty
    )
  )
}

.mark_concrete_path <- function(path) {
  igraph::set_edge_attr(
    path,
    "is_virtual",
    value = rep(FALSE, igraph::ecount(path))
  )
}

.combine_virtual_fallback_paths <- function(paths, from_key) {
  edges <- lapply(
    paths,
    \(path) igraph::as_data_frame(path, what = "edges")
  )
  edges <- edges[vapply(edges, nrow, integer(1)) > 0L]
  if (length(edges) == 0L) {
    return(
      igraph::make_empty_graph(n = 1L, directed = TRUE) |>
        igraph::set_vertex_attr("name", value = from_key)
    )
  }

  edges <- do.call(rbind, edges)
  edges <- edges[order(edges$step), ]
  edge_keys <- paste(
    edges$from,
    edges$to,
    edges$enzyme,
    edges$is_virtual,
    sep = "\r"
  )
  edges <- edges[!duplicated(edge_keys), ]
  igraph::graph_from_data_frame(edges, directed = TRUE)
}

.abort_virtual_fallback <- function(max_steps, max_virtual_steps) {
  cli::cli_abort(c(
    "No synthesis path found for 1 target(s) within {.val {max_steps}} steps.",
    "i" = "Virtual fallback was limited to {.val {max_virtual_steps}} step(s)."
  ))
}
