# Backward biosynthesis search for virtual enzymes -------------------------

.decide_virtual_starting_glycan <- function(glycan) {
  graph <- glyrepr::get_structure_graphs(glycan)
  core <- .virtual_n_glycan_core(glycan)
  core_matches <- .virtual_core_matches(glycan, core)

  if (length(core_matches) > 0L) {
    start_graph <- igraph::induced_subgraph(graph, core_matches[[1]])
  } else {
    root <- which(igraph::degree(graph, mode = "in") == 0L)
    start_graph <- igraph::induced_subgraph(graph, root)
  }

  glyrepr::as_glycan_structure(start_graph)
}

.virtual_n_glycan_core <- function(glycan) {
  structure_level <- .glycan_structure_level(glycan)
  if (identical(structure_level, "basic")) {
    return(glyrepr::n_glycan_core(
      linkage = FALSE,
      mono_type = "generic"
    ))
  }
  if (identical(structure_level, "topological")) {
    return(glyrepr::n_glycan_core(linkage = FALSE))
  }
  glyrepr::n_glycan_core()
}

.virtual_core_matches <- function(glycan, core) {
  structure_level <- .glycan_structure_level(glycan)
  mode <- if (identical(structure_level, "partial")) {
    "lenient"
  } else {
    "strict"
  }

  tryCatch(
    glymotif::match_motif(
      glycan,
      core,
      alignment = "core",
      mode = mode
    )[[1]],
    error = function(e) list()
  )
}

.prepare_virtual_start <- function(from, to) {
  from_level <- .glycan_structure_level(from)
  to_level <- .glycan_structure_level(to)
  if (
    to_level %in%
      c("topological", "basic") &&
      .structure_level_rank(from_level) > .structure_level_rank(to_level)
  ) {
    return(.reduce_structure_level(from, to_level))
  }
  from
}

.perform_virtual_synthesis <- function(
  from_g,
  to_gs,
  max_steps,
  filter = NULL
) {
  from_key <- as.character(from_g)[[1]]
  endpoint_keys <- character(length(to_gs))
  all_edges <- list()
  found <- logical(length(to_gs))

  for (i in seq_along(to_gs)) {
    result <- .virtual_trim_target(
      from_g,
      to_gs[i],
      max_steps,
      filter
    )
    found[[i]] <- result$found
    endpoint_keys[[i]] <- result$endpoint_key
    if (length(result$all_edges) > 0L) {
      all_edges <- c(all_edges, result$all_edges)
    }
  }

  if (!all(found)) {
    cli::cli_abort(
      "No synthesis path found for {sum(!found)} target(s) within {.val {max_steps}} steps."
    )
  }

  search_result <- list(
    found_keys = endpoint_keys,
    all_edges = all_edges,
    parent = rlang::env(),
    parent_enzyme = rlang::env(),
    parent_step = rlang::env()
  )
  path <- build_synthesis_result_graph(
    search_result,
    from_key,
    unique(endpoint_keys),
    max_steps
  )

  reachable <- igraph::subcomponent(
    path,
    which(igraph::V(path)$name == from_key),
    mode = "out"
  )
  reachable_keys <- igraph::V(path)$name[reachable]
  if (!all(endpoint_keys %in% reachable_keys)) {
    cli::cli_abort(
      "No synthesis path found for {sum(!endpoint_keys %in% reachable_keys)} target(s) within {.val {max_steps}} steps."
    )
  }

  path
}

.virtual_trim_target <- function(from_g, to_g, max_steps, filter) {
  from_key <- as.character(from_g)[[1]]
  to_key <- as.character(to_g)[[1]]
  from_size <- .virtual_glycan_size(from_g)
  to_size <- .virtual_glycan_size(to_g)
  steps_needed <- to_size - from_size

  if (
    steps_needed < 0L ||
      steps_needed > max_steps ||
      !.virtual_contains_start(to_g, from_g)
  ) {
    return(list(
      found = FALSE,
      endpoint_key = NA_character_,
      all_edges = list()
    ))
  }
  if (steps_needed == 0L) {
    return(list(
      found = TRUE,
      endpoint_key = from_key,
      all_edges = list()
    ))
  }

  network_level <- .glycan_structure_level(to_g)
  queue <- list(to_g)
  queue_keys <- to_key
  visited <- rlang::env()
  rlang::env_poke(visited, to_key, TRUE)
  all_edges <- list()
  found <- FALSE

  for (depth in seq_len(steps_needed)) {
    candidates <- .virtual_trim_frontier(
      queue,
      queue_keys,
      from_g,
      from_size,
      network_level
    )
    if (length(candidates) == 0L) {
      break
    }

    candidates <- .filter_virtual_candidates(candidates, filter)
    if (length(candidates) == 0L) {
      break
    }

    new_queue <- list()
    new_queue_keys <- character()
    for (candidate in candidates) {
      all_edges[[length(all_edges) + 1L]] <- list(
        from = candidate$key,
        to = candidate$product_key,
        enzyme = candidate$enzyme,
        step = candidate$step
      )

      if (identical(candidate$key, from_key)) {
        found <- TRUE
      } else if (!rlang::env_has(visited, candidate$key)) {
        rlang::env_poke(visited, candidate$key, TRUE)
        new_queue[[length(new_queue) + 1L]] <- candidate$glycan
        new_queue_keys[[length(new_queue_keys) + 1L]] <- candidate$key
      }
    }

    queue <- new_queue
    queue_keys <- new_queue_keys
  }

  list(
    found = found,
    endpoint_key = if (found) to_key else NA_character_,
    all_edges = all_edges
  )
}

.virtual_trim_frontier <- function(
  queue,
  queue_keys,
  from_g,
  from_size,
  network_level
) {
  candidates <- list()
  from_key <- as.character(from_g)[[1]]

  for (i in seq_along(queue)) {
    current <- queue[[i]]
    current_key <- queue_keys[[i]]
    graph <- glyrepr::get_structure_graphs(current)
    root <- which(igraph::degree(graph, mode = "in") == 0L)
    leaves <- which(igraph::degree(graph, mode = "out") == 0L)
    leaves <- setdiff(leaves, root)

    for (leaf in leaves) {
      incoming <- igraph::incident(graph, leaf, mode = "in")
      enzyme <- .virtual_enzyme_name(
        mono = igraph::vertex_attr(graph, "mono", index = leaf),
        linkage = igraph::edge_attr(graph, "linkage", index = incoming),
        structure_level = network_level
      )
      precursor_graph <- igraph::delete_vertices(graph, leaf)
      precursor <- glyrepr::as_glycan_structure(precursor_graph)
      precursor_size <- .virtual_glycan_size(precursor)

      if (
        precursor_size < from_size ||
          !.virtual_contains_start(precursor, from_g)
      ) {
        next
      }

      precursor_key <- as.character(precursor)[[1]]
      if (precursor_size == from_size) {
        precursor <- from_g
        precursor_key <- from_key
      }

      candidates[[length(candidates) + 1L]] <- list(
        glycan = precursor,
        key = precursor_key,
        product_key = current_key,
        enzyme = enzyme,
        step = igraph::vcount(graph) - from_size
      )
    }
  }

  candidates
}

.filter_virtual_candidates <- function(candidates, filter) {
  if (is.null(filter) || length(candidates) == 0L) {
    return(candidates)
  }

  glycans <- do.call(c, purrr::map(candidates, "glycan"))
  keep <- filter(glycans)
  checkmate::assert_logical(
    keep,
    len = length(candidates),
    any.missing = FALSE
  )
  candidates[keep]
}

.virtual_contains_start <- function(glycan, from) {
  glycan_level <- .glycan_structure_level(glycan)
  from_level <- .glycan_structure_level(from)
  mode <- if (
    identical(glycan_level, "intact") &&
      identical(from_level, "intact")
  ) {
    "strict"
  } else {
    "lenient"
  }

  tryCatch(
    isTRUE(glymotif::have_motif(
      glycan,
      from,
      alignment = "core",
      mode = mode
    )),
    error = function(e) FALSE
  )
}

.virtual_enzyme_name <- function(mono, linkage, structure_level) {
  if (!identical(structure_level, "intact")) {
    return(paste0(mono, "T"))
  }

  anomer <- substr(linkage, 1L, 1L)
  acceptor_position <- sub(".*-", "", linkage)
  paste0(anomer, acceptor_position, mono, "T")
}

.virtual_glycan_size <- function(glycan) {
  igraph::vcount(glyrepr::get_structure_graphs(glycan))
}

.validate_virtual_enzymes <- function(enzymes) {
  if (!is.null(enzymes)) {
    cli::cli_abort(c(
      "{.arg enzymes} must be {.code NULL} when {.code method = \"virtual\"}.",
      "i" = "Virtual-enzyme tracing does not use known enzyme rules."
    ))
  }
  invisible(NULL)
}
