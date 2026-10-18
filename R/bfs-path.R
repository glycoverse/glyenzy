#' Breadth-First Search for Glycan Synthesis Paths
#'
#' Core BFS algorithms for finding synthesis paths between glycan structures.
#' These functions provide the algorithmic foundation for both path_biosynthesis()
#' and trace_biosynthesis().

#' Perform BFS search for synthesis paths between glycan structures
#'
#' This function instantiates the BFS R6 engine and returns its results.
#'
#' @param from_g Starting glycan structure (single glyrepr::glycan_structure)
#' @param to_gs Target glycan structures (glyrepr::glycan_structure vector)
#' @param enzymes List of `glyenzy_enzyme` objects to use
#' @param max_steps Maximum number of synthesis steps to explore
#' @param filter Optional function to filter glycan structures at each step
#' @param from_key Optional pre-computed string key for starting glycan
#' @param to_keys Optional pre-computed string keys for target glycans
#' @param structure_level Structure level used for generated products
#' @param target_match Target matching strategy
#' @param allow_partial Whether to return the explored search state when one or
#'   more targets are not reached.
#'
#' @returns List with search results including found paths and exploration data
#'
#' @examples
#' \dontrun{
#' from_g <- glyparse::parse_iupac_condensed("Gal(b1-4)GlcNAc(b1-")
#' to_gs <- c(
#'   glyparse::parse_iupac_condensed("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"),
#'   glyparse::parse_iupac_condensed("Fuc(a1-3)Gal(b1-4)GlcNAc(b1-")
#' )
#' result <- bfs_synthesis_search(from_g, to_gs, "ST6GAL1", 3)
#' }
#'
#' @noRd
bfs_synthesis_search <- function(
  from_g,
  to_gs,
  enzymes,
  max_steps,
  filter = NULL,
  from_key = NULL,
  to_keys = NULL,
  structure_level = "intact",
  target_match = c("key", "whole"),
  allow_partial = FALSE
) {
  target_match <- match.arg(target_match)
  checkmate::assert_flag(allow_partial)
  engine <- BfsSynthesisSearch$new(
    from_g = from_g,
    to_gs = to_gs,
    enzymes = enzymes,
    max_steps = max_steps,
    filter = filter,
    from_key = from_key,
    to_keys = to_keys,
    structure_level = structure_level,
    target_match = target_match,
    allow_partial = allow_partial
  )

  engine$run()
}

# Encode the canonical colored tree topology used by standard glycan actions.
# Vertex names are intentionally excluded because they are igraph identifiers
# rather than glycan semantics.
.bfs_graph_signature <- function(graph) {
  edges <- igraph::as_edgelist(graph, names = FALSE)
  incoming_linkages <- rep(NA_character_, igraph::vcount(graph))
  if (length(edges) > 0L) {
    incoming_linkages[edges[, 2]] <- igraph::edge_attr(graph, "linkage")
  }
  vertex_labels <- paste0(
    .bfs_signature_values(igraph::vertex_attr(graph, "mono")),
    .bfs_signature_values(igraph::vertex_attr(graph, "sub")),
    .bfs_signature_values(incoming_linkages)
  )
  label_levels <- sort(unique(vertex_labels))
  colors <- match(vertex_labels, label_levels)
  labeling <- igraph::canonical_permutation(
    graph,
    colors = colors
  )$labeling

  canonical_edges <- if (length(edges) == 0L) {
    edges
  } else {
    matrix(labeling[edges], ncol = 2L)
  }
  if (nrow(canonical_edges) > 1L) {
    canonical_edges <- canonical_edges[
      order(canonical_edges[, 1], canonical_edges[, 2]),
      ,
      drop = FALSE
    ]
  }
  edge_indices <- if (length(canonical_edges) == 0L) {
    ""
  } else {
    paste(as.integer(t(canonical_edges)), collapse = ",")
  }

  paste0(
    as.integer(igraph::is_directed(graph)),
    ";",
    igraph::vcount(graph),
    ";",
    igraph::ecount(graph),
    ";",
    edge_indices,
    ";",
    .bfs_signature_strings(label_levels),
    ";",
    paste(colors[order(labeling)], collapse = ",")
  )
}

.bfs_signature_strings <- function(values) {
  paste0(.bfs_signature_values(values), collapse = "")
}

.bfs_signature_values <- function(values) {
  values <- enc2utf8(as.character(values))
  missing <- is.na(values)
  values[missing] <- ""
  lengths <- nchar(values, type = "bytes")
  prefixes <- ifelse(missing, "N:", paste0(lengths, ":"))
  paste0(prefixes, values)
}

# Put the targets most likely to contain an intermediate first so scalar
# matching can stop early. Larger and more substituted structures dominate
# the common multi-target case where targets are related biosynthetic states.
.order_bfs_pruning_targets <- function(target_graphs) {
  if (length(target_graphs) < 2L) {
    return(target_graphs)
  }

  sizes <- vapply(target_graphs, igraph::vcount, numeric(1))
  substituent_counts <- vapply(
    target_graphs,
    function(target_graph) {
      sum(lengths(lapply(
        igraph::vertex_attr(target_graph, "sub"),
        .substituent_tokens
      )))
    },
    integer(1)
  )
  target_graphs[order(-sizes, -substituent_counts)]
}

#' BFS synthesis search as an R6 engine
#'
#' Stores search inputs and R extension callbacks for the native BFS engine.
#' The C++ engine owns expansion, pruning, deduplication and queue traversal.
#'
#' R6 fields expose the mutable state used across the BFS expansion steps. This
#' keeps the algorithm logic cohesive and avoids error-prone argument plumbing.
#'
#' @noRd
BfsSynthesisSearch <- R6::R6Class(
  "BfsSynthesisSearch",
  public = list(
    from_g = NULL,
    to_gs = NULL,
    enzymes = NULL,
    max_steps = NULL,
    filter = NULL,
    from_key = NULL,
    to_keys = NULL,
    structure_level = NULL,
    target_match = NULL,
    allow_partial = NULL,
    remaining_targets_map = NULL,
    queue = NULL,
    queue_graphs = NULL,
    queue_keys = NULL,
    visited = NULL,
    parent = NULL,
    parent_enzyme = NULL,
    parent_step = NULL,
    all_edges = NULL,
    found_keys_storage = NULL,
    found_tail = 0L,
    step = 0L,
    native_stats = NULL,

    initialize = function(
      from_g,
      to_gs,
      enzymes,
      max_steps,
      filter = NULL,
      from_key = NULL,
      to_keys = NULL,
      structure_level = "intact",
      target_match = "key",
      allow_partial = FALSE
    ) {
      self$from_g <- from_g
      self$to_gs <- to_gs
      self$enzymes <- enzymes
      self$max_steps <- max_steps
      self$filter <- filter
      self$structure_level <- .validate_structure_level(structure_level)
      self$target_match <- match.arg(target_match, c("key", "whole"))
      self$allow_partial <- allow_partial
      self$from_key <- if (is.null(from_key)) {
        as.character(from_g)[1]
      } else {
        from_key
      }
      self$to_keys <- if (is.null(to_keys)) as.character(to_gs) else to_keys

      self$remaining_targets_map <- fastmap::fastmap()
      for (target_key in self$to_keys) {
        self$remaining_targets_map$set(target_key, TRUE)
      }

      self$queue <- list(self$from_g)
      self$queue_keys <- c(self$from_key)
      self$visited <- rlang::env()
      rlang::env_poke(self$visited, self$from_key, TRUE)
      self$parent <- rlang::env()
      self$parent_enzyme <- rlang::env()
      self$parent_step <- rlang::env()
      self$all_edges <- list()

      self$found_keys_storage <- character(max(1L, length(self$to_keys)))
      self$found_tail <- 0L
      self$step <- 0L

      private$source_graph <- glyrepr::get_structure_graphs(self$from_g)
      self$queue_graphs <- list(private$source_graph)
      private$target_graphs <- glyrepr::get_structure_graphs(
        self$to_gs,
        return_list = TRUE
      )
      private$target_match_mode <- .glymotif_mode(self$to_gs)
      unique_target_graphs <- private$target_graphs[
        !duplicated(self$to_keys)
      ]
      private$pruning_target_graphs <- .order_bfs_pruning_targets(
        unique_target_graphs
      )
      private$product_match_mode <- .glymotif_mode(self$from_g)
      private$rule_plan <- .bfs_search_rule_plan(self$enzymes)
      private$product_cache <- fastmap::fastmap()
      private$n_core_graph <- glyrepr::get_structure_graphs(
        .n_glycan_starting_glycan("virtual")
      )
      private$pre_mgat2_graph <- glyrepr::get_structure_graphs(
        glyrepr::as_glycan_structure(
          "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
        )
      )
    },

    run = function() {
      result <- .bfs_native_run(self, private)
      if (length(result$missing_target_keys) && !self$allow_partial) {
        missing_targets <- result$missing_target_keys
        cli::cli_abort(
          "No synthesis path found for {length(missing_targets)} target(s) within {.val {self$max_steps}} steps."
        )
      }
      result
    }
  ),

  private = list(
    source_graph = NULL,
    target_graphs = NULL,
    pruning_target_graphs = NULL,
    target_match_mode = NULL,
    product_match_mode = NULL,
    rule_plan = NULL,
    product_cache = NULL,
    n_core_graph = NULL,
    pre_mgat2_graph = NULL,

    # Prune graphs after only the root-position normalization required by
    # glymotif, before paying for canonical branch ordering and IUPAC keys.
    prepare_graph_products = function(
      product_graphs,
      product_mode = private$product_match_mode
    ) {
      if (length(product_graphs) == 0L) {
        return(list(
          products = NULL,
          graphs = list(),
          keys = character()
        ))
      }

      prepared_products <- purrr::map(
        product_graphs,
        private$prepare_graph_product,
        product_mode = product_mode
      )
      keep <- vapply(prepared_products, `[[`, logical(1), "keep")
      prepared_products <- prepared_products[keep]
      if (length(prepared_products) == 0L) {
        return(list(
          products = NULL,
          graphs = list(),
          keys = character()
        ))
      }

      product_graphs <- lapply(prepared_products, `[[`, "graph")
      product_keys <- vapply(
        prepared_products,
        `[[`,
        character(1),
        "key"
      )
      unique_products <- !duplicated(product_keys)
      list(
        products = NULL,
        graphs = product_graphs[unique_products],
        keys = unname(product_keys[unique_products])
      )
    },

    # Cache standard graph products before target pruning and canonicalization.
    prepare_graph_product = function(product_graph, product_mode) {
      signature <- paste0(
        product_mode,
        "|",
        .bfs_graph_signature(product_graph)
      )
      if (
        !is.null(private$product_cache) &&
          private$product_cache$has(signature)
      ) {
        return(private$product_cache$get(signature))
      }

      product_graph <- .move_glycan_root_last(product_graph)
      keep <- private$is_promising_product(
        product_graph,
        product_mode = product_mode
      )
      if (keep) {
        keyed <- .canonicalize_valid_glycan_graphs(list(product_graph))
        prepared_product <- list(
          keep = TRUE,
          graph = keyed$graphs[[1]],
          key = unname(keyed$keys[[1]])
        )
      } else {
        prepared_product <- list(keep = FALSE)
      }

      if (!is.null(private$product_cache)) {
        private$product_cache$set(signature, prepared_product)
      }
      prepared_product
    },

    # Preserve structure materialization for filters and custom S3 enzymes.
    prepare_products = function(products) {
      if (length(products) == 0L) {
        return(list(
          products = products,
          graphs = list(),
          keys = character()
        ))
      }

      product_match_mode <- .glymotif_mode(products)
      product_graphs <- glyrepr::get_structure_graphs(
        products,
        return_list = TRUE
      )
      keep <- purrr::map_lgl(
        product_graphs,
        private$is_promising_product,
        product_mode = product_match_mode
      )

      list(
        products = products[keep],
        graphs = product_graphs[keep],
        keys = unname(as.character(products[keep]))
      )
    },

    is_promising_product = function(product_graph, product_mode) {
      .is_promising_intermediate_graph(
        product_graph,
        target_graphs = private$pruning_target_graphs,
        product_mode = product_mode,
        target_mode = private$target_match_mode,
        n_core_graph = private$n_core_graph,
        pre_mgat2_graph = private$pre_mgat2_graph
      )
    }
  )
)

#' Decide if the glycan is ready for MGAT2 or has been synthesized by MGAT2
#'
#' For N-glycans, mgat2_ready() returns TRUE if the α1-3/6 mannoses
#' on the α1-6 arm of the core Man have already been trimmed (i.e. the structure
#' is past the MAN2A1/2 trimming stage and is ready for / has passed MGAT2).
#' For non-N-glycans, the returned value is not biologically meaningful and is
#' only used downstream together with .is_n_glycan().
#'
#' @param glycans A glyrepr::glycan_structure vector.
#'
#' @returns A logical vector of the same length as `glycans`.
#' @noRd
mgat2_ready <- function(glycans) {
  !.have_motif_substituent_subset(
    glycans,
    "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
}

#' Check if any product is a promising intermediate
#'
#' A glycan is not a promising intermediate if it isn't a substructure (core-aligned)
#' of any of the target glycans.
#' N-glycans not ready for MGAT2 are always considered as promising intermediates.
#'
#' This pruning is crucial for the performance of the BFS search.
#' By removing unpromising intermediates, the search space is reduced significantly.
#'
#' This pruning is based on the assumption that,
#' apart from N-glycans before mannose trimming,
#' all other enzymatic steps in the biosynthesis of a glycan are
#' catalyzed by glycosyltransferases.
#'
#' @param products A glyrepr::glycan_structure vector.
#' @param target_glycans A glyrepr::glycan_structure vector.
#' @returns A logical vector of the same length as `products`.
#' @noRd
is_promising_intermediate <- function(products, target_glycans) {
  product_graphs <- glyrepr::get_structure_graphs(
    products,
    return_list = TRUE
  )
  target_graphs <- glyrepr::get_structure_graphs(
    target_glycans,
    return_list = TRUE
  )
  n_core_graph <- glyrepr::get_structure_graphs(
    .n_glycan_starting_glycan("virtual")
  )
  pre_mgat2_graph <- glyrepr::get_structure_graphs(
    glyrepr::as_glycan_structure(
      "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    )
  )
  purrr::map_lgl(
    product_graphs,
    .is_promising_intermediate_graph,
    target_graphs = target_graphs,
    product_mode = .glymotif_mode(products),
    target_mode = .glymotif_mode(target_glycans),
    n_core_graph = n_core_graph,
    pre_mgat2_graph = pre_mgat2_graph
  )
}

# Check one canonical product graph using graphs prepared for a BFS search.
.is_promising_intermediate_graph <- function(
  product_graph,
  target_graphs,
  product_mode,
  target_mode,
  n_core_graph,
  pre_mgat2_graph
) {
  if (
    .target_contains_intermediate_graph(
      product_graph,
      target_graphs,
      target_mode
    )
  ) {
    return(TRUE)
  }

  is_n_glycan <- length(.g_match_motif_substituent_subset(
    product_graph,
    n_core_graph,
    mode = product_mode
  )) >
    0L
  if (!is_n_glycan) {
    return(FALSE)
  }

  is_pre_mgat2 <- length(.g_match_motif_substituent_subset(
    product_graph,
    pre_mgat2_graph,
    mode = product_mode
  )) >
    0L
  if (!is_pre_mgat2) {
    return(FALSE)
  }

  # Before MGAT2, extra terminal Glc and Man residues may still be removed by
  # the supported N-glycan hydrolases. Remove only those reversible residues,
  # then require every irreversible decoration to fit at least one target.
  trimmed_graph <- .trim_pruning_n_glycan_graph(product_graph)
  .target_contains_intermediate_graph(
    trimmed_graph,
    target_graphs,
    target_mode
  )
}

.target_contains_intermediate_graph <- function(
  product_graph,
  target_graphs,
  target_mode
) {
  purrr::some(
    target_graphs,
    function(target_graph) {
      length(.g_match_motif_substituent_subset(
        target_graph,
        product_graph,
        alignment = "core",
        mode = target_mode
      )) >
        0L
    }
  )
}

.trim_pruning_n_glycan_graph <- function(graph) {
  repeat {
    monosaccharides <- igraph::vertex_attr(graph, "mono")
    substituents <- igraph::vertex_attr(graph, "sub")
    removable <- which(
      igraph::degree(graph, mode = "out") == 0L &
        monosaccharides %in% c("Glc", "Man", "Hex") &
        lengths(lapply(substituents, .substituent_tokens)) == 0L
    )
    if (length(removable) == 0L) {
      break
    }
    graph <- igraph::delete_vertices(graph, removable)
  }
  igraph::V(graph)$name <- as.character(seq_len(igraph::vcount(graph)))
  graph
}

# Find every vertex that can reach at least one target with one reverse
# traversal through a temporary multi-target sink.
.bfs_reverse_reachable <- function(graph, target_vertices) {
  target_vertices <- unique(as.integer(target_vertices))
  if (length(target_vertices) == 0L) {
    return(integer())
  }

  sink <- igraph::vcount(graph) + 1L
  graph_with_sink <- igraph::add_vertices(graph, 1L)
  sink_edges <- as.vector(rbind(
    target_vertices,
    rep.int(sink, length(target_vertices))
  ))
  graph_with_sink <- igraph::add_edges(graph_with_sink, sink_edges)
  reachable <- igraph::subcomponent(
    graph_with_sink,
    sink,
    mode = "in"
  )
  setdiff(as.integer(reachable), sink)
}

#' Build result graph from BFS search results
#'
#' Constructs an igraph object representing synthesis paths from BFS search results.
#' Uses bidirectional reachability to prune dead-end branches and keep only
#' paths that connect the source to target(s).
#'
#' @param search_result Results from BFS search containing edges and found keys
#' @param from_key Starting glycan key
#' @param to_keys Target glycan keys
#' @param max_steps Maximum steps for error reporting
#'
#' @returns igraph object representing synthesis path(s)
#' @noRd
build_synthesis_result_graph <- function(
  search_result,
  from_key,
  to_keys,
  max_steps
) {
  if (length(search_result$all_edges) == 0L) {
    # Return single-node graph for trivial case
    vertices <- tibble::tibble(name = from_key)
    return(igraph::graph_from_data_frame(
      tibble::tibble(
        from = character(0),
        to = character(0),
        enzyme = character(0),
        step = integer(0)
      ),
      directed = TRUE,
      vertices = vertices
    ))
  }

  edges_df <- tibble::tibble(
    from = vapply(search_result$all_edges, `[[`, character(1), "from"),
    to = vapply(search_result$all_edges, `[[`, character(1), "to"),
    enzyme = vapply(search_result$all_edges, `[[`, character(1), "enzyme"),
    step = unlist(
      lapply(search_result$all_edges, `[[`, "step"),
      use.names = FALSE
    )
  )
  edges_df <- unique(edges_df)

  g_all <- igraph::graph_from_data_frame(edges_df, directed = TRUE)

  # Check if all targets are present in the graph
  missing_targets <- setdiff(to_keys, igraph::V(g_all)$name)
  if (length(missing_targets) > 0L) {
    cli::cli_abort(
      "No synthesis path found for {length(missing_targets)} target(s) within {.val {max_steps}} steps."
    )
  }

  # ALGORITHM: Multi-Target Dead-End Pruning via Bidirectional Reachability
  # =======================================================================
  # Problem: BFS exploration creates many "dead-end" branches that don't lead to targets.
  # These should be excluded from the "all paths" result to show only valid synthesis routes.
  #
  # Solution: Use graph reachability to identify vertices that are both:
  # 1. Reachable from the starting glycan (forward reachability)
  # 2. Can reach at least one target glycan (backward reachability from any target)
  #
  # Algorithm:
  # 1. Compute forward reachable set: vertices reachable from 'from' via out-edges
  # 2. Connect every target to a temporary sink
  # 3. Compute one backward reachable set from that sink
  # 4. Keep intersection: vertices that lie on at least one valid from→(any target) path
  # 5. Return induced subgraph containing only these vertices and their connecting edges
  #
  # Complexity: O(V + E + number_of_targets)
  # Result: All dead-end branches are automatically pruned

  vid_from <- which(igraph::V(g_all)$name == from_key)
  vid_to_list <- vapply(
    to_keys,
    function(key) {
      which(igraph::V(g_all)$name == key)
    },
    integer(1)
  )

  # Forward reachability: vertices reachable from starting glycan
  reach_from <- igraph::subcomponent(g_all, vid_from, mode = "out")

  # Backward reachability: vertices that can reach any target glycan
  reach_to_any <- .bfs_reverse_reachable(g_all, vid_to_list)

  # Keep only vertices that are on valid synthesis paths to any target
  keep <- intersect(reach_from, reach_to_any)
  if (length(keep) == 0L) {
    cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
  }

  igraph::induced_subgraph(g_all, vids = keep)
}
