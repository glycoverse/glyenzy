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
#' Encapsulates breadth-first search state for synthesizing glycans into an R6
#' object. This replaces the earlier functional implementation that required
#' passing numerous state-tracking arguments between helper functions.
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
      private$rule_plan <- .prepare_bfs_rule_plan(self$enzymes)
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
      initial_targets <- private$matched_target_keys(
        self$from_key,
        private$source_graph
      )
      if (length(initial_targets) > 0L) {
        private$record_found_targets(
          rep(self$from_key, length(initial_targets)),
          initial_targets
        )
      }

      if (self$remaining_targets_map$size() == 0L) {
        return(list(
          found_keys = self$found_keys_storage[seq_len(self$found_tail)],
          all_edges = list(),
          parent = rlang::env(),
          parent_enzyme = rlang::env(),
          parent_step = rlang::env(),
          missing_target_keys = character()
        ))
      }

      while (
        length(self$queue_keys) > 0L &&
          self$step < self$max_steps &&
          self$remaining_targets_map$size() > 0L
      ) {
        self$step <- self$step + 1L

        found_keys <- private$expand_frontier()
        if (length(found_keys$endpoint_keys) > 0L) {
          private$record_found_targets(
            found_keys$endpoint_keys,
            found_keys$target_keys
          )
        }
      }

      if (
        self$remaining_targets_map$size() > 0L &&
          !self$allow_partial
      ) {
        missing_targets <- self$remaining_targets_map$keys()
        cli::cli_abort(
          "No synthesis path found for {length(missing_targets)} target(s) within {.val {self$max_steps}} steps."
        )
      }

      list(
        found_keys = if (self$found_tail == 0L) {
          character()
        } else {
          self$found_keys_storage[seq_len(self$found_tail)]
        },
        all_edges = self$all_edges,
        parent = self$parent,
        parent_enzyme = self$parent_enzyme,
        parent_step = self$parent_step,
        missing_target_keys = self$remaining_targets_map$keys()
      )
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

    # Expand entire BFS frontier for current step.
    # Applies every enzyme to all glycans in the queue, collects successors,
    # and accumulates any target hits produced during this level.
    expand_frontier = function() {
      can_batch <- is.null(self$filter) &&
        all(vapply(
          self$enzymes,
          .can_batch_bfs_enzyme,
          logical(1)
        ))
      if (can_batch) {
        return(private$expand_frontier_batched())
      }
      private$expand_frontier_scalar()
    },

    # Generate a whole frontier from shared rule jobs, then replay each
    # glycan-enzyme cell in the original order.
    expand_frontier_batched = function() {
      frontier_keys <- self$queue_keys
      frontier_graphs <- self$queue_graphs
      if (length(frontier_graphs) != length(frontier_keys)) {
        frontier_graphs <- purrr::map(
          self$queue,
          glyrepr::get_structure_graphs
        )
      }
      new_queue_graphs <- list()
      new_queue_keys <- character(0)
      found_endpoint_keys <- character(0)
      found_target_keys <- character(0)

      chunk_size <- 32L
      chunk_starts <- seq.int(1L, length(frontier_graphs), by = chunk_size)
      for (chunk_start in chunk_starts) {
        chunk_idx <- seq.int(
          chunk_start,
          min(chunk_start + chunk_size - 1L, length(frontier_graphs))
        )
        chunk_graphs <- frontier_graphs[chunk_idx]
        chunk_keys <- frontier_keys[chunk_idx]
        rule_results <- .apply_bfs_rule_frontier(
          chunk_graphs,
          rep(private$product_match_mode, length(chunk_graphs)),
          private$rule_plan,
          structure_level = self$structure_level,
          glycan_keys = chunk_keys
        )
        plan_products <- private$prepare_plan_products(
          rule_results,
          length(chunk_graphs)
        )

        for (i in seq_along(chunk_graphs)) {
          curr_key <- chunk_keys[[i]]
          for (enzyme_idx in seq_along(self$enzymes)) {
            ez <- self$enzymes[[enzyme_idx]]
            plan_id <- private$rule_plan$enzyme_plan_ids[[enzyme_idx]]
            expansion_result <- private$integrate_products(
              curr_key,
              ez,
              plan_products[[plan_id]][[i]]
            )

            if (length(expansion_result$new_graphs) > 0L) {
              start_idx <- length(new_queue_graphs)
              for (j in seq_along(expansion_result$new_graphs)) {
                new_queue_graphs[[
                  start_idx + j
                ]] <- expansion_result$new_graphs[[
                  j
                ]]
                new_queue_keys[start_idx + j] <- expansion_result$new_keys[[j]]
              }
            }

            if (length(expansion_result$found_endpoint_keys) > 0L) {
              found_endpoint_keys <- c(
                found_endpoint_keys,
                expansion_result$found_endpoint_keys
              )
              found_target_keys <- c(
                found_target_keys,
                expansion_result$found_target_keys
              )
            }
          }
        }
      }

      self$queue <- list()
      self$queue_graphs <- new_queue_graphs
      self$queue_keys <- new_queue_keys

      list(endpoint_keys = found_endpoint_keys, target_keys = found_target_keys)
    },

    # Retain scalar expansion when a user filter can have observable state.
    expand_frontier_scalar = function() {
      frontier <- self$queue
      frontier_keys <- self$queue_keys
      frontier_graphs <- self$queue_graphs
      if (length(frontier_graphs) != length(frontier_keys)) {
        frontier_graphs <- purrr::map(
          frontier,
          glyrepr::get_structure_graphs
        )
      }

      new_queue <- list()
      new_queue_graphs <- list()
      new_queue_keys <- character(0)
      found_endpoint_keys <- character(0)
      found_target_keys <- character(0)

      for (i in seq_along(frontier)) {
        curr_g <- frontier[[i]]
        curr_key <- frontier_keys[[i]]
        curr_graph <- frontier_graphs[[i]]

        for (enzyme_idx in seq_along(self$enzymes)) {
          ez <- self$enzymes[[enzyme_idx]]
          expansion_result <- private$expand_single(
            curr_g,
            curr_graph,
            curr_key,
            ez,
            private$rule_plan$prepared_rules[[enzyme_idx]],
            .glymotif_mode(curr_g)
          )

          if (length(expansion_result$new_graphs) > 0L) {
            new_structures <- expansion_result$new_structures
            if (length(new_structures) == 0L) {
              graph_lookup <- expansion_result$new_graphs
              names(graph_lookup) <- expansion_result$new_keys
              products <- glyrepr::new_glycan_structure(
                expansion_result$new_keys,
                graph_lookup
              )
              new_structures <- lapply(
                seq_along(products),
                function(j) products[j]
              )
            }
            start_idx <- length(new_queue)
            for (j in seq_along(new_structures)) {
              new_queue[[start_idx + j]] <- new_structures[[j]]
              new_queue_graphs[[start_idx + j]] <- expansion_result$new_graphs[[
                j
              ]]
              new_queue_keys[start_idx + j] <- expansion_result$new_keys[[j]]
            }
          }

          if (length(expansion_result$found_endpoint_keys) > 0L) {
            found_endpoint_keys <- c(
              found_endpoint_keys,
              expansion_result$found_endpoint_keys
            )
            found_target_keys <- c(
              found_target_keys,
              expansion_result$found_target_keys
            )
          }
        }
      }

      self$queue <- new_queue
      self$queue_graphs <- new_queue_graphs
      self$queue_keys <- new_queue_keys

      list(endpoint_keys = found_endpoint_keys, target_keys = found_target_keys)
    },

    # Expand one glycan by applying a single enzyme.
    # Generates candidate products, filters them, registers edges, updates
    # parent bookkeeping, and returns the new frontier entries plus targets hit.
    expand_single = function(
      curr_g,
      curr_graph,
      curr_key,
      ez,
      prepared_rules,
      rule_match_mode
    ) {
      if (.uses_standard_graph_action(ez)) {
        product_graphs <- .apply_enzyme_prepared_graphs(
          curr_graph,
          ez,
          prepared_rules,
          structure_level = self$structure_level,
          mode = rule_match_mode
        )
        prepared_products <- private$prepare_graph_products(
          product_graphs,
          product_mode = rule_match_mode
        )
      } else {
        products <- .apply_enzyme(
          curr_g,
          ez,
          structure_level = self$structure_level
        )[[1]]
        prepared_products <- private$prepare_products(products)
      }
      private$integrate_products(curr_key, ez, prepared_products)
    },

    # Prepare each shared rule result once, then replay unique enzyme plans.
    prepare_plan_products = function(rule_results, frontier_size) {
      prepared_rule_results <- purrr::map(
        rule_results,
        function(rule_result) {
          purrr::map(rule_result, private$prepare_graph_products)
        }
      )

      purrr::map(
        private$rule_plan$enzyme_plans,
        function(rule_ids) {
          lapply(
            seq_len(frontier_size),
            function(frontier_idx) {
              prepared_products <- lapply(
                rule_ids,
                function(rule_id) {
                  prepared_rule_results[[rule_id]][[frontier_idx]]
                }
              )
              private$combine_prepared_graph_products(prepared_products)
            }
          )
        }
      )
    },

    # Combine prepared rule cells in rule order and retain the first graph for
    # every canonical product key, matching glycan-vector unique() semantics.
    combine_prepared_graph_products = function(prepared_products) {
      if (length(prepared_products) == 0L) {
        return(list(
          products = NULL,
          graphs = list(),
          keys = character()
        ))
      }

      product_graphs <- unlist(
        lapply(prepared_products, `[[`, "graphs"),
        recursive = FALSE
      )
      product_keys <- unlist(
        lapply(prepared_products, `[[`, "keys"),
        use.names = FALSE
      )
      unique_products <- !duplicated(product_keys)
      list(
        products = NULL,
        graphs = product_graphs[unique_products],
        keys = product_keys[unique_products]
      )
    },

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
    },

    # Apply the user filter and replay BFS bookkeeping in scalar order.
    integrate_products = function(curr_key, ez, prepared_products) {
      products <- prepared_products$products
      product_graphs <- prepared_products$graphs
      prod_keys <- prepared_products$keys
      if (length(prod_keys) == 0L) {
        return(private$empty_expansion())
      }

      if (!is.null(self$filter)) {
        if (is.null(products)) {
          graph_lookup <- product_graphs
          names(graph_lookup) <- prod_keys
          products <- glyrepr::new_glycan_structure(
            prod_keys,
            graph_lookup
          )
        }
        keep <- self$filter(products)
        checkmate::assert_logical(
          keep,
          len = length(products),
          any.missing = FALSE
        )
        products <- products[keep]
        product_graphs <- product_graphs[keep]
        prod_keys <- prod_keys[keep]
        if (length(products) == 0L) {
          return(private$empty_expansion())
        }
      }

      new_structures <- list()
      new_graphs <- list()
      new_keys <- character(0)
      found_endpoint_keys <- character(0)
      found_target_keys <- character(0)

      for (j in seq_along(prod_keys)) {
        pk <- prod_keys[[j]]

        self$all_edges[[length(self$all_edges) + 1L]] <- list(
          from = curr_key,
          to = pk,
          enzyme = ez$name,
          step = self$step
        )

        if (!rlang::env_has(self$visited, pk)) {
          rlang::env_poke(self$visited, pk, TRUE)
          rlang::env_poke(self$parent, pk, curr_key)
          rlang::env_poke(self$parent_enzyme, pk, ez$name)
          rlang::env_poke(self$parent_step, pk, self$step)

          new_graphs[[length(new_graphs) + 1L]] <- product_graphs[[j]]
          if (!is.null(products)) {
            new_structures[[length(new_structures) + 1L]] <- products[j]
          }
          new_keys[length(new_keys) + 1L] <- pk
        }

        matched_targets <- private$matched_target_keys(
          pk,
          product_graphs[[j]]
        )
        if (length(matched_targets) > 0L) {
          found_endpoint_keys <- c(
            found_endpoint_keys,
            rep(pk, length(matched_targets))
          )
          found_target_keys <- c(found_target_keys, matched_targets)
        }
      }

      list(
        new_structures = new_structures,
        new_graphs = new_graphs,
        new_keys = new_keys,
        found_endpoint_keys = found_endpoint_keys,
        found_target_keys = found_target_keys
      )
    },

    empty_expansion = function() {
      list(
        new_structures = list(),
        new_graphs = list(),
        new_keys = character(0),
        found_endpoint_keys = character(0),
        found_target_keys = character(0)
      )
    },

    # Integrate newly discovered targets into tracking buffers and remaining set.
    # Resizes storage if needed and removes matched targets from the active goal set.
    record_found_targets = function(endpoint_keys, target_keys) {
      if (length(endpoint_keys) == 0L || length(target_keys) == 0L) {
        return(invisible(NULL))
      }

      remaining <- vapply(
        target_keys,
        self$remaining_targets_map$has,
        logical(1)
      )
      endpoint_keys <- endpoint_keys[remaining]
      target_keys <- target_keys[remaining]
      if (length(endpoint_keys) == 0L) {
        return(invisible(NULL))
      }

      required <- self$found_tail + length(endpoint_keys)
      if (required > length(self$found_keys_storage)) {
        new_len <- max(required, max(1L, length(self$found_keys_storage) * 2L))
        length(self$found_keys_storage) <- new_len
      }

      idx <- seq.int(self$found_tail + 1L, required)
      self$found_keys_storage[idx] <- endpoint_keys
      self$found_tail <- required

      for (target_key in target_keys) {
        if (self$remaining_targets_map$has(target_key)) {
          self$remaining_targets_map$remove(target_key)
        }
      }
      invisible(NULL)
    },

    matched_target_keys = function(key, glycan_graph) {
      if (self$remaining_targets_map$size() == 0L) {
        return(character())
      }

      if (identical(self$target_match, "key")) {
        if (self$remaining_targets_map$has(key)) {
          return(key)
        }
        return(character())
      }

      remaining_target_keys <- self$remaining_targets_map$keys()
      target_idx <- match(remaining_target_keys, self$to_keys)
      target_graphs <- private$target_graphs[target_idx]
      target_matches <- purrr::map_lgl(
        target_graphs,
        function(target_graph) {
          glymotif::.g_have_motif(
            glycan_graph,
            target_graph,
            alignment = "whole",
            mode = "lenient"
          )
        }
      )
      remaining_target_keys[target_matches]
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

  edges_df <- do.call(
    rbind,
    purrr::map(
      search_result$all_edges,
      ~ tibble::tibble(
        from = .x$from,
        to = .x$to,
        enzyme = .x$enzyme,
        step = .x$step
      )
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
