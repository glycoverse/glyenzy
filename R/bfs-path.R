#' Breadth-First Search for Glycan Synthesis Paths
#'
#' Core BFS algorithms for finding synthesis paths between glycan structures.
#' These functions provide the algorithmic foundation for both find_synthesis_path()
#' and rebuild_biosynthesis().

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
  to_keys = NULL
) {
  engine <- BfsSynthesisSearch$new(
    from_g = from_g,
    to_gs = to_gs,
    enzymes = enzymes,
    max_steps = max_steps,
    filter = filter,
    from_key = from_key,
    to_keys = to_keys
  )

  engine$run()
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
    target_key_map = NULL,
    remaining_targets_map = NULL,
    queue = NULL,
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
      to_keys = NULL
    ) {
      self$from_g <- from_g
      self$to_gs <- to_gs
      self$enzymes <- enzymes
      self$max_steps <- max_steps
      self$filter <- filter
      self$from_key <- if (is.null(from_key)) as.character(from_g)[1] else from_key
      self$to_keys <- if (is.null(to_keys)) as.character(to_gs) else to_keys

      self$target_key_map <- fastmap::fastmap()
      self$remaining_targets_map <- fastmap::fastmap()
      for (target_key in self$to_keys) {
        self$target_key_map$set(target_key, TRUE)
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
    },

    run = function() {
      trivial_targets <- self$to_keys[self$to_keys == self$from_key]
      if (length(trivial_targets) > 0L) {
        return(list(
          found_keys = trivial_targets,
          all_edges = list(),
          parent = rlang::env(),
          parent_enzyme = rlang::env(),
          parent_step = rlang::env()
        ))
      }

      while (
        length(self$queue) > 0L &&
          self$step < self$max_steps &&
          self$remaining_targets_map$size() > 0L
      ) {
        self$step <- self$step + 1L

        found_keys <- private$expand_frontier()
        if (length(found_keys) > 0L) {
          private$record_found_keys(found_keys)
        }
      }

      if (self$remaining_targets_map$size() > 0L) {
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
        parent_step = self$parent_step
      )
    }
  ),

  private = list(
    # Expand entire BFS frontier for current step.
    # Applies every enzyme to all glycans in the queue, collects successors,
    # and accumulates any target hits produced during this level.
    expand_frontier = function() {
      frontier <- self$queue
      frontier_keys <- self$queue_keys

      new_queue <- list()
      new_queue_keys <- character(0)
      found_keys <- character(0)

      for (i in seq_along(frontier)) {
        curr_g <- frontier[[i]]
        curr_key <- frontier_keys[[i]]

        for (ez in self$enzymes) {
          expansion_result <- private$expand_single(curr_g, curr_key, ez)

          if (length(expansion_result$new_structures) > 0L) {
            start_idx <- length(new_queue)
            for (j in seq_along(expansion_result$new_structures)) {
              new_queue[[start_idx + j]] <- expansion_result$new_structures[[j]]
              new_queue_keys[start_idx + j] <- expansion_result$new_keys[[j]]
            }
          }

          if (length(expansion_result$found_keys) > 0L) {
            found_keys <- c(found_keys, expansion_result$found_keys)
          }
        }
      }

      self$queue <- new_queue
      self$queue_keys <- new_queue_keys

      found_keys
    },

    # Expand one glycan by applying a single enzyme.
    # Generates candidate products, filters them, registers edges, updates
    # parent bookkeeping, and returns the new frontier entries plus targets hit.
    expand_single = function(curr_g, curr_key, ez) {
      zero_products <- function() {
        list(
          new_structures = list(),
          new_keys = character(0),
          found_keys = character(0)
        )
      }

      products <- suppressMessages(glyenzy::apply_enzyme(curr_g, ez))
      if (length(products) == 0L) return(zero_products())

      keep <- is_promising_intermediate(products, self$to_gs)
      products <- products[keep]
      if (length(products) == 0L) return(zero_products())

      if (!is.null(self$filter)) {
        keep <- self$filter(products)
        checkmate::assert_logical(keep, len = length(products), any.missing = FALSE)
        products <- products[keep]
        if (length(products) == 0L) return(zero_products())
      }

      prod_keys <- as.character(products)
      new_structures <- list()
      new_keys <- character(0)
      found_keys <- character(0)

      for (j in seq_along(products)) {
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

          new_structures[[length(new_structures) + 1L]] <- products[j]
          new_keys[length(new_keys) + 1L] <- pk
        }

        if (self$target_key_map$has(pk)) {
          found_keys <- c(found_keys, pk)
        }
      }

      list(
        new_structures = new_structures,
        new_keys = new_keys,
        found_keys = found_keys
      )
    },

    # Integrate newly discovered targets into tracking buffers and remaining set.
    # Resizes storage if needed and removes matched targets from the active goal set.
    record_found_keys = function(keys) {
      required <- self$found_tail + length(keys)
      if (required > length(self$found_keys_storage)) {
        new_len <- max(required, max(1L, length(self$found_keys_storage) * 2L))
        length(self$found_keys_storage) <- new_len
      }

      idx <- seq.int(self$found_tail + 1L, required)
      self$found_keys_storage[idx] <- keys
      self$found_tail <- required

      for (found_key in keys) {
        if (self$remaining_targets_map$has(found_key)) {
          self$remaining_targets_map$remove(found_key)
        }
      }
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
  !glymotif::have_motif(glycans, "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
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
  res_mat <- glymotif::have_motifs(target_glycans, products, alignment = "core")
  res <- colSums(res_mat) > 0L
  res[.is_n_glycan(products) & !mgat2_ready(products)] <- TRUE
  res
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
build_synthesis_result_graph <- function(search_result, from_key, to_keys, max_steps) {
  if (length(search_result$all_edges) == 0L) {
    # Return single-node graph for trivial case
    vertices <- tibble::tibble(name = from_key)
    return(igraph::graph_from_data_frame(
      tibble::tibble(from = character(0), to = character(0), enzyme = character(0), step = integer(0)),
      directed = TRUE, vertices = vertices
    ))
  }

  edges_df <- do.call(rbind, purrr::map(search_result$all_edges, ~ tibble::tibble(
    from = .x$from, to = .x$to, enzyme = .x$enzyme, step = .x$step
  )))
  edges_df <- unique(edges_df)

  g_all <- igraph::graph_from_data_frame(edges_df, directed = TRUE)

  # Check if all targets are present in the graph
  missing_targets <- setdiff(to_keys, igraph::V(g_all)$name)
  if (length(missing_targets) > 0L) {
    cli::cli_abort("No synthesis path found for {length(missing_targets)} target(s) within {.val {max_steps}} steps.")
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
  # 2. For each target, compute backward reachable set: vertices that can reach this target via in-edges
  # 3. Union all backward reachable sets to get vertices that can reach any target
  # 4. Keep intersection: vertices that lie on at least one valid from→(any target) path
  # 5. Return induced subgraph containing only these vertices and their connecting edges
  #
  # Complexity: O(V + E) × (1 + number_of_targets) for reachability computations
  # Result: All dead-end branches are automatically pruned

  vid_from <- which(igraph::V(g_all)$name == from_key)
  vid_to_list <- sapply(to_keys, function(key) which(igraph::V(g_all)$name == key))

  # Forward reachability: vertices reachable from starting glycan
  reach_from <- igraph::subcomponent(g_all, vid_from, mode = "out")

  # Backward reachability: vertices that can reach any target glycan
  reach_to_any <- integer(0)
  for (vid_to in vid_to_list) {
    reach_to_current <- igraph::subcomponent(g_all, vid_to, mode = "in")
    reach_to_any <- union(reach_to_any, reach_to_current)
  }

  # Keep only vertices that are on valid synthesis paths to any target
  keep <- intersect(reach_from, reach_to_any)
  if (length(keep) == 0L) {
    cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
  }

  igraph::induced_subgraph(g_all, vids = keep)
}
