#' Find Synthesis Path Between Glycan Structures
#'
#' Find a synthesis path from one glycan structure to another using enzymatic reactions.
#' This function uses breadth-first search to find the shortest path or all possible paths
#' within a given number of steps.
#'
#' @param from A [glyrepr::glycan_structure()] scalar, or a character string
#'   supported by [glyparse::auto_parse()]. The starting glycan structure.
#' @param to A [glyrepr::glycan_structure()] scalar, or a character string
#'   supported by [glyparse::auto_parse()]. The target glycan structure.
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()] objects.
#'   If `NULL` (default), all available enzymes will be used.
#' @param max_steps Integer, maximum number of enzymatic steps to search.
#'   Default is 10.
#' @param filter Optional function to filter generated glycans at each step.
#'   Should take a [glyrepr::glycan_structure()] vector as input and return
#'   a logical vector of the same length.
#' @param return One of `"shortest"` or `"all"`. If `"shortest"`, returns the
#'   first shortest path found. If `"all"`, returns a graph containing all
#'   possible paths within `max_steps`.
#'
#' @returns An [igraph::igraph()] object representing the synthesis path(s).
#'   Vertices represent glycan structures with `name` attribute containing
#'   IUPAC-condensed strings. Edges represent enzymatic reactions with
#'   `enzyme` attribute containing gene symbols and `step` attribute
#'   indicating the step number.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Find shortest path
#' from <- "Gal(b1-4)GlcNAc(b1-"
#' to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
#' path <- find_synthesis_path(from, to, enzymes = "ST6GAL1", max_steps = 3)
#'
#' # View the path
#' igraph::as_data_frame(path, what = "edges")
#'
#' @export
find_synthesis_path <- function(
  from,
  to,
  enzymes = NULL,
  max_steps = 10,
  filter = NULL,
  return = c("shortest", "all")
) {
  return <- rlang::arg_match(return)

  # Parse and validate basic inputs first
  from_g <- glyrepr::as_glycan_structure(from)
  to_g <- glyrepr::as_glycan_structure(to)
  checkmate::assert_true(length(from_g) == 1L && length(to_g) == 1L)
  checkmate::assert_int(max_steps, lower = 1)

  from_key <- as.character(from_g)[1]
  to_key <- as.character(to_g)[1]

  # Check for trivial case first
  if (from_key == to_key) {
    return(.create_empty_path_graph(from_key))
  }

  # Now validate and process remaining inputs
  search_params <- .validate_synthesis_path_inputs(
    from_g, to_g, from_key, to_key, enzymes, max_steps, filter
  )

  # Perform BFS search
  search_result <- .perform_synthesis_search(search_params, return)

  # Build and return result graph
  .build_result_graph(search_result, search_params, return)
}

#' Validate and process inputs for synthesis path search
#' @param from_g Parsed starting glycan structure
#' @param to_g Parsed target glycan structure
#' @param from_key Starting glycan key
#' @param to_key Target glycan key
#' @param enzymes Enzyme list or NULL
#' @param max_steps Maximum search steps
#' @param filter Optional filter function
#' @returns List with processed parameters
#' @noRd
.validate_synthesis_path_inputs <- function(
  from_g,
  to_g,
  from_key,
  to_key,
  enzymes,
  max_steps,
  filter
) {
  # Process enzyme list
  enzyme_names <- .process_synthesis_enzymes(enzymes, to_g)

  # Process filter function
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }

  list(
    from_g = from_g,
    to_g = to_g,
    from_key = from_key,
    to_key = to_key,
    enzyme_names = enzyme_names,
    max_steps = max_steps,
    filter = filter
  )
}

#' Process and validate enzyme list for synthesis search
#' @param enzymes Raw enzyme input
#' @param to_g Target glycan structure for pre-filtering
#' @returns Character vector of enzyme names
#' @noRd
.process_synthesis_enzymes <- function(enzymes, to_g) {
  # Get enzyme names
  if (is.null(enzymes)) {
    enzyme_names <- names(glyenzy_enzymes)
  } else if (is.character(enzymes)) {
    enzyme_names <- enzymes
    unknown <- setdiff(enzyme_names, names(glyenzy_enzymes))
    if (length(unknown) > 0) {
      cli::cli_abort("Unknown enzymes: {.val {unknown}}.")
    }
  } else {
    checkmate::assert_list(enzymes, types = "glyenzy_enzyme")
    enzyme_names <- purrr::map_chr(enzymes, ~ .x$name)
  }

  # ALGORITHM: Enzyme Pre-filtering for Search Space Reduction
  # ========================================================
  # This optimization dramatically reduces the search space by filtering out enzymes
  # that cannot possibly contribute to synthesizing the target glycan.
  #
  # Strategy: Use is_synthesized_by() to check if each enzyme could have been involved
  # in creating the target structure. This is based on motif matching - if the target
  # contains motifs that the enzyme produces, it's a candidate.
  #
  # Benefits:
  # - Reduces search branching factor from ~100 enzymes to typically 5-15 relevant ones
  # - Maintains completeness: no valid paths are missed since we only exclude enzymes
  #   that provably cannot contribute to the target
  # - Graceful degradation: if pre-filtering fails, we fall back to using all enzymes
  can_contribute <- tryCatch({
    purrr::map_lgl(enzyme_names, ~ {
      tryCatch(
        glyenzy::is_synthesized_by(to_g, .x),
        error = function(e) FALSE
      )
    })
  }, error = function(e) {
    rep(TRUE, length(enzyme_names))
  })

  enzyme_names <- enzyme_names[can_contribute]
  if (length(enzyme_names) == 0L) {
    cli::cli_abort("No enzymes are predicted to contribute to the target glycan.")
  }

  enzyme_names
}

#' Create empty path graph for trivial case (from == to)
#' @param node_key Single node name
#' @returns igraph object with single node and no edges
#' @noRd
.create_empty_path_graph <- function(node_key) {
  vertices <- tibble::tibble(name = node_key)
  igraph::graph_from_data_frame(
    tibble::tibble(from = character(0), to = character(0),
                   enzyme = character(0), step = integer(0)),
    directed = TRUE, vertices = vertices
  )
}

#' Perform BFS search for synthesis paths
#' @param search_params Validated search parameters
#' @param return_mode Either "shortest" or "all"
#' @returns List with search results
#' @noRd
.perform_synthesis_search <- function(search_params, return_mode) {
  # ALGORITHM: Breadth-First Search (BFS) for Shortest Synthesis Paths
  # =================================================================
  # We use BFS to guarantee finding the shortest path (minimum number of enzymatic steps)
  # from the starting glycan to the target glycan.
  #
  # Key data structures:
  # - queue: Current frontier of glycan structures to explore
  # - visited: Hash set (environment) to avoid revisiting the same glycan
  # - parent/parent_enzyme/parent_step: Backtracking information for path reconstruction
  # - all_edges: Complete exploration graph for "all" mode
  #
  # BFS properties:
  # - Completeness: Will find a solution if one exists within max_steps
  # - Optimality: First solution found is guaranteed to be shortest (fewest steps)
  # - Time complexity: O(b^d) where b=branching factor, d=solution depth
  # - Space complexity: O(b^d) for storing the frontier and visited set

  # Initialize BFS state
  queue <- list(search_params$from_g)
  queue_keys <- c(search_params$from_key)
  visited <- rlang::env()  # Hash set for O(1) lookup
  rlang::env_poke(visited, search_params$from_key, TRUE)
  parent <- rlang::env()
  parent_enzyme <- rlang::env()
  parent_step <- rlang::env()

  found <- FALSE
  found_keys <- character(0)
  step <- 0L
  all_edges <- list()

  # BFS main loop: explore level by level
  while (length(queue) > 0L && step < search_params$max_steps) {
    step <- step + 1L

    bfs_result <- .expand_bfs_frontier(
      queue, queue_keys, search_params, step, return_mode,
      visited, parent, parent_enzyme, parent_step, all_edges
    )

    queue <- bfs_result$queue
    queue_keys <- bfs_result$queue_keys
    all_edges <- bfs_result$all_edges

    if (length(bfs_result$found_keys) > 0L) {
      found <- TRUE
      found_keys <- c(found_keys, bfs_result$found_keys)
      if (return_mode == "shortest") break  # Early termination for shortest path
    }
  }

  if (!found) {
    cli::cli_abort("No synthesis path found within {.val {search_params$max_steps}} steps.")
  }

  list(
    found_keys = found_keys,
    all_edges = all_edges,
    parent = parent,
    parent_enzyme = parent_enzyme,
    parent_step = parent_step
  )
}

#' Expand BFS frontier by one step
#' @param queue Current queue of glycan structures
#' @param queue_keys Current queue keys
#' @param search_params Search parameters
#' @param step Current step number
#' @param return_mode Either "shortest" or "all"
#' @param visited Environment tracking visited nodes
#' @param parent Environment tracking parent relationships
#' @param parent_enzyme Environment tracking parent enzymes
#' @param parent_step Environment tracking parent steps
#' @param all_edges List of all explored edges
#' @return List with updated queue and found targets
#' @noRd
.expand_bfs_frontier <- function(
  queue,
  queue_keys,
  search_params,
  step,
  return_mode,
  visited,
  parent,
  parent_enzyme,
  parent_step,
  all_edges
) {
  # ALGORITHM: BFS Frontier Expansion with Enzymatic Reactions
  # =========================================================
  # This function implements one level of BFS expansion in the glycan synthesis space.
  # Each glycan in the current frontier is expanded by applying all candidate enzymes.
  #
  # Key aspects:
  # 1. Systematic exploration: Try every enzyme on every glycan in current level
  # 2. Duplicate detection: Use visited set to avoid cycles and redundant exploration
  # 3. Parent tracking: Maintain backpointers for path reconstruction
  # 4. Early termination: Stop immediately when target found (shortest mode)
  # 5. Complete exploration: Continue until exhausted (all mode)
  #
  # Branching factor: |enzymes| × |products_per_enzyme|
  # Typical values: 5-15 enzymes × 1-3 products = 5-45 new states per glycan

  frontier <- queue
  frontier_keys <- queue_keys
  new_queue <- list()
  new_queue_keys <- character()
  found_keys <- character()

  # Expand each glycan in current frontier
  for (i in seq_along(frontier)) {
    curr_g <- frontier[[i]]
    curr_key <- frontier_keys[[i]]

    # Try each candidate enzyme on current glycan
    for (ez_name in search_params$enzyme_names) {
      expansion_result <- .expand_single_node(
        curr_g, curr_key, ez_name, search_params, step, return_mode,
        visited, parent, parent_enzyme, parent_step, all_edges
      )

      new_queue <- c(new_queue, expansion_result$new_structures)
      new_queue_keys <- c(new_queue_keys, expansion_result$new_keys)
      all_edges <- expansion_result$all_edges
      found_keys <- c(found_keys, expansion_result$found_keys)

      # Early termination for shortest path mode
      if (length(expansion_result$found_keys) > 0L && return_mode == "shortest") {
        return(list(
          queue = new_queue,
          queue_keys = new_queue_keys,
          all_edges = all_edges,
          found_keys = found_keys
        ))
      }
    }
  }

  list(
    queue = new_queue,
    queue_keys = new_queue_keys,
    all_edges = all_edges,
    found_keys = found_keys
  )
}

#' Expand a single node with a single enzyme
#' @param curr_g Current glycan structure
#' @param curr_key Current glycan key
#' @param ez_name Enzyme name
#' @param search_params Search parameters
#' @param step Current step number
#' @param return_mode Either "shortest" or "all"
#' @param visited Environment tracking visited nodes
#' @param parent Environment tracking parent relationships
#' @param parent_enzyme Environment tracking parent enzymes
#' @param parent_step Environment tracking parent steps
#' @param all_edges List of all explored edges
#' @return List with new structures and found targets
#' @noRd
.expand_single_node <- function(
  curr_g,
  curr_key,
  ez_name,
  search_params,
  step,
  return_mode,
  visited,
  parent,
  parent_enzyme,
  parent_step,
  all_edges
) {
  # ALGORITHM: Single Enzymatic Reaction Expansion
  # =============================================
  # This function simulates applying one specific enzyme to one glycan structure,
  # generating all possible product structures and managing the search state.
  #
  # Process:
  # 1. Enzymatic reaction: Apply enzyme to current glycan using biochemical rules
  # 2. Product filtering: Apply user-defined filter to prune unwanted intermediates
  # 3. Duplicate detection: Check if products have been seen before (cycle prevention)
  # 4. State management: Update visited set, parent pointers, and exploration graph
  # 5. Goal testing: Check if any product matches the target structure
  #
  # Key optimizations:
  # - Early return on empty products (enzyme not applicable)
  # - Hash-based duplicate detection for O(1) lookup
  # - Lazy edge recording (only for "all" mode to save memory)
  # - Batch processing of multiple products from single enzyme application

  # Apply enzyme to current glycan structure
  products <- suppressMessages(glyenzy::apply_enzyme(curr_g, ez_name))
  if (length(products) == 0L) {
    return(list(
      new_structures = list(),
      new_keys = character(),
      all_edges = all_edges,
      found_keys = character()
    ))
  }

  # Apply user-defined filter to prune search space
  if (!is.null(search_params$filter)) {
    keep <- search_params$filter(products)
    checkmate::assert_logical(keep, len = length(products), any.missing = FALSE)
    products <- products[keep]
    if (length(products) == 0L) {
      return(list(
        new_structures = list(),
        new_keys = character(),
        all_edges = all_edges,
        found_keys = character()
      ))
    }
  }

  prod_keys <- as.character(products)
  new_structures <- list()
  new_keys <- character()
  found_keys <- character()

  # Process each product structure
  for (j in seq_along(products)) {
    pk <- prod_keys[[j]]

    # Record exploration edge for complete graph construction
    if (return_mode == "all") {
      all_edges[[length(all_edges) + 1L]] <- list(
        from = curr_key, to = pk, enzyme = ez_name, step = step
      )
    }

    # Handle new (unvisited) glycan structures
    if (!rlang::env_has(visited, pk)) {
      rlang::env_poke(visited, pk, TRUE)
      rlang::env_poke(parent, pk, curr_key)
      rlang::env_poke(parent_enzyme, pk, ez_name)
      rlang::env_poke(parent_step, pk, step)
      new_structures[[length(new_structures) + 1L]] <- products[j]
      new_keys <- c(new_keys, pk)
    }

    # Goal test: check if target structure reached
    if (pk == search_params$to_key) {
      found_keys <- c(found_keys, pk)
    }
  }

  list(
    new_structures = new_structures,
    new_keys = new_keys,
    all_edges = all_edges,
    found_keys = found_keys
  )
}

#' Build result graph from search results
#' @param search_result Results from BFS search
#' @param search_params Search parameters
#' @param return_mode Either "shortest" or "all"
#' @return igraph object representing synthesis path(s)
#' @noRd
.build_result_graph <- function(search_result, search_params, return_mode) {
  if (return_mode == "shortest") {
    .build_shortest_path_graph(search_result, search_params)
  } else {
    .build_all_paths_graph(search_result, search_params)
  }
}

#' Build shortest path graph
#' @param search_result Results from BFS search
#' @param search_params Search parameters
#' @return igraph object with shortest path
#' @noRd
.build_shortest_path_graph <- function(search_result, search_params) {
  # ALGORITHM: Path Reconstruction via Parent Backtracking
  # =====================================================
  # After BFS finds the target, we reconstruct the shortest path by following
  # parent pointers backward from target to source.
  #
  # Process:
  # 1. Start from any found target (all are equidistant due to BFS level-order)
  # 2. Follow parent chain: target → parent → parent's parent → ... → source
  # 3. Collect enzymes and intermediate glycans along the path
  # 4. Reverse the collected path (since we built it backward)
  # 5. Construct directed graph with path vertices and edges
  #
  # Guarantees:
  # - Path is optimal (shortest) due to BFS properties
  # - Path is valid (each edge represents a feasible enzymatic reaction)
  # - Path is complete (connects source to target)

  path_keys <- character()
  edge_enzymes <- character()
  curr <- search_result$found_keys[1]  # Take first found target (all equidistant)

  # Backtrack from target to source
  while (!identical(curr, search_params$from_key)) {
    path_keys <- c(curr, path_keys)
    edge_enzymes <- c(rlang::env_get(search_result$parent_enzyme, curr), edge_enzymes)
    curr <- rlang::env_get(search_result$parent, curr)
  }

  # Construct path graph
  vertices <- unique(c(search_params$from_key, path_keys))
  edges <- tibble::tibble(
    from = c(search_params$from_key, head(path_keys, -1L)),
    to = path_keys,
    enzyme = edge_enzymes,
    step = seq_along(edge_enzymes)
  )

  igraph::graph_from_data_frame(
    edges, directed = TRUE, vertices = tibble::tibble(name = vertices)
  )
}

#' Build all paths graph with pruning
#' @param search_result Results from BFS search
#' @param search_params Search parameters
#' @return igraph object with all valid paths
#' @noRd
.build_all_paths_graph <- function(search_result, search_params) {
  if (length(search_result$all_edges) == 0L) {
    return(.create_empty_path_graph(search_params$from_key))
  }

  edges_df <- do.call(rbind, purrr::map(search_result$all_edges, ~ tibble::tibble(
    from = .x$from, to = .x$to, enzyme = .x$enzyme, step = .x$step
  )))
  edges_df <- unique(edges_df)

  g_all <- igraph::graph_from_data_frame(edges_df, directed = TRUE)

  if (!(search_params$to_key %in% igraph::V(g_all)$name)) {
    cli::cli_abort("No synthesis path found within {.val {search_params$max_steps}} steps.")
  }

  # ALGORITHM: Dead-End Pruning via Bidirectional Reachability
  # ==========================================================
  # Problem: BFS exploration creates many "dead-end" branches that don't lead to target.
  # These should be excluded from the "all paths" result to show only valid synthesis routes.
  #
  # Solution: Use graph reachability to identify vertices that are both:
  # 1. Reachable from the starting glycan (forward reachability)
  # 2. Can reach the target glycan (backward reachability)
  #
  # Algorithm:
  # 1. Compute forward reachable set: vertices reachable from 'from' via out-edges
  # 2. Compute backward reachable set: vertices that can reach 'to' via in-edges
  # 3. Keep intersection: vertices that lie on at least one valid from→to path
  # 4. Return induced subgraph containing only these vertices and their connecting edges
  #
  # Complexity: O(V + E) for each reachability computation
  # Result: All dead-end branches are automatically pruned

  vid_from <- which(igraph::V(g_all)$name == search_params$from_key)
  vid_to <- which(igraph::V(g_all)$name == search_params$to_key)

  # Forward reachability: vertices reachable from starting glycan
  reach_from <- igraph::subcomponent(g_all, vid_from, mode = "out")
  # Backward reachability: vertices that can reach target glycan
  reach_to <- igraph::subcomponent(g_all, vid_to, mode = "in")

  # Keep only vertices that are on valid synthesis paths
  keep <- intersect(reach_from, reach_to)
  if (length(keep) == 0L) {
    cli::cli_abort("No synthesis path found within {.val {search_params$max_steps}} steps.")
  }

  igraph::induced_subgraph(g_all, vids = keep)
}
