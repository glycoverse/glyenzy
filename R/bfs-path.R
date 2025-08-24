#' Breadth-First Search for Glycan Synthesis Paths
#'
#' Core BFS algorithms for finding synthesis paths between glycan structures.
#' These functions provide the algorithmic foundation for both find_synthesis_path()
#' and rebuild_biosynthesis().

#' Perform BFS search for synthesis paths between glycan structures
#'
#' This function implements a breadth-first search algorithm to find synthesis
#' paths from a starting glycan to one or more target glycans using enzymatic
#' reactions.
#'
#' @param from_g Starting glycan structure (single glyrepr::glycan_structure)
#' @param to_g Target glycan structure (single glyrepr::glycan_structure)
#' @param enzyme_names Character vector of enzyme names to use
#' @param max_steps Maximum number of synthesis steps to explore
#' @param filter Optional function to filter glycan structures at each step
#' @param from_key Optional pre-computed string key for starting glycan
#' @param to_key Optional pre-computed string key for target glycan
#'
#' @returns List with search results including found paths and exploration data
#'
#' @examples
#' \dontrun{
#' from_g <- glyparse::parse_iupac_condensed("Gal(b1-4)GlcNAc(b1-")
#' to_g <- glyparse::parse_iupac_condensed("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-")
#' result <- bfs_synthesis_search(from_g, to_g, "ST6GAL1", 3)
#' }
#'
#' @noRd
bfs_synthesis_search <- function(
  from_g,
  to_g,
  enzyme_names,
  max_steps,
  filter = NULL,
  from_key = NULL,
  to_key = NULL
) {
  # Compute keys if not provided
  if (is.null(from_key)) from_key <- as.character(from_g)[1]
  if (is.null(to_key)) to_key <- as.character(to_g)[1]

  # Check for trivial case
  if (from_key == to_key) {
    return(list(
      found_keys = to_key,
      all_edges = list(),
      parent = rlang::env(),
      parent_enzyme = rlang::env(),
      parent_step = rlang::env()
    ))
  }

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
  queue <- list(from_g)
  queue_keys <- c(from_key)
  visited <- rlang::env()  # Hash set for O(1) lookup
  rlang::env_poke(visited, from_key, TRUE)
  parent <- rlang::env()
  parent_enzyme <- rlang::env()
  parent_step <- rlang::env()

  found <- FALSE
  found_keys <- character(0)
  step <- 0L
  all_edges <- list()

  # BFS main loop: explore level by level
  while (length(queue) > 0L && step < max_steps) {
    step <- step + 1L

    bfs_result <- .expand_bfs_frontier_core(
      queue, queue_keys, to_key, enzyme_names, filter, step,
      visited, parent, parent_enzyme, parent_step, all_edges
    )

    queue <- bfs_result$queue
    queue_keys <- bfs_result$queue_keys
    all_edges <- bfs_result$all_edges

    if (length(bfs_result$found_keys) > 0L) {
      found <- TRUE
      found_keys <- c(found_keys, bfs_result$found_keys)
    }
  }

  if (!found) {
    cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
  }

  list(
    found_keys = found_keys,
    all_edges = all_edges,
    parent = parent,
    parent_enzyme = parent_enzyme,
    parent_step = parent_step
  )
}

#' Expand BFS frontier by one step (core algorithm)
#'
#' @param queue Current queue of glycan structures
#' @param queue_keys Current queue keys
#' @param to_key Target glycan key for goal testing
#' @param enzyme_names Available enzyme names
#' @param filter Optional filter function
#' @param step Current step number
#' @param visited Environment tracking visited nodes
#' @param parent Environment tracking parent relationships
#' @param parent_enzyme Environment tracking parent enzymes
#' @param parent_step Environment tracking parent steps
#' @param all_edges List of all explored edges
#' @return List with updated queue and found targets
#' @noRd
.expand_bfs_frontier_core <- function(
  queue,
  queue_keys,
  to_key,
  enzyme_names,
  filter,
  step,
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
    for (ez_name in enzyme_names) {
      expansion_result <- .expand_single_node_core(
        curr_g, curr_key, ez_name, to_key, filter, step,
        visited, parent, parent_enzyme, parent_step, all_edges
      )

      new_queue <- c(new_queue, expansion_result$new_structures)
      new_queue_keys <- c(new_queue_keys, expansion_result$new_keys)
      all_edges <- expansion_result$all_edges
      found_keys <- c(found_keys, expansion_result$found_keys)
    }
  }

  list(
    queue = new_queue,
    queue_keys = new_queue_keys,
    all_edges = all_edges,
    found_keys = found_keys
  )
}

#' Expand a single node with a single enzyme (core algorithm)
#'
#' @param curr_g Current glycan structure
#' @param curr_key Current glycan key
#' @param ez_name Enzyme name
#' @param to_key Target glycan key for goal testing
#' @param filter Optional filter function
#' @param step Current step number
#' @param visited Environment tracking visited nodes
#' @param parent Environment tracking parent relationships
#' @param parent_enzyme Environment tracking parent enzymes
#' @param parent_step Environment tracking parent steps
#' @param all_edges List of all explored edges
#' @return List with new structures and found targets
#' @noRd
.expand_single_node_core <- function(
  curr_g,
  curr_key,
  ez_name,
  to_key,
  filter,
  step,
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
  if (!is.null(filter)) {
    keep <- filter(products)
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
    all_edges[[length(all_edges) + 1L]] <- list(
      from = curr_key, to = pk, enzyme = ez_name, step = step
    )

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
    if (pk == to_key) {
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

#' Build result graph from BFS search results
#'
#' Constructs an igraph object representing synthesis paths from BFS search results.
#' Uses bidirectional reachability to prune dead-end branches and keep only
#' paths that connect the source to target.
#'
#' @param search_result Results from BFS search containing edges and found keys
#' @param from_key Starting glycan key
#' @param to_key Target glycan key
#' @param max_steps Maximum steps for error reporting
#'
#' @returns igraph object representing synthesis path(s)
#' @noRd
build_synthesis_result_graph <- function(search_result, from_key, to_key, max_steps) {
  if (length(search_result$all_edges) == 0L) {
    # Return single-node graph for trivial case
    vertices <- tibble::tibble(name = from_key)
    return(igraph::graph_from_data_frame(
      tibble::tibble(from = character(0), to = character(0),
                     enzyme = character(0), step = integer(0)),
      directed = TRUE, vertices = vertices
    ))
  }

  edges_df <- do.call(rbind, purrr::map(search_result$all_edges, ~ tibble::tibble(
    from = .x$from, to = .x$to, enzyme = .x$enzyme, step = .x$step
  )))
  edges_df <- unique(edges_df)

  g_all <- igraph::graph_from_data_frame(edges_df, directed = TRUE)

  if (!(to_key %in% igraph::V(g_all)$name)) {
    cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
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

  vid_from <- which(igraph::V(g_all)$name == from_key)
  vid_to <- which(igraph::V(g_all)$name == to_key)

  # Forward reachability: vertices reachable from starting glycan
  reach_from <- igraph::subcomponent(g_all, vid_from, mode = "out")
  # Backward reachability: vertices that can reach target glycan
  reach_to <- igraph::subcomponent(g_all, vid_to, mode = "in")

  # Keep only vertices that are on valid synthesis paths
  keep <- intersect(reach_from, reach_to)
  if (length(keep) == 0L) {
    cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
  }

  igraph::induced_subgraph(g_all, vids = keep)
}
