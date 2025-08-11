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
#' @return An [igraph::igraph()] object representing the synthesis path(s).
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
find_synthesis_path <- function(from, to, enzymes = NULL, max_steps = 10,
                                filter = NULL, return = c("shortest", "all")) {
  return <- rlang::arg_match(return)
  
  # Parse and validate inputs
  from_g <- glyrepr::as_glycan_structure(from)
  to_g <- glyrepr::as_glycan_structure(to)
  checkmate::assert_true(length(from_g) == 1L && length(to_g) == 1L)
  checkmate::assert_int(max_steps, lower = 1)
  
  # Get enzyme list
  if (is.null(enzymes)) {
    # Use all available enzymes
    enzyme_names <- names(glyenzy_enzymes)
  } else {
    # Validate provided enzymes
    if (is.character(enzymes)) {
      enzyme_names <- enzymes
      # Check if all enzymes exist
      unknown <- setdiff(enzyme_names, names(glyenzy_enzymes))
      if (length(unknown) > 0) {
        cli::cli_abort("Unknown enzymes: {.val {unknown}}.")
      }
    } else {
      # List of enzyme objects
      checkmate::assert_list(enzymes, types = "glyenzy_enzyme")
      enzyme_names <- purrr::map_chr(enzymes, ~ .x$name)
    }
  }
  
  from_key <- as.character(from_g)[1]
  to_key <- as.character(to_g)[1]

  # Check if from and to are the same
  if (from_key == to_key) {
    # Return a single-node graph
    vertices <- tibble::tibble(name = from_key)
    return(igraph::graph_from_data_frame(
      tibble::tibble(from = character(0), to = character(0),
                     enzyme = character(0), step = integer(0)),
      directed = TRUE, vertices = vertices
    ))
  }

  # Pre-filter enzymes using is_synthesized_by
  can_contribute <- tryCatch({
    purrr::map_lgl(enzyme_names, ~ {
      tryCatch(
        glyenzy::is_synthesized_by(to_g, .x),
        error = function(e) FALSE
      )
    })
  }, error = function(e) {
    # If pre-filtering fails, keep all enzymes
    rep(TRUE, length(enzyme_names))
  })

  enzyme_names <- enzyme_names[can_contribute]
  if (length(enzyme_names) == 0L) {
    cli::cli_abort("No enzymes are predicted to contribute to the target glycan.")
  }
  
  # BFS state
  queue <- list(from_g)
  queue_keys <- c(from_key)
  visited <- rlang::env()
  rlang::env_poke(visited, from_key, TRUE)
  parent <- rlang::env()          # child_key -> parent_key
  parent_enzyme <- rlang::env()   # child_key -> enzyme_name
  parent_step <- rlang::env()     # child_key -> step number
  
  found <- FALSE
  found_keys <- character(0)
  step <- 0L
  
  # For "all" mode, we need to track all edges
  all_edges <- list()
  
  while (length(queue) > 0L && step < max_steps) {
    step <- step + 1L
    
    # Current frontier
    frontier <- queue
    frontier_keys <- queue_keys
    queue <- list()
    queue_keys <- character()
    
    # Expand each node by all enzymes
    for (i in seq_along(frontier)) {
      curr_g <- frontier[[i]]
      curr_key <- frontier_keys[[i]]
      
      # Try each enzyme
      for (ez_name in enzyme_names) {
        products <- suppressMessages(glyenzy::apply_enzyme(curr_g, ez_name))
        if (length(products) == 0L) next
        
        # Optional pruning
        if (!is.null(filter)) {
          keep <- filter(products)
          checkmate::assert_logical(keep, len = length(products), any.missing = FALSE)
          products <- products[keep]
          if (length(products) == 0L) next
        }
        
        prod_keys <- as.character(products)
        
        # For each product, register and enqueue if new
        for (j in seq_along(products)) {
          pk <- prod_keys[[j]]
          
          # For "all" mode, always record the edge
          if (return == "all") {
            all_edges[[length(all_edges) + 1L]] <- list(
              from = curr_key, to = pk, enzyme = ez_name, step = step
            )
          }
          
          # Check if this is a new node
          if (!rlang::env_has(visited, pk)) {
            rlang::env_poke(visited, pk, TRUE)
            rlang::env_poke(parent, pk, curr_key)
            rlang::env_poke(parent_enzyme, pk, ez_name)
            rlang::env_poke(parent_step, pk, step)
            queue[[length(queue) + 1L]] <- products[j]
            queue_keys <- c(queue_keys, pk)
          }
          
          # Goal test
          if (pk == to_key) {
            found <- TRUE
            found_keys <- c(found_keys, pk)
            if (return == "shortest") {
              break  # Found shortest path, exit immediately
            }
          }
        }
        if (found && return == "shortest") break
      } # enzymes
      if (found && return == "shortest") break
    } # frontier
    
    if (found && return == "shortest") break
  } # while
  
  if (!found) {
    cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
  }
  
  if (return == "shortest") {
    # Reconstruct shortest path
    path_keys <- character()
    edge_enzymes <- character()
    curr <- found_keys[1]  # Take the first found path
    
    while (!identical(curr, from_key)) {
      path_keys <- c(curr, path_keys)
      edge_enzymes <- c(rlang::env_get(parent_enzyme, curr), edge_enzymes)
      curr <- rlang::env_get(parent, curr)
    }
    
    vertices <- unique(c(from_key, path_keys))
    edges <- tibble::tibble(
      from = c(from_key, head(path_keys, -1L)),
      to = path_keys,
      enzyme = edge_enzymes,
      step = seq_along(edge_enzymes)
    )
    
    igraph::graph_from_data_frame(edges, directed = TRUE, 
                                  vertices = tibble::tibble(name = vertices))
  } else {
    # Return all valid paths as the union subgraph of all from->to paths
    # Build explored graph from all_edges, then prune to only vertices/edges
    # that are both reachable from `from` and can reach `to`.
    if (length(all_edges) == 0L) {
      vertices <- tibble::tibble(name = from_key)
      return(igraph::graph_from_data_frame(
        tibble::tibble(from = character(0), to = character(0),
                       enzyme = character(0), step = integer(0)),
        directed = TRUE, vertices = vertices
      ))
    }

    edges_df <- do.call(rbind, purrr::map(all_edges, ~ tibble::tibble(
      from = .x$from, to = .x$to, enzyme = .x$enzyme, step = .x$step
    )))
    edges_df <- unique(edges_df)

    g_all <- igraph::graph_from_data_frame(edges_df, directed = TRUE)

    # If `to` never seen in explored graph, no valid path exists
    if (!(to_key %in% igraph::V(g_all)$name)) {
      cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
    }

    vid_from <- which(igraph::V(g_all)$name == from_key)
    vid_to   <- which(igraph::V(g_all)$name == to_key)

    # Vertices reachable from `from`
    reach_from <- igraph::subcomponent(g_all, vid_from, mode = "out")
    # Vertices that can reach `to` in original graph = vertices in the in-component of `to`
    reach_to   <- igraph::subcomponent(g_all, vid_to, mode = "in")

    keep <- intersect(reach_from, reach_to)
    if (length(keep) == 0L) {
      cli::cli_abort("No synthesis path found within {.val {max_steps}} steps.")
    }

    igraph::induced_subgraph(g_all, vids = keep)
  }
}
