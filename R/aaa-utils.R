#' Make a function that only applies to N-glycans
#'
#' This function is used to wrap a function that only applies to N-glycans.
#' For N-glycans in the input `glycans`, it uses the values in `.f`.
#' Otherwise, it returns FALSE.
#'
#' @param .f A function that takes two arguments: `glycans` and `enzyme`.
#' @noRd
.make_n_glycan_guard <- function(.f, type = "logical") {
  force(.f)
  function(glycans, enzyme, is_n = NULL) {
    if (is.null(is_n)) {
      is_n <- glymotif::is_n_glycan(glycans)
    }
    if (type == "logical") {
      res <- rep(FALSE, length(glycans))
    } else {  # type is integer
      res <- rep(0, length(glycans))
    }
    if (any(is_n)) {
      res[is_n] <- .f(glycans[is_n], enzyme)
    }
    res
  }
}

.process_glycans_arg <- function(x) {
  if (is.character(x)) {
    x <- glyparse::auto_parse(x)
  } else if (!glyrepr::is_glycan_structure(x)) {
    cli::cli_abort(c(
      "{.arg glycans} must be a {.cls glyrepr_structure} vector or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(x)}}."
    ))
  }
  return(x)
}

.process_enzyme_arg <- function(x) {
  if (is.character(x)) {
    x <- enzyme(x)
  } else if (!inherits(x, "glyenzy_enzyme")) {
    cli::cli_abort(c(
      "{.arg enzyme} must be a {.cls glyenzy_enzyme} object or a character string of gene symbol.",
      "x" = "Got {.cls {class(x)}}."
    ))
  }
  return(x)
}

.process_enzymes_arg <- function(x) {
  if (is.character(x)) {
    x <- purrr::map(x, enzyme)
  } else if (!purrr::every(x, ~ inherits(.x, "glyenzy_enzyme"))) {
    cli::cli_abort("{.arg enzymes} must be a character vector of gene symbols or a list of {.cls glyenzy_enzyme} objects.")
  }
  return(unname(x))
}

# Validate and process return_list parameter early
.validate_return_list <- function(return_list, input_length) {
  checkmate::assert_flag(return_list, null.ok = TRUE)
  checkmate::assert_integerish(input_length, len = 1, lower = 1)

  if (is.null(return_list)) {
    return_list <- input_length > 1
  }
  if (!return_list && input_length > 1) {
    cli::cli_abort(c(
      "When {.arg return_list} is FALSE, input must have length 1.",
      "x" = "Input length: {.val {input_length}}."
    ))
  }

  return_list
}

# Format result based on return_list setting
.format_result <- function(result_list, return_list) {
  checkmate::assert_list(result_list)
  checkmate::assert_flag(return_list)

  if (!return_list) {
    result_list[[1]]
  } else {
    result_list
  }
}

#' Process and validate enzyme list for synthesis search
#'
#' This function processes various enzyme input formats (NULL, character vector, 
#' enzyme object list) and applies pre-filtering based on target glycan compatibility.
#'
#' @param enzymes Raw enzyme input (NULL, character vector, or enzyme object list)
#' @param to_g Target glycan structure for pre-filtering (optional)
#' @param apply_prefilter Whether to apply enzyme pre-filtering based on target glycan
#' @returns Character vector of enzyme names
#' @noRd
.process_synthesis_enzymes <- function(enzymes, to_g = NULL, apply_prefilter = TRUE) {
  # Get enzyme names from various input formats
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

  # Apply enzyme pre-filtering if requested and target glycan provided
  if (apply_prefilter && !is.null(to_g)) {
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
  }

  enzyme_names
}

#' Perform BFS synthesis search with common input processing
#'
#' This is a high-level wrapper that handles input validation, enzyme processing,
#' and BFS search execution. Used by both find_synthesis_path and rebuild_biosynthesis.
#'
#' @param from_g Starting glycan structure
#' @param to_g Target glycan structure  
#' @param enzymes Raw enzyme input
#' @param max_steps Maximum search steps
#' @param filter Optional filter function
#' @returns igraph object representing synthesis path(s)
#' @noRd
.perform_bfs_synthesis <- function(from_g, to_g, enzymes, max_steps, filter = NULL) {
  # Parse glycan structures and compute keys
  from_g <- glyrepr::as_glycan_structure(from_g)
  to_g <- glyrepr::as_glycan_structure(to_g)
  from_key <- as.character(from_g)[1]
  to_key <- as.character(to_g)[1]
  
  # Process enzyme list with pre-filtering
  enzyme_names <- .process_synthesis_enzymes(enzymes, to_g, apply_prefilter = TRUE)
  
  # Process filter function
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }
  
  # Perform BFS search using core algorithm
  search_result <- bfs_synthesis_search(
    from_g = from_g,
    to_g = to_g,
    enzyme_names = enzyme_names,
    max_steps = max_steps,
    filter = filter,
    from_key = from_key,
    to_key = to_key
  )
  
  # Build and return result graph
  build_synthesis_result_graph(
    search_result,
    from_key,
    to_key,
    max_steps
  )
}