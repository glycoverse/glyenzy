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
      is_n <- .is_n_glycan(glycans)
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

.is_n_glycan <- function(x) {
  glymotif::have_motif(x, "N-Glycan core basic")
}

.process_glycan_arg <- function(x) {
  if (length(x) == 1L) {
    .process_glycans_arg(x)
  } else {
    cli::cli_abort(c(
      "{.arg x} must have length 1.",
      "x" = "Got {.val {length(x)}}."
    ))
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

  is_concrete <- .is_concrete_glycan(x)
  if (!all(is_concrete)) {
    cli::cli_abort(c(
      "All glycans must have concrete monosaccharides (e.g. Gal, GlcNAc, etc.).",
      "x" = "These glycans are not concrete: {.val {unique(x[!is_concrete])}}."
    ))
  }

  has_intact_linkages <- .has_intact_linkages(x)
  if (!all(has_intact_linkages)) {
    cli::cli_abort(c(
      "All linkages must be intact (no `?`).",
      "x" = "These glycans have unknown linkages: {.val {unique(x[!has_intact_linkages])}}."
    ))
  }

  has_substituents <- .has_substituents(x)
  if (any(has_substituents)) {
    cli::cli_abort(c(
      "Glycans with substituents are not supported.",
      "x" = "These glycans have substituents: {.val {unique(x[has_substituents])}}.",
      "i" = "Use {.fn glyrepr::remove_substituents} to get clean glycans."
    ))
  }

  x
}

.is_concrete_glycan <- function(x) {
  glyrepr::get_mono_type(x) == "concrete"
}

.has_intact_linkages <- function(x) {
  all_linkages_intact <- function(graph) {
    (all(!stringr::str_detect(igraph::E(graph)$linkage, stringr::fixed("?"))) &&
     !stringr::str_detect(graph$anomer, stringr::fixed("?")))
  }
  glyrepr::smap_lgl(x, all_linkages_intact)
}

.has_substituents <- function(x) {
  has_sub_single <- function(graph) {
    purrr::some(igraph::V(graph)$sub, ~ .x != "")
  }
  glyrepr::smap_lgl(x, has_sub_single)
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

#' Process and validate enzyme list for synthesis search
#'
#' This function processes various enzyme input formats (NULL, character vector, 
#' enzyme object list) and applies pre-filtering based on target glycan compatibility.
#'
#' @param enzymes Raw enzyme input (NULL, character vector, or enzyme object list)
#' @param glycans Target glycan structures for pre-filtering (optional)
#' @param apply_prefilter Whether to apply enzyme pre-filtering based on target glycans
#' @returns A list of `glyenzy_enzyme` objects.
#' @noRd
.process_enzymes_arg <- function(enzymes, glycans = NULL, apply_prefilter = TRUE) {
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

  # Apply enzyme pre-filtering if requested and target glycans provided
  if (apply_prefilter && !is.null(glycans)) {
    can_contribute <- tryCatch({
      purrr::map_lgl(enzyme_names, ~ {
        # Check if enzyme can contribute to any target glycan
        any(purrr::map_lgl(glycans, function(target) {
          tryCatch(
            glyenzy::is_synthesized_by(target, .x),
            error = function(e) FALSE
          )
        }))
      })
    }, error = function(e) {
      rep(TRUE, length(enzyme_names))
    })

    enzyme_names <- enzyme_names[can_contribute]
    if (length(enzyme_names) == 0L) {
      cli::cli_abort("No enzymes are predicted to contribute to any target glycan.")
    }
  }

  unname(glyenzy_enzymes[enzyme_names])
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

#' Perform BFS synthesis search with common input processing
#'
#' This is a high-level wrapper that handles input validation, enzyme processing,
#' and BFS search execution. Used by both find_synthesis_path and rebuild_biosynthesis.
#'
#' @param from_g Starting glycan structure
#' @param to_gs Target glycan structures
#' @param enzymes List of `glyenzy_enzyme` objects to use
#' @param max_steps Maximum search steps
#' @param filter Optional filter function
#' @returns igraph object representing synthesis path(s)
#' @noRd
.perform_bfs_synthesis <- function(from_g, to_gs, enzymes, max_steps, filter = NULL) {
  # Parse glycan structures and compute keys
  from_key <- as.character(from_g)[1]
  to_keys <- as.character(to_gs)

  # Perform BFS search using core algorithm
  search_result <- bfs_synthesis_search(
    from_g = from_g,
    to_gs = to_gs,
    enzymes = enzymes,
    max_steps = max_steps,
    filter = filter,
    from_key = from_key,
    to_keys = to_keys
  )
  
  # Build and return result graph
  build_synthesis_result_graph(
    search_result,
    from_key,
    to_keys,
    max_steps
  )
}