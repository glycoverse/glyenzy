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
    } else {
      # type is integer
      res <- rep(0, length(glycans))
    }
    if (any(is_n)) {
      res[is_n] <- .f(glycans[is_n], enzyme)
    }
    res
  }
}

.is_n_glycan <- function(x) {
  glymotif::have_motif(x, "GlcNAc(b1-4)GlcNAc(b1-")
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

#' Get product alignment for final-product enzyme matching
#'
#' @param rule A `glyenzy_enzyme_rule` object.
#'
#' @returns An alignment string.
#' @noRd
.product_alignment <- function(rule) {
  if (rule$acceptor_alignment %in% c("whole", "core")) {
    return("core")
  }
  "substructure"
}

#' Process and validate enzyme list for synthesis search
#'
#' This function processes various enzyme input formats (NULL, character vector,
#' enzyme object list) and applies pre-filtering based on target glycan compatibility.
#'
#' @param enzymes Raw enzyme input (NULL, character vector, or enzyme object list).
#' @param glycans Target glycan structures for pre-filtering (optional).
#' @param apply_prefilter Whether to apply enzyme pre-filtering based on target glycans.
#' @returns A list of `glyenzy_enzyme` objects.
#' @noRd
.process_enzymes_arg <- function(
  enzymes,
  glycans = NULL,
  apply_prefilter = TRUE
) {
  enzyme_names <- .enzyme_names_from_arg(enzymes)
  enzyme_names <- .prefilter_enzyme_names(
    enzyme_names,
    glycans,
    apply_prefilter
  )

  unname(glyenzy_enzymes[enzyme_names])
}

#' Get enzyme names from supported enzyme input forms
#'
#' @param enzymes Raw enzyme input (`NULL`, character vector, or enzyme object
#'   list).
#' @returns A character vector of enzyme names.
#' @noRd
.enzyme_names_from_arg <- function(enzymes) {
  if (is.null(enzymes)) {
    to_keep <- !purrr::map_lgl(
      glyenzy_enzymes,
      ~ .is_starter_gt(.x) || .is_npre_gt(.x)
    )
    return(names(glyenzy_enzymes)[to_keep])
  }

  if (is.character(enzymes)) {
    .check_known_enzyme_names(enzymes)
    return(enzymes)
  }

  .enzyme_names_from_list(enzymes)
}

#' Validate enzyme names against the built-in enzyme table
#'
#' @param enzyme_names A character vector of enzyme names.
#' @returns The input `enzyme_names`, invisibly.
#' @noRd
.check_known_enzyme_names <- function(enzyme_names) {
  unknown <- setdiff(enzyme_names, names(glyenzy_enzymes))
  if (length(unknown) > 0) {
    cli::cli_abort("Unknown enzymes: {.val {unknown}}.")
  }

  invisible(enzyme_names)
}

#' Extract names from a list of enzyme objects
#'
#' @param enzymes A list of `glyenzy_enzyme` objects.
#' @returns A character vector of enzyme names.
#' @noRd
.enzyme_names_from_list <- function(enzymes) {
  checkmate::assert_list(enzymes, types = "glyenzy_enzyme")
  purrr::map_chr(enzymes, "name")
}

#' Optionally prefilter enzyme names against target glycans
#'
#' @param enzyme_names A character vector of enzyme names.
#' @param glycans Target glycan structures for pre-filtering.
#' @param apply_prefilter Whether to apply enzyme pre-filtering based on target
#'   glycans.
#' @returns A character vector of enzyme names.
#' @noRd
.prefilter_enzyme_names <- function(
  enzyme_names,
  glycans,
  apply_prefilter
) {
  if (!apply_prefilter || is.null(glycans)) {
    return(enzyme_names)
  }

  can_contribute <- .can_enzymes_contribute(enzyme_names, glycans)
  enzyme_names <- enzyme_names[can_contribute]
  if (length(enzyme_names) == 0L) {
    cli::cli_abort(
      "No enzymes are predicted to contribute to any target glycan."
    )
  }

  enzyme_names
}

#' Check whether enzymes can contribute to any target glycan
#'
#' @param enzyme_names A character vector of enzyme names.
#' @param glycans Target glycan structures.
#' @returns A logical vector with one value per enzyme name.
#' @noRd
.can_enzymes_contribute <- function(enzyme_names, glycans) {
  tryCatch(
    purrr::map_lgl(
      enzyme_names,
      ~ any(purrr::map_lgl(
        glycans,
        .enzyme_contributes_to_target,
        enzyme_name = .x
      ))
    ),
    error = function(e) {
      rep(TRUE, length(enzyme_names))
    }
  )
}

#' Check whether one enzyme can contribute to one target glycan
#'
#' @param target A target glycan structure.
#' @param enzyme_name A single enzyme name.
#' @returns `TRUE` when the enzyme can synthesize `target`; otherwise `FALSE`.
#' @noRd
.enzyme_contributes_to_target <- function(target, enzyme_name) {
  tryCatch(
    glyenzy::have_enzyme(target, enzyme_name),
    error = function(e) FALSE
  )
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

#' Extract enzyme names from one trace-derived synthesis graph
#'
#' @param glycan A length-one `glyrepr_structure` vector.
#'
#' @returns A character vector of enzyme names from graph edges.
#' @noRd
.trace_enzyme_edges_single <- function(glycan) {
  path <- trace_biosynthesis(glycan)
  igraph::E(path)$enzyme
}

#' Extract trace-derived enzyme names for each glycan independently
#'
#' @param glycans A `glyrepr_structure` vector.
#'
#' @returns A list of character vectors.
#' @noRd
.trace_enzyme_edges <- function(glycans) {
  purrr::map(glycans, .trace_enzyme_edges_single)
}

#' Perform BFS synthesis search with common input processing
#'
#' This is a high-level wrapper that handles input validation, enzyme processing,
#' and BFS search execution. Used by both path_biosynthesis and trace_biosynthesis.
#'
#' @param from_g Starting glycan structure
#' @param to_gs Target glycan structures
#' @param enzymes List of `glyenzy_enzyme` objects to use
#' @param max_steps Maximum search steps
#' @param filter Optional filter function
#' @returns igraph object representing synthesis path(s)
#' @noRd
.perform_bfs_synthesis <- function(
  from_g,
  to_gs,
  enzymes,
  max_steps,
  filter = NULL
) {
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
