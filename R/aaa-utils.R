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
  .have_motif(x, "GlcNAc(b1-4)GlcNAc(b1-")
}

.process_glycan_arg <- function(x, allow_generic = FALSE) {
  if (length(x) == 1L) {
    .process_glycans_arg(x, allow_generic = allow_generic)
  } else {
    cli::cli_abort(c(
      "{.arg x} must have length 1.",
      "x" = "Got {.val {length(x)}}."
    ))
  }
}

.process_glycans_arg <- function(x, allow_generic = FALSE) {
  if (is.character(x)) {
    x <- glyparse::auto_parse(x)
  } else if (!glyrepr::is_glycan_structure(x)) {
    cli::cli_abort(c(
      "{.arg glycans} must be a {.cls glyrepr_structure} vector or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(x)}}."
    ))
  }

  is_concrete <- .is_concrete_glycan(x)
  if (!allow_generic && !all(is_concrete)) {
    cli::cli_abort(c(
      "All glycans must have concrete monosaccharides (e.g. Gal, GlcNAc, etc.).",
      "x" = "These glycans are not concrete: {.val {unique(x[!is_concrete])}}."
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

  .warn_non_intact_glycans(x)
  x
}

.is_concrete_glycan <- function(x) {
  glyrepr::get_mono_type(x) == "concrete"
}

.has_substituents <- function(x) {
  has_sub_single <- function(graph) {
    purrr::some(igraph::V(graph)$sub, ~ .x != "")
  }
  glyrepr::smap_lgl(x, has_sub_single)
}

#' Warn when glycan structures require lenient motif matching
#'
#' @param x A `glyrepr_structure` vector.
#' @returns `x`, invisibly.
#' @noRd
.warn_non_intact_glycans <- function(x) {
  structure_level <- .glycan_structure_level(x)
  if (.is_intact_structure_level(structure_level)) {
    return(invisible(x))
  }

  cli::cli_warn(c(
    "Using lenient motif matching for non-intact glycan structures.",
    "i" = "Detected structure level: {.val {structure_level}}.",
    "i" = "Results may be less reliable when exact linkages are missing or ambiguous."
  ))
  invisible(x)
}

#' Get the structure level of a glycan vector
#'
#' @param x A `glyrepr_structure` vector.
#' @returns A scalar character structure level, or `NA_character_`.
#' @noRd
.glycan_structure_level <- function(x) {
  suppressWarnings(glyrepr::get_structure_level(x))
}

#' Check whether a structure level is intact
#'
#' @param x A scalar character structure level.
#' @returns A logical scalar.
#' @noRd
.is_intact_structure_level <- function(x) {
  is.na(x) || identical(x, "intact")
}

#' Select the motif matching mode for a glycan vector
#'
#' @param glycans A `glyrepr_structure` vector.
#' @returns `"strict"` for intact structures and `"lenient"` otherwise.
#' @noRd
.glymotif_mode <- function(glycans) {
  if (.is_intact_structure_level(.glycan_structure_level(glycans))) {
    return("strict")
  }
  "lenient"
}

#' Validate a requested output structure level
#'
#' @param structure_level Requested structure level.
#' @returns A scalar character structure level.
#' @noRd
.validate_structure_level <- function(structure_level) {
  checkmate::assert_choice(structure_level, c("intact", "topological", "basic"))
  structure_level
}

#' Reduce glycan structures to a requested output level
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param structure_level Requested structure level.
#' @returns A `glyrepr_structure` vector.
#' @noRd
.reduce_structure_level <- function(glycans, structure_level) {
  structure_level <- .validate_structure_level(structure_level)
  if (identical(structure_level, "intact")) {
    return(glycans)
  }

  suppressWarnings(glyrepr::reduce_structure_level(glycans, structure_level))
}

#' Match a motif using the glycan-appropriate glymotif mode
#'
#' @inheritParams glymotif::match_motif
#' @noRd
.match_motif <- function(
  glycans,
  motif,
  ...,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  glymotif::match_motif(
    glycans,
    motif,
    ...,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = .glymotif_mode(glycans)
  )
}

#' Match motifs using the glycan-appropriate glymotif mode
#'
#' @inheritParams glymotif::match_motifs
#' @noRd
.match_motifs <- function(
  glycans,
  motifs,
  ...,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  glymotif::match_motifs(
    glycans,
    motifs,
    ...,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = .glymotif_mode(glycans)
  )
}

#' Check for a motif using the glycan-appropriate glymotif mode
#'
#' @inheritParams glymotif::have_motif
#' @noRd
.have_motif <- function(
  glycans,
  motif,
  ...,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  glymotif::have_motif(
    glycans,
    motif,
    ...,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = .glymotif_mode(glycans)
  )
}

#' Check for motifs using the glycan-appropriate glymotif mode
#'
#' @inheritParams glymotif::have_motifs
#' @noRd
.have_motifs <- function(
  glycans,
  motifs,
  ...,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  glymotif::have_motifs(
    glycans,
    motifs,
    ...,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = .glymotif_mode(glycans)
  )
}

#' Count motif matches using the glycan-appropriate glymotif mode
#'
#' @inheritParams glymotif::count_motif
#' @noRd
.count_motif <- function(
  glycans,
  motif,
  ...,
  alignment = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  glymotif::count_motif(
    glycans,
    motif,
    ...,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = .glymotif_mode(glycans)
  )
}

#' Count motif matches using the glycan-appropriate glymotif mode
#'
#' @inheritParams glymotif::count_motifs
#' @noRd
.count_motifs <- function(
  glycans,
  motifs,
  ...,
  alignments = NULL,
  ignore_linkages = FALSE,
  strict_sub = TRUE,
  match_degree = NULL
) {
  glymotif::count_motifs(
    glycans,
    motifs,
    ...,
    alignments = alignments,
    ignore_linkages = ignore_linkages,
    strict_sub = strict_sub,
    match_degree = match_degree,
    mode = .glymotif_mode(glycans)
  )
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
  enzymes <- .enzymes_from_arg(enzymes)
  enzymes <- .prefilter_enzymes(
    enzymes,
    glycans,
    apply_prefilter
  )

  unname(enzymes)
}

#' Get enzymes from supported input forms
#'
#' @param enzymes Raw enzyme input (`NULL`, character vector, or enzyme object
#'   list).
#' @returns A named list of `glyenzy_enzyme` objects.
#' @noRd
.enzymes_from_arg <- function(enzymes) {
  if (is.null(enzymes)) {
    to_keep <- !purrr::map_lgl(
      glyenzy_enzymes,
      ~ .is_starter_gt(.x) || .is_npre_gt(.x)
    )
    return(glyenzy_enzymes[to_keep])
  }

  if (is.character(enzymes)) {
    .check_known_enzyme_names(enzymes)
    return(glyenzy_enzymes[enzymes])
  }

  .enzymes_from_list(enzymes)
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

#' Normalize a list of enzyme objects
#'
#' @param enzymes A list of `glyenzy_enzyme` objects.
#' @returns A named list of `glyenzy_enzyme` objects.
#' @noRd
.enzymes_from_list <- function(enzymes) {
  checkmate::assert_list(enzymes, types = "glyenzy_enzyme")
  names(enzymes) <- purrr::map_chr(enzymes, "name")
  enzymes
}

#' Optionally prefilter enzymes against target glycans
#'
#' @param enzymes A named list of `glyenzy_enzyme` objects.
#' @param glycans Target glycan structures for pre-filtering.
#' @param apply_prefilter Whether to apply enzyme pre-filtering based on target
#'   glycans.
#' @returns A named list of `glyenzy_enzyme` objects.
#' @noRd
.prefilter_enzymes <- function(
  enzymes,
  glycans,
  apply_prefilter
) {
  if (!apply_prefilter || is.null(glycans)) {
    return(enzymes)
  }

  can_contribute <- .can_enzymes_contribute(enzymes, glycans)
  enzymes <- enzymes[can_contribute]
  if (length(enzymes) == 0L) {
    cli::cli_abort(
      "No enzymes are predicted to contribute to any target glycan."
    )
  }

  enzymes
}

#' Check whether enzymes can contribute to any target glycan
#'
#' @param enzymes A named list of `glyenzy_enzyme` objects.
#' @param glycans Target glycan structures.
#' @returns A logical vector with one value per enzyme.
#' @noRd
.can_enzymes_contribute <- function(enzymes, glycans) {
  tryCatch(
    purrr::map_lgl(
      enzymes,
      ~ any(purrr::map_lgl(
        glycans,
        .enzyme_contributes_to_target,
        enzyme = .x
      ))
    ),
    error = function(e) {
      rep(TRUE, length(enzymes))
    }
  )
}

#' Check whether one enzyme can contribute to one target glycan
#'
#' @param target A target glycan structure.
#' @param enzyme A `glyenzy_enzyme` object.
#' @returns `TRUE` when the enzyme can synthesize `target`; otherwise `FALSE`.
#' @noRd
.enzyme_contributes_to_target <- function(target, enzyme) {
  tryCatch(
    .have_enzyme_motif(target, enzyme),
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
#' @param enzymes A list of `glyenzy_enzyme` objects, or `NULL` to use defaults.
#'
#' @returns A character vector of enzyme names from graph edges.
#' @noRd
.trace_enzyme_edges_single <- function(glycan, enzymes = NULL) {
  path <- trace_biosynthesis(glycan, enzymes = enzymes)
  igraph::E(path)$enzyme
}

#' Extract trace-derived enzyme names for each glycan independently
#'
#' @param glycans A `glyrepr_structure` vector.
#'
#' @param enzymes A list of `glyenzy_enzyme` objects, or `NULL` to use defaults.
#'
#' @returns A list of character vectors.
#' @noRd
.trace_enzyme_edges <- function(glycans, enzymes = NULL) {
  purrr::map(glycans, .trace_enzyme_edges_single, enzymes = enzymes)
}

#' Merge one enzyme into the default trace enzyme set
#'
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @returns A list of `glyenzy_enzyme` objects.
#' @noRd
.trace_enzymes_with <- function(enzyme) {
  to_keep <- !purrr::map_lgl(
    glyenzy_enzymes,
    ~ .is_starter_gt(.x) || .is_npre_gt(.x)
  )
  enzymes <- glyenzy_enzymes[to_keep]
  enzymes[[enzyme$name]] <- enzyme
  unname(enzymes)
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
  target_structure_level <- .glycan_structure_level(to_gs)
  search_structure_level <- .bfs_search_structure_level(target_structure_level)
  target_match <- .bfs_target_match(target_structure_level)

  from_g <- .reduce_structure_level(from_g, search_structure_level)
  to_gs <- .reduce_structure_level(to_gs, search_structure_level)

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
    to_keys = to_keys,
    structure_level = search_structure_level,
    target_match = target_match
  )

  # Build and return result graph
  build_synthesis_result_graph(
    search_result,
    from_key,
    unique(search_result$found_keys),
    max_steps
  )
}

#' Choose the BFS product structure level from target structures
#'
#' @param target_structure_level Target glycan structure level.
#' @returns A structure level accepted by `apply_enzyme()`.
#' @noRd
.bfs_search_structure_level <- function(target_structure_level) {
  if (identical(target_structure_level, "topological")) {
    return("topological")
  }
  if (identical(target_structure_level, "basic")) {
    return("basic")
  }
  "intact"
}

#' Choose the BFS target matching strategy from target structures
#'
#' @param target_structure_level Target glycan structure level.
#' @returns A target matching strategy.
#' @noRd
.bfs_target_match <- function(target_structure_level) {
  if (identical(target_structure_level, "partial")) {
    return("whole")
  }
  "key"
}
