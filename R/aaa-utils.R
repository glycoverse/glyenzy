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
  .have_motif_substituent_subset(
    x,
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
}

.process_glycan_arg <- function(
  x,
  allow_generic = FALSE,
  warn_non_intact = TRUE
) {
  if (length(x) == 1L) {
    .process_glycans_arg(
      x,
      allow_generic = allow_generic,
      warn_non_intact = warn_non_intact
    )
  } else {
    cli::cli_abort(c(
      "{.arg x} must have length 1.",
      "x" = "Got {.val {length(x)}}."
    ))
  }
}

.process_glycans_arg <- function(
  x,
  allow_generic = FALSE,
  warn_non_intact = TRUE
) {
  if (is.character(x)) {
    x <- glyparse::auto_parse(x)
  } else if (!glyrepr::is_glycan_structure(x)) {
    cli::cli_abort(c(
      "{.arg glycans} must be a {.cls glyrepr_structure} vector or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(x)}}."
    ))
  }

  mono_type <- glyrepr::get_mono_type(x)
  is_concrete <- .is_concrete_glycan(x)
  non_concrete <- !is_concrete & !is.na(mono_type)
  if (!allow_generic && any(non_concrete)) {
    cli::cli_abort(c(
      "All glycans must have concrete monosaccharides (e.g. Gal, GlcNAc, etc.).",
      "x" = "These glycans are not concrete: {.val {unique(x[non_concrete])}}."
    ))
  }

  has_unsupported_substituents <- .has_unsupported_substituents(x)
  has_unsupported_substituents[is.na(has_unsupported_substituents)] <- FALSE
  if (any(has_unsupported_substituents)) {
    cli::cli_abort(c(
      "Only sulfate substituents are supported.",
      "x" = "These glycans have unsupported substituents: {.val {unique(x[has_unsupported_substituents])}}.",
      "i" = "Use {.fn glyrepr::remove_substituents} to get clean glycans."
    ))
  }

  if (warn_non_intact) {
    .warn_non_intact_glycans(x)
  }
  x
}

#' Process glycans supplied to a biosynthesis tracing function
#'
#' @param x Glycan structures or parseable glycan strings.
#' @param call The calling environment for validation errors.
#' @returns A validated `glyrepr_structure` vector.
#' @noRd
.process_biosynthesis_glycans_arg <- function(
  x,
  call = rlang::caller_env()
) {
  x <- .process_glycans_arg(
    x,
    allow_generic = TRUE,
    warn_non_intact = FALSE
  )
  .validate_biosynthesis_glycans(x, call = call)
  .warn_non_intact_glycans(x)
  x
}

#' Process a biosynthesis path's starting and target glycans
#'
#' @param from Starting glycan structure or parseable glycan string.
#' @param to Target glycan structure or parseable glycan string.
#' @param call The calling environment for validation errors.
#' @returns A list containing validated `from` and `to` structures.
#' @noRd
.process_biosynthesis_path_args <- function(
  from,
  to,
  call = rlang::caller_env()
) {
  from <- .process_glycan_arg(
    from,
    allow_generic = TRUE,
    warn_non_intact = FALSE
  )
  to <- .process_glycan_arg(
    to,
    allow_generic = TRUE,
    warn_non_intact = FALSE
  )
  glycans <- c(from, to)
  .validate_biosynthesis_glycans(glycans, call = call)
  .warn_non_intact_glycans(glycans)
  list(from = from, to = to)
}

#' Validate glycans used in biosynthesis tracing
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param call The calling environment for validation errors.
#' @returns `glycans`, invisibly.
#' @noRd
.validate_biosynthesis_glycans <- function(
  glycans,
  call = rlang::caller_env()
) {
  mono_types <- unique(glyrepr::get_mono_type(glycans))
  structure_levels <- unique(.glycan_structure_levels(glycans))
  valid_mono_type <- length(mono_types) == 1L &&
    !is.na(mono_types) &&
    mono_types %in% c("concrete", "generic")
  valid_structure_level <- length(structure_levels) == 1L &&
    !is.na(structure_levels) &&
    structure_levels %in% c("intact", "topological")

  if (valid_mono_type && valid_structure_level) {
    return(invisible(glycans))
  }

  displayed_mono_types <- ifelse(is.na(mono_types), "NA", mono_types)
  displayed_structure_levels <- ifelse(
    is.na(structure_levels),
    "NA",
    structure_levels
  )
  cli::cli_abort(
    c(
      "Biosynthesis glycans must share one monosaccharide type and one structure level.",
      "x" = "Detected monosaccharide type(s): {.val {displayed_mono_types}}.",
      "x" = "Detected structure level(s): {.val {displayed_structure_levels}}.",
      "i" = "Supported monosaccharide types are {.val concrete} and {.val generic}; mixed-residue and missing glycans are not supported.",
      "i" = "Supported structure levels are {.val intact} and {.val topological}; partial structures are not supported.",
      "i" = "Use {.fn glyrepr::convert_to_generic} or {.fn glyrepr::remove_linkages} to standardize inputs."
    ),
    call = call
  )
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

#' Split a normalized substituent string into tokens
#'
#' @param x A scalar substituent string.
#' @returns A character vector of substituent tokens.
#' @noRd
.substituent_tokens <- function(x) {
  if (length(x) == 0L || is.na(x) || identical(x, "")) {
    return(character())
  }
  tokens <- strsplit(x, ",", fixed = TRUE)[[1]]
  tokens[tokens != ""]
}

#' Identify substituent tokens supported by glyenzy
#'
#' @param x A character vector of substituent tokens.
#' @returns A logical vector.
#' @noRd
.supported_substituents <- function(x) {
  grepl("^(\\?|[0-9]+)S$", x)
}

#' Detect glycans containing unsupported substituents
#'
#' @param x A `glyrepr_structure` vector.
#' @returns A logical vector.
#' @noRd
.has_unsupported_substituents <- function(x) {
  has_unsupported_single <- function(graph) {
    tokens <- unlist(
      purrr::map(igraph::V(graph)$sub, .substituent_tokens),
      use.names = FALSE
    )
    length(tokens) > 0L && !all(.supported_substituents(tokens))
  }
  glyrepr::smap_lgl(x, has_unsupported_single)
}

#' Check whether required substituents are contained in available ones
#'
#' Each required token must match a distinct available token. Unknown positions
#' are compatible with concrete positions only in lenient mode.
#'
#' @param available_sub A scalar substituent string on the target residue.
#' @param required_sub A scalar substituent string on the motif residue.
#' @param mode Matching mode.
#' @returns A logical scalar.
#' @noRd
.substituent_tokens_contained <- function(
  available_sub,
  required_sub,
  mode = c("strict", "lenient")
) {
  mode <- match.arg(mode)
  available <- .substituent_tokens(available_sub)
  required <- .substituent_tokens(required_sub)
  if (length(required) == 0L) {
    return(TRUE)
  }
  if (length(required) > length(available)) {
    return(FALSE)
  }

  token_matches <- function(available_token, required_token) {
    if (identical(available_token, required_token)) {
      return(TRUE)
    }
    if (!identical(mode, "lenient")) {
      return(FALSE)
    }

    available_parts <- regmatches(
      available_token,
      regexec("^(\\?|[0-9]+)(.+)$", available_token)
    )[[1]]
    required_parts <- regmatches(
      required_token,
      regexec("^(\\?|[0-9]+)(.+)$", required_token)
    )[[1]]
    if (length(available_parts) == 0L || length(required_parts) == 0L) {
      return(FALSE)
    }
    identical(available_parts[[3]], required_parts[[3]]) &&
      (identical(available_parts[[2]], "?") ||
        identical(required_parts[[2]], "?"))
  }

  candidates <- purrr::map(
    required,
    ~ which(purrr::map_lgl(available, token_matches, required_token = .x))
  )
  if (any(lengths(candidates) == 0L)) {
    return(FALSE)
  }
  candidates <- candidates[order(lengths(candidates))]
  assign_next <- function(i, used) {
    if (i > length(candidates)) {
      return(TRUE)
    }
    available_idx <- setdiff(candidates[[i]], used)
    any(purrr::map_lgl(
      available_idx,
      ~ assign_next(i + 1L, c(used, .x))
    ))
  }
  assign_next(1L, integer())
}

#' Match a motif graph with substituent-subset semantics
#'
#' Topology and residues are matched by glymotif after substituent attributes
#' are blanked. Candidate mappings are then filtered so every motif
#' substituent occurs on the corresponding glycan residue.
#'
#' @param glycan_graph Target glycan graph.
#' @param motif_graph Motif graph.
#' @param alignment Motif alignment.
#' @param ignore_linkages Whether to ignore linkage attributes.
#' @param match_degree Optional degree matching mode.
#' @param mode Matching mode.
#' @returns A list of integer node mappings.
#' @noRd
.g_match_motif_substituent_subset <- function(
  glycan_graph,
  motif_graph,
  alignment = "substructure",
  ignore_linkages = FALSE,
  match_degree = NULL,
  mode = c("strict", "lenient")
) {
  mode <- match.arg(mode)
  glycan_base <- igraph::set_vertex_attr(
    glycan_graph,
    "sub",
    value = rep("", igraph::vcount(glycan_graph))
  )
  motif_base <- igraph::set_vertex_attr(
    motif_graph,
    "sub",
    value = rep("", igraph::vcount(motif_graph))
  )
  matches <- glymotif::.g_match_motif(
    glycan_base,
    motif_base,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    strict_sub = TRUE,
    match_degree = match_degree,
    mode = mode
  )
  if (length(matches) == 0L) {
    return(matches)
  }

  glycan_subs <- igraph::vertex_attr(glycan_graph, "sub")
  motif_subs <- igraph::vertex_attr(motif_graph, "sub")
  keep <- purrr::map_lgl(matches, function(mapping) {
    purrr::every(seq_along(mapping), function(i) {
      .substituent_tokens_contained(
        glycan_subs[[mapping[[i]]]],
        motif_subs[[i]],
        mode = mode
      )
    })
  })
  matches[keep]
}

#' Match a motif using substituent-subset semantics
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param motif A scalar `glyrepr_structure` or parseable motif string.
#' @param alignment Motif alignment.
#' @param ignore_linkages Whether to ignore linkage attributes.
#' @param match_degree Optional degree matching mode.
#' @returns A nested list of node mappings.
#' @noRd
.match_motif_substituent_subset <- function(
  glycans,
  motif,
  alignment = "substructure",
  ignore_linkages = FALSE,
  match_degree = NULL
) {
  if (is.character(motif)) {
    motif <- glyparse::auto_parse(motif)
  }
  checkmate::assert_class(motif, "glyrepr_structure")
  checkmate::assert_true(length(motif) == 1L)
  mode <- .glymotif_mode(glycans)
  glycan_graphs <- glyrepr::get_structure_graphs(
    glycans,
    return_list = TRUE
  )
  motif_graph <- glyrepr::get_structure_graphs(
    motif,
    return_list = FALSE
  )
  purrr::map(
    glycan_graphs,
    .g_match_motif_substituent_subset,
    motif_graph = motif_graph,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    match_degree = match_degree,
    mode = mode
  )
}

#' Match motifs using substituent-subset semantics
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param motifs A `glyrepr_structure` vector.
#' @param alignments Motif alignments.
#' @returns A list indexed by motif, then glycan.
#' @noRd
.match_motifs_substituent_subset <- function(
  glycans,
  motifs,
  alignments = "substructure"
) {
  if (is.character(motifs)) {
    motifs <- glyparse::auto_parse(motifs)
  }
  alignments <- rep_len(alignments, length(motifs))
  res <- purrr::map2(
    seq_along(motifs),
    alignments,
    ~ .match_motif_substituent_subset(
      glycans,
      motifs[.x],
      alignment = .y
    )
  )
  names(res) <- names(motifs)
  res
}

#' Check for a motif using substituent-subset semantics
#'
#' @inheritParams .match_motif_substituent_subset
#' @returns A logical vector.
#' @noRd
.have_motif_substituent_subset <- function(
  glycans,
  motif,
  alignment = "substructure",
  ignore_linkages = FALSE,
  match_degree = NULL
) {
  lengths(.match_motif_substituent_subset(
    glycans,
    motif,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    match_degree = match_degree
  )) >
    0L
}

#' Check for motifs using substituent-subset semantics
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param motifs A `glyrepr_structure` vector.
#' @param alignments Motif alignments.
#' @returns A logical matrix with glycans in rows and motifs in columns.
#' @noRd
.have_motifs_substituent_subset <- function(
  glycans,
  motifs,
  alignments = "substructure"
) {
  matches <- .match_motifs_substituent_subset(glycans, motifs, alignments)
  if (length(matches) == 0L) {
    return(matrix(
      logical(),
      nrow = length(glycans),
      ncol = 0L
    ))
  }
  res <- do.call(
    cbind,
    purrr::map(matches, ~ lengths(.x) > 0L)
  )
  colnames(res) <- names(motifs)
  rownames(res) <- names(glycans)
  res
}

#' Count a motif using substituent-subset semantics
#'
#' @inheritParams .match_motif_substituent_subset
#' @returns An integer vector.
#' @noRd
.count_motif_substituent_subset <- function(
  glycans,
  motif,
  alignment = "substructure",
  ignore_linkages = FALSE,
  match_degree = NULL
) {
  lengths(.match_motif_substituent_subset(
    glycans,
    motif,
    alignment = alignment,
    ignore_linkages = ignore_linkages,
    match_degree = match_degree
  ))
}

#' Count motifs using substituent-subset semantics
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param motifs A `glyrepr_structure` vector.
#' @param alignments Motif alignments.
#' @returns An integer matrix with glycans in rows and motifs in columns.
#' @noRd
.count_motifs_substituent_subset <- function(
  glycans,
  motifs,
  alignments = "substructure"
) {
  matches <- .match_motifs_substituent_subset(glycans, motifs, alignments)
  if (length(matches) == 0L) {
    return(matrix(
      integer(),
      nrow = length(glycans),
      ncol = 0L
    ))
  }
  res <- do.call(cbind, purrr::map(matches, lengths))
  colnames(res) <- names(motifs)
  rownames(res) <- names(glycans)
  res
}

#' Evaluate OR-valued rule requirements
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param rule A `glyenzy_enzyme_rule`.
#' @returns A logical vector.
#' @noRd
.rule_requirements_met <- function(glycans, rule) {
  if (length(rule$requires) == 0L) {
    return(rep(TRUE, length(glycans)))
  }
  requirement_matches <- purrr::map(
    rule$requires,
    function(requirement) {
      .have_motif_substituent_subset(
        glycans,
        requirement$motif,
        alignment = requirement$alignment
      )
    }
  )
  Reduce(`|`, requirement_matches)
}

#' Check whether a sulfate position is available on a residue
#'
#' @param graph A glycan graph.
#' @param vertex A node index.
#' @param sulfate A sulfate token such as `"6S"`.
#' @returns A logical scalar.
#' @noRd
.sulfate_position_available <- function(graph, vertex, sulfate) {
  position <- sub("S$", "", sulfate)
  if (identical(position, "?")) {
    return(TRUE)
  }

  existing_substituents <- .substituent_tokens(
    igraph::vertex_attr(graph, "sub", index = vertex)
  )
  existing_positions <- sub("[^0-9?].*$", "", existing_substituents)
  if (position %in% existing_positions) {
    return(FALSE)
  }

  incoming_edges <- igraph::incident(graph, vertex, mode = "in")
  incoming_linkages <- igraph::edge_attr(
    graph,
    "linkage",
    index = incoming_edges
  )
  donor_positions <- sub("-.*$", "", incoming_linkages)
  donor_positions <- sub("^[ab?]", "", donor_positions)
  donor_positions <- donor_positions[
    grepl("^[0-9]+$", donor_positions)
  ]
  if (length(incoming_edges) == 0L) {
    root_anomer <- igraph::graph_attr(graph, "anomer")
    if (length(root_anomer) == 1L) {
      root_position <- sub("^[ab?]", "", root_anomer)
      if (grepl("^[0-9]+$", root_position)) {
        donor_positions <- c(donor_positions, root_position)
      }
    }
  }
  if (position %in% donor_positions) {
    return(FALSE)
  }

  outgoing_edges <- igraph::incident(graph, vertex, mode = "out")
  existing_linkages <- igraph::edge_attr(
    graph,
    "linkage",
    index = outgoing_edges
  )
  linkage_positions <- sub(".*-", "", existing_linkages)
  linkage_positions <- linkage_positions[
    linkage_positions != "?" &
      !grepl("/", linkage_positions, fixed = TRUE)
  ]
  !position %in% linkage_positions
}

#' Warn when glycan structures require lenient motif matching
#'
#' @param x A `glyrepr_structure` vector.
#' @returns `x`, invisibly.
#' @noRd
.warn_non_intact_glycans <- function(x) {
  structure_levels <- .glycan_structure_levels(x)
  mono_types <- glyrepr::get_mono_type(x)
  non_intact <- !is.na(structure_levels) & structure_levels != "intact"
  non_concrete <- !is.na(mono_types) & mono_types != "concrete"
  if (!any(non_intact | non_concrete)) {
    return(invisible(x))
  }

  detected_levels <- unique(structure_levels[non_intact])
  detected_types <- unique(mono_types[non_concrete])
  details <- character()
  if (length(detected_levels) > 0L) {
    details <- c(
      details,
      "i" = "Detected structure level(s): {.val {detected_levels}}."
    )
  }
  if (length(detected_types) > 0L) {
    details <- c(
      details,
      "i" = "Detected monosaccharide type(s): {.val {detected_types}}."
    )
  }

  cli::cli_warn(c(
    "Using lenient motif matching for non-concrete or non-intact glycan structures.",
    details,
    "i" = "Results may be less reliable when exact residue identities, linkages, or anomers are missing or ambiguous."
  ))
  invisible(x)
}

#' Get element-wise structure levels for a glycan vector
#'
#' @param x A `glyrepr_structure` vector.
#' @returns A character vector with one level per structure.
#' @noRd
.glycan_structure_levels <- function(x) {
  suppressWarnings(glyrepr::get_structure_level(x))
}

#' Summarize structure levels for algorithms that require one level
#'
#' @param x A `glyrepr_structure` vector.
#' @returns A scalar structure level, or `NA_character_`.
#' @noRd
.glycan_structure_level <- function(x) {
  levels <- unique(.glycan_structure_levels(x))
  levels <- levels[!is.na(levels)]
  if (length(levels) == 0L) {
    return(NA_character_)
  }
  if (length(levels) > 1L) {
    cli::cli_abort("Glycans must have the same structure level.")
  }
  levels[[1]]
}

#' Get the shared monosaccharide type of a glycan vector
#'
#' @param x A `glyrepr_structure` vector.
#' @returns A scalar monosaccharide type, or `NA_character_`.
#' @noRd
.glycan_mono_type <- function(x) {
  mono_types <- unique(glyrepr::get_mono_type(x))
  mono_types <- mono_types[!is.na(mono_types)]
  if (length(mono_types) == 0L) {
    return(NA_character_)
  }
  if (length(mono_types) > 1L) {
    cli::cli_abort("Glycans must have the same monosaccharide type.")
  }
  mono_types[[1]]
}

#' Check whether a structure level is intact
#'
#' @param x A scalar character structure level.
#' @returns A logical scalar.
#' @noRd
.is_intact_structure_level <- function(x) {
  x <- x[!is.na(x)]
  length(x) == 0L || all(x == "intact")
}

#' Select the motif matching mode for a glycan vector
#'
#' @param glycans A `glyrepr_structure` vector.
#' @returns `"strict"` for intact structures and `"lenient"` otherwise.
#' @noRd
.glymotif_mode <- function(glycans) {
  structure_levels <- .glycan_structure_levels(glycans)
  mono_types <- glyrepr::get_mono_type(glycans)
  if (
    .is_intact_structure_level(structure_levels) &&
      length(mono_types) > 0L &&
      all(!is.na(mono_types) & mono_types == "concrete")
  ) {
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
  checkmate::assert_choice(structure_level, c("intact", "topological"))
  structure_level
}

#' Validate that output structure level does not reduce input resolution
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param structure_level Requested structure level.
#' @returns `structure_level`, invisibly.
#' @noRd
.validate_output_structure_level <- function(glycans, structure_level) {
  input_structure_levels <- .glycan_structure_levels(glycans)
  input_structure_levels <- input_structure_levels[
    !is.na(input_structure_levels)
  ]
  if (length(input_structure_levels) == 0L) {
    return(invisible(structure_level))
  }

  if (
    any(
      .structure_level_rank(structure_level) <
        .structure_level_rank(input_structure_levels)
    )
  ) {
    input_structure_level <- unique(input_structure_levels)
    cli::cli_abort(c(
      "{.arg structure_level} cannot be lower than the input glycan structure level.",
      "x" = "Input level(s): {.val {input_structure_level}}; requested: {.val {structure_level}}."
    ))
  }

  invisible(structure_level)
}

#' Rank glycan structure levels by information resolution
#'
#' @param structure_level A scalar structure level.
#' @returns An integer rank where larger values preserve more information.
#' @noRd
.structure_level_rank <- function(structure_level) {
  ranks <- c(
    topological = 1L,
    partial = 2L,
    intact = 3L
  )
  unname(ranks[structure_level])
}

#' Remove linkages from glycan structures for topological searches
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param structure_level Requested structure level.
#' @returns A `glyrepr_structure` vector.
#' @noRd
.remove_linkages_for_level <- function(glycans, structure_level) {
  structure_level <- .validate_structure_level(structure_level)
  if (identical(structure_level, "intact")) {
    return(glycans)
  }

  glyrepr::remove_linkages(glycans)
}

# Reduce trusted product graphs before constructing a glycan structure vector.
.prepare_valid_glycan_graphs <- function(graphs, structure_level) {
  structure_level <- .validate_structure_level(structure_level)
  if (identical(structure_level, "intact")) {
    return(graphs)
  }

  purrr::map(graphs, glyrepr::remove_linkages)
}

# Put the reducing-end root in the position expected by glymotif's core
# matcher without paying for full IUPAC branch ordering. The graph is trusted
# to be a valid out-tree, so it has exactly one in-degree-zero vertex.
.move_glycan_root_last <- function(graph) {
  n_vertices <- igraph::vcount(graph)
  roots <- which(igraph::degree(graph, mode = "in") == 0L)
  if (length(roots) != 1L) {
    glyrepr::validate_glycan_graph(graph)
    cli::cli_abort("Internal error: a glycan graph must have exactly one root.")
  }
  root <- roots[[1]]
  if (root == n_vertices) {
    return(graph)
  }

  other_vertices <- setdiff(seq_len(n_vertices), root)
  permutation <- integer(n_vertices)
  permutation[other_vertices] <- seq_along(other_vertices)
  permutation[root] <- n_vertices
  igraph::permute(graph, permutation)
}

# Canonicalize trusted graphs and generate their authoritative IUPAC keys.
.canonicalize_valid_glycan_graphs <- function(graphs) {
  input_names <- names(graphs)
  graphs <- unname(graphs)
  graphs <- purrr::map(graphs, glyrepr::canonicalize_glycan_graph)
  keys <- purrr::map_chr(graphs, glyrepr::graph_to_iupac)
  names(keys) <- input_names
  list(keys = keys, graphs = graphs)
}

# Assemble valid, mutually compatible graphs without repeating semantic
# validation. Graph canonicalization, key generation, and lookup deduplication
# remain required representation steps.
.new_glycan_structure_from_valid_graphs <- function(graphs) {
  if (length(graphs) == 0L) {
    return(glyrepr::new_glycan_structure())
  }

  keyed <- .canonicalize_valid_glycan_graphs(graphs)
  unique_graphs <- !duplicated(unname(keyed$keys))
  graph_lookup <- keyed$graphs[unique_graphs]
  names(graph_lookup) <- unname(keyed$keys[unique_graphs])

  glyrepr::new_glycan_structure(keyed$keys, graph_lookup)
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
      .enzyme_contributes_to_targets,
      targets = glycans
    ),
    error = function(e) {
      rep(TRUE, length(enzymes))
    }
  )
}

# Check all targets together so vectorized motif preparation is shared.
.enzyme_contributes_to_targets <- function(enzyme, targets) {
  vectorized <- tryCatch(
    .have_enzyme_motif(targets, enzyme),
    error = function(e) NULL
  )
  if (is.logical(vectorized) && length(vectorized) == length(targets)) {
    return(any(vectorized))
  }

  any(purrr::map_lgl(
    targets,
    .enzyme_contributes_to_target,
    enzyme = enzyme
  ))
}

#' Check whether one enzyme can contribute to one target glycan
#'
#' @param target A target glycan structure.
#' @param enzyme A `glyenzy_enzyme` object.
#' @returns `TRUE` when the enzyme can synthesize `target`; otherwise `FALSE`.
#' @noRd
.enzyme_contributes_to_target <- function(target, enzyme) {
  if (!.enzyme_glycan_type_mask(target, enzyme)) {
    return(FALSE)
  }
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
#' @param max_virtual_steps Maximum number of target-directed virtual steps.
#' @returns igraph object representing synthesis path(s)
#' @noRd
.perform_bfs_synthesis <- function(
  from_g,
  to_gs,
  enzymes,
  max_steps,
  filter = NULL,
  max_virtual_steps = 0L
) {
  target_structure_level <- .glycan_structure_level(to_gs)
  search_structure_level <- .bfs_search_structure_level(target_structure_level)
  target_match <- .bfs_target_match(to_gs)

  from_g <- .remove_linkages_for_level(from_g, search_structure_level)
  to_gs <- .remove_linkages_for_level(to_gs, search_structure_level)

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
    target_match = target_match,
    allow_partial = max_virtual_steps > 0L
  )

  if (length(search_result$missing_target_keys) > 0L) {
    return(.perform_virtual_fallback_synthesis(
      from_g = from_g,
      to_gs = to_gs,
      enzymes = enzymes,
      max_steps = max_steps,
      filter = filter,
      max_virtual_steps = max_virtual_steps,
      structure_level = search_structure_level,
      target_match = target_match,
      strict_result = search_result
    ))
  }

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
  "intact"
}

#' Choose the BFS target matching strategy from target structures
#'
#' @param target_glycans Target glycan structures.
#' @returns A target matching strategy.
#' @noRd
.bfs_target_match <- function(target_glycans) {
  if (!identical(.glycan_mono_type(target_glycans), "concrete")) {
    return("whole")
  }
  "key"
}

# Return only targets that are not already satisfied by the starting glycan.
.bfs_targets_requiring_synthesis <- function(from_g, to_gs) {
  target_structure_level <- .glycan_structure_level(to_gs)
  search_structure_level <- .bfs_search_structure_level(target_structure_level)
  target_match <- .bfs_target_match(to_gs)

  from_g <- .remove_linkages_for_level(from_g, search_structure_level)
  to_gs <- .remove_linkages_for_level(to_gs, search_structure_level)

  if (identical(target_match, "key")) {
    matched <- as.character(to_gs) == as.character(from_g)[1]
  } else {
    from_graph <- glyrepr::get_structure_graphs(from_g)
    target_graphs <- glyrepr::get_structure_graphs(
      to_gs,
      return_list = TRUE
    )
    matched <- purrr::map_lgl(
      target_graphs,
      function(target_graph) {
        glymotif::.g_have_motif(
          from_graph,
          target_graph,
          alignment = "whole",
          mode = "lenient"
        )
      }
    )
  }

  to_gs[!matched]
}
