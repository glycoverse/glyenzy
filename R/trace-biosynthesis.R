#' Trace the Biosynthetic Path of Glycans
#'
#' Reconstruct biosynthetic pathways for one or more glycans. The default
#' method uses known enzyme rules in a forward, multi-target breadth-first
#' search. To infer structure-driven paths without enzyme specificity, use
#' [trace_biosynthesis_virtual()].
#'
#' @inheritSection have_enzyme Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()] vector, or a character vector
#'   of strings supported by [glyparse::auto_parse()]. Can also be a single glycan.
#'   If multiple glycans are provided, the starting structure will be decided by the first glycan.
#'   Therefore, please make sure `glycans` are not of mixed glycan types.
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()]
#'   objects. If `NULL` (default), all available enzymes will be used.
#' @param max_steps Integer, maximum number of enzymatic steps to search.
#'   Default is 20.
#' @param filter Optional function to filter generated glycans at each step.
#'   Should take a [glyrepr::glycan_structure()] vector as input and return
#'   a logical vector of the same length. It filters generated products.
#'
#' @returns An [igraph::igraph()] object representing the synthesis path(s).
#'   Vertices represent glycan structures, with IUPAC-condensed strings in the
#'   `name` attribute. Every edge has a `step` attribute indicating the forward
#'   synthesis step and an `enzyme` attribute containing its gene symbol.
#'   Multiple enzymes catalysing the same substrate-to-product transition are
#'   represented by parallel edges.
#'
#'   For multiple targets, the graph includes all synthesis paths needed to
#'   reach every target glycan.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Rebuild the biosynthetic pathway of a single glycan
#' glycan <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' path <- trace_biosynthesis(glycan, max_steps = 20)
#'
#' # Rebuild pathways for multiple glycans
#' glycans <- c(
#'   "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-",
#'   "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' )
#' path <- trace_biosynthesis(glycans, max_steps = 20)
#'
#' # View the path
#' igraph::as_data_frame(path, what = "edges")
#'
#' @export
trace_biosynthesis <- function(
  glycans,
  enzymes = NULL,
  max_steps = 20,
  filter = NULL
) {
  # Parse and validate basic inputs first
  glycans <- .process_glycans_arg(glycans, allow_generic = TRUE)
  checkmate::assert_int(max_steps, lower = 1)
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }

  enzymes <- .process_enzymes_arg(
    enzymes,
    glycans = glycans,
    apply_prefilter = TRUE
  )

  # Find all possible paths using unified BFS logic
  starting_glycan <- .decide_starting_glycan(glycans[1])
  .perform_bfs_synthesis(starting_glycan, glycans, enzymes, max_steps, filter)
}

.decide_starting_glycan <- function(glycan) {
  if (.can_reliably_detect_n_glycan(glycan) && .is_n_glycan(glycan)) {
    start <- .n_glycan_starting_glycan()
  } else if (.have_motif(glycan, "GalNAc(a1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("GalNAc(a1-")
  } else if (.have_motif(glycan, "GlcNAc(b1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("GlcNAc(b1-")
  } else if (.have_motif(glycan, "Man(a1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("Man(a1-")
  } else if (.have_motif(glycan, "Fuc(a1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("Fuc(a1-")
  } else if (.have_motif(glycan, "Glc(b1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("Glc(b1-")
  } else {
    cli::cli_abort(c(
      "Cannot decide the starting point for the given glycan(s).",
      "i" = "Currenly, only N-glycans and O-glycans are supported."
    ))
  }
  start
}

.decide_virtual_starting_glycan <- function(glycan) {
  graph <- glyrepr::get_structure_graphs(glycan)
  core <- .n_glycan_starting_glycan(
    "virtual",
    .glycan_structure_level(glycan)
  )
  core_matches <- .virtual_core_matches(glycan, core)

  if (length(core_matches) > 0L) {
    start_graph <- igraph::induced_subgraph(graph, core_matches[[1]])
  } else {
    root <- which(igraph::degree(graph, mode = "in") == 0L)
    start_graph <- igraph::induced_subgraph(graph, root)
  }

  .new_glycan_structure_from_valid_graphs(list(start_graph))
}

.n_glycan_starting_glycan <- function(
  method = c("enzymatic", "virtual"),
  structure_level = "intact"
) {
  method <- match.arg(method)
  checkmate::assert_choice(
    structure_level,
    c("intact", "partial", "topological", "basic")
  )

  if (identical(method, "enzymatic")) {
    return(glyparse::parse_iupac_condensed(
      "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ))
  }

  if (identical(structure_level, "basic")) {
    return(glyrepr::n_glycan_core(
      linkage = FALSE,
      mono_type = "generic"
    ))
  }
  if (identical(structure_level, "topological")) {
    return(glyrepr::n_glycan_core(linkage = FALSE))
  }
  glyrepr::n_glycan_core()
}

.virtual_core_matches <- function(glycan, core) {
  structure_level <- .glycan_structure_level(glycan)
  mode <- if (identical(structure_level, "partial")) {
    "lenient"
  } else {
    "strict"
  }

  tryCatch(
    glymotif::match_motif(
      glycan,
      core,
      alignment = "core",
      mode = mode
    )[[1]],
    error = function(e) list()
  )
}

.can_reliably_detect_n_glycan <- function(glycan) {
  !identical(.glycan_structure_level(glycan), "basic")
}
