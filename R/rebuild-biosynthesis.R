#' Rebuild the Biosynthetic Path of Glycans
#'
#' Reconstruct the biosynthetic pathway for one or more glycans using enzymatic reactions.
#' This function uses a multi-target breadth-first search to find all feasible pathways
#' that can synthesize all the target glycans.
#'
#' For N-glycans, the starting structure is assumed to be "Glc(3)Man(9)GlcNAc(2)",
#' the N-glycan precursor transfered to Asn by OST.
#' For O-glycans, the starting structure is assumed to be "GalNAc(a1-".
#'
#' @inheritSection is_synthesized_by Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()] vector, or a character vector
#'   of strings supported by [glyparse::auto_parse()]. Can also be a single glycan.
#'   If multiple glycans are provided, the starting structure will be decided by the first glycan.
#'   Therefore, please make sure `glycans` are not a mixed vector of N- and O-glycans.
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()] objects.
#'   If `NULL` (default), all available enzymes will be used.
#' @param max_steps Integer, maximum number of enzymatic steps to search.
#'   Default is 20.
#' @param filter Optional function to filter generated glycans at each step.
#'   Should take a [glyrepr::glycan_structure()] vector as input and return
#'   a logical vector of the same length.
#'   It will be applied to all the generated glycans at each BFS step for pruning.
#'
#' @returns An [igraph::igraph()] object representing the synthesis path(s).
#'   Vertices represent glycan structures with `name` attribute containing
#'   IUPAC-condensed strings. Edges represent enzymatic reactions with
#'   `enzyme` attribute containing gene symbols and `step` attribute
#'   indicating the step number. For multiple targets, the graph includes
#'   all synthesis paths needed to reach every target glycan.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Rebuild the biosynthetic pathway of a single glycan
#' glycan <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' path <- rebuild_biosynthesis(glycan, max_steps = 20)
#'
#' # Rebuild pathways for multiple glycans
#' glycans <- c(
#'   "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-",
#'   "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' )
#' path <- rebuild_biosynthesis(glycans, max_steps = 20)
#'
#' # View the path
#' igraph::as_data_frame(path, what = "edges")
#'
#' @export
rebuild_biosynthesis <- function(
  glycans,
  enzymes = NULL,
  max_steps = 20,
  filter = NULL
) {
  # Parse and validate basic inputs first
  glycans <- .process_glycans_arg(glycans)
  enzymes <- .process_enzymes_arg(enzymes, glycans, apply_prefilter = TRUE)
  checkmate::assert_int(max_steps, lower = 1)
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }

  # Find all possible paths using unified BFS logic
  starting_glycan <- .decide_starting_glycan(glycans[1])  # Use first glycan to decide starting point
  .perform_bfs_synthesis(starting_glycan, glycans, enzymes, max_steps, filter)
}

.decide_starting_glycan <- function(glycan) {
  if (glymotif::is_n_glycan(glycan)) {
    start <- glyparse::parse_iupac_condensed("Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
  } else {
    start <- glyparse::parse_iupac_condensed("GalNAc(a1-")
  }
  return(start)
}