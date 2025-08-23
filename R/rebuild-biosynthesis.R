#' Rebuild the Biosynthetic Path of a Glycan
#'
#' Reconstruct the biosynthetic pathway of a glycan from a given starting structure
#' using enzymatic reactions.
#' This function uses a combination of depth-first search and heuristics to find a likely pathway,
#' or all feasible pathways.
#'
#' For N-glycans, the starting structure is assumed to be "Glc(3)Man(9)GlcNAc(2)",
#' the N-glycan precursor transfered to Asn by OST.
#' For O-glycans, the starting structure is assumed to be "GalNAc(a1-".
#'
#' @param glycan A [glyrepr::glycan_structure()] scalar, or a character string
#'   supported by [glyparse::auto_parse()].
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()] objects.
#'   If `NULL` (default), all available enzymes will be used.
#' @param max_steps Integer, maximum number of enzymatic steps to search.
#'   Default is 20.
#' @param filter Optional function to filter generated glycans at each step.
#'   Should take a [glyrepr::glycan_structure()] vector as input and return
#'   a logical vector of the same length.
#' @param return One of `"one"` or `"all"`. If `"one"`, returns the
#'   first shortest path found. If `"all"`, returns a graph containing all
#'   possible paths within `max_steps`. Default is `"one"`.
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
#' # Rebuild the biosynthetic pathway of a glycan
#' glycan <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' path <- rebuild_biosynthesis(glycan, max_steps = 20)
#'
#' # View the path
#' igraph::as_data_frame(path, what = "edges")
#'
#' @export
rebuild_biosynthesis <- function(
  glycan,
  enzymes = NULL,
  max_steps = 20,
  filter = NULL,
  return = c("one", "all")
) {
  return <- rlang::arg_match(return)

  # Parse and validate basic inputs first
  glycan <- glyrepr::as_glycan_structure(glycan)
  checkmate::assert_true(length(glycan) == 1L)
  checkmate::assert_int(max_steps, lower = 1)

  # Find all possible paths
  starting_glycan <- .decide_starting_glycan(glycan)
  find_synthesis_path(starting_glycan, glycan, enzymes, max_steps, filter, return)
}

.decide_starting_glycan <- function(glycan) {
  if (glymotif::is_n_glycan(glycan)) {
    start <- glyparse::parse_iupac_condensed("Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
  } else {
    start <- glyparse::parse_iupac_condensed("GalNAc(a1-")
  }
  return(start)
}