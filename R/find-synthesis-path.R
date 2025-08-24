#' Find Synthesis Path Between Glycan Structures
#'
#' Find a synthesis path from one glycan structure to another using enzymatic reactions.
#' This function uses breadth-first search to find the shortest path or all possible paths
#' within a given number of steps.
#'
#' @inheritSection rebuild_biosynthesis Important notes
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
#'   It will be applied to all the generated glycans at each BFS step for pruning.
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
#' # Find shortest path
#' from <- "Gal(b1-4)GlcNAc(b1-"
#' to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
#' path <- find_synthesis_path(from, to, enzymes = "ST6GAL1", max_steps = 3)
#'
#' # View the path
#' igraph::as_data_frame(path, what = "edges")
#'
#' @export
find_synthesis_path <- function(
  from,
  to,
  enzymes = NULL,
  max_steps = 10,
  filter = NULL
) {
  # Parse and validate basic inputs first
  from <- .process_glycan_arg(from)
  to <- .process_glycan_arg(to)
  enzymes <- .process_enzymes_arg(enzymes, apply_prefilter = FALSE)
  checkmate::assert_int(max_steps, lower = 1)
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }
  # Perform BFS search using unified logic
  .perform_bfs_synthesis(from, to, enzymes, max_steps, filter)
}