#' Find a Biosynthesis Path Between Glycan Structures
#'
#' Find biosynthetic paths from one glycan structure to another. The default
#' method uses known enzyme rules in a forward breadth-first search. The
#' virtual-enzyme method trims `to` backward to `from` and returns every
#' possible residue-addition order.
#'
#' @inheritSection trace_biosynthesis Important notes
#' @inheritSection trace_biosynthesis Virtual enzymes
#'
#' @param from A [glyrepr::glycan_structure()] scalar, or a character string
#'   supported by [glyparse::auto_parse()]. The starting glycan structure.
#' @param to A [glyrepr::glycan_structure()] scalar, or a character string
#'   supported by [glyparse::auto_parse()]. The target glycan structure.
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()] objects.
#'   If `NULL` (default), all available enzymes will be used.
#'   Must be `NULL` when `method = "virtual"`.
#' @param max_steps Integer, maximum number of enzymatic steps to search.
#'   Default is 10.
#' @param filter Optional function to filter generated glycans at each step.
#'   Should take a [glyrepr::glycan_structure()] vector as input and return
#'   a logical vector of the same length.
#'   For `method = "enzymatic"`, it filters generated products. For
#'   `method = "virtual"`, it filters generated precursors during backward
#'   trimming.
#' @param method Biosynthesis inference method. `"enzymatic"` (default) uses
#'   known enzyme rules in a forward search. `"virtual"` uses virtual enzymes
#'   in a backward search that removes terminal residues from `to` until `from`
#'   is reached.
#'
#' @returns An [igraph::igraph()] object representing the synthesis path(s).
#'   Vertices represent glycan structures with `name` attribute containing
#'   IUPAC-condensed strings. Edges represent enzymatic reactions with
#'   `enzyme` attribute containing gene symbols or virtual-enzyme names and
#'   `step` attribute indicating the forward synthesis step.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Find shortest path
#' from <- "Gal(b1-4)GlcNAc(b1-"
#' to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
#' path <- path_biosynthesis(from, to, enzymes = "ST6GAL1", max_steps = 3)
#'
#' # Ignore known enzyme specificity and infer additions by trimming backward
#' virtual_path <- path_biosynthesis(
#'   "Gal(b1-3)GalNAc(a1-",
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
#'   method = "virtual"
#' )
#'
#' # View the path
#' igraph::as_data_frame(path, what = "edges")
#'
#' @export
path_biosynthesis <- function(
  from,
  to,
  enzymes = NULL,
  max_steps = 10,
  filter = NULL,
  method = c("enzymatic", "virtual")
) {
  # Parse and validate basic inputs first
  method <- match.arg(method)
  from <- .process_glycan_arg(from, allow_generic = TRUE)
  to <- .process_glycan_arg(to, allow_generic = TRUE)
  checkmate::assert_int(max_steps, lower = 1)
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }

  if (identical(method, "virtual")) {
    .validate_virtual_enzymes(enzymes)
    from <- .prepare_virtual_start(from, to)
    return(.perform_virtual_synthesis(from, to, max_steps, filter))
  }

  enzymes <- .process_enzymes_arg(enzymes, apply_prefilter = FALSE)
  # Perform BFS search using unified logic
  .perform_bfs_synthesis(from, to, enzymes, max_steps, filter)
}
