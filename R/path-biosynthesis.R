#' Find a Biosynthesis Path Between Glycan Structures
#'
#' Find biosynthetic paths from one glycan structure to another. The default
#' method uses known enzyme rules in a forward breadth-first search. To infer
#' structure-driven paths without enzyme specificity, use
#' [path_biosynthesis_virtual()].
#'
#' @inheritSection trace_biosynthesis Important notes
#' @inheritSection trace_biosynthesis Input compatibility
#' @inheritSection trace_biosynthesis Virtual fallback
#'
#' @param from A [glyrepr::glycan_structure()] scalar, or a character string
#'   supported by [glyparse::auto_parse()]. The starting glycan structure.
#' @param to A [glyrepr::glycan_structure()] scalar, or a character string
#'   supported by [glyparse::auto_parse()]. The target glycan structure.
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()]
#'   objects. If `NULL` (default), all available enzymes will be used.
#' @param max_steps Integer, maximum number of enzymatic steps to search.
#'   Default is 10.
#' @param max_virtual_steps Integer, maximum number of target-directed virtual
#'   enzyme steps allowed when no fully enzymatic path exists.
#'   Default is `0L`, which disables virtual fallback.
#'   See the "Virtual fallback" section for more details.
#' @param filter Optional function to filter generated glycans at each step.
#'   Should take a [glyrepr::glycan_structure()] vector as input and return
#'   a logical vector of the same length. It filters generated products.
#'
#' @returns A `glyenzy_biosynthesis_network` object inheriting from
#'   [igraph::igraph()] and representing the synthesis path(s). Vertices
#'   represent glycan structures, with IUPAC-condensed strings in the `name`
#'   attribute and a logical `target` attribute indicating the target glycan.
#'   At most one directed edge connects each substrate and product. Every edge
#'   has an integer `step`, a logical `is_virtual`, a scalar display label in
#'   `enzyme`, and a list-valued `enzymes` attribute containing concrete enzyme
#'   candidates. For known reactions, `enzyme` combines the candidates with
#'   `" / "`; for virtual reactions, it contains the structural virtual-enzyme
#'   name and `enzymes` is empty unless concrete candidates were annotated.
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
  max_virtual_steps = 0L
) {
  # Parse and validate inputs first
  glycans <- .process_biosynthesis_path_args(from, to)
  from <- glycans$from
  to <- glycans$to
  checkmate::assert_int(max_steps, lower = 1)
  checkmate::assert_int(max_virtual_steps, lower = 0)
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }

  enzymes <- .process_enzymes_arg(enzymes, apply_prefilter = FALSE)
  # Perform BFS search using unified logic
  .perform_bfs_synthesis(
    from,
    to,
    enzymes,
    max_steps,
    filter,
    max_virtual_steps
  ) |>
    .finalize_biosynthesis_network(
      to,
      match = .bfs_target_match(to)
    )
}
