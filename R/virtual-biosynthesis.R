#' Trace a Virtual Biosynthetic Path of Glycans
#'
#' Reconstruct every structure-driven biosynthetic path for one or more
#' glycans by trimming terminal residues and sulfate groups backward. Unlike
#' [trace_biosynthesis()], this does not require known enzyme rules.
#'
#' @inheritSection have_enzyme Important notes
#' @inheritSection trace_biosynthesis Input compatibility
#' @inheritParams trace_biosynthesis
#'
#' @section Virtual enzymes:
#' Each edge is named for the residue added by that step. Intact glycans include
#' the linkage anomer and acceptor position, so a beta-1,4-linked GlcNAc is
#' labeled `"b4GlcNAcT"`. Topological glycans omit linkage
#' information and use `"GlcNAcT"`. Generic topological glycans use
#' their preserved generic residue names, such as `"HexNAcT"`.
#'
#' Sulfation is represented as its own atomic transition. Sulfate additions at
#' positions 3 and 6 use `"3SulfoT"` and `"6SulfoT"`; an unknown or other
#' position uses `"?SulfoT"`. A sulfated terminal residue is therefore
#' desulfated before the residue itself can be trimmed.
#'
#' Virtual tracing starts N-glycans at the N-glycan core and all other glycans
#' at their reducing-end root residue. Sulfates are removed from these
#' automatically selected starts. In [path_biosynthesis_virtual()], the
#' explicit `from` glycan is always the virtual starting structure, including
#' any sulfate groups it contains; those sulfates must also occur in `to`.
#' These networks represent structural possibilities rather than biological
#' feasibility.
#'
#' Generic structures do not retain glycan-class metadata. A generic structure
#' matching the N-glycan-core topology is therefore assumed to be an N-glycan;
#' use [path_biosynthesis_virtual()] with an explicit `from` when that topology
#' belongs to another glycan class.
#'
#' @param enzymes A character vector of concrete gene symbols or abstract
#'   activity names, or a list of [enzyme()] objects. Concrete and abstract
#'   enzymes cannot be mixed. Used only when `annotate_enzymes` is `TRUE`; if
#'   `NULL`, all available concrete enzymes are considered.
#' @param annotate_enzymes Whether to annotate each virtual transition with
#'   concrete enzymes whose rules can perform it. Defaults to `FALSE`.
#'
#' @returns A `glyenzy_virtual_biosynthesis_network` object inheriting from
#' `glyenzy_biosynthesis_network` and [igraph::igraph()]. Vertices contain
#' IUPAC-condensed strings in `name` and a logical `target` attribute indicating
#' whether each vertex is a target glycan. At most one directed edge connects
#' each substrate and product. Edges have a forward `step`, `is_virtual = TRUE`,
#' a virtual-enzyme `enzyme` label, and a list-valued `enzymes` attribute. When
#' `annotate_enzymes` is `TRUE`, `enzymes` contains every candidate concrete
#' enzyme for each transition; otherwise each element is empty.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' virtual_path <- trace_biosynthesis_virtual(
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
#' )
#'
#' annotated_path <- trace_biosynthesis_virtual(
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
#'   annotate_enzymes = TRUE
#' )
#'
#' @export
trace_biosynthesis_virtual <- function(
  glycans,
  enzymes = NULL,
  annotate_enzymes = FALSE
) {
  glycans <- .process_biosynthesis_glycans_arg(glycans)
  checkmate::assert_flag(annotate_enzymes)

  if (annotate_enzymes) {
    enzymes <- .process_enzymes_arg(enzymes, apply_prefilter = FALSE)
  } else {
    .validate_virtual_enzymes(enzymes)
  }

  path <- .perform_virtual_synthesis(
    .decide_virtual_starting_glycan(glycans[1]),
    glycans
  )
  if (annotate_enzymes) {
    path <- .amplify_virtual_edges(path, enzymes)
  }
  .finalize_biosynthesis_network(
    path,
    glycans,
    match = .bfs_target_match(glycans),
    virtual = TRUE
  )
}

#' Find a Virtual Biosynthesis Path Between Glycan Structures
#'
#' Infer every structure-driven biosynthetic path from `from` to `to` by
#' trimming terminal residues and sulfate groups from `to` backward to `from`.
#' Unlike [path_biosynthesis()], this does not require known enzyme rules.
#'
#' @inheritSection trace_biosynthesis_virtual Virtual enzymes
#' @inheritSection trace_biosynthesis Input compatibility
#' @inheritParams path_biosynthesis
#'
#' @param enzymes A character vector of concrete gene symbols or abstract
#'   activity names, or a list of [enzyme()] objects. Concrete and abstract
#'   enzymes cannot be mixed. Used only when `annotate_enzymes` is `TRUE`; if
#'   `NULL`, all available concrete enzymes are considered.
#' @param annotate_enzymes Whether to annotate each virtual transition with
#'   concrete enzymes whose rules can perform it. Defaults to `FALSE`.
#'
#' @returns A `glyenzy_virtual_biosynthesis_network` object inheriting from
#' `glyenzy_biosynthesis_network` and [igraph::igraph()]. Vertices contain
#' IUPAC-condensed strings in `name` and a logical `target` attribute indicating
#' the target glycan. At most one directed edge connects each substrate and
#' product. Edges have a forward `step`, `is_virtual = TRUE`, a virtual-enzyme
#' `enzyme` label, and a list-valued `enzymes` attribute. When
#' `annotate_enzymes` is `TRUE`, `enzymes` contains every candidate concrete
#' enzyme for each transition; otherwise each element is empty.
#'
#' @examples
#' virtual_path <- path_biosynthesis_virtual(
#'   "Gal(b1-3)GalNAc(a1-",
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
#' )
#'
#' annotated_path <- path_biosynthesis_virtual(
#'   "Gal(b1-3)GalNAc(a1-",
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
#'   annotate_enzymes = TRUE
#' )
#'
#' @export
path_biosynthesis_virtual <- function(
  from,
  to,
  enzymes = NULL,
  annotate_enzymes = FALSE
) {
  glycans <- .process_biosynthesis_path_args(from, to)
  from <- glycans$from
  to <- glycans$to
  checkmate::assert_flag(annotate_enzymes)

  if (annotate_enzymes) {
    enzymes <- .process_enzymes_arg(enzymes, apply_prefilter = FALSE)
  } else {
    .validate_virtual_enzymes(enzymes)
  }

  path <- .perform_virtual_synthesis(
    from,
    to
  )
  if (annotate_enzymes) {
    path <- .amplify_virtual_edges(path, enzymes)
  }
  .finalize_biosynthesis_network(
    path,
    to,
    match = .bfs_target_match(to),
    virtual = TRUE
  )
}
