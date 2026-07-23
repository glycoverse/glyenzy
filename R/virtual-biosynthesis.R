#' Trace a Virtual Biosynthetic Path of Glycans
#'
#' Reconstruct every structure-driven biosynthetic path for one or more
#' glycans by trimming terminal residues and sulfate groups backward. Unlike
#' [trace_biosynthesis()], this does not require known enzyme rules.
#'
#' @inheritSection have_enzyme Important notes
#' @inheritParams trace_biosynthesis
#'
#' @section Virtual enzymes:
#' Each edge is named for the residue added by that step. Intact glycans include
#' the linkage anomer and acceptor position, so a beta-1,4-linked GlcNAc is
#' labeled `"b4GlcNAcT"`. Partial and topological glycans omit linkage
#' information and use `"GlcNAcT"`; basic glycans use the generic residue name,
#' such as `"HexNAcT"`.
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
#' Basic structures do not retain glycan-class metadata. A basic structure
#' matching the generic N-glycan-core topology is therefore assumed to be an
#' N-glycan; use [path_biosynthesis_virtual()] with an explicit `from` when
#' that topology belongs to another glycan class.
#'
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()]
#'   objects. Used only when `annotate_enzymes` is `TRUE`; if `NULL`, all
#'   available enzymes are considered.
#' @param annotate_enzymes Whether to annotate each virtual transition with
#'   concrete enzymes whose rules can perform it. Defaults to `FALSE`.
#'
#' @returns An [igraph::igraph()] object representing the synthesis path(s).
#' Vertices contain IUPAC-condensed strings in `name`; edges have a forward
#' `step` and virtual-enzyme `enzyme` attribute. When `annotate_enzymes` is
#' `TRUE`, `concrete_enzymes` is a list of character vectors containing every
#' candidate concrete enzyme for each transition.
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
  glycans <- .process_glycans_arg(glycans, allow_generic = TRUE)
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
  path
}

#' Find a Virtual Biosynthesis Path Between Glycan Structures
#'
#' Infer every structure-driven biosynthetic path from `from` to `to` by
#' trimming terminal residues and sulfate groups from `to` backward to `from`.
#' Unlike [path_biosynthesis()], this does not require known enzyme rules.
#'
#' @inheritSection trace_biosynthesis_virtual Virtual enzymes
#' @inheritParams path_biosynthesis
#'
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()]
#'   objects. Used only when `annotate_enzymes` is `TRUE`; if `NULL`, all
#'   available enzymes are considered.
#' @param annotate_enzymes Whether to annotate each virtual transition with
#'   concrete enzymes whose rules can perform it. Defaults to `FALSE`.
#'
#' @returns An [igraph::igraph()] object representing the synthesis path(s).
#' Vertices contain IUPAC-condensed strings in `name`; edges have a forward
#' `step` and virtual-enzyme `enzyme` attribute. When `annotate_enzymes` is
#' `TRUE`, `concrete_enzymes` is a list of character vectors containing every
#' candidate concrete enzyme for each transition.
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
  from <- .process_glycan_arg(from, allow_generic = TRUE)
  to <- .process_glycan_arg(to, allow_generic = TRUE)
  checkmate::assert_flag(annotate_enzymes)

  if (annotate_enzymes) {
    enzymes <- .process_enzymes_arg(enzymes, apply_prefilter = FALSE)
  } else {
    .validate_virtual_enzymes(enzymes)
  }

  path <- .perform_virtual_synthesis(
    .prepare_virtual_start(from, to),
    to
  )
  if (annotate_enzymes) {
    path <- .amplify_virtual_edges(path, enzymes)
  }
  path
}
