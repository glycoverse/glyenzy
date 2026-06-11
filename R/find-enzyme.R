#' Identify Potentially Involved Enzymes
#'
#' This function returns all possible isoenzymes associated with the biosynthetic
#' steps of the input glycan.
#' Note that this function ignores the residues in glycans
#' that cannot be matched to any enzyme rules.
#'
#' @inheritSection have_enzyme Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#' @param return_list If `NULL` (default),
#'   return a list of character vectors when `glycans` has length greater than 1,
#'   and a single character vector when `glycans` has length 1.
#'   Set to `TRUE` to always return a list.
#'   This can be useful when you are working programmatically with unknown input length.
#'   Note that when `return_list = FALSE` and `length(glycans) > 1`,
#'   an error will be thrown.
#' @param method Method used to infer enzyme involvement.
#'   `"motif"` checks product motifs directly in each glycan.
#'   `"path"` extracts enzymes from [trace_biosynthesis()] results, which is
#'   more accurate but slower.
#'
#' @return A character vector or a list of character vectors (see `return_list` parameter),
#'   each containing the names of enzymes involved in the biosynthesis of the corresponding glycan.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Use `glycan_structure()`
#' glycans <- auto_parse(c(
#'   "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
#'   "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' ))
#' find_enzyme(glycans)
#'
#' # Or use characters directly
#' find_enzyme("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
#'
#' # Use reconstructed biosynthesis paths
#' find_enzyme(glycans, method = "path")
#'
#' @export
find_enzyme <- function(
  glycans,
  return_list = NULL,
  method = c("motif", "path")
) {
  method <- match.arg(method)
  glycans <- .process_glycans_arg(glycans)
  return_list <- .validate_return_list(return_list, length(glycans))
  res <- switch(
    method,
    motif = .find_enzyme_motif(glycans),
    path = .find_enzyme_path(glycans)
  )

  .format_result(res, return_list)
}

#' Identify enzymes using final-glycan motif matching
#'
#' @param glycans A `glyrepr_structure` vector.
#'
#' @returns A list of character vectors.
#' @noRd
.find_enzyme_motif <- function(glycans) {
  # Compute is_n once for all enzymes to avoid repeated computation
  is_n <- .is_n_glycan(glycans)
  masks <- purrr::map(glyenzy_enzymes, ~ .safe_have_enzyme(glycans, .x, is_n))
  mast_mat <- do.call(cbind, masks)
  res <- purrr::map(
    seq_along(glycans),
    ~ names(glyenzy_enzymes)[mast_mat[.x, ]]
  )
  res
}

#' Identify trace-derived enzymes for each glycan
#'
#' @param glycans A `glyrepr_structure` vector.
#'
#' @returns A list of character vectors.
#' @noRd
.find_enzyme_path <- function(glycans) {
  edges <- .trace_enzyme_edges(glycans)
  purrr::map(edges, unique)
}

# Like `.have_enzyme_motif()`, but returns FALSE instead of throwing error.
.safe_have_enzyme <- function(glycans, enzyme, is_n) {
  tryCatch(
    .have_enzyme_motif(glycans, enzyme, is_n),
    error = function(e) rep(FALSE, length(glycans))
  )
}
