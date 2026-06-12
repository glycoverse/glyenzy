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
  masks <- purrr::map(glyenzy_enzymes, ~ .safe_have_enzyme(glycans, .x))
  mask_mat <- do.call(cbind, masks)
  res <- purrr::map(
    seq_along(glycans),
    ~ names(glyenzy_enzymes)[mask_mat[.x, ]]
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
  purrr::map(glycans, .find_enzyme_path_single)
}

#' Identify trace-derived enzymes for one glycan
#'
#' @param glycan A length-one `glyrepr_structure` vector.
#'
#' @returns A character vector of enzyme names.
#' @noRd
.find_enzyme_path_single <- function(glycan) {
  npre_enzymes <- .find_npre_enzymes(glycan)
  traced_enzymes <- .trace_enzyme_edges_single(glycan)
  unique(c(traced_enzymes, npre_enzymes))
}

#' Identify N-glycan precursor enzymes for one glycan
#'
#' @param glycan A length-one `glyrepr_structure` vector.
#'
#' @returns A character vector of enzyme names.
#' @noRd
.find_npre_enzymes <- function(glycan) {
  npre_enzymes <- purrr::keep(glyenzy_enzymes, .is_npre_gt)
  masks <- purrr::map_lgl(
    npre_enzymes,
    ~ .have_enzyme_motif.glyenzy_npre_gt_enzyme(glycan, .x)
  )
  names(npre_enzymes)[masks]
}

# Like `.have_enzyme_motif()`, but returns FALSE instead of throwing error.
.safe_have_enzyme <- function(glycans, enzyme) {
  tryCatch(
    .have_enzyme_motif(glycans, enzyme),
    error = function(e) rep(FALSE, length(glycans))
  )
}
