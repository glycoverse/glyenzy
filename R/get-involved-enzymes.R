#' Identify Potentially Involved Enzymes
#'
#' This function returns all possible isoenzymes associated with the biosynthetic
#' steps of the input glycan.
#'
#' @inheritSection is_synthesized_by Important notes
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
#' get_involved_enzymes(glycans)
#'
#' # Or use characters directly
#' get_involved_enzymes("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
#'
#' @export
get_involved_enzymes <- function(glycans, return_list = NULL) {
  glycans <- .process_glycans_arg(glycans)
  return_list <- .validate_return_list(return_list, length(glycans))

  # Compute is_n once for all enzymes to avoid repeated computation
  is_n <- .is_n_glycan(glycans)
  masks <- purrr::map(glyenzy_enzymes, ~ .safe_is_synthesized_by(glycans, .x, is_n))
  mast_mat <- do.call(cbind, masks)
  res <- purrr::map(seq_along(glycans), ~ names(glyenzy_enzymes)[mast_mat[.x, ]])
  .format_result(res, return_list)
}

# Like `.is_synthesized_by()`, but returns FALSE instead of throwing error.
.safe_is_synthesized_by <- function(glycans, enzyme, is_n) {
  tryCatch(
    .is_synthesized_by(glycans, enzyme, is_n),
    error = function(e) rep(FALSE, length(glycans))
  )
}