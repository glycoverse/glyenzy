#' Identify Potentially Involved Enzymes
#'
#' This function returns all possible isoenzymes associated with the biosynthetic
#' steps of the input glycan.
#'
#' @inheritSection is_synthesized_by Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#'
#' @return A list of character vectors, each containing the names of enzymes
#'   involved in the biosynthesis of the corresponding glycan.
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
get_involved_enzymes <- function(glycans) {
  if (is.character(glycans)) {
    glycans <- glyparse::auto_parse(glycans)
  } else if (!glyrepr::is_glycan_structure(glycans)) {
    cli::cli_abort(c(
      "{.arg glycans} must be a {.cls glyrepr_structure} vector or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(glycans)}}."
    ))
  }

  # Compute is_n once for all enzymes to avoid repeated computation
  is_n <- glymotif::is_n_glycan(glycans)
  masks <- purrr::map(glyenzy_enzymes, ~ .safe_is_synthesized_by(glycans, .x, is_n))
  mast_mat <- do.call(cbind, masks)
  purrr::map(seq_along(glycans), ~ names(glyenzy_enzymes)[mast_mat[.x, ]])
}

# Like `.is_synthesized_by()`, but returns FALSE instead of throwing error.
.safe_is_synthesized_by <- function(glycans, enzyme, is_n) {
  tryCatch(
    .is_synthesized_by(glycans, enzyme, is_n),
    error = function(e) rep(FALSE, length(glycans))
  )
}