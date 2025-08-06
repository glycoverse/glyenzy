#' Is a Glycan Synthesized by an Enzyme?
#'
#' @description
#' Glycans are synthesized by a series of enzymatic reactions.
#' This function checks if an enzyme is involved in the biosynthesis of a glycan.
#'
#' Only works for glycans with "concrete" residues (e.g. "Glc", "GalNAc"),
#' not "generic" residues (e.g. "Hex", "HexNAc").
#'
#' @details
#' The implementation is simple: for each rule of the enzyme,
#' check if the product motif exists in the glycan.
#' If any rule is satisfied, return TRUE.
#'
#' For N-glycans, we have additional logic to handle special cases.
#' The products of ALG enzymes and MGAT1 are further trimmed by other exoglycosidases.
#' Therefore, checking the product doesn't reflect the enzyme's involvement.
#' We checks special motif markers for these enzymes.
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#' @param enzyme An [enzyme()] or a gene symbol.
#'
#' @return A logical vector of the same length as `glycans`.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Use `glycan_structure()` and `enzyme()`
#' glycan <- auto_parse("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-")
#' is_synthesized_by(glycan, enzyme("ST6GAL1"))
#'
#' # Or use characters directly
#' is_synthesized_by("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-", "ST6GAL1")
#'
#' # Vectorized input
#' glycans <- c(
#'   "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
#'   "Gal(b1-4)GlcNAc(b1-"
#' )
#' is_synthesized_by(glycans, "ST6GAL1")
#'
#' @export
is_synthesized_by <- function(glycans, enzyme) {
  # Process `glycans` argument
  if (is.character(glycans)) {
    glycans <- glyparse::auto_parse(glycans)
  } else if (!glyrepr::is_glycan_structure(glycans)) {
    cli::cli_abort(c(
      "Input glycans must be a {.cls glyrepr_structure} vector or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(glycans)}}."
    ))
  }
  # Process `enzyme` argument
  if (is.character(enzyme)) {
    enzyme <- glyenzy::enzyme(enzyme)
  } else if (!inherits(enzyme, "glyenzy_enzyme")) {
    cli::cli_abort(c(
      "Input enzyme must be a {.cls glyenzy_enzyme} object or a character string of gene symbol.",
      "x" = "Got {.cls {class(enzyme)}}."
    ))
  }

  .is_synthesized_by(glycans, enzyme)
}

#' Is a Glycan Synthesized by an Enzyme? (Internal)
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @noRd
.is_synthesized_by <- function(glycans, enzyme) {
  dplyr::if_else(
    glymotif::is_n_glycan(glycans),
    .is_synthesized_by_n_glycan(glycans, enzyme),
    .is_synthesized_by_default(glycans, enzyme)
  )
}

#' Is a N-Glycan Synthesized by an Enzyme?
#'
#' N-glycans are special because of its elongation-trmming-elongation biosynthesis machenism.
#' Firstly, Glc(3)Man(9)GlcNAc(2) is synthezised by ALG enzymes.
#' Then, it is trimmed to Man(5)GlcNAc(2).
#' Finally, it is elongated again to form various N-glycans.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @noRd
.is_synthesized_by_n_glycan <- function(glycans, enzyme) {
  # Special case for ALG family
  # These enzymes involve in the the biosynthesis process of all N-glycans.
  if (stringr::str_starts(enzyme$name, "ALG")) {
    return(TRUE)
  }

  # Special case for MGAT1
  # After MGAT1 adds the b1-2 GlcNAc, two mannoses are trimmed.
  # Therefore, checking the product doesn't reflect its involvement.
  if (enzyme$name == "MGAT1") {
    return(glymotif::have_motif(glycans, "GlcNAc(b1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"))
  }

  # Default case
  .is_synthesized_by_default(glycans, enzyme)
}

#' Is a Glycan Synthesized by an Enzyme? (Default Case)
#'
#' Check all rules of the `enzyme`.
#' Returns TRUE if any rule is satisfied.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @noRd
.is_synthesized_by_default <- function(glycans, enzyme) {
  purrr::some(enzyme$rules, ~ glymotif::have_motif(glycans, .x$product, "substructure"))
}