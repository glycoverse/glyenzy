#' Does an Enzyme Involve in the Biosynthesis of a Glycan?
#'
#' @description
#' Glycans are synthesized by a series of enzymatic reactions.
#' This function checks if an enzyme is involved in the biosynthesis of
#' a vector of glycan structures.
#'
#' This function only works for glycans with "concrete" residues (e.g. "Glc", "GalNAc"),
#' not "generic" residues (e.g. "Hex", "HexNAc").
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#' @param enzyme An [enzyme()] or a gene symbol.
#'
#' @return A logical vector.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Use `glycan_structure()` and `enzyme()`
#' glycan <- auto_parse("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-")
#' involve(glycan, enzyme("ST6GAL1"))
#'
#' # Or use characters directly
#' involve("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-", "ST6GAL1")
#'
#' @export
involve <- function(glycans, enzyme) {
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
  # Check if enzyme is involved
  mat <- glymotif::have_motifs(glycans, enzyme$markers$motifs, alignments = enzyme$markers$alignments)
  unname(rowSums(mat) > 0)
}
