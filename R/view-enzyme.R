#' View Residues Added or Modified by an Enzyme
#'
#' Visualize where an enzyme adds or modifies residues in a glycan structure.
#'
#' `view_enzyme()` matches one `enzyme` against one `glycan` with the same
#' matching rules used by [match_enzyme()], then draws the glycan with the
#' matched residues highlighted.
#'
#' @inheritSection have_enzyme Important notes
#'
#' @param glycan A [glyrepr::glycan_structure()], or a glycan structure string
#'   supported by [glyparse::auto_parse()].
#' @param enzyme A glycosyltransferase or sulfotransferase [enzyme()], or a gene
#'   symbol for one. Glycoside hydrolases are not supported.
#'
#' @returns
#' A `ggplot` object returned by [glydraw::draw_cartoon()]. If no match is found,
#' the glycan is drawn without highlighted residues and a cli alert is emitted.
#'
#' @seealso [match_enzyme()], [glydraw::draw_cartoon()]
#'
#' @examples
#' glycan <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
#'
#' \dontrun{
#' view_enzyme(glycan, "ST3GAL3")
#' }
#' @export
view_enzyme <- function(glycan, enzyme) {
  if (length(glycan) != 1L) {
    cli::cli_abort(c(
      "Only one glycan can be visualized at a time.",
      "x" = "Got {.val {length(glycan)}} glycans."
    ))
  }

  glycan <- .process_glycans_arg(glycan)
  idx <- match_enzyme(glycan, enzyme)[[1]]

  if (length(idx) == 0L) {
    cli::cli_alert_danger(
      "No residues added or modified by the enzyme were found in the glycan."
    )
    glydraw::draw_cartoon(glycan, highlight = integer(0))
  } else {
    glydraw::draw_cartoon(glycan, highlight = idx)
  }
}
