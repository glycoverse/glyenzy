# TODO: Add support for starter GTs.

#' Match Residues Added by an Enzyme
#'
#' This function finds residues in glycans that match the product motifs of a
#' glycosyltransferase and returns their node indices.
#'
#' @inheritSection have_enzyme Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()] vector.
#' @param enzyme A glycosyltransferase [enzyme()] or a gene symbol for one.
#'   Glycoside hydrolases are not supported.
#'
#' @return A list of integer vectors with the same length as `glycans`.
#'   Each integer vector contains node indices for residues added by `enzyme`
#'   in the corresponding glycan.
#'
#' @examples
#' glycan <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
#' match_enzyme(glycan, "ST3GAL3")
#'
#' @export
match_enzyme <- function(glycans, enzyme) {
  glycan_names <- names(glycans)
  if (!glyrepr::is_glycan_structure(glycans)) {
    cli::cli_abort(
      "{.arg glycans} must be a {.cls glyrepr_structure} vector."
    )
  }
  glycans <- .process_glycans_arg(glycans)
  enzyme <- .process_enzyme_arg(enzyme)
  .validate_match_enzyme_type(enzyme)

  res <- .match_enzyme(glycans, enzyme)
  if (!is.null(glycan_names)) {
    names(res) <- glycan_names
  }
  res
}

#' Match residues added by a glycosyltransferase
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @returns A list of integer vectors with the same length as `glycans`.
#' @noRd
.match_enzyme <- function(glycans, enzyme) {
  if (length(enzyme$rules) == 0L) {
    return(rep(list(integer()), length(glycans)))
  }

  rule_res <- purrr::map(enzyme$rules, ~ .match_enzyme_rule(glycans, .x))
  res <- purrr::pmap(rule_res, c)
  purrr::map(res, .unique_match_indices)
}

#' Normalize matched node indices
#'
#' @param x A combined integer vector, or `NULL` when no rule matched.
#'
#' @returns An integer vector.
#' @noRd
.unique_match_indices <- function(x) {
  if (is.null(x)) {
    return(integer())
  }
  unique(x)
}

#' Match residues added by one enzyme rule
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param rule A `glyenzy_enzyme_rule` object.
#'
#' @returns A list of integer vectors with the same length as `glycans`.
#' @noRd
.match_enzyme_rule <- function(glycans, rule) {
  product_alignment <- .product_alignment(rule)
  product_matches <- glymotif::match_motif(
    glycans,
    rule$product,
    alignment = product_alignment
  )
  acceptor_matches <- .match_product_acceptor_motifs(
    glycans,
    rule,
    product_alignment
  )

  purrr::map2(
    product_matches,
    acceptor_matches,
    ~ .match_enzyme_rule_single(.x, .y, rule)
  )
}

#' Match acceptor motifs inside final product glycans
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param rule A `glyenzy_enzyme_rule` object.
#' @param product_alignment The alignment mode used for final product matching.
#'
#' @returns A nested list of acceptor match indices.
#' @noRd
.match_product_acceptor_motifs <- function(
  glycans,
  rule,
  product_alignment
) {
  res <- glymotif::match_motif(
    glycans,
    rule$acceptor,
    alignment = product_alignment
  )
  if (length(rule$rejects) > 0) {
    rej_match_res <- glymotif::match_motifs(
      glycans,
      rule$rejects,
      product_alignment
    )
    res <- .reject_matches(res, rej_match_res)
  }
  res
}

#' Match residues added by one enzyme rule in one glycan
#'
#' @param product_matches Product motif matches in one glycan.
#' @param acceptor_matches Acceptor motif matches in one glycan.
#' @param rule A `glyenzy_enzyme_rule` object.
#'
#' @returns An integer vector of matched node indices.
#' @noRd
.match_enzyme_rule_single <- function(product_matches, acceptor_matches, rule) {
  added_indices <- purrr::map(
    product_matches,
    .match_product_acceptor,
    acceptor_matches = acceptor_matches,
    product_idx = rule$product_idx
  )
  unlist(added_indices, use.names = FALSE)
}

#' Match a product match to its corresponding acceptor match
#'
#' @param product_match An integer vector mapping product nodes to glycan nodes.
#' @param acceptor_matches Acceptor motif matches in the same glycan.
#' @param product_idx The product motif index of the added residue.
#'
#' @returns A length-one integer vector if a corresponding acceptor match exists,
#'   or `integer(0)` otherwise.
#' @noRd
.match_product_acceptor <- function(
  product_match,
  acceptor_matches,
  product_idx
) {
  acceptor_idx <- product_match[-product_idx]
  is_matching_acceptor <- purrr::map_lgl(
    acceptor_matches,
    ~ setequal(.x, acceptor_idx)
  )
  if (!any(is_matching_acceptor)) {
    return(integer())
  }

  setdiff(product_match, acceptor_matches[[which(is_matching_acceptor)[[1]]]])
}

#' Validate that an enzyme can be used by `match_enzyme()`
#'
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @noRd
.validate_match_enzyme_type <- function(enzyme) {
  if (enzyme$type != "GT") {
    cli::cli_abort(
      "{.fn match_enzyme} only supports glycosyltransferases."
    )
  }
  invisible(enzyme)
}
