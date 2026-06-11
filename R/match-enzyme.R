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
#' @param method Method used to decide whether the enzyme is involved.
#'   `"motif"` matches product motifs directly in each glycan.
#'   `"path"` matches substrates and products from [trace_biosynthesis()]
#'   results back to each glycan, which is more accurate but slower.
#'
#' @return A list of integer vectors with the same length as `glycans`.
#'   Each integer vector contains node indices for residues added by `enzyme`
#'   in the corresponding glycan.
#'
#' @examples
#' glycan <- glyrepr::as_glycan_structure("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
#' match_enzyme(glycan, "ST3GAL3")
#' match_enzyme(glycan, "ST3GAL3", method = "path")
#'
#' @export
match_enzyme <- function(glycans, enzyme, method = c("motif", "path")) {
  method <- match.arg(method)
  glycan_names <- names(glycans)
  if (!glyrepr::is_glycan_structure(glycans)) {
    cli::cli_abort(
      "{.arg glycans} must be a {.cls glyrepr_structure} vector."
    )
  }
  glycans <- .process_glycans_arg(glycans)
  enzyme <- .process_enzyme_arg(enzyme)
  .validate_match_enzyme_type(enzyme)

  res <- switch(
    method,
    motif = .match_enzyme_motif(glycans, enzyme),
    path = .match_enzyme_path(glycans, enzyme)
  )
  if (!is.null(glycan_names)) {
    names(res) <- glycan_names
  }
  res
}

#' Match residues using final-glycan motif matching
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @returns A list of integer vectors with the same length as `glycans`.
#' @noRd
.match_enzyme_motif <- function(glycans, enzyme) {
  .match_enzyme(glycans, enzyme)
}

#' Match residues only when the enzyme appears in traced paths
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @returns A list of integer vectors with the same length as `glycans`.
#' @noRd
.match_enzyme_path <- function(glycans, enzyme) {
  purrr::map(
    glycans,
    ~ .match_enzyme_path_single(.x, enzyme$name)
  )
}

#' Match residues added by traced enzyme edges in one glycan
#'
#' @param glycan A length-one `glyrepr_structure` vector.
#' @param enzyme_name An enzyme name.
#'
#' @returns An integer vector of residue indices.
#' @noRd
.match_enzyme_path_single <- function(glycan, enzyme_name) {
  path <- trace_biosynthesis(glycan)
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[edges$enzyme == enzyme_name, , drop = FALSE]
  if (nrow(edges) == 0L) {
    return(integer())
  }

  edge_matches <- purrr::map2(
    edges$from,
    edges$to,
    ~ .match_traced_edge_added_residues(glycan, .x, .y)
  )
  .unique_match_indices(unlist(edge_matches, use.names = FALSE))
}

#' Match residues added between one traced substrate-product pair
#'
#' @param glycan A length-one `glyrepr_structure` vector.
#' @param substrate A traced substrate glycan string.
#' @param product A traced product glycan string.
#'
#' @returns An integer vector of residue indices.
#' @noRd
.match_traced_edge_added_residues <- function(glycan, substrate, product) {
  substrate <- glyparse::auto_parse(substrate)
  product <- glyparse::auto_parse(product)
  substrate_matches <- glymotif::match_motif(
    glycan,
    substrate,
    alignment = "substructure"
  )[[1]]
  product_matches <- glymotif::match_motif(
    glycan,
    product,
    alignment = "substructure"
  )[[1]]

  added_residues <- purrr::map(
    product_matches,
    .match_traced_product_substrate,
    substrate_matches = substrate_matches
  )
  unlist(added_residues, use.names = FALSE)
}

#' Match one traced product occurrence to substrate occurrences
#'
#' @param product_match A product motif match in the final glycan.
#' @param substrate_matches Substrate motif matches in the final glycan.
#'
#' @returns An integer vector of residue indices.
#' @noRd
.match_traced_product_substrate <- function(product_match, substrate_matches) {
  is_matching_substrate <- purrr::map_lgl(
    substrate_matches,
    ~ all(.x %in% product_match)
  )
  if (!any(is_matching_substrate)) {
    return(integer())
  }

  substrate_match <- substrate_matches[[which(is_matching_substrate)[[1]]]]
  setdiff(product_match, substrate_match)
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
  UseMethod(".match_enzyme", enzyme)
}

.match_enzyme.glyenzy_starter_gt_enzyme <- function(glycans, enzyme) {
  if (length(enzyme$rules) == 0L) {
    return(rep(list(integer()), length(glycans)))
  }

  have_product <- .have_enzyme_by_type.glyenzy_gt_enzyme(glycans, enzyme)
  n_mono <- glyrepr::count_mono(glycans)
  purrr::map2(
    have_product,
    n_mono,
    ~ if (.x) .y else integer()
  )
}

.match_enzyme.glyenzy_gt_enzyme <- function(glycans, enzyme) {
  if (length(enzyme$rules) == 0L) {
    return(rep(list(integer()), length(glycans)))
  }

  rule_res <- purrr::map(enzyme$rules, ~ .match_enzyme_rule(glycans, .x))
  res <- purrr::pmap(rule_res, c)
  purrr::map(res, .unique_match_indices)
}

.match_enzyme.glyenzy_enzyme <- function(glycans, enzyme) {
  if (.is_starter_gt(enzyme)) {
    return(.match_enzyme.glyenzy_starter_gt_enzyme(glycans, enzyme))
  }

  switch(
    enzyme$type,
    GT = .match_enzyme.glyenzy_gt_enzyme(glycans, enzyme),
    GH = cli::cli_abort(
      "{.fn match_enzyme} only supports glycosyltransferases."
    ),
    cli::cli_abort("Unsupported enzyme type: {enzyme$type}")
  )
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
  sort(unique(x))
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
  if (!inherits(enzyme, "glyenzy_gt_enzyme") && enzyme$type != "GT") {
    cli::cli_abort(
      "{.fn match_enzyme} only supports glycosyltransferases."
    )
  }
  invisible(enzyme)
}
