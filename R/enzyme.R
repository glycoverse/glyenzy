#' Create a Glyco-Enzyme
#'
#' An glyco-enzyme is a biological catalyst that processes one glycan into another.
#' Two types of enzymes are supported: glycosyltransferases (GT) and exoglycosidases (GD).
#' Glycosyltransferases add a monosaccharide to the terminal of a glycan,
#' while exoglycosidases remove a monosaccharide from the terminal of a glycan.
#'
#' @details
#' `acceptor` and `product` must follow these rules:
#'
#' 1. `acceptor` and `product` must be a single glycan structure.
#' 2. For GTs, `product` must have exactly one more residue that `acceptor`,
#'    and the extra residue must locate at the terminal.
#' 3. For GDs, `acceptor` must have exactly one more residue that `product`,
#'    and the extra residue must locate at the terminal.
#'
#' @param name The name of the enzyme.
#' @param acceptor A `glyrepr_structure` object or a glycan structure string representing the acceptor.
#'   The acceptor is the substrate that the enzyme acts on.
#' @param product A `glyrepr_structure` object or a glycan structure string representing the product.
#'   The product is the glycan after the enzyme acts.
#' @param type The type of the enzyme, "GT" for glycosyltransferase or "GD" for exoglycosidase.
#' @param species The species of the enzyme, e.g. "human" or "mouse".
#'
#' @return A `glyenzy_enzyme` object.
#' @export
create_enzyme <- function(name, acceptor, product, type, species) {
  acceptor <- .convert_to_glycan_structure(acceptor, "acceptor")
  product <- .convert_to_glycan_structure(product, "product")
  validate_enzyme(new_enzyme(name, acceptor, product, type, species))
}

#' Create a new enzyme object
#'
#' @param name The name of the enzyme.
#' @param acceptor A `glyrepr_structure` object representing the acceptor.
#' @param product A `glyrepr_structure` object representing the product.
#' @param type The type of the enzyme, "GT" for glycosyltransferase or "GD" for exoglycosidase.
#' @param species The species of the enzyme, e.g. "human" or "mouse".
#'
#' @return A `glyenzy_enzyme` object.
#' @noRd
new_enzyme <- function(name, acceptor, product, type, species) {
  structure(
    list(name = name, acceptor = acceptor, product = product, type = type, species = species),
    class = "glyenzy_enzyme"
  )
}

.convert_to_glycan_structure <- function(x, name) {
  if (is.character(x)) {
    glyparse::auto_parse(x)
  } else if (glyrepr::is_glycan_structure(x)) {
    x
  } else {
    cli::cli_abort(c(
      "{.arg {name}} must be either a character string or a glyrepr_structure object.",
      "x" = "Got class: {.cls {class(x)}}."
    ))
  }
}

validate_enzyme <- function(x) {
  checkmate::assert_class(x, "glyenzy_enzyme")

  # Check length of acceptor and product
  if (length(x$acceptor) != 1) {
    cli::cli_abort("The {.arg acceptor} must be a single structure.")
  }
  if (length(x$product) != 1) {
    cli::cli_abort("The {.arg product} must be a single structure.")
  }

  # Based on "type", check acceptor and product
  if (x$type == "GT") {
    .assert_product_acceptor(x$acceptor, x$product, "acceptor", "product")
  } else {
    .assert_product_acceptor(x$product, x$acceptor, "product", "acceptor")
  }

  invisible(x)
}

.assert_product_acceptor <- function(smaller, larger, smaller_name, larger_name) {
  smaller_graph <- glyrepr::get_structure_graphs(smaller)
  larger_graph <- glyrepr::get_structure_graphs(larger)
  match_res <- glymotif::match_motif(larger, smaller)[[1]]  # only one glycan, so `[[1]]`
  # `smaller` is a substructure of `larger`
  if (length(match_res) != 1) {
    cli::cli_abort("{.arg {smaller_name}} must be a substructure of {.arg {larger_name}}.")
  }
  match_res <- match_res[[1]]  # only one match, so `[[1]]`
  # `larger` has only one more residue than `smaller`
  if (igraph::vcount(larger_graph) - igraph::vcount(smaller_graph) != 1) {
    cli::cli_abort("{.arg {larger_name}} must have exactly one more residue than {.arg {smaller_name}}.")
  }
  # the extra residue in `larger` has an out-degree of 0
  extra_residue_i <- setdiff(1:igraph::vcount(larger_graph), match_res)
  if (igraph::degree(larger_graph, extra_residue_i, mode = "out") != 0) {
    cli::cli_abort("The extra residue in {.arg {larger_name}} must be at the terminal.")
  }
}

#' @export
print.glyenzy_enzyme <- function(x, ...) {
  cli::cli_h2("Glyco-enzyme: {.field {x$name}}")

  type_desc <- if (x$type == "GT") "Glycosyltransferase" else "Exoglycosidase"
  cli::cli_alert_info("Type: {.val {type_desc}} ({.field {x$type}})")
  cli::cli_alert_info("Species: {.field {x$species}}")

  cli::cli_h3("Reaction")
  acceptor_str <- as.character(x$acceptor)
  product_str <- as.character(x$product)

  if (x$type == "GT") {
    cli::cli_text("{.field {acceptor_str}} + donor -> {.field {product_str}}")
  } else {
    cli::cli_text("{.field {acceptor_str}} -> {.field {product_str}} + released")
  }
}