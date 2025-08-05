#' Create a new enzyme object
#'
#' @param name The name of the enzyme.
#' @param rules A list of `glyenzy_enzyme_rule` objects.
#' @param rejects A `glyenzy_motif_set` object representing the motifs that the enzyme rejects.
#' @param markers A `glyenzy_motif_set` object representing the motifs that
#'   indicate the enzyme is involved in the biosynthesis of a glycan.
#' @param type The type of the enzyme, "GT" for glycosyltransferase or "GD" for exoglycosidase.
#' @param species The species of the enzyme, e.g. "human" or "mouse".
#'
#' @return A `glyenzy_enzyme` object.
#' @noRd
new_enzyme <- function(name, rules, rejects, markers, type, species) {
  checkmate::assert_string(name)
  checkmate::assert_list(rules, types = "glyenzy_enzyme_rule")
  checkmate::assert_class(rejects, "glyenzy_motif_set")
  checkmate::assert_class(markers, "glyenzy_motif_set")
  checkmate::assert_choice(type, c("GT", "GD"))
  checkmate::assert_string(species)

  structure(
    list(name = name, rules = rules, rejects = rejects, markers = markers, type = type, species = species),
    class = "glyenzy_enzyme"
  )
}

#' Validate a `glyenzy_enzyme` object
#'
#' This function validate `rules`, `rejects`, and `markers`
#' by calling their own validation functions.
#'
#' @param x A `glyenzy_enzyme` object.
#' @noRd
validate_enzyme <- function(x) {
  checkmate::assert_class(x, "glyenzy_enzyme")
  if (purrr::some(x$rules, ~ .x$type != x$type)) {
    cli::cli_abort(c(
      "All rules must have the same type as the enzyme.",
      "i" = "Enzyme type: {.val {x$type}}.",
      "x" = "Got rule types: {.val {purrr::map_chr(x$rules, ~ .x$type)}}."
    ))
  }
  purrr::walk(x$rules, validate_enzyme_rule)
  validate_motif_set(x$rejects)
  validate_motif_set(x$markers)
  invisible(x)
}

new_enzyme_rule <- function(acceptor, product, acceptor_alignment, type) {
  checkmate::assert_class(acceptor, "glyrepr_structure")
  checkmate::assert_class(product, "glyrepr_structure")
  checkmate::assert_choice(acceptor_alignment, c("substructure", "core", "terminal", "whole"))
  checkmate::assert_choice(type, c("GT", "GD"))

  structure(
    list(acceptor = acceptor, product = product, acceptor_alignment = acceptor_alignment, type = type),
    class = "glyenzy_enzyme_rule"
  )
}

validate_enzyme_rule <- function(x) {
  checkmate::assert_class(x, "glyenzy_enzyme_rule")
  if (length(x$acceptor) != 1) {
    cli::cli_abort("The {.arg acceptor} must be a single structure.")
  }
  if (length(x$product) != 1) {
    cli::cli_abort("The {.arg product} must be a single structure.")
  }

  # Check acceptor and product
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
  match_res <- glymotif::match_motif(larger, smaller, alignment = "core")[[1]]  # only one glycan, so `[[1]]`
  # `smaller` is a substructure of `larger`
  if (length(match_res) != 1) {
    cli::cli_abort("{.arg {smaller_name}} must be a substructure of {.arg {larger_name}}.")
  }
  match_res <- match_res[[1]]  # only one match, so `[[1]]`
  # `larger` has only one more residue than `smaller`
  if (igraph::vcount(larger_graph) - igraph::vcount(smaller_graph) != 1) {
    cli::cli_abort("{.arg {larger_name}} must have exactly one more residue than {.arg {smaller_name}}.")
  }
}

new_motif_set <- function(motifs, alignments) {
  checkmate::assert_class(motifs, "glyrepr_structure")
  checkmate::assert_subset(alignments, c("substructure", "core", "terminal", "whole"))

  structure(
    list(motifs = motifs, alignments = alignments),
    class = "glyenzy_motif_set"
  )
}

validate_motif_set <- function(x) {
  checkmate::assert_class(x, "glyenzy_motif_set")
  if (length(x$motifs) != length(x$alignments)) {
    cli::cli_abort(c(
      "The length of {.arg motifs} must be equal to the length of {.arg alignments}.",
      "i" = "Length of {.arg motifs}: {.val {length(x$motifs)}}.",
      "i" = "Length of {.arg alignments}: {.val {length(x$alignments)}}."
    ))
  }

  invisible(x)
}

#' Print method for glyenzy_enzyme objects
#'
#' @param x A `glyenzy_enzyme` object.
#' @param ... Additional arguments passed to print methods.
#' @export
print.glyenzy_enzyme <- function(x, ...) {
  cli::cli_h1("Enzyme: {.field {x$name}}")

  # Basic information
  cli::cli_alert_info("Type: {.val {x$type}} ({.emph {if (x$type == 'GT') 'Glycosyltransferase' else 'Exoglycosidase'}})")
  cli::cli_alert_info("Species: {.val {x$species}}")

  # Rules section
  cli::cli_h2("Rules ({.val {length(x$rules)}})")
  if (length(x$rules) > 0) {
    purrr::iwalk(x$rules, function(rule, i) {
      cli::cli_alert("Rule {.val {i}}: {.field {rule$acceptor_alignment}} alignment")
      cli::cli_text("  Acceptor: {.val {as.character(rule$acceptor)}}")
      cli::cli_text("  Product:  {.val {as.character(rule$product)}}")
    })
  } else {
    cli::cli_text("  {.emph No rules defined}")
  }

  # Rejects section
  cli::cli_h2("Rejects ({.val {length(x$rejects$motifs)}})")
  if (length(x$rejects$motifs) > 0) {
    purrr::iwalk(x$rejects$motifs, function(motif, i) {
      cli::cli_text("  {.val {i}}: {.val {as.character(motif)}} ({.field {x$rejects$alignments[i]}})")
    })
  } else {
    cli::cli_text("  {.emph No reject motifs defined}")
  }

  # Markers section
  cli::cli_h2("Markers ({.val {length(x$markers$motifs)}})")
  if (length(x$markers$motifs) > 0) {
    purrr::iwalk(x$markers$motifs, function(motif, i) {
      cli::cli_text("  {.val {i}}: {.val {as.character(motif)}} ({.field {x$markers$alignments[i]}})")
    })
  } else {
    cli::cli_text("  {.emph No marker motifs defined}")
  }

  invisible(x)
}