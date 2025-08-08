#' Determine Whether a Glycan Is Synthesized by a Given Enzyme
#'
#' @description
#' Glycans are produced through a series of enzymatic reactions.
#' This function checks whether a specific enzyme participates
#' in the biosynthesis of a given glycan (or glycans).
#'
#' @details
#'
#' # Important notes
#'
#' ## Inclusiveness
#'
#' The algorithm takes an intentionally inclusive approach,
#' assuming that all possible isoenzymes capable of catalyzing
#' a given reaction may be involved.
#' Therefore, results should be interpreted with caution.
#'
#' For example, in humans, detection of the motif "Neu5Ac(a2-3)Gal(b1-" will return
#' both "ST3GAL3" and "ST3GAL4". In reality, only one of them might be active, depending
#' on factors such as tissue specificity.
#'
#' ## Only "concrete" glycans
#'
#' The function only works for glycans containing **concrete** residues
#' (e.g., `"Glc"`, `"GalNAc"`), and not for glycans with **generic**
#' residues (e.g., `"Hex"`, `"HexNAc"`).
#'
#' ## Incomplete glycan structures
#'
#' If the glycan structure is incomplete or partially degraded,
#' the result may be misleading.
#'
#' # Algorithm
#'
#' The basic approach is straightforward: for each reaction rule
#' associated with the enzyme, the function checks whether the
#' corresponding product motif appears in the glycan.
#' If any rule matches, the function returns `TRUE`.
#'
#' For N-glycans, additional logic is applied to handle special cases.
#' Products of **ALG** enzymes and **MGAT1** are often further trimmed
#' by exoglycosidases, meaning that the final glycan product may no longer
#' contain the original motif.
#' In these cases, the function instead looks for specific motif markers
#' to determine enzyme involvement.
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
      "{.arg glycans} must be a {.cls glyrepr_structure} vector or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(glycans)}}."
    ))
  }
  # Process `enzyme` argument
  if (is.character(enzyme)) {
    enzyme <- glyenzy::enzyme(enzyme)
  } else if (!inherits(enzyme, "glyenzy_enzyme")) {
    cli::cli_abort(c(
      "{.arg enzyme} must be a {.cls glyenzy_enzyme} object or a character string of gene symbol.",
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
#' @param is_n A logical vector indicating which elements of `glycans` are N-glycans.
#'   If `NULL`, it will be computed by [glymotif::is_n_glycan()].
#' @noRd
.is_synthesized_by <- function(glycans, enzyme, is_n = NULL) {
  .f <- switch(enzyme$name,
    ALG1 = , ALG2 = , ALG3 = , ALG6 = , ALG8 = , ALG9 = , ALG10 = ,
    ALG11 = , ALG12 = , ALG13 = , ALG14 = , DPAGT1 = .is_synthesized_by_alg,
    MGAT1 = .is_synthesized_by_mgat1,
    MOGS = , MAN1B1 = .is_synthesized_by_mogs_man1b1,
    MAN1A1 = , MAN1A2 = , MAN1C1 = .is_synthesized_by_man123,
    MAN2A1 = , MAN2A2 = .is_synthesized_by_man2a12,
    GANAB = .is_synthesized_by_ganab,
    .is_synthesized_by_default
  )
  .f(glycans, enzyme, is_n)
}

#' Make a function that only applies to N-glycans
#'
#' This function is used to wrap a function that only applies to N-glycans.
#' For N-glycans in the input `glycans`, it uses the values in `.f`.
#' Otherwise, it returns FALSE.
#'
#' @param .f A function that takes two arguments: `glycans` and `enzyme`.
#' @noRd
.make_n_glycan_guard <- function(.f) {
  force(.f)
  function(glycans, enzyme, is_n = NULL) {
    if (is.null(is_n)) {
      is_n <- glymotif::is_n_glycan(glycans)
    }
    res <- rep(FALSE, length(glycans))
    if (any(is_n)) {
      res[is_n] <- .f(glycans[is_n], enzyme)
    }
    res
  }
}

# Special case for ALG family and DPAGT1
# These enzymes involve in the the biosynthesis process of all N-glycans.
.is_synthesized_by_alg <- function(glycans, enzyme) {
  rep(TRUE, length(glycans))
}
.is_synthesized_by_alg <- .make_n_glycan_guard(.is_synthesized_by_alg)

# Special case for MGAT1
# After MGAT1 adds the b1-2 GlcNAc, two mannoses are trimmed.
# Therefore, checking the product doesn't reflect its involvement.
.is_synthesized_by_mgat1 <- function(glycans, enzyme) {
  glymotif::have_motif(glycans, "GlcNAc(b1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
}
.is_synthesized_by_mgat1 <- .make_n_glycan_guard(.is_synthesized_by_mgat1)

# Special case for MOGS and MAN1B1
# These exoglycosidases catalyze only one trimming step.
# If the acceptor motif is not found, it is synthesized by the enzyme.
.is_synthesized_by_mogs_man1b1 <- function(glycans, enzyme) {
  !glymotif::have_motif(glycans, enzyme$rules[[1]]$acceptor)
}
.is_synthesized_by_mogs_man1b1 <- .make_n_glycan_guard(.is_synthesized_by_mogs_man1b1)

# Special case for MAN1A1, MAN1A2, and MAN1C1
# These exoglycosidases catalyze Man(a1-2) trimming.
.is_synthesized_by_man123 <- function(glycans, enzyme) {
  glycan_type <- glymotif::n_glycan_type(glycans)
  n_man <- glyrepr::count_mono(glycans, "Man")
  dplyr::case_when(
    glycan_type %in% c("complex", "hybrid", "paucimannose") ~ TRUE,
    n_man == 9 ~ FALSE,
    n_man <= 7 ~ TRUE,
    # 8 mannoses
    TRUE ~ glymotif::have_motif(
      glycans,
      "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      alignment = "core"
    )
  )
}
.is_synthesized_by_man123 <- .make_n_glycan_guard(.is_synthesized_by_man123)

# Special case for MAN2A1 and MAN2A2
# These exoglycosidases catalyze Man(a1-3) and Man(a1-6) trimming
# after MGAT1 adds the b1-2 GlcNAc.
.is_synthesized_by_man2a12 <- function(glycans, enzyme) {
  !glymotif::have_motif(
    glycans,
    "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    alignment = "core"
  )
}
.is_synthesized_by_man2a12 <- .make_n_glycan_guard(.is_synthesized_by_man2a12)

# Special case for GANAB
# This enzyme removes two Glcs.
.is_synthesized_by_ganab <- function(glycans, enzyme) {
  !glymotif::have_motif(
    glycans,
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    alignment = "core"
  )
}
.is_synthesized_by_ganab <- .make_n_glycan_guard(.is_synthesized_by_ganab)

#' Is a Glycan Synthesized by an Enzyme? (Default Case)
#'
#' Check all rules of the `enzyme`.
#' Returns TRUE if any rule is satisfied.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @param ... Ignored.
#' @noRd
.is_synthesized_by_default <- function(glycans, enzyme, ...) {
  if (enzyme$type == "GD") {
    cli::cli_abort("Exoglycosidases except a few involved in N-glycan biosynthesis are not supported yet.")
  }
  products <- do.call(c, purrr::map(enzyme$rules, ~ .x$product))
  have_products_mat <- glymotif::have_motifs(glycans, products)
  unname(rowSums(have_products_mat) > 0)
}