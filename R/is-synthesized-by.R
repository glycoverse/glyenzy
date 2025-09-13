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
#' Here are some important notes for all functions in the `glyenzy` package.
#'
#' ## Applicability
#'
#' All algorithms and enzyme information in glyenzy are applicable only to humans,
#' and specifically to N-glycans and O-GalNAc glycans.
#' Results may be inaccurate for other types of glycans (e.g., GAGs, glycolipids)
#' or for glycans in other species (e.g., plants, insects).
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
#' ## Substituents
#'
#' Subtituents (e.g. sulfation, phosphorylation) is not supported yet,
#' and the algorithms might fail for glycans with subtituents.
#' If your glycans contains substituents,
#' use [glyrepr::remove_substituents()] to get clean glycans.
#'
#' ## Incomplete glycan structures
#'
#' If the glycan structure is incomplete or partially degraded,
#' the result may be misleading.
#'
#' ## Starting points
#'
#' Throughout `glyenzy`, the starting glycan is the Glc(3)Man(9)GlcNAc(2) precursor for N-glycans,
#' and GalNAc(a1- for O-glycans.
#' This means that enzymes involved in N-glycan precursor biosynthesis, mainly ALGs,
#' and OST, which transfered the precursor to Asn, are not considered here.
#' Similarly, GALNTs for O-glycans are not considered.
#'
#' # Algorithm
#'
#' The basic approach is straightforward: for each reaction rule
#' associated with the enzyme, the function checks whether the
#' corresponding product motif appears in the glycan.
#' If any rule matches, the function returns `TRUE`.
#'
#' For N-glycans, additional logic is applied to handle special cases.
#' Products of **MGAT1** are often further trimmed
#' by glycoside hydrolases, meaning that the final glycan product may no longer
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
  glycans <- .process_glycans_arg(glycans)
  enzyme <- .process_enzyme_arg(enzyme)
  .is_synthesized_by(glycans, enzyme)
}

#' Is a Glycan Synthesized by an Enzyme? (Internal)
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @param is_n A logical vector indicating which elements of `glycans` are N-glycans.
#'   If `NULL`, it will be computed.
#' @noRd
.is_synthesized_by <- function(glycans, enzyme, is_n = NULL) {
  .f <- switch(enzyme$name,
    MGAT1 = .is_synthesized_by_mgat1,
    MOGS = .is_synthesized_by_mogs,
    MAN1B1 = .is_synthesized_by_man1b1,
    MAN1A1 = , MAN1A2 = , MAN1C1 = .is_synthesized_by_man123,
    MAN2A1 = , MAN2A2 = .is_synthesized_by_man2a12,
    GANAB = .is_synthesized_by_ganab,
    .is_synthesized_by_default
  )
  .f(glycans, enzyme, is_n)
}



# Special case for MGAT1
# After MGAT1 adds the b1-2 GlcNAc, two mannoses are trimmed.
# Therefore, checking the product doesn't reflect its involvement.
.is_synthesized_by_mgat1 <- function(glycans, enzyme) {
  glymotif::have_motif(glycans, "GlcNAc(b1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
}
.is_synthesized_by_mgat1 <- .make_n_glycan_guard(.is_synthesized_by_mgat1)

# Special case for MOGS
# This glycoside hydrolase catalyzes only one trimming step.
# If the acceptor motif is not found, it is synthesized by the enzyme.
.is_synthesized_by_mogs <- function(glycans, enzyme) {
  !glymotif::have_motif(glycans, enzyme$rules[[1]]$acceptor)
}
.is_synthesized_by_mogs <- .make_n_glycan_guard(.is_synthesized_by_mogs)

# Special case for MAN1B1
# This glycoside hydrolase catalyzes only one trimming step from Man(9)GlcNAc(2) to Man(8)GlcNAc(2).
# One special case is Man(5)GlcNAc(2), which can and cannot be synthesized by MAN1B1.
# We just assume it is synthesized by MAN1B1, as this is the major route.
# Please check Fig. 115.2 of Handbook of Glycosyltransferases and Related Genes for details.
.is_synthesized_by_man1b1 <- function(glycans, enzyme) {
  !glymotif::have_motif(
    glycans,
    "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    alignment = "core"
  )
}
.is_synthesized_by_man1b1 <- .make_n_glycan_guard(.is_synthesized_by_man1b1)

# Special case for MAN1A1, MAN1A2, and MAN1C1
# These glycoside hydrolases catalyze Man(a1-2) trimming.
.is_synthesized_by_man123 <- function(glycans, enzyme) {
  motifs <- c(
    "a_branch" = "Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "b_branch" = "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "c_branch" = "Man(a1-2)Man(a1-6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  have_motifs_mat <- glymotif::have_motifs(glycans, motifs, alignment = "core")
  n_man <- glyrepr::count_mono(glycans, "Man")
  dplyr::case_when(
    rowSums(have_motifs_mat) == 0L ~ TRUE,
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
# These glycoside hydrolases catalyze Man(a1-3) and Man(a1-6) trimming
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
  fn <- switch(
    enzyme$type,
    GH = .is_synthesized_by_gh,
    GT = .is_synthesized_by_gt,
    cli::cli_abort("Unsupported enzyme type: {enzyme$type}")
  )
  fn(glycans, enzyme)
}

.is_synthesized_by_gh <- function(glycans, enzyme) {
  cli::cli_abort("Glycoside hydrolases except a few involved in N-glycan biosynthesis are not supported yet.")
}

.is_synthesized_by_gt <- function(glycans, enzyme) {
  products <- do.call(c, purrr::map(enzyme$rules, ~ .x$product))
  acceptor_alignments <- purrr::map_chr(enzyme$rules, ~ .x$acceptor_alignment)
  product_alignments <- dplyr::if_else(
    acceptor_alignments %in% c("whole", "core"),
    "core", "substructure"
  )
  have_products_mat <- glymotif::have_motifs(glycans, products, product_alignments)
  unname(rowSums(have_products_mat) > 0)
}