#' Count Enzyme Involvement
#'
#' Count how many times an enzyme is involved in the biosynthesis of a glycan.
#'
#' @inheritSection have_enzyme Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#' @param enzyme An [enzyme()] or a gene symbol.
#'
#' @return An integer vector of the same length as `glycans`.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Use `glycan_structure()` and `enzyme()`
#' glycan <- auto_parse("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-")
#' count_enzyme(glycan, enzyme("ST6GAL1"))
#'
#' # Or use characters directly
#' count_enzyme("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-", "ST6GAL1")
#'
#' # Vectorized input
#' glycans <- c(
#'   "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
#'   "Gal(b1-4)GlcNAc(b1-"
#' )
#' count_enzyme(glycans, "ST6GAL1")
#'
#' @export
count_enzyme <- function(glycans, enzyme) {
  glycans <- .process_glycans_arg(glycans)
  enzyme <- .process_enzyme_arg(enzyme)
  .count_enzyme(glycans, enzyme)
}

#' Count Enzyme Steps (Internal)
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @param is_n A logical vector indicating which elements of `glycans` are N-glycans.
#'   If `NULL`, it will be computed.
#' @noRd
.count_enzyme <- function(glycans, enzyme, is_n = NULL) {
  .f <- switch(
    enzyme$name,
    MGAT1 = .count_enzyme_mgat1,
    MOGS = .count_enzyme_mogs,
    MAN1B1 = .count_enzyme_man1b1,
    MAN1A1 = ,
    MAN1A2 = .count_enzyme_man12,
    MAN1C1 = .count_enzyme_man3,
    MAN2A1 = ,
    MAN2A2 = .count_enzyme_man2a12,
    GANAB = .count_enzyme_ganab,
    .count_enzyme_default
  )
  .f(glycans, enzyme, is_n)
}

# Here we use `.have_enzyme` functions to handle N-glycans.
# See the corresponding functions in `have_enzyme.R` for details.

.count_enzyme_mgat1 <- function(glycans, enzyme, is_n) {
  1L * .have_enzyme_mgat1(glycans, enzyme, is_n)
}

.count_enzyme_mogs <- function(glycans, enzyme, is_n) {
  1L * .have_enzyme_mogs(glycans, enzyme, is_n)
}

.count_enzyme_man1b1 <- function(glycans, enzyme, is_n) {
  1L * .have_enzyme_man1b1(glycans, enzyme, is_n)
}

# Special case for MAN1A1 and MAN1A2
.count_enzyme_man12 <- function(glycans, enzyme) {
  motifs <- c(
    "a_branch" = "Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "b_branch" = "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "c_branch" = "Man(a1-2)Man(a1-6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  have_motifs_mat <- glymotif::have_motifs(glycans, motifs, alignment = "core")
  special_glycan <- "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  is_special <- glymotif::have_motif(
    glycans,
    special_glycan,
    alignment = "whole"
  )
  dplyr::case_when(
    rowSums(have_motifs_mat) == 0L ~ 3L,
    is_special ~ 1L,
    have_motifs_mat[, "b_branch"] ~ 9L - glyrepr::count_mono(glycans, "Man"),
    TRUE ~ 8L - glyrepr::count_mono(glycans, "Man")
  )
}
.count_enzyme_man12 <- .make_n_glycan_guard(
  .count_enzyme_man12,
  type = "integer"
)

.count_enzyme_man3 <- function(glycans, enzyme) {
  motifs <- c(
    "a_branch" = "Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "b_branch" = "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "c_branch" = "Man(a1-2)Man(a1-6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  have_motifs_mat <- glymotif::have_motifs(glycans, motifs, alignment = "core")
  dplyr::case_when(
    rowSums(have_motifs_mat) == 0L ~ 3L,
    have_motifs_mat[, "b_branch"] ~ 9L - glyrepr::count_mono(glycans, "Man"),
    TRUE ~ 8L - glyrepr::count_mono(glycans, "Man")
  )
}
.count_enzyme_man3 <- .make_n_glycan_guard(.count_enzyme_man3, type = "integer")

.count_enzyme_man2a12 <- function(glycans, enzyme) {
  branch_motif <- "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  2L - glymotif::count_motif(glycans, branch_motif, alignment = "core")
}
.count_enzyme_man2a12 <- .make_n_glycan_guard(
  .count_enzyme_man2a12,
  type = "integer"
)

.count_enzyme_ganab <- function(glycans, enzyme) {
  man9_motif <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  res <- dplyr::if_else(
    glymotif::have_motif(glycans, man9_motif, alignment = "core"),
    2L - glymotif::count_motif(glycans, "Glc(a1-"),
    2L
  )
  res[res < 0] <- 0L
  res
}
.count_enzyme_ganab <- .make_n_glycan_guard(
  .count_enzyme_ganab,
  type = "integer"
)

.count_enzyme_default <- function(glycans, enzyme, ...) {
  .count_enzyme_by_type(glycans, enzyme)
}

#' Count enzyme involvement using type-level behavior
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @returns An integer vector.
#' @noRd
.count_enzyme_by_type <- function(glycans, enzyme) {
  UseMethod(".count_enzyme_by_type", enzyme)
}

.count_enzyme_by_type.glyenzy_gh_enzyme <- function(glycans, enzyme) {
  cli::cli_abort(
    "Glycoside hydrolases except a few involved in N-glycan biosynthesis are not supported yet."
  )
}

.count_enzyme_by_type.glyenzy_gt_enzyme <- function(glycans, enzyme) {
  products <- do.call(c, purrr::map(enzyme$rules, ~ .x$product))
  product_alignments <- purrr::map_chr(enzyme$rules, .product_alignment)
  count_products_mat <- glymotif::count_motifs(
    glycans,
    products,
    alignment = product_alignments
  )
  unname(rowSums(count_products_mat))
}

.count_enzyme_by_type.glyenzy_enzyme <- function(glycans, enzyme) {
  switch(
    enzyme$type,
    GH = .count_enzyme_by_type.glyenzy_gh_enzyme(glycans, enzyme),
    GT = .count_enzyme_by_type.glyenzy_gt_enzyme(glycans, enzyme),
    cli::cli_abort("Unsupported enzyme type: {enzyme$type}")
  )
}

.count_enzyme_gh <- .count_enzyme_by_type.glyenzy_gh_enzyme
.count_enzyme_gt <- .count_enzyme_by_type.glyenzy_gt_enzyme
