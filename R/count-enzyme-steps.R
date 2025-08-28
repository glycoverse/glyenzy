#' Count Enzyme Steps
#'
#' Count how many times an enzyme is involved in the biosynthesis of a glycan.
#'
#' @inheritSection is_synthesized_by Important notes
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
#' count_enzyme_steps(glycan, enzyme("ST6GAL1"))
#'
#' # Or use characters directly
#' count_enzyme_steps("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-", "ST6GAL1")
#'
#' # Vectorized input
#' glycans <- c(
#'   "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
#'   "Gal(b1-4)GlcNAc(b1-"
#' )
#' count_enzyme_steps(glycans, "ST6GAL1")
#'
#' @export
count_enzyme_steps <- function(glycans, enzyme) {
  glycans <- .process_glycans_arg(glycans)
  enzyme <- .process_enzyme_arg(enzyme)
  .count_enzyme_steps(glycans, enzyme)
}

#' Count Enzyme Steps (Internal)
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @param is_n A logical vector indicating which elements of `glycans` are N-glycans.
#'   If `NULL`, it will be computed by [glymotif::is_n_glycan()].
#' @noRd
.count_enzyme_steps <- function(glycans, enzyme, is_n = NULL) {
  .f <- switch(enzyme$name,
    MGAT1 = .count_enzyme_steps_mgat1,
    MOGS = .count_enzyme_steps_mogs,
    MAN1B1 = .count_enzyme_steps_man1b1,
    MAN1A1 = ,
    MAN1A2 = .count_enzyme_steps_man12,
    MAN1C1 = .count_enzyme_steps_man3,
    MAN2A1 = ,
    MAN2A2 = .count_enzyme_steps_man2a12,
    GANAB = .count_enzyme_steps_ganab,
    .count_enzyme_steps_default
  )
  .f(glycans, enzyme, is_n)
}

# Here we use `.is_synthesized_by` functions to handle N-glycans.
# See the corresponding functions in `is_synthesized_by.R` for details.



.count_enzyme_steps_mgat1 <- function(glycans, enzyme, is_n) {
  1L * .is_synthesized_by_mgat1(glycans, enzyme, is_n)
}

.count_enzyme_steps_mogs <- function(glycans, enzyme, is_n) {
  1L * .is_synthesized_by_mogs(glycans, enzyme, is_n)
}

.count_enzyme_steps_man1b1 <- function(glycans, enzyme, is_n) {
  1L * .is_synthesized_by_man1b1(glycans, enzyme, is_n)
}

# Special case for MAN1A1 and MAN1A2
.count_enzyme_steps_man12 <- function(glycans, enzyme) {
  b_branch_motif <- "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  special_glycan <- "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  dplyr::case_when(
    glymotif::n_glycan_type(glycans) != "highmannose" ~ 3L,
    glymotif::have_motif(glycans, special_glycan, alignment = "whole") ~ 1L,
    glymotif::have_motif(glycans, b_branch_motif, alignment = "core") ~ 9L - glyrepr::count_mono(glycans, "Man"),
    TRUE ~ 8L - glyrepr::count_mono(glycans, "Man")
  )
}
.count_enzyme_steps_man12 <- .make_n_glycan_guard(.count_enzyme_steps_man12, type = "integer")

.count_enzyme_steps_man3 <- function(glycans, enzyme) {
  b_branch_motif <- "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  dplyr::case_when(
    glymotif::n_glycan_type(glycans) != "highmannose" ~ 3L,
    glymotif::have_motif(glycans, b_branch_motif, alignment = "core") ~ 9L - glyrepr::count_mono(glycans, "Man"),
    TRUE ~ 8L - glyrepr::count_mono(glycans, "Man")
  )
}
.count_enzyme_steps_man3 <- .make_n_glycan_guard(.count_enzyme_steps_man3, type = "integer")

.count_enzyme_steps_man2a12 <- function(glycans, enzyme) {
  branch_motif <- "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  2L - glymotif::count_motif(glycans, branch_motif, alignment = "core")
}
.count_enzyme_steps_man2a12 <- .make_n_glycan_guard(.count_enzyme_steps_man2a12, type = "integer")

.count_enzyme_steps_ganab <- function(glycans, enzyme) {
  man9_motif <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  res <- dplyr::if_else(
    glymotif::have_motif(glycans, man9_motif, alignment = "core"),
    2L - glymotif::count_motif(glycans, "Glc(a1-"),
    2L
  )
  res[res < 0] <- 0L
  res
}
.count_enzyme_steps_ganab <- .make_n_glycan_guard(.count_enzyme_steps_ganab, type = "integer")

.count_enzyme_steps_default <- function(glycans, enzyme, ...) {
  fn <- switch(
    enzyme$type,
    GH = .count_enzyme_steps_gh,
    GT = .count_enzyme_steps_gt,
    cli::cli_abort("Unsupported enzyme type: {enzyme$type}")
  )
  fn(glycans, enzyme)
}

.count_enzyme_steps_gh <- function(glycans, enzyme) {
  cli::cli_abort("Glycoside hydrolases except a few involved in N-glycan biosynthesis are not supported yet.")
}

.count_enzyme_steps_gt <- function(glycans, enzyme) {
  products <- do.call(c, purrr::map(enzyme$rules, ~ .x$product))
  count_products_mat <- glymotif::count_motifs(glycans, products)
  unname(rowSums(count_products_mat))
}