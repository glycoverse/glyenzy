#' Count Enzyme Involvement
#'
#' Count how many times an enzyme is involved in the biosynthesis of a glycan.
#'
#' @inheritSection have_enzyme Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#' @param enzyme An [enzyme()] or a gene symbol.
#' @param method Method used to count enzyme involvement.
#'   `"motif"` counts product motifs directly in each glycan.
#'   `"path"` counts enzyme-labeled edges in [trace_biosynthesis()] results,
#'   which is more accurate but slower.
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
#' # Use reconstructed biosynthesis paths
#' count_enzyme(glycans, "ST6GAL1", method = "path")
#'
#' @export
count_enzyme <- function(glycans, enzyme, method = c("motif", "path")) {
  method <- match.arg(method)
  glycans <- .process_glycans_arg(glycans)
  enzyme <- .process_enzyme_arg(enzyme)
  switch(
    method,
    motif = .count_enzyme_motif(glycans, enzyme),
    path = .count_enzyme_path(glycans, enzyme)
  )
}

#' Count Enzyme Steps (Internal)
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @noRd
.count_enzyme_motif <- function(glycans, enzyme) {
  UseMethod(".count_enzyme_motif", enzyme)
}

.count_enzyme_motif.glyenzy_npre_gt_enzyme <- function(glycans, enzyme) {
  n_steps <- length(enzyme$rules)
  res <- rep(0L, length(glycans))
  res[.is_n_glycan(glycans)] <- n_steps
  res
}

.count_enzyme_motif.glyenzy_gh_enzyme <- function(glycans, enzyme) {
  .f <- switch(
    enzyme$name,
    MOGS = .count_enzyme_mogs,
    MAN1B1 = .count_enzyme_man1b1,
    MAN1A1 = ,
    MAN1A2 = .count_enzyme_man12,
    MAN1C1 = .count_enzyme_man3,
    MAN2A1 = ,
    MAN2A2 = .count_enzyme_man2a12,
    GANAB = .count_enzyme_ganab,
    NULL
  )
  if (is.null(.f)) {
    cli::cli_abort(
      "Glycoside hydrolases except a few involved in N-glycan biosynthesis are not supported yet."
    )
  }
  .f(glycans, enzyme)
}

.count_enzyme_motif.glyenzy_gt_enzyme <- function(glycans, enzyme) {
  .f <- switch(
    enzyme$name,
    MGAT1 = .count_enzyme_mgat1,
    .count_enzyme_motif_gt_default
  )
  .f(glycans, enzyme)
}

.count_enzyme_motif.glyenzy_st_enzyme <- function(glycans, enzyme) {
  lengths(.match_st_enzyme_events(glycans, enzyme))
}

#' Count trace-derived enzyme edges for each glycan
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @returns An integer vector.
#' @noRd
.count_enzyme_path <- function(glycans, enzyme) {
  UseMethod(".count_enzyme_path", enzyme)
}

.count_enzyme_path.glyenzy_npre_gt_enzyme <- function(glycans, enzyme) {
  .count_enzyme_motif.glyenzy_npre_gt_enzyme(glycans, enzyme)
}

.count_enzyme_path.glyenzy_st_enzyme <- function(glycans, enzyme) {
  lengths(.match_st_enzyme_path_events(glycans, enzyme))
}

.count_enzyme_path.glyenzy_enzyme <- function(glycans, enzyme) {
  edges <- .trace_enzyme_edges(glycans, enzymes = .trace_enzymes_with(enzyme))
  unname(purrr::map_int(edges, ~ sum(.x == enzyme$name)))
}

# Here we use `.have_enzyme_motif` functions to handle N-glycans.
# See the corresponding functions in `have_enzyme.R` for details.

.count_enzyme_mgat1 <- function(glycans, enzyme) {
  1L * .have_enzyme_mgat1(glycans, enzyme)
}

.count_enzyme_mogs <- function(glycans, enzyme) {
  1L * .have_enzyme_mogs(glycans, enzyme)
}

.count_enzyme_man1b1 <- function(glycans, enzyme) {
  1L * .have_enzyme_man1b1(glycans, enzyme)
}

# Special case for MAN1A1 and MAN1A2
.count_enzyme_man12 <- function(glycans, enzyme) {
  motifs <- c(
    "a_branch" = "Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "b_branch" = "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "c_branch" = "Man(a1-2)Man(a1-6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  have_motifs_mat <- .have_motifs_substituent_subset(
    glycans,
    motifs,
    alignments = "core"
  )
  special_glycan <- "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  is_special <- .have_motif_substituent_subset(
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
  have_motifs_mat <- .have_motifs_substituent_subset(
    glycans,
    motifs,
    alignments = "core"
  )
  dplyr::case_when(
    rowSums(have_motifs_mat) == 0L ~ 3L,
    have_motifs_mat[, "b_branch"] ~ 9L - glyrepr::count_mono(glycans, "Man"),
    TRUE ~ 8L - glyrepr::count_mono(glycans, "Man")
  )
}
.count_enzyme_man3 <- .make_n_glycan_guard(.count_enzyme_man3, type = "integer")

.count_enzyme_man2a12 <- function(glycans, enzyme) {
  branch_motif <- "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  2L -
    .count_motif_substituent_subset(
      glycans,
      branch_motif,
      alignment = "core"
    )
}
.count_enzyme_man2a12 <- .make_n_glycan_guard(
  .count_enzyme_man2a12,
  type = "integer"
)

.count_enzyme_ganab <- function(glycans, enzyme) {
  man9_motif <- "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  res <- dplyr::if_else(
    .have_motif_substituent_subset(
      glycans,
      man9_motif,
      alignment = "core"
    ),
    2L - .count_motif_substituent_subset(glycans, "Glc(a1-"),
    2L
  )
  res[res < 0] <- 0L
  res
}
.count_enzyme_ganab <- .make_n_glycan_guard(
  .count_enzyme_ganab,
  type = "integer"
)

.count_enzyme_motif_gt_default <- function(glycans, enzyme) {
  if (length(enzyme$rules) == 0L) {
    return(integer(length(glycans)))
  }
  counts <- purrr::map(
    enzyme$rules,
    function(rule) {
      rule_counts <- .count_motif_substituent_subset(
        glycans,
        rule$product,
        alignment = .product_alignment(rule)
      )
      requirements_met <- .rule_requirements_met(glycans, rule)
      rule_counts[!requirements_met] <- 0L
      rule_counts
    }
  )
  count_products_mat <- do.call(cbind, counts)
  unname(rowSums(count_products_mat))
}
