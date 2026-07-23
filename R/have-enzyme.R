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
#' Known-enzyme algorithms and enzyme information in glyenzy are applicable only
#' to humans, and specifically to N-glycans and O-glycans.
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
#' ## Concrete glycans by default
#'
#' Most functions only work for glycans containing **concrete** residues
#' (e.g., `"Glc"`, `"GalNAc"`), and not for glycans with **generic**
#' residues (e.g., `"Hex"`, `"HexNAc"`).
#' Reduced-level inputs with generic residues are supported where explicitly
#' documented, such as `apply_enzyme(structure_level = "basic")`,
#' `trace_biosynthesis()`, and `path_biosynthesis()`.
#'
#' ## Substituents
#'
#' Sulfate substituents are supported. Other substituents, such as
#' phosphorylation and methylation, are not supported. Use
#' [glyrepr::remove_substituents()] when unsupported substituents are present.
#'
#' ## Incomplete glycan structures
#'
#' If the glycan structure is incomplete or partially degraded,
#' the result may be misleading.
#' Glycans with a [glyrepr::get_structure_level()] other than `"intact"`
#' are matched with the lenient motif matching mode in glymotif,
#' and a warning is raised because enzyme predictions may be less reliable.
#'
#' ## Starting points
#'
#' For known-enzyme path inference:
#'
#' - For N-glycans, the starting structure is assumed to be "Glc(3)Man(9)GlcNAc(2)",
#'   the N-glycan precursor transferred to Asn by OST.
#' - For O-GalNAc glycans, the starting structure is assumed to be "GalNAc(a1-".
#' - For O-GlcNAc glycans, the starting structure is assumed to be "GlcNAc(b1-".
#' - For O-Man glycans, the starting structure is assumed to be "Man(a1-".
#' - For O-Fuc glycans, the starting structure is assumed to be "Fuc(a1-".
#' - For O-Glc glycans, the starting structure is assumed to be "Glc(b1-".
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
#' @param method Method used to infer enzyme involvement.
#'   `"motif"` checks product motifs directly in each glycan.
#'   `"path"` extracts enzymes from [trace_biosynthesis()] results, which is
#'   more accurate but slower.
#'
#' @return A logical vector of the same length as `glycans`.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Use `glycan_structure()` and `enzyme()`
#' glycan <- auto_parse("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-")
#' have_enzyme(glycan, enzyme("ST6GAL1"))
#'
#' # Or use characters directly
#' have_enzyme("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-", "ST6GAL1")
#'
#' # Vectorized input
#' glycans <- c(
#'   "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
#'   "Gal(b1-4)GlcNAc(b1-"
#' )
#' have_enzyme(glycans, "ST6GAL1")
#'
#' # Use reconstructed biosynthesis paths
#' have_enzyme(glycans, "ST6GAL1", method = "path")
#'
#' @export
have_enzyme <- function(glycans, enzyme, method = c("motif", "path")) {
  method <- match.arg(method)
  glycans <- .process_glycans_arg(glycans)
  enzyme <- .process_enzyme_arg(enzyme)
  switch(
    method,
    motif = .have_enzyme_motif(glycans, enzyme),
    path = .have_enzyme_path(glycans, enzyme)
  )
}

#' Is a Glycan Synthesized by an Enzyme? (Internal)
#'
#' This function skips argument validation and conversion.
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @noRd
.have_enzyme_motif <- function(glycans, enzyme) {
  UseMethod(".have_enzyme_motif", enzyme)
}

.have_enzyme_motif.glyenzy_npre_gt_enzyme <- function(glycans, enzyme) {
  .is_n_glycan(glycans)
}

.have_enzyme_motif.glyenzy_gh_enzyme <- function(glycans, enzyme) {
  .f <- switch(
    enzyme$name,
    MOGS = .have_enzyme_mogs,
    MAN1B1 = .have_enzyme_man1b1,
    MAN1A1 = ,
    MAN1A2 = ,
    MAN1C1 = .have_enzyme_man123,
    MAN2A1 = ,
    MAN2A2 = .have_enzyme_man2a12,
    GANAB = .have_enzyme_ganab,
    NULL
  )
  if (is.null(.f)) {
    cli::cli_abort(
      "Glycoside hydrolases except a few involved in N-glycan biosynthesis are not supported yet."
    )
  }
  .f(glycans, enzyme)
}

.have_enzyme_motif.glyenzy_gt_enzyme <- function(glycans, enzyme) {
  .f <- switch(
    enzyme$name,
    MGAT1 = .have_enzyme_mgat1,
    .have_enzyme_motif_gt_default
  )
  .f(glycans, enzyme)
}

.have_enzyme_motif.glyenzy_st_enzyme <- function(glycans, enzyme) {
  .have_enzyme_motif_product_default(glycans, enzyme)
}

.have_enzyme_motif_gt_default <- function(glycans, enzyme) {
  .have_enzyme_motif_product_default(glycans, enzyme)
}

#' Check product rules with sulfate-subset and requirement semantics
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#' @returns A logical vector.
#' @noRd
.have_enzyme_motif_product_default <- function(glycans, enzyme) {
  if (length(enzyme$rules) == 0L) {
    return(rep(FALSE, length(glycans)))
  }
  have_products <- purrr::map(
    enzyme$rules,
    function(rule) {
      .have_motif_substituent_subset(
        glycans,
        rule$product,
        alignment = .product_alignment(rule)
      ) &
        .rule_requirements_met(glycans, rule)
    }
  )
  have_products_mat <- do.call(cbind, have_products)
  unname(rowSums(have_products_mat) > 0)
}

# Special case for MGAT1
# After MGAT1 adds the b1-2 GlcNAc, two mannoses are trimmed.
# Therefore, checking the product doesn't reflect its involvement.
.have_enzyme_mgat1 <- function(glycans, enzyme) {
  .have_motif_substituent_subset(
    glycans,
    "GlcNAc(b1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
}
.have_enzyme_mgat1 <- .make_n_glycan_guard(.have_enzyme_mgat1)

# Special case for MOGS
# This glycoside hydrolase catalyzes only one trimming step.
# If the acceptor motif is not found, it is synthesized by the enzyme.
.have_enzyme_mogs <- function(glycans, enzyme) {
  !.have_motif_substituent_subset(
    glycans,
    enzyme$rules[[1]]$acceptor
  )
}
.have_enzyme_mogs <- .make_n_glycan_guard(.have_enzyme_mogs)

# Special case for MAN1B1
# This glycoside hydrolase catalyzes only one trimming step from Man(9)GlcNAc(2) to Man(8)GlcNAc(2).
# One special case is Man(5)GlcNAc(2), which can and cannot be synthesized by MAN1B1.
# We just assume it is synthesized by MAN1B1, as this is the major route.
# Please check Fig. 115.2 of Handbook of Glycosyltransferases and Related Genes for details.
.have_enzyme_man1b1 <- function(glycans, enzyme) {
  !.have_motif_substituent_subset(
    glycans,
    "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    alignment = "core"
  )
}
.have_enzyme_man1b1 <- .make_n_glycan_guard(.have_enzyme_man1b1)

# Special case for MAN1A1, MAN1A2, and MAN1C1
# These glycoside hydrolases catalyze Man(a1-2) trimming.
.have_enzyme_man123 <- function(glycans, enzyme) {
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
  n_man <- glyrepr::count_mono(glycans, "Man")
  dplyr::case_when(
    rowSums(have_motifs_mat) == 0L ~ TRUE,
    n_man == 9 ~ FALSE,
    n_man <= 7 ~ TRUE,
    # 8 mannoses
    TRUE ~ .have_motif_substituent_subset(
      glycans,
      "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      alignment = "core"
    )
  )
}
.have_enzyme_man123 <- .make_n_glycan_guard(.have_enzyme_man123)

# Special case for MAN2A1 and MAN2A2
# These glycoside hydrolases catalyze Man(a1-3) and Man(a1-6) trimming
# after MGAT1 adds the b1-2 GlcNAc.
.have_enzyme_man2a12 <- function(glycans, enzyme) {
  !.have_motif_substituent_subset(
    glycans,
    "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    alignment = "core"
  )
}
.have_enzyme_man2a12 <- .make_n_glycan_guard(.have_enzyme_man2a12)

# Special case for GANAB
# This enzyme removes two Glcs.
.have_enzyme_ganab <- function(glycans, enzyme) {
  !.have_motif_substituent_subset(
    glycans,
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    alignment = "core"
  )
}
.have_enzyme_ganab <- .make_n_glycan_guard(.have_enzyme_ganab)

#' Determine whether each glycan has an enzyme in its traced path
#'
#' @param glycans A `glyrepr_structure` vector.
#' @param enzyme A `glyenzy_enzyme` object.
#'
#' @returns A logical vector.
#' @noRd
.have_enzyme_path <- function(glycans, enzyme) {
  UseMethod(".have_enzyme_path", enzyme)
}

.have_enzyme_path.glyenzy_npre_gt_enzyme <- function(glycans, enzyme) {
  .have_enzyme_motif.glyenzy_npre_gt_enzyme(glycans, enzyme)
}

.have_enzyme_path.glyenzy_enzyme <- function(glycans, enzyme) {
  edges <- .trace_enzyme_edges(glycans, enzymes = .trace_enzymes_with(enzyme))
  unname(purrr::map_lgl(edges, ~ enzyme$name %in% .x))
}
