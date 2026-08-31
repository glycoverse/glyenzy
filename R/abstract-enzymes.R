#' Abstract enzymes for N-glycan biosynthesis
#'
#' Return a separate collection of simplified N-glycan reaction rules. These
#' enzymes represent reaction activities rather than individual genes. They
#' are never included in [db_enzymes()], [enzyme()], [enzymes_from_rnaseq()],
#' or the default enzyme sets of other functions.
#'
#' @param return_str If `FALSE` (default), return a named list of enzyme
#'   objects. If `TRUE`, return their names. These names are labels, not gene
#'   symbols accepted by [enzyme()] or other character enzyme arguments.
#'
#' @returns A named list of `glyenzy_abstract_enzyme` objects, or a character
#'   vector when `return_str = TRUE`. Each object also inherits from
#'   `glyenzy_gt_enzyme` or `glyenzy_gh_enzyme`, and `glyenzy_enzyme`.
#'
#' @details
#' The collection contains ManI, GnTI, ManII, GnTII, a6FucT, GnTIV, GnTV,
#' b4GalT, iGnT, a3SiaT, and GlcH: 11 enzymes with 14 reaction rules.
#' All have `glycan_type = "N"` and `species = "unspecified"`. The additional
#' `localization` field records `"cis"`, `"medial"`, `"trans"`, or `"ER"`;
#' it is descriptive and does not impose compartment order during tracing.
#'
#' The rules use ordinary, site-specific `rejects` motifs. In networks
#' generated from the default precursor, ManII removes the alpha1-3 mannose
#' before the alpha1-6 mannose. Its GlcNAc prerequisite must be terminal,
#' whereas the corresponding a6FucT prerequisite can be extended.
#'
#' iGnT rejects the first LacNAc extension on both the beta1-2 and beta1-4
#' antennae of the alpha1-3 arm, without excluding sites on another arm.
#' Because only iGnT initiates LacNAc extension in this collection, that arm
#' remains unextended in networks generated from the default precursor.
#' The alpha1-6 arm can undergo repeated extensions.
#'
#' These fixed terminal motifs do not enforce arbitrary ancestor paths or
#' node degrees on external inputs. A pre-extended alpha1-3 arm can escape
#' the iGnT rejects; an extended extra branch can escape the ManII reject.
#' Ambiguous linkages use the existing motif-matching semantics, including
#' potentially compatible reject matches. Product-motif statistics do not
#' apply substrate rejects and do not establish that a product is reachable.
#'
#' GlcH combines the three glucose-removal steps represented by MOGS and
#' GANAB in the ordinary database. It removes one residue per step from the
#' complete Glc3Man9GlcNAc2 precursor, retaining Glc2 and Glc1 intermediates.
#' Thus the collection can trace from the usual [trace_biosynthesis()]
#' starting structure without adding ordinary enzymes or changing the start.
#'
#' Pass objects explicitly, for example `abstract_enzymes()[["ManI"]]` to
#' [apply_enzyme()], or the whole list to [trace_biosynthesis()]. Abstract
#' reactions are explicit enzyme rules and have `is_virtual = FALSE` in
#' biosynthesis networks. Existing restrictions on GH involvement statistics
#' and residue matching still apply; this collection does not add GH support
#' to those interfaces.
#'
#' @references
#' The simplified Golgi rules follow the user-supplied reaction table based
#' on the N-glycosylation model of Spahn et al. (2016), with branch notation
#' described by Liang et al. (2020). The supplied table takes precedence over
#' differences in published versions. See the bundled rule-data attribution.
#'
#' Spahn et al. (2016). A Markov chain model for N-linked protein
#' glycosylation. \doi{10.1016/j.ymben.2015.10.007}.
#'
#' Liang et al. (2020). A Markov model of glycosylation elucidates isozyme
#' specificity and glycosyltransferase interactions for glycoengineering.
#' \doi{10.1016/j.crbiot.2020.01.001}.
#'
#' @examples
#' abstract_enzymes(return_str = TRUE)
#' abstract_enzymes()[["ManII"]]
#'
#' man5 <- paste0(
#'   "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]",
#'   "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' )
#' apply_enzyme(man5, abstract_enzymes()[["GnTI"]])
#' trace_biosynthesis(man5, enzymes = abstract_enzymes(), max_virtual_steps = 0L)
#'
#' @export
abstract_enzymes <- function(return_str = FALSE) {
  checkmate::assert_flag(return_str)
  if (return_str) names(glyenzy_abstract_enzymes) else glyenzy_abstract_enzymes
}
