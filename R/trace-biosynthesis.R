#' Trace the Biosynthetic Path of Glycans
#'
#' Reconstruct biosynthetic pathways for one or more glycans. The default
#' method uses known enzyme rules in a forward, multi-target breadth-first
#' search. The virtual and hybrid methods instead trim targets backward to infer
#' every possible residue-addition order. Hybrid tracing additionally annotates
#' each virtual transition with concrete enzymes whose rules can perform it.
#'
#' @inheritSection have_enzyme Important notes
#'
#' @section Virtual enzymes:
#' With `method = "virtual"`, each edge is named for the residue added by that
#' step. Intact glycans include the linkage anomer and acceptor position, so a
#' beta-1,4-linked GlcNAc is labeled `"b4GlcNAcT"`. Partial and topological
#' glycans omit linkage information and use `"GlcNAcT"`; basic glycans use the
#' generic residue name, such as `"HexNAcT"`.
#'
#' Virtual tracing starts N-glycans at the N-glycan core and all other glycans
#' at their reducing-end root residue. In [path_biosynthesis()], the explicit
#' `from` glycan is always the virtual starting structure. These networks do not
#' apply organism-specific substrate rules and represent structural
#' possibilities rather than biological feasibility.
#'
#' With `method = "hybrid"`, the same virtual network is returned and every edge
#' gains a list-valued `concrete_enzymes` attribute. It contains all candidate
#' enzymes whose rules can transform that edge's substrate into its product.
#' The vector is empty when no candidate rule supports the transition. The
#' scalar `enzyme` attribute continues to contain the virtual-enzyme name.
#'
#' Basic structures do not retain glycan-class metadata. A basic structure
#' matching the generic N-glycan-core topology is therefore assumed to be an
#' N-glycan; use [path_biosynthesis()] with an explicit `from` when that
#' topology belongs to another glycan class.
#'
#' @param glycans A [glyrepr::glycan_structure()] vector, or a character vector
#'   of strings supported by [glyparse::auto_parse()]. Can also be a single glycan.
#'   If multiple glycans are provided, the starting structure will be decided by the first glycan.
#'   Therefore, please make sure `glycans` are not of mixed glycan types.
#' @param enzymes A character vector of gene symbols, or a list of [enzyme()]
#'   objects. If `NULL` (default), all available enzymes will be used. With
#'   `method = "hybrid"`, this selects the candidate enzymes used to annotate
#'   virtual transitions. Must be `NULL` when `method = "virtual"`.
#' @param max_steps Integer, maximum number of enzymatic steps to search.
#'   Default is 20.
#' @param filter Optional function to filter generated glycans at each step.
#'   Should take a [glyrepr::glycan_structure()] vector as input and return
#'   a logical vector of the same length.
#'   For `method = "enzymatic"`, it filters generated products. For
#'   `method = "virtual"` or `method = "hybrid"`, it filters generated
#'   precursors during backward trimming.
#' @param method Biosynthesis inference method. `"enzymatic"` (default) uses
#'   known enzyme rules in a forward search. `"virtual"` uses virtual enzymes
#'   in a backward search that removes terminal residues from the targets.
#'   `"hybrid"` builds the same virtual network and annotates each transition
#'   with concrete enzyme candidates.
#'
#' @returns An [igraph::igraph()] object representing the synthesis path(s).
#'   Vertices represent glycan structures, with IUPAC-condensed strings in the
#'   `name` attribute. Every edge has a `step` attribute indicating the forward
#'   synthesis step. Other edge attributes depend on `method`:
#'
#'   - `method = "enzymatic"`: returns a multi-edge graph. Each edge represents
#'     one concrete enzyme, whose gene symbol is stored in `enzyme`. When
#'     multiple enzymes catalyze the same substrate-to-product transition, they
#'     are represented by parallel edges.
#'   - `method = "virtual"`: has one edge per unique substrate-to-product
#'     transition. `enzyme` contains the virtual-enzyme name inferred from the
#'     residue and linkage added in that transition.
#'   - `method = "hybrid"`: has the same single-edge topology as `"virtual"`.
#'     `enzyme` contains the virtual-enzyme name, and `concrete_enzymes` is a
#'     list of character vectors containing every candidate concrete enzyme.
#'     The vector is empty when no candidate rule supports that transition.
#'
#'   For multiple targets, the graph includes all synthesis paths needed to
#'   reach every target glycan.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Rebuild the biosynthetic pathway of a single glycan
#' glycan <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' path <- trace_biosynthesis(glycan, max_steps = 20)
#'
#' # Rebuild pathways for multiple glycans
#' glycans <- c(
#'   "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-",
#'   "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
#' )
#' path <- trace_biosynthesis(glycans, max_steps = 20)
#'
#' # Build an enzyme-agnostic network using virtual enzymes
#' virtual_path <- trace_biosynthesis(
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
#'   method = "virtual"
#' )
#'
#' # Annotate virtual transitions with exact rule-matched enzyme candidates
#' hybrid_path <- trace_biosynthesis(
#'   "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
#'   method = "hybrid"
#' )
#'
#' # View the path
#' igraph::as_data_frame(path, what = "edges")
#'
#' @export
trace_biosynthesis <- function(
  glycans,
  enzymes = NULL,
  max_steps = 20,
  filter = NULL,
  method = c("enzymatic", "virtual", "hybrid")
) {
  # Parse and validate basic inputs first
  method <- match.arg(method)
  glycans <- .process_glycans_arg(glycans, allow_generic = TRUE)
  checkmate::assert_int(max_steps, lower = 1)
  if (!is.null(filter)) {
    filter <- rlang::as_function(filter)
  }

  if (method %in% c("virtual", "hybrid")) {
    concrete_enzymes <- NULL
    if (identical(method, "virtual")) {
      .validate_virtual_enzymes(enzymes)
    } else {
      concrete_enzymes <- .process_enzymes_arg(
        enzymes,
        apply_prefilter = FALSE
      )
    }
    starting_glycan <- .decide_virtual_starting_glycan(glycans[1])
    path <- .perform_virtual_synthesis(
      starting_glycan,
      glycans,
      max_steps,
      filter
    )
    if (identical(method, "hybrid")) {
      path <- .amplify_virtual_edges(path, concrete_enzymes)
    }
    return(path)
  }

  enzymes <- .process_enzymes_arg(
    enzymes,
    glycans = glycans,
    apply_prefilter = TRUE
  )

  # Find all possible paths using unified BFS logic
  starting_glycan <- .decide_starting_glycan(glycans[1]) # Use first glycan to decide starting point
  .perform_bfs_synthesis(starting_glycan, glycans, enzymes, max_steps, filter)
}

.decide_starting_glycan <- function(glycan) {
  if (.can_reliably_detect_n_glycan(glycan) && .is_n_glycan(glycan)) {
    start <- glyparse::parse_iupac_condensed(
      "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    )
  } else if (.have_motif(glycan, "GalNAc(a1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("GalNAc(a1-")
  } else if (.have_motif(glycan, "GlcNAc(b1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("GlcNAc(b1-")
  } else if (.have_motif(glycan, "Man(a1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("Man(a1-")
  } else if (.have_motif(glycan, "Fuc(a1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("Fuc(a1-")
  } else if (.have_motif(glycan, "Glc(b1-", alignment = "core")) {
    start <- glyparse::parse_iupac_condensed("Glc(b1-")
  } else {
    cli::cli_abort(c(
      "Cannot decide the starting point for the given glycan(s).",
      "i" = "Currenly, only N-glycans and O-glycans are supported."
    ))
  }
  start
}

.can_reliably_detect_n_glycan <- function(glycan) {
  !identical(.glycan_structure_level(glycan), "basic")
}
