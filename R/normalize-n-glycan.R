#' Normalize Topological N-Glycans
#'
#' Normalize selected Gal and Neu5Ac distributions in topological N-glycans.
#' For the supported cases, glycans with the same numbers of antennae, Gal,
#' and Neu5Ac are mapped to a common representative topology.
#'
#' @details
#' # Why does this function exist?
#'
#' [trace_biosynthesis()] and related functions support topological glycans.
#' However, N-glycan annotations can assign residues to specific antennae even
#' when the available data do not resolve those positions. These assignments
#' can introduce unsupported distinctions into a biosynthesis network.
#'
#' For example, consider the two glycans below.
#'
#' Glycan 1:
#' ```
#' Neu5Ac - Gal - GlcNAc
#'                       \
#' Neu5Ac - Gal - GlcNAc - Man
#'                            \
#'                             Man - GlcNAc - GlcNAc -
#'                            /
#'          Gal - GlcNAc - Man
#' ```
#'
#' Glycan 2:
#' ```
#' Neu5Ac - Gal - GlcNAc
#'                       \
#'          Gal - GlcNAc - Man
#'                            \
#'                             Man - GlcNAc - GlcNAc -
#'                            /
#' Neu5Ac - Gal - GlcNAc - Man
#' ```
#'
#' Both are tri-antennary glycans with two sialic acids. When mass spectrometry
#' data do not resolve which antennae carry the sialic acids, either topology
#' may be used to annotate the same unresolved structure. Nevertheless,
#' tracing these annotations can produce different biosynthetic paths and
#' introduce redundant branches when multiple glycans are traced together.
#'
#' Glycan 3 below illustrates the problem. It has one sialic acid and could
#' serve as a precursor to a disialylated glycan. Of the two topologies above,
#' only Glycan 2 can be obtained from Glycan 3 by adding a single sialic acid.
#' Tracing Glycan 1 and Glycan 3 together therefore requires separate paths
#' rather than connecting them through a single sialylation step. When the
#' sialic acid positions are unresolved, this distinction reflects the
#' annotation rather than evidence for separate biosynthetic routes.
#'
#' Glycan 3:
#' ```
#'          Gal - GlcNAc
#'                       \
#'          Gal - GlcNAc - Man
#'                            \
#'                             Man - GlcNAc - GlcNAc -
#'                            /
#' Neu5Ac - Gal - GlcNAc - Man
#' ```
#'
#' This function addresses the problem through a normalization convention.
#' It applies rules to reposition Gal and Neu5Ac. In this example, mapping
#' Glycan 1 to Glycan 2 allows the monosialylated Glycan 3 to connect to the
#' disialylated representative through a single sialylation step. This
#' convention does not imply that distinct isomers are biologically equivalent.
#'
#' # Algorithm
#'
#' The algorithm repositions only Gal and Neu5Ac residues.
#' To explain the algorithm, we name the branches for tri- and tetra-antennary
#' glycans.
#'
#' Tri-antennary:
#'
#' ```
#' A: GlcNAc
#'           \
#' B: GlcNAc - Man
#'                 \
#'                  Man - GlcNAc - GlcNAc -
#'                 /
#' C: GlcNAc - Man
#' ```
#'
#' Tetra-antennary:
#'
#' ```
#' A: GlcNAc
#'           \
#' B: GlcNAc - Man
#'                 \
#'                  Man - GlcNAc - GlcNAc -
#'                 /
#' C: GlcNAc - Man
#'           /
#' D: GlcNAc
#' ```
#'
#' First, Gal residues are repositioned in the following cases:
#' - Three antennae, one Gal: positioned on C.
#' - Three antennae, two Gal: positioned on A and C.
#' - Four antennae, two Gal: positioned on A and C.
#'
#' Next, Neu5Ac residues are repositioned in the following cases:
#' - Three antennae, two Gal, one Neu5Ac: positioned on the C Gal.
#' - Three antennae, three Gal, one Neu5Ac: positioned on the C Gal.
#' - Three antennae, three Gal, two Neu5Ac: positioned on the A and C Gal.
#' - Four antennae, three Gal, one Neu5Ac: with Gal on A, B, and C,
#'   Neu5Ac is positioned on the C Gal.
#' - Four antennae, three Gal, two Neu5Ac: with Gal on A, B, and C,
#'   Neu5Ac residues are positioned on the A and C Gal.
#' - Four antennae, four Gal, two Neu5Ac: positioned on the A and C Gal.
#'
#' Other distributions are left unchanged. Branch labels describe topology,
#' not linkage positions or drawing order. For three antennae, C is on the
#' Man with one antenna; A and B are interchangeable. For four antennae,
#' A/B and C/D are the two pairs. With three Gal, A/B is the pair with two
#' Gal, and C is the occupied antenna of the other pair.
#'
#' Only concrete, topological structures are accepted. Use
#' [glyrepr::remove_linkages()] explicitly before calling this function on
#' intact or partial structures. Missing values are preserved.
#'
#' Only simple GlcNAc, Gal-GlcNAc, and Neu5Ac-Gal-GlcNAc antennae are
#' normalized. Non-N-glycans, other antenna counts, and structures with
#' extended or modified antennae, other sialic acids, or floating parts or
#' substituents are returned unchanged. Core fucose on the reducing-end
#' GlcNAc and a bisecting GlcNAc on the central Man are preserved.
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#'
#' @returns A [glyrepr::glycan_structure()] vector of normalized glycans.
#' @examples
#' x <- glyparse::auto_parse(paste0(
#'   "Gal(b1-4)GlcNAc(b1-2)[GlcNAc(b1-6)]Man(a1-6)",
#'   "[GlcNAc(b1-2)Man(a1-3)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' ))
#' normalize_n_glycan(glyrepr::remove_linkages(x))
#'
#' @export
normalize_n_glycan <- function(glycans) {
  input_names <- names(glycans)
  if (is.character(glycans)) {
    glycans <- glyparse::auto_parse(glycans)
  } else if (!glyrepr::is_glycan_structure(glycans)) {
    cli::cli_abort(
      "{.arg glycans} must be a glycan structure vector or a character vector."
    )
  }
  present <- !is.na(glycans)
  if (any(glyrepr::get_mono_type(glycans)[present] != "concrete")) {
    cli::cli_abort("{.arg glycans} must use concrete monosaccharides.")
  }
  if (any(.glycan_structure_levels(glycans)[present] != "topological")) {
    cli::cli_abort(c(
      "{.arg glycans} must contain only topological structures.",
      "i" = "Use {.fn glyrepr::remove_linkages} on intact or partial structures first."
    ))
  }
  result <- glyrepr::smap_structure(glycans, .normalize_n_glycan_graph)
  attr(result, "names") <- input_names
  result
}

# Recognize the entire supported scaffold before moving any residues.
.normalization_n_arms <- function(graph) {
  if (
    length(igraph::graph_attr(graph, "floating_parts")) > 0L ||
      length(igraph::graph_attr(graph, "floating_substituents")) > 0L ||
      any(nzchar(igraph::vertex_attr(graph, "sub")))
  ) {
    return(NULL)
  }
  mono <- igraph::vertex_attr(graph, "mono")
  children <- igraph::adjacent_vertices(graph, igraph::V(graph), mode = "out")
  children <- lapply(children, as.integer)
  root <- which(igraph::degree(graph, mode = "in") == 0L)
  if (length(root) != 1L || mono[root] != "GlcNAc") {
    return(NULL)
  }
  root_children <- children[[root]]
  inner <- root_children[mono[root_children] == "GlcNAc"]
  fuc <- root_children[mono[root_children] == "Fuc"]
  if (
    length(inner) != 1L ||
      length(fuc) > 1L ||
      length(root_children) != 1L + length(fuc) ||
      any(lengths(children[fuc]) != 0L)
  ) {
    return(NULL)
  }
  central <- children[[inner]]
  if (length(central) != 1L || mono[central] != "Man") {
    return(NULL)
  }
  central_children <- children[[central]]
  sides <- central_children[mono[central_children] == "Man"]
  bisect <- central_children[mono[central_children] == "GlcNAc"]
  if (
    length(sides) != 2L ||
      length(bisect) > 1L ||
      length(central_children) != 2L + length(bisect) ||
      any(lengths(children[bisect]) != 0L)
  ) {
    return(NULL)
  }
  arms <- children[sides]
  sizes <- lengths(arms)
  # Three antennae split 2 + 1 between the two Man residues; four split 2 + 2.
  if (!all(sort(sizes) == c(1L, 2L)) && !all(sizes == 2L)) {
    return(NULL)
  }
  tips <- unlist(arms, use.names = FALSE)
  if (any(mono[tips] != "GlcNAc")) {
    return(NULL)
  }
  for (tip in tips) {
    gal <- children[[tip]]
    if (length(gal) == 0L) {
      next
    }
    if (length(gal) != 1L || mono[gal] != "Gal") {
      return(NULL)
    }
    sia <- children[[gal]]
    if (length(sia) == 0L) {
      next
    }
    if (
      length(sia) != 1L ||
        mono[sia] != "Neu5Ac" ||
        length(children[[sia]]) != 0L
    ) {
      return(NULL)
    }
  }
  # Put the double-antenna (or double-Gal) side first. Within each side,
  # occupied antennae come first; remaining ties are topologically symmetric.
  arms <- lapply(arms, function(ids) {
    ids[order(lengths(children[ids]), decreasing = TRUE)]
  })
  counts <- vapply(arms, function(ids) sum(lengths(children[ids])), integer(1))
  arms[order(sizes, counts, decreasing = TRUE)]
}

.normalize_n_glycan_graph <- function(graph) {
  arms <- .normalization_n_arms(graph)
  if (is.null(arms)) {
    return(graph)
  }

  # Ordered tips are A, B, C (three antennae) or A, B, C, D (four).
  # Index 3 is C; indices 1 and 3 are A and C, one on each Man side.
  tips <- unlist(arms, use.names = FALSE)
  # Edge columns 1 and 2 hold parent and child vertex IDs, respectively.
  edges <- igraph::as_edgelist(graph, names = FALSE)
  original_edges <- edges
  mono <- igraph::vertex_attr(graph, "mono")
  gal <- which(mono == "Gal")
  sia <- which(mono == "Neu5Ac")
  # Rule counts: n = antennae, g = Gal residues, s = Neu5Ac residues.
  n <- length(tips)
  g <- length(gal)
  s <- length(sia)

  gal_targets <- NULL
  # Three antennae with one Gal: place it on the single-antenna side (C).
  if (n == 3L && g == 1L) {
    gal_targets <- tips[3L]
  }
  # Either antenna count with two Gal: place one on each Man side (A, C).
  if ((n == 3L && g == 2L) || (n == 4L && g == 2L)) {
    gal_targets <- tips[c(1L, 3L)]
  }
  if (!is.null(gal_targets)) {
    edges <- .normalization_reparent(edges, gal, gal_targets)
  }

  sia_targets <- NULL
  # One Neu5Ac goes to C with 2 or 3 Gal on three antennae, or 3 Gal on four.
  if (
    (n == 3L && g %in% c(2L, 3L) && s == 1L) ||
      (n == 4L && g == 3L && s == 1L)
  ) {
    sia_targets <- tips[3L]
  }
  # Two Neu5Ac go to A and C with 3 Gal on three antennae, or 3 or 4 on four.
  if (
    (n == 3L && g == 3L && s == 2L) ||
      (n == 4L && g %in% c(3L, 4L) && s == 2L)
  ) {
    sia_targets <- tips[c(1L, 3L)]
  }
  if (!is.null(sia_targets)) {
    gal_parents <- edges[match(gal, edges[, 2L]), 1L]
    edges <- .normalization_reparent(
      edges,
      sia,
      gal[match(sia_targets, gal_parents)]
    )
  }
  if (identical(edges, original_edges)) {
    return(graph)
  }

  edge_attrs <- igraph::edge_attr(graph)
  graph <- igraph::delete_edges(graph, igraph::E(graph))
  igraph::add_edges(graph, as.vector(t(edges)), attr = edge_attrs)
}

# Keep already occupied target sites, moving only the remaining residues.
.normalization_reparent <- function(edges, residues, targets) {
  rows <- match(residues, edges[, 2L])
  parents <- edges[rows, 1L]
  move <- !parents %in% targets
  edges[rows[move], 1L] <- setdiff(targets, parents)
  edges
}
