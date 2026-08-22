# Determine broad glycan classes from reducing-end structure

.glycan_types <- function(glycans) {
  checkmate::assert_class(glycans, "glyrepr_structure")

  if (length(glycans) == 0L) {
    return(character())
  }

  graphs <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)
  purrr::map2_chr(graphs, as.character(glycans), .glycan_type_graph)
}

.glycan_type_graph <- function(graph, structure = NULL) {
  if (is.null(graph) || igraph::vcount(graph) == 0L) {
    return(NA_character_)
  }

  mono <- igraph::vertex_attr(graph, "mono")
  if (length(mono) == 0L || anyNA(mono)) {
    return(NA_character_)
  }

  if (.graph_has_glycan_type_core(graph, "GlcNAc(b1-4)GlcNAc(?1-")) {
    return("N")
  }

  if (
    .graph_has_glycan_type_core(graph, "Xyl(??-?)Glc(?1-") ||
      .graph_has_glycan_type_core(graph, "Glc(a1-2)Gal(?1-")
  ) {
    return("O")
  }

  roots <- which(igraph::degree(graph, mode = "in") == 0L)
  if (length(roots) != 1L) {
    return(NA_character_)
  }

  root_mono <- mono[[roots[[1]]]]
  if (
    identical(root_mono, "Man") &&
      ((!is.null(structure) && grepl("Man\\(b[0-9?]+-$", structure)) ||
        (is.null(structure) &&
          (.graph_has_glycan_type_core(graph, "Man(a1-3)Man(?1-") ||
            .graph_has_glycan_type_core(graph, "Man(a1-6)Man(?1-"))))
  ) {
    return(NA_character_)
  }
  if (root_mono %in% c("GalNAc", "Man", "Fuc", "GlcNAc", "Xyl")) {
    return("O")
  }
  if (root_mono %in% c("Glc", "Gal")) {
    return("lipid/free")
  }

  NA_character_
}

.graph_has_glycan_type_core <- function(graph, motif) {
  motif_graph <- glyrepr::get_structure_graphs(
    glyparse::auto_parse(motif),
    return_list = FALSE
  )
  length(.g_match_motif_substituent_subset(
    graph,
    motif_graph,
    alignment = "core",
    mode = "lenient"
  )) >
    0L
}

.enzyme_glycan_type_mask <- function(glycans, enzyme) {
  if (is.null(enzyme$glycan_type)) {
    return(rep(TRUE, length(glycans)))
  }

  types <- .glycan_types(glycans)
  vapply(
    types,
    .glycan_type_is_compatible,
    logical(1),
    supported_types = enzyme$glycan_type
  )
}

.enzyme_supports_glycan_graph <- function(graph, enzyme) {
  if (is.null(enzyme$glycan_type)) {
    return(TRUE)
  }
  type <- .glycan_type_graph(graph)
  .glycan_type_is_compatible(type, enzyme$glycan_type)
}

.glycan_type_is_compatible <- function(type, supported_types) {
  if (is.na(type)) {
    return(TRUE)
  }
  if (identical(type, "lipid/free")) {
    return(any(c("lipid", "free") %in% supported_types))
  }
  type %in% supported_types
}

.validate_glycan_type <- function(glycan_type) {
  if (is.null(glycan_type)) {
    return(NULL)
  }

  checkmate::assert_character(glycan_type, min.len = 1L, any.missing = FALSE)
  checkmate::assert_subset(glycan_type, c("N", "O", "lipid", "free"))
  unique(intersect(c("N", "O", "lipid", "free"), glycan_type))
}
