# Compact records are safe for fully specified concrete trees. Reduced-level
# and metadata-bearing structures retain the graph-level matching semantics.
.virtual_native_records <- function(glycans) {
  if (
    !identical(.glycan_mono_type(glycans), "concrete") ||
      !identical(.glycan_structure_level(glycans), "intact")
  ) {
    return(NULL)
  }
  graphs <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)
  records <- lapply(graphs, .bfs_native_graph)
  if (any(vapply(records, is.null, logical(1)))) {
    return(NULL)
  }
  stats::setNames(records, as.character(glycans))
}

.virtual_native_synthesis <- function(from, to) {
  from_records <- .virtual_native_records(from)
  to_records <- .virtual_native_records(to)
  if (is.null(from_records) || is.null(to_records)) {
    return(NULL)
  }
  pack <- function(records) {
    lapply(seq_along(records), function(i) {
      list(
        record = records[[i]],
        graph = NULL,
        structure = NULL,
        key = names(records)[[i]]
      )
    })
  }
  cpp_virtual_trim(
    pack(from_records)[[1L]],
    pack(to_records),
    .bfs_native_dictionary(),
    base::order,
    identical(Sys.getlocale("LC_COLLATE"), "C")
  )
}

.virtual_enzymes_arg <- function(enzymes, glycans) {
  if (is.null(enzymes)) {
    builtin <- .bfs_default_enzymes(glycans)
    if (!is.null(builtin)) {
      return(unname(builtin))
    }
  }
  .process_enzymes_arg(enzymes, apply_prefilter = FALSE)
}

.amplify_virtual_edges <- function(path, enzymes, records = NULL) {
  if (
    igraph::ecount(path) > 0L &&
      length(enzymes) > 0L &&
      all(vapply(enzymes, .bfs_is_builtin_enzyme, logical(1)))
  ) {
    keys <- igraph::V(path)$name
    if (is.null(records)) {
      records <- .virtual_native_records(glyparse::auto_parse(keys))
    }
    if (!is.null(records)) {
      plan <- .bfs_search_rule_plan(enzymes)
      plan$enzymes <- enzymes
      annotations <- cpp_virtual_annotate(
        records[keys],
        matrix(as.integer(igraph::as_edgelist(path, names = FALSE)), ncol = 2L),
        plan,
        glyenzy_bfs_structures$type_motifs,
        .bfs_native_dictionary()
      )
      if (!is.null(annotations)) {
        return(igraph::set_edge_attr(
          path,
          "concrete_enzymes",
          value = annotations
        ))
      }
    }
  }
  .amplify_virtual_edges_r(path, lapply(enzymes, .bfs_resolve_enzyme))
}
