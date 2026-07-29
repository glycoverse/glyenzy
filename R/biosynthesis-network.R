.new_biosynthesis_network <- function(graph, virtual = FALSE) {
  checkmate::assert_class(graph, "igraph")
  checkmate::assert_flag(virtual)

  network_classes <- c(
    if (virtual) "glyenzy_virtual_biosynthesis_network",
    "glyenzy_biosynthesis_network"
  )
  class(graph) <- unique(c(network_classes, class(graph)))
  graph
}

.set_biosynthesis_targets <- function(
  graph,
  targets,
  match = c("key", "whole")
) {
  match <- match.arg(match)
  vertex_names <- igraph::vertex_attr(graph, "name")

  if (identical(match, "key")) {
    is_target <- vertex_names %in% as.character(targets)
  } else {
    vertex_graphs <- glyrepr::get_structure_graphs(
      glyparse::auto_parse(vertex_names),
      return_list = TRUE
    )
    target_graphs <- glyrepr::get_structure_graphs(
      targets,
      return_list = TRUE
    )
    is_target <- purrr::map_lgl(
      vertex_graphs,
      function(vertex_graph) {
        purrr::some(
          target_graphs,
          function(target_graph) {
            glymotif::.g_have_motif(
              vertex_graph,
              target_graph,
              alignment = "whole",
              mode = "lenient"
            )
          }
        )
      }
    )
  }

  igraph::set_vertex_attr(graph, "target", value = is_target)
}
