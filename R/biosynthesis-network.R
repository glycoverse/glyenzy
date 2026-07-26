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
