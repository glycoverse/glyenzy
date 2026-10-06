library(glyenzy)
ns <- asNamespace("glyenzy")
internal <- function(name) get(name, ns)
base_dir <- "data-raw/bfs-native"
Rcpp::sourceCpp(file.path(base_dir, "native.cpp"))
compact <- function(g) {
  stopifnot(
    all(igraph::V(g)$sub == ""),
    !isTRUE(igraph::graph_attr(g, "alditol"))
  )
  stopifnot(is.null(igraph::graph_attr(g, "floating_parts")))
  list(
    n = as.integer(igraph::vcount(g)),
    edges = matrix(as.integer(igraph::as_edgelist(g, names = FALSE)), ncol = 2),
    attributes = igraph::graph_attr(g),
    vertices = igraph::vertex_attr(g),
    edge_attributes = igraph::edge_attr(g)
  )
}
graph <- function(x) glyrepr::get_structure_graphs(x)
prepare <- function(from, to, enzymes, steps) {
  from <- glyrepr::as_glycan_structure(from)
  to <- glyrepr::as_glycan_structure(to)
  stopifnot(
    internal(".glymotif_mode")(from) == "strict",
    internal(".glymotif_mode")(to) == "strict"
  )
  ez <- lapply(enzymes, enzyme)
  rules <- lapply(ez, function(e) {
    stopifnot(
      identical(class(e), c("glyenzy_gt_enzyme", "glyenzy_enzyme")) ||
        identical(class(e), c("glyenzy_gh_enzyme", "glyenzy_enzyme"))
    )
    stopifnot(
      is.null(e$glycan_type) ||
        internal(".glycan_type_is_compatible")(
          internal(".glycan_type_graph")(graph(from)),
          e$glycan_type
        )
    )
    lapply(e$rules, function(r) {
      list(
        acceptor = compact(graph(r$acceptor)),
        rejects = lapply(
          glyrepr::get_structure_graphs(r$rejects, return_list = TRUE),
          compact
        ),
        requires = lapply(r$requires, function(q) {
          list(motif = compact(graph(q$motif)), alignment = q$alignment)
        }),
        alignment = r$acceptor_alignment,
        site = r$acceptor_idx,
        gt = inherits(e, "glyenzy_gt_enzyme"),
        mono = if (is.null(r$new_residue)) "" else r$new_residue,
        link = if (is.null(r$new_linkage)) "" else r$new_linkage
      )
    })
  })
  monos <- glyrepr::available_monosaccharides()
  list(
    from = from,
    to = to,
    enzymes = ez,
    names = enzymes,
    steps = steps,
    native = list(
      source = compact(graph(from)),
      targets = lapply(
        glyrepr::get_structure_graphs(to, return_list = TRUE),
        compact
      ),
      enzyme_rules = rules,
      dictionary = data.frame(
        concrete = monos,
        generic = glyrepr::convert_to_generic(monos)
      ),
      ncore = compact(graph(internal(".n_glycan_starting_glycan")("virtual"))),
      pre = compact(graph(glyrepr::as_glycan_structure(
        "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
      ))),
      max_steps = steps
    )
  )
}
decode <- function(z, w) {
  keys <- vapply(
    z$graphs,
    function(x) {
      g <- igraph::make_empty_graph(x$n, directed = TRUE)
      g <- igraph::add_edges(g, as.vector(t(x$edges)))
      for (a in names(x$attributes)) {
        g <- igraph::set_graph_attr(g, a, x$attributes[[a]])
      }
      for (a in names(x$vertices)) {
        g <- igraph::set_vertex_attr(g, a, value = x$vertices[[a]])
      }
      g <- igraph::set_vertex_attr(
        g,
        "name",
        value = as.character(seq_len(x$n))
      )
      g <- igraph::set_edge_attr(
        g,
        "linkage",
        value = x$edge_attributes$linkage
      )
      glyrepr::graph_to_iupac(glyrepr::canonicalize_glycan_graph(g))
    },
    character(1)
  )
  list(
    edges = data.frame(
      from = keys[z$from + 1L],
      to = keys[z$to + 1L],
      enzyme = w$names[z$enzyme + 1L],
      step = z$step
    ),
    found = keys[z$found + 1L],
    parent = setNames(keys[z$parent[-1] + 1L], keys[-1]),
    parent_enzyme = setNames(w$names[z$parent_enzyme[-1] + 1L], keys[-1]),
    parent_step = setNames(z$parent_step[-1], keys[-1]),
    missing = z$missing
  )
}
reference <- function(w) {
  internal("bfs_synthesis_search")(
    w$from,
    w$to,
    w$enzymes,
    w$steps,
    allow_partial = TRUE
  )
}
normalize_ref <- function(z) {
  named <- function(e) {
    x <- as.list(e)
    unlist(x, use.names = TRUE)
  }
  list(
    edges = do.call(rbind, lapply(z$all_edges, as.data.frame)),
    found = z$found_keys,
    parent = named(z$parent),
    parent_enzyme = named(z$parent_enzyme),
    parent_step = named(z$parent_step),
    missing = length(z$missing_target_keys)
  )
}
parity <- function(a, b) {
  sorted <- function(x) x[order(names(x))]
  edge_sort <- function(x) {
    x <- x[do.call(order, x), , drop = FALSE]
    rownames(x) <- NULL
    x
  }
  c(
    edges = identical(edge_sort(a$edges), edge_sort(b$edges)),
    edge_order = identical(a$edges, b$edges),
    found = identical(sort(a$found), sort(b$found)),
    found_order = identical(a$found, b$found),
    parent = identical(sorted(a$parent), sorted(b$parent)),
    parent_enzyme = identical(sorted(a$parent_enzyme), sorted(b$parent_enzyme)),
    parent_step = identical(sorted(a$parent_step), sorted(b$parent_step)),
    missing = identical(as.integer(a$missing), as.integer(b$missing))
  )
}
