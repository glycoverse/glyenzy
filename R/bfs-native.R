# Public igraph accessors are used only at the R/native boundary. Objects with
# additional attributes or unresolved metadata retain their R graph operations;
# the native engine still owns their search and bookkeeping.
.bfs_native_graph <- function(graph) {
  ga <- igraph::graph_attr(graph)
  va <- igraph::vertex_attr(graph)
  ea <- igraph::edge_attr(graph)
  attributes_present <- any(vapply(
    list(ga$anomer, ga$alditol, va$mono, va$sub, ea$linkage),
    function(x) !is.null(attributes(x)),
    logical(1)
  ))
  if (
    attributes_present ||
      length(setdiff(names(ga), c("anomer", "alditol"))) ||
      length(setdiff(names(va), c("name", "mono", "sub"))) ||
      length(setdiff(names(ea), "linkage")) ||
      anyNA(c(va$mono, va$sub, ea$linkage, ga$anomer)) ||
      any(!grepl("^[ab?][12?]-([1-9](/[1-9])*|[?])$", ea$linkage))
  ) {
    return(NULL)
  }
  list(
    n = as.integer(igraph::vcount(graph)),
    edges = matrix(
      as.integer(igraph::as_edgelist(graph, names = FALSE)),
      ncol = 2L
    ),
    attributes = ga,
    vertices = va,
    edge_attributes = ea
  )
}

.bfs_native_restore <- function(record) {
  graph <- igraph::make_empty_graph(record$n, directed = TRUE)
  graph <- igraph::add_edges(graph, as.vector(t(record$edges)))
  igraph::graph_attr(graph) <- record$attributes
  igraph::vertex_attr(graph) <- c(
    list(name = as.character(seq_len(record$n))),
    record$vertices[c("mono", "sub")]
  )
  igraph::edge_attr(graph) <- record$edge_attributes
  graph
}

.bfs_native_pack <- function(graph, key, structure = NULL) {
  list(
    record = .bfs_native_graph(graph),
    graph = graph,
    key = key,
    structure = structure
  )
}

.bfs_native_run <- function(self, private) {
  plan <- private$rule_plan
  rules <- lapply(plan$rules, function(job) {
    prepared <- job$prepared_rule
    list(
      acceptor = .bfs_native_graph(prepared$acceptor),
      rejects = lapply(prepared$rejects, .bfs_native_graph),
      requires = lapply(prepared$requires, function(x) {
        list(motif = .bfs_native_graph(x$motif), alignment = x$alignment)
      }),
      alignment = job$rule$acceptor_alignment,
      site = job$rule$acceptor_idx,
      action = if (inherits(job$enzyme, "glyenzy_st_enzyme")) {
        "ST"
      } else if (inherits(job$enzyme, "glyenzy_gh_enzyme")) {
        "GH"
      } else {
        "GT"
      },
      mono = job$rule$new_residue,
      linkage = job$rule$new_linkage,
      sulfate = job$rule$new_substituent,
      product_sub = job$rule$product_sub
    )
  })
  enzymes <- lapply(seq_along(self$enzymes), function(i) {
    e <- self$enzymes[[i]]
    list(
      name = e$name,
      standard = .uses_standard_graph_action(e),
      inert = .can_batch_bfs_enzyme(e) && !.uses_standard_graph_action(e),
      rules = plan$enzyme_rule_ids[[i]],
      types = e$glycan_type
    )
  })
  get_graph <- function(node) {
    if (!is.null(node$graph)) node$graph else .bfs_native_restore(node$record)
  }
  get_structure <- function(node) {
    if (!is.null(node$structure)) {
      return(node$structure)
    }
    graph <- get_graph(node)
    glyrepr::new_glycan_structure(
      node$key,
      stats::setNames(list(graph), node$key)
    )
  }
  # This callback is for user S3 actions and attribute/floating graph operations,
  # never an alternative R BFS. Preserve each complete enzyme-cell call.
  expand <- function(node, enzyme_idx) {
    e <- self$enzymes[[enzyme_idx]]
    graph <- get_graph(node)
    structure <- get_structure(node)
    if (.uses_standard_graph_action(e)) {
      raw <- .apply_enzyme_prepared_graphs(
        graph,
        e,
        plan$prepared_rules[[enzyme_idx]],
        self$structure_level,
        .glymotif_mode(structure)
      )
      p <- private$prepare_graph_products(raw, .glymotif_mode(structure))
    } else {
      products <- .apply_enzyme(
        structure,
        e,
        structure_level = self$structure_level
      )[[1]]
      p <- private$prepare_products(products)
    }
    list(
      nodes = lapply(seq_along(p$keys), function(j) {
        .bfs_native_pack(
          p$graphs[[j]],
          p$keys[[j]],
          if (is.null(p$products)) NULL else p$products[j]
        )
      }),
      products = p$products
    )
  }
  filter <- function(nodes, products) {
    if (is.null(products)) {
      keys <- vapply(nodes, `[[`, character(1), "key")
      graphs <- lapply(nodes, get_graph)
      products <- glyrepr::new_glycan_structure(
        keys,
        stats::setNames(graphs, keys)
      )
    }
    keep <- self$filter(products)
    checkmate::assert_logical(keep, len = length(products), any.missing = FALSE)
    keep
  }
  prune <- function(node, mode) {
    private$is_promising_product(.move_glycan_root_last(get_graph(node)), mode)
  }
  match_target <- function(node, target_idx) {
    glymotif::.g_have_motif(
      get_graph(node),
      private$target_graphs[[target_idx]],
      alignment = "whole",
      mode = "lenient"
    )
  }
  frontier_graphs <- self$queue_graphs
  if (length(frontier_graphs) != length(self$queue_keys)) {
    frontier_graphs <- lapply(self$queue, glyrepr::get_structure_graphs)
  }
  frontier <- lapply(seq_along(self$queue_keys), function(i) {
    .bfs_native_pack(
      frontier_graphs[[i]],
      self$queue_keys[[i]],
      if (length(self$queue) >= i) self$queue[[i]] else NULL
    )
  })
  monos <- glyrepr::available_monosaccharides()
  type_motifs <- c(
    "GlcNAc(b1-4)GlcNAc(?1-",
    "Xyl(??-?)Glc(?1-",
    "Glc(a1-2)Gal(?1-",
    "Man(a1-3)Man(?1-",
    "Man(a1-6)Man(?1-"
  )
  config <- list(
    source = .bfs_native_pack(private$source_graph, self$from_key, self$from_g),
    frontier = frontier,
    targets = lapply(private$target_graphs, .bfs_native_graph),
    target_keys = self$to_keys,
    remaining = self$remaining_targets_map$keys(),
    rules = rules,
    enzymes = enzymes,
    max_steps = self$max_steps,
    step = self$step,
    filter = !is.null(self$filter),
    scalar = !is.null(self$filter) ||
      !all(vapply(self$enzymes, .can_batch_bfs_enzyme, logical(1))),
    topological = self$structure_level == "topological",
    whole = self$target_match == "whole",
    product_lenient = private$product_match_mode == "lenient",
    target_lenient = private$target_match_mode == "lenient",
    ncore = .bfs_native_graph(private$n_core_graph),
    pre = .bfs_native_graph(private$pre_mgat2_graph),
    type_motifs = lapply(
      glyrepr::get_structure_graphs(
        glyrepr::as_glycan_structure(type_motifs),
        return_list = TRUE
      ),
      .bfs_native_graph
    ),
    dictionary = data.frame(
      concrete = monos,
      generic = glyrepr::convert_to_generic(monos)
    ),
    byte_order = identical(Sys.getlocale("LC_COLLATE"), "C"),
    visited = self$visited,
    parent = self$parent,
    parent_enzyme = self$parent_enzyme,
    parent_step = self$parent_step
  )
  out <- cpp_bfs_search(
    config,
    list(
      expand = expand,
      filter = filter,
      prune = prune,
      match_target = match_target,
      order = base::order,
      mode = function(node) .glymotif_mode(get_structure(node)) == "lenient"
    )
  )
  self$step <- out$step
  self$all_edges <- c(self$all_edges, out$all_edges)
  self$native_stats <- out$stats
  self$queue_keys <- vapply(out$frontier, `[[`, character(1), "key")
  self$queue_graphs <- lapply(out$frontier, get_graph)
  self$queue <- if (config$scalar) {
    lapply(out$frontier, get_structure)
  } else {
    list()
  }
  self$found_keys_storage <- c(
    self$found_keys_storage[seq_len(self$found_tail)],
    out$found_keys
  )
  self$found_tail <- length(self$found_keys_storage)
  for (key in setdiff(
    self$remaining_targets_map$keys(),
    out$missing_target_keys
  )) {
    self$remaining_targets_map$remove(key)
  }
  list(
    found_keys = self$found_keys_storage,
    all_edges = self$all_edges,
    parent = self$parent,
    parent_enzyme = self$parent_enzyme,
    parent_step = self$parent_step,
    missing_target_keys = out$missing_target_keys
  )
}
