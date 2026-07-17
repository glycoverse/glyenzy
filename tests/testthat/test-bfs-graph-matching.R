test_that("prepared graph matching preserves site-specific rejects", {
  glycan <- glyrepr::as_glycan_structure(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  rule <- new_enzyme_rule(
    acceptor = glyrepr::as_glycan_structure("Gal(b1-4)GlcNAc(b1-"),
    product = glyrepr::as_glycan_structure(
      "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
    ),
    acceptor_alignment = "terminal",
    rejects = glyrepr::as_glycan_structure(
      "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-"
    )
  )
  graph <- glyrepr::get_structure_graphs(glycan)

  expected <- .match_rule(glycan, rule)[[1]]
  actual <- .match_rule_graph(
    graph,
    rule,
    .prepare_rule_graphs(rule),
    mode = "strict"
  )

  expect_equal(actual, expected)
  expect_length(actual, 1L)
})

test_that("prepared graph matching preserves lenient reduced-level matches", {
  rule <- enzyme("C1GALT1")$rules[[1]]
  intact <- glyrepr::as_glycan_structure("GalNAc(a1-")
  glycans <- list(
    glyrepr::remove_linkages(intact),
    glyrepr::reduce_structure_level(intact, "basic")
  )

  for (glycan in glycans) {
    expected <- .match_rule(glycan, rule)[[1]]
    actual <- .match_rule_graph(
      glyrepr::get_structure_graphs(glycan),
      rule,
      .prepare_rule_graphs(rule),
      mode = "lenient"
    )

    expect_equal(actual, expected)
  }
})

test_that("prepared graph matching preserves duplicate-reject errors", {
  enzyme <- make_enzyme(
    name = "TEST_DUPLICATE_REJECTS",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "GalNAc(a1-",
      acceptor_alignment = "core",
      rejects = c(
        "Fuc(a1-2)GalNAc(a1-",
        "Fuc(a1-2)GalNAc(a1-"
      ),
      product = "GlcNAc(b1-3)GalNAc(a1-"
    ))
  )

  expect_snapshot(
    trace_biosynthesis(
      "GlcNAc(b1-3)GalNAc(a1-",
      enzymes = list(enzyme),
      max_steps = 1
    ),
    error = TRUE
  )
})

test_that("prepared graph pruning matches glycan-level pruning", {
  products <- glyrepr::as_glycan_structure(c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-3)GalNAc(a1-",
    "GlcNAc(b1-3)GalNAc(a1-"
  ))
  targets <- glyrepr::as_glycan_structure(c(
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",
    "GlcNAc(b1-3)GalNAc(a1-"
  ))
  product_graphs <- glyrepr::get_structure_graphs(products, return_list = TRUE)
  target_graphs <- glyrepr::get_structure_graphs(targets, return_list = TRUE)
  n_core_graph <- glyrepr::get_structure_graphs(
    glyrepr::as_glycan_structure(
      "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    )
  )
  pre_mgat2_graph <- glyrepr::get_structure_graphs(
    glyrepr::as_glycan_structure(
      "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    )
  )

  expected <- is_promising_intermediate(products, targets)
  expect_equal(expected, c(TRUE, FALSE, TRUE, TRUE))
  actual <- purrr::map_lgl(
    product_graphs,
    .is_promising_intermediate_graph,
    target_graphs = target_graphs,
    product_mode = "strict",
    target_mode = "strict",
    n_core_graph = n_core_graph,
    pre_mgat2_graph = pre_mgat2_graph
  )

  expect_equal(actual, expected)
})

test_that("BFS rule plans share only equivalent standard rules", {
  enzymes <- list(
    enzyme("MAN2A1"),
    enzyme("MAN2A2"),
    enzyme("B4GALT2"),
    enzyme("B4GALT3")
  )
  plan <- .prepare_bfs_rule_plan(enzymes)

  expect_length(plan$rules, 2L)
  expect_true(all(vapply(enzymes, .can_batch_bfs_enzyme, logical(1))))
  expect_equal(plan$enzyme_rule_ids[[1]], plan$enzyme_rule_ids[[2]])
  expect_equal(plan$enzyme_rule_ids[[3]], plan$enzyme_rule_ids[[4]])
  expect_equal(plan$enzyme_plan_ids, c(1L, 1L, 2L, 2L))

  custom <- enzyme("B4GALT3")
  class(custom) <- c("custom_gt_enzyme", class(custom))
  private_plan <- .prepare_bfs_rule_plan(list(enzyme("B4GALT2"), custom))

  expect_length(private_plan$rules, 2L)
  expect_false(.can_batch_bfs_enzyme(custom))
  expect_equal(private_plan$enzyme_plan_ids, c(1L, 2L))

  custom_starter <- enzyme("DPAGT1")
  class(custom_starter) <- c("custom_starter", class(custom_starter))
  expect_false(.can_batch_bfs_enzyme(custom_starter))
})

test_that("batched rule jobs preserve scalar products for every cell", {
  frontier <- glyrepr::as_glycan_structure(c(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  enzymes <- list(
    enzyme("MAN1A1"),
    enzyme("MAN1A2"),
    enzyme("MAN2A1"),
    enzyme("B4GALT1"),
    enzyme("ST3GAL4")
  )
  plan <- .prepare_bfs_rule_plan(enzymes)
  graphs <- glyrepr::get_structure_graphs(frontier, return_list = TRUE)
  modes <- vapply(
    seq_along(frontier),
    function(i) .glymotif_mode(frontier[i]),
    character(1)
  )
  batched <- .apply_bfs_rule_frontier(
    graphs,
    modes,
    plan,
    structure_level = "intact"
  )

  for (frontier_idx in seq_along(frontier)) {
    for (enzyme_idx in seq_along(enzymes)) {
      expected <- .apply_enzyme_prepared(
        graphs[[frontier_idx]],
        enzymes[[enzyme_idx]],
        plan$prepared_rules[[enzyme_idx]],
        structure_level = "intact",
        mode = modes[[frontier_idx]]
      )
      actual <- .collect_bfs_rule_products(
        batched,
        plan$enzyme_rule_ids[[enzyme_idx]],
        frontier_idx
      )

      expect_equal(unname(as.character(actual)), unname(as.character(expected)))
    }
  }
})

test_that("batched rule jobs preserve reduced-level products", {
  intact <- glyrepr::as_glycan_structure("GalNAc(a1-")
  glycans <- list(
    topological = glyrepr::remove_linkages(intact),
    basic = glyrepr::reduce_structure_level(intact, "basic")
  )
  enzymes <- list(enzyme("C1GALT1"))
  plan <- .prepare_bfs_rule_plan(enzymes)

  for (structure_level in names(glycans)) {
    glycan <- glycans[[structure_level]]
    graph <- glyrepr::get_structure_graphs(glycan)
    batched <- .apply_bfs_rule_frontier(
      list(graph),
      "lenient",
      plan,
      structure_level = structure_level
    )
    actual <- .collect_bfs_rule_products(
      batched,
      plan$enzyme_rule_ids[[1]],
      1L
    )
    expected <- .apply_enzyme_prepared(
      graph,
      enzymes[[1]],
      plan$prepared_rules[[1]],
      structure_level = structure_level,
      mode = "lenient"
    )

    expect_equal(unname(as.character(actual)), unname(as.character(expected)))
  }
})

test_that("shared rule jobs retain distinct enzyme edges", {
  make_core1_enzyme <- function(name) {
    make_enzyme(
      name = name,
      type = "GT",
      species = "human",
      rules = list(list(
        acceptor = "GalNAc(a1-",
        acceptor_alignment = "core",
        rejects = NULL,
        product = "Gal(b1-3)GalNAc(a1-"
      ))
    )
  }
  path <- trace_biosynthesis(
    "Gal(b1-3)GalNAc(a1-",
    enzymes = list(make_core1_enzyme("E1"), make_core1_enzyme("E2")),
    max_steps = 1
  )
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(edges$enzyme, c("E1", "E2"))
  expect_equal(edges$to, rep("Gal(b1-3)GalNAc(a1-", 2L))
})

test_that("batched expansion rejects occupied acceptor positions", {
  make_branch_enzyme <- function(name, product) {
    make_enzyme(
      name = name,
      type = "GT",
      species = "human",
      rules = list(list(
        acceptor = "GalNAc(a1-",
        acceptor_alignment = "core",
        rejects = NULL,
        product = product
      ))
    )
  }
  target <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  path <- trace_biosynthesis(
    target,
    enzymes = list(
      make_branch_enzyme("C3", "Gal(b1-3)GalNAc(a1-"),
      make_branch_enzyme("C6", "GlcNAc(b1-6)GalNAc(a1-")
    ),
    max_steps = 2
  )
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(edges$to[edges$step == 2L], rep(target, 2L))
  expect_equal(edges$enzyme[edges$step == 2L], c("C6", "C3"))
})

test_that("batched expansion processes bounded chunks in frontier order", {
  from <- glyrepr::as_glycan_structure("GalNAc(a1-")
  target <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  engine <- BfsSynthesisSearch$new(
    from_g = from,
    to_gs = target,
    enzymes = list(enzyme("C1GALT1")),
    max_steps = 1L
  )
  frontier_size <- 33L
  frontier_keys <- paste0("source-", seq_len(frontier_size))
  engine$queue <- rep(list(from), frontier_size)
  engine$queue_keys <- frontier_keys

  result <- engine$run()
  edge_sources <- vapply(
    result$all_edges,
    function(edge) edge$from,
    character(1)
  )

  expect_equal(edge_sources, frontier_keys)
})

test_that("custom enzyme actions retain scalar frontier order", {
  action_order <- character()
  action_method <- function(
    enzyme,
    graph,
    indices_to_act_on,
    rule,
    structure_level = "intact"
  ) {
    action_order <<- c(action_order, enzyme$name)
    .apply_rule_action.glyenzy_gt_enzyme(
      enzyme,
      graph,
      indices_to_act_on,
      rule,
      structure_level
    )
  }
  rlang::local_bindings(
    .apply_rule_action.test_stateful_gt_enzyme = action_method,
    .env = globalenv()
  )
  make_custom_enzyme <- function(name) {
    custom <- make_enzyme(
      name = name,
      type = "GT",
      species = "human",
      rules = list(list(
        acceptor = "GalNAc(a1-",
        acceptor_alignment = "core",
        rejects = NULL,
        product = "Gal(b1-3)GalNAc(a1-"
      ))
    )
    class(custom) <- c("test_stateful_gt_enzyme", class(custom))
    custom
  }
  from <- glyrepr::as_glycan_structure("GalNAc(a1-")
  target <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  engine <- BfsSynthesisSearch$new(
    from_g = from,
    to_gs = target,
    enzymes = list(make_custom_enzyme("E1"), make_custom_enzyme("E2")),
    max_steps = 1L
  )
  engine$queue <- rep(list(from), 2L)
  engine$queue_keys <- c("source-1", "source-2")

  result <- engine$run()

  expect_equal(action_order, c("E1", "E2", "E1", "E2"))
  expect_length(result$all_edges, 4L)
})

test_that("custom enzyme-level actions retain S3 dispatch", {
  calls <- 0L
  action_method <- function(glycans, enzyme, structure_level = "intact") {
    calls <<- calls + 1L
    rep(
      list(glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")),
      length(glycans)
    )
  }
  rlang::local_bindings(
    .apply_enzyme.test_custom_enzyme = action_method,
    .env = globalenv()
  )
  custom <- make_enzyme(
    name = "CUSTOM",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "GalNAc(a1-",
      acceptor_alignment = "core",
      rejects = NULL,
      product = "Gal(b1-3)GalNAc(a1-"
    ))
  )
  class(custom) <- c("test_custom_enzyme", class(custom))

  path <- trace_biosynthesis(
    "Gal(b1-3)GalNAc(a1-",
    enzymes = list(custom),
    max_steps = 1
  )
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(calls, 1L)
  expect_equal(edges$to, "Gal(b1-3)GalNAc(a1-")
})

test_that("BFS batching does not group enzymes by name", {
  make_named_enzyme <- function(product) {
    make_enzyme(
      name = "DUP",
      type = "GT",
      species = "human",
      rules = list(list(
        acceptor = "GalNAc(a1-",
        acceptor_alignment = "core",
        rejects = NULL,
        product = product
      ))
    )
  }
  targets <- c(
    "Gal(b1-3)GalNAc(a1-",
    "GlcNAc(b1-3)GalNAc(a1-"
  )
  path <- trace_biosynthesis(
    targets,
    enzymes = list(
      make_named_enzyme(targets[[1]]),
      make_named_enzyme(targets[[2]])
    ),
    max_steps = 1
  )
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(edges$enzyme, c("DUP", "DUP"))
  expect_equal(edges$to, targets)
})

test_that("stateful filters retain scalar enzyme calls", {
  make_core1_enzyme <- function(name) {
    make_enzyme(
      name = name,
      type = "GT",
      species = "human",
      rules = list(list(
        acceptor = "GalNAc(a1-",
        acceptor_alignment = "core",
        rejects = NULL,
        product = "Gal(b1-3)GalNAc(a1-"
      ))
    )
  }
  calls <- list()
  filter <- function(products) {
    calls[[length(calls) + 1L]] <<- as.character(products)
    rep(length(calls) == 2L, length(products))
  }
  path <- trace_biosynthesis(
    "Gal(b1-3)GalNAc(a1-",
    enzymes = list(make_core1_enzyme("E1"), make_core1_enzyme("E2")),
    max_steps = 1,
    filter = filter
  )
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(calls, rep(list("Gal(b1-3)GalNAc(a1-"), 2L))
  expect_equal(edges$enzyme, "E2")
})
