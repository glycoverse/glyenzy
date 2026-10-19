test_that("shared rules preserve selected enzyme order and duplicate labels", {
  enzymes <- lapply(
    c("B4GALT3", "MAN2A2", "B4GALT2", "MAN2A1", "B4GALT3"),
    enzyme
  )
  actual <- .bfs_search_rule_plan(enzymes)
  expected <- .prepare_bfs_rule_plan(enzymes, prepare_strict = FALSE)
  expect_identical(actual$enzyme_rule_ids, expected$enzyme_rule_ids)
  expect_null(actual$prepared_rules)
  expect_equal(actual$native_rules, lapply(expected$rules, .bfs_native_rule))
})

test_that("modified rules and custom dispatch bypass shared rule preparation", {
  builtin <- enzyme("C1GALT1")
  .bfs_search_rule_plan(list(builtin))
  edited <- builtin
  edited$rules <- enzyme("B3GNT6")$rules
  custom <- builtin
  class(custom) <- c("custom_gt_enzyme", class(custom))
  for (enzyme in list(edited, custom)) {
    expect_equal(
      .bfs_search_rule_plan(list(enzyme)),
      .prepare_bfs_rule_plan(list(enzyme), prepare_strict = FALSE)
    )
  }
})

test_that("bundled plans stay synchronized with the enzyme database", {
  expect_identical(
    .build_bfs_prefilter_plan(glyenzy_enzymes),
    glyenzy_bfs_prefilter
  )
  expect_identical(
    .build_bfs_builtin_rule_plan(glyenzy_enzymes),
    glyenzy_bfs_rules
  )
  actual <- .build_bfs_structures()
  expected <- glyenzy_bfs_structures
  expect_identical(actual$type_motifs, expected$type_motifs)
  for (name in setdiff(names(actual), "type_motifs")) {
    expect_identical(
      as.character(actual[[name]]),
      as.character(expected[[name]])
    )
    expect_identical(
      .bfs_native_graph(glyrepr::get_structure_graphs(actual[[name]])),
      .bfs_native_graph(glyrepr::get_structure_graphs(expected[[name]]))
    )
  }
})

test_that("default tracing does not load or prepare full enzyme graphs", {
  local_mocked_bindings(
    glyenzy_enzymes = NULL,
    .prepare_rule_graphs = function(...) stop("Unexpected runtime preparation"),
    .build_bfs_prefilter_plan = function(...) stop("Unexpected plan build"),
    .build_bfs_builtin_rule_plan = function(...) stop("Unexpected plan build")
  )
  graph <- trace_biosynthesis("Gal(b1-3)GalNAc(a1-")
  expect_identical(igraph::edge_attr(graph, "enzymes"), list("C1GALT1"))
})

test_that("compact enzyme descriptors retain full-enzyme search results", {
  targets <- c(
    "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
  )
  enzymes <- .enzymes_from_arg(NULL)
  for (target in targets) {
    actual <- trace_biosynthesis(target)
    expected <- trace_biosynthesis(target, enzymes = enzymes)
    expect_identical(
      igraph::as_data_frame(actual, what = "both"),
      igraph::as_data_frame(expected, what = "both")
    )
  }
  target <- glyparse::auto_parse(targets[[2]])
  refs <- .bfs_default_enzymes(target)
  expected <- .can_enzymes_contribute(enzymes, target)
  local_mocked_bindings(cpp_bfs_have_motifs = function(...) stop("Unavailable"))
  expect_identical(.can_enzymes_contribute(refs, target), expected)
})

test_that("compact enzymes resolve correctly for graph-metadata callbacks", {
  from <- glyparse::auto_parse("GalNAc(a1-")
  graph <- igraph::set_graph_attr(
    glyrepr::get_structure_graphs(from),
    "custom",
    TRUE
  )
  from <- glyrepr::new_glycan_structure(
    as.character(from),
    stats::setNames(list(graph), as.character(from))
  )
  target <- glyparse::auto_parse("Gal(b1-3)GalNAc(a1-")
  refs <- .bfs_default_enzymes(target)["C1GALT1"]
  expected <- .perform_bfs_synthesis(from, target, list(enzyme("C1GALT1")), 1L)
  actual <- .perform_bfs_synthesis(from, target, unname(refs), 1L)
  expect_identical(
    igraph::as_data_frame(actual, what = "both"),
    igraph::as_data_frame(expected, what = "both")
  )
})
