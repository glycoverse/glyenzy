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
