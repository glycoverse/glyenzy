test_that("virtual enzymes trace an intact glycan backward", {
  target <- "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"

  path <- trace_biosynthesis(target, method = "virtual")
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "b4GlcNAcT"))
  expect_equal(edges$step, 1:2)
  expect_equal(edges$from[[1]], "GalNAc(a1-")
  expect_equal(edges$to[[2]], target)
})

test_that("virtual enzymes include every order for independent branches", {
  target <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"

  path <- trace_biosynthesis(target, method = "virtual")
  vertices <- igraph::as_data_frame(path, what = "vertices")
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_setequal(
    vertices$name,
    c(
      "GalNAc(a1-",
      "Gal(b1-3)GalNAc(a1-",
      "GlcNAc(b1-6)GalNAc(a1-",
      target
    )
  )
  expect_equal(nrow(edges), 4L)
  expect_setequal(edges$enzyme, c("b3GalT", "b6GlcNAcT"))

  substrate_sizes <- purrr::map_int(
    glyparse::auto_parse(edges$from),
    ~ igraph::vcount(glyrepr::get_structure_graphs(.x))
  )
  product_sizes <- purrr::map_int(
    glyparse::auto_parse(edges$to),
    ~ igraph::vcount(glyrepr::get_structure_graphs(.x))
  )
  expect_equal(product_sizes - substrate_sizes, rep(1L, nrow(edges)))
  expect_equal(edges$step, product_sizes - 1L)
  expect_identical(igraph::is_dag(path), TRUE)
})

test_that("virtual tracing protects the N-glycan core", {
  target <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  path <- trace_biosynthesis(target, method = "virtual")
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(nrow(edges), 1L)
  expect_equal(edges$from, as.character(glyrepr::n_glycan_core()))
  expect_equal(edges$to, target)
  expect_equal(edges$enzyme, "b2GlcNAcT")
  expect_equal(edges$step, 1L)
  expect_equal(
    glymotif::have_motif(
      glyparse::auto_parse(igraph::V(path)$name),
      glyrepr::n_glycan_core(),
      alignment = "core"
    ),
    rep(TRUE, igraph::vcount(path))
  )
})

test_that("virtual enzyme names follow the target structure level", {
  topological <- glyrepr::reduce_structure_level(
    glyrepr::o_glycan_core_2(),
    "topological"
  )
  basic <- glyrepr::reduce_structure_level(
    glyrepr::o_glycan_core_2(),
    "basic"
  )
  partial <- glyparse::auto_parse(
    "Gal(b1-3)[GlcNAc(b1-?)]GalNAc(a1-"
  )

  topological_path <- suppressWarnings(
    trace_biosynthesis(topological, method = "virtual")
  )
  partial_path <- suppressWarnings(
    trace_biosynthesis(partial, method = "virtual")
  )
  basic_path <- suppressWarnings(
    trace_biosynthesis(basic, method = "virtual")
  )

  expect_setequal(
    igraph::E(topological_path)$enzyme,
    c("GalT", "GlcNAcT")
  )
  expect_setequal(
    igraph::E(partial_path)$enzyme,
    c("GalT", "GlcNAcT")
  )
  expect_setequal(
    igraph::E(basic_path)$enzyme,
    c("HexT", "HexNAcT")
  )
})

test_that("virtual tracing uses the root for non-N-glycans", {
  target <- "Gal(b1-4)Xyl(b1-"

  path <- trace_biosynthesis(target, method = "virtual")
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(edges$from, "Xyl(b1-")
  expect_equal(edges$to, target)
  expect_equal(edges$enzyme, "b4GalT")
})

test_that("virtual tracing combines multiple targets from one root", {
  targets <- c(
    "GlcNAc(b1-4)Gal(b1-",
    "Fuc(a1-2)Gal(b1-"
  )

  path <- trace_biosynthesis(targets, method = "virtual")
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(igraph::vcount(path), 3L)
  expect_equal(igraph::ecount(path), 2L)
  expect_setequal(edges$enzyme, c("b4GlcNAcT", "a2FucT"))
  expect_equal(unique(edges$from), "Gal(b1-")
  expect_equal(edges$step, rep(1L, 2L))
})

test_that("virtual N-glycan labels follow reduced structure levels", {
  intact <- glyparse::auto_parse(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  partial <- glyparse::auto_parse(
    "GlcNAc(b1-?)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  topological <- glyrepr::reduce_structure_level(intact, "topological")
  basic <- glyrepr::reduce_structure_level(intact, "basic")

  paths <- list(
    intact = trace_biosynthesis(intact, method = "virtual"),
    partial = suppressWarnings(
      trace_biosynthesis(partial, method = "virtual")
    ),
    topological = suppressWarnings(
      trace_biosynthesis(topological, method = "virtual")
    ),
    basic = suppressWarnings(
      trace_biosynthesis(basic, method = "virtual")
    )
  )

  expect_equal(igraph::E(paths$intact)$enzyme, "b2GlcNAcT")
  expect_equal(igraph::E(paths$partial)$enzyme, "GlcNAcT")
  expect_equal(igraph::E(paths$topological)$enzyme, "GlcNAcT")
  expect_equal(igraph::E(paths$basic)$enzyme, "HexNAcT")
})

test_that("virtual path tracing trims to the requested starting glycan", {
  from <- "Gal(b1-3)GalNAc(a1-"
  to <- "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"

  path <- path_biosynthesis(from, to, method = "virtual")
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(nrow(edges), 1L)
  expect_equal(edges$from, from)
  expect_equal(edges$to, to)
  expect_equal(edges$enzyme, "b4GlcNAcT")
  expect_equal(edges$step, 1L)
})

test_that("virtual paths map partial precursors to the requested start", {
  from <- "Gal(b1-3)GalNAc(a1-"
  to <- "Gal(b1-?)[GlcNAc(b1-6)]GalNAc(a1-"

  path <- suppressWarnings(path_biosynthesis(
    from,
    to,
    max_steps = 1,
    method = "virtual"
  ))
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(nrow(edges), 1L)
  expect_equal(edges$from, from)
  expect_equal(edges$to, as.character(glyparse::auto_parse(to)))
  expect_equal(edges$enzyme, "GlcNAcT")
})

test_that("virtual tracing shares a root across partial targets", {
  targets <- glyparse::auto_parse(c(
    "GalNAc(a1-",
    "Gal(b1-3)GalNAc(?1-"
  ))

  path <- suppressWarnings(trace_biosynthesis(
    targets,
    method = "virtual"
  ))
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(igraph::vcount(path), 2L)
  expect_equal(edges$from, "GalNAc(a1-")
  expect_equal(edges$to, as.character(targets[[2]]))
  expect_equal(edges$enzyme, "GalT")
})

test_that("virtual paths handle a trivial starting target", {
  glycan <- "Gal(b1-3)GalNAc(a1-"

  path <- path_biosynthesis(glycan, glycan, method = "virtual")

  expect_equal(igraph::vcount(path), 1L)
  expect_equal(igraph::ecount(path), 0L)
  expect_equal(igraph::V(path)$name, glycan)
})

test_that("virtual tracing preserves the enzymatic default", {
  target <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  enzymes <- "ST6GAL1"

  default <- path_biosynthesis(
    "Gal(b1-4)GlcNAc(b1-",
    target,
    enzymes = enzymes
  )
  explicit <- path_biosynthesis(
    "Gal(b1-4)GlcNAc(b1-",
    target,
    enzymes = enzymes,
    method = "enzymatic"
  )

  expect_equal(
    igraph::as_data_frame(default, what = "edges"),
    igraph::as_data_frame(explicit, what = "edges")
  )
})

test_that("virtual tracing validates incompatible options", {
  expect_snapshot(
    error = TRUE,
    trace_biosynthesis(
      "Gal(b1-3)GalNAc(a1-",
      enzymes = "C1GALT1",
      method = "virtual"
    )
  )
})

test_that("virtual tracing reports unreachable paths", {
  target <- "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"

  expect_snapshot(
    error = TRUE,
    trace_biosynthesis(target, method = "virtual", max_steps = 1)
  )
  expect_snapshot(
    error = TRUE,
    trace_biosynthesis(
      target,
      method = "virtual",
      filter = function(glycans) length(glycans) == 0L
    )
  )
  expect_snapshot(
    error = TRUE,
    path_biosynthesis(
      "Man(a1-3)GlcNAc(b1-",
      target,
      method = "virtual"
    )
  )
})
