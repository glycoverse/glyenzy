test_that("virtual enzymes trace an intact glycan backward", {
  target <- "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"

  path <- trace_biosynthesis_virtual(target)
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "b4GlcNAcT"))
  expect_equal(edges$step, 1:2)
  expect_equal(edges$from[[1]], "GalNAc(a1-")
  expect_equal(edges$to[[2]], target)
  expect_null(igraph::edge_attr(path, "concrete_enzymes"))
})

test_that("virtual enzymes add sulfate groups as atomic actions", {
  target <- "Gal3S6S(b1-3)GalNAc(a1-"

  path <- trace_biosynthesis_virtual(target)
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(nrow(edges), 5L)
  expect_setequal(
    edges$enzyme,
    c("b3GalT", "3SulfoT", "6SulfoT")
  )
  expect_setequal(edges$step, 1:3)

  substrate_sizes <- vapply(
    seq_len(nrow(edges)),
    \(i) .virtual_glycan_size(glyparse::auto_parse(edges$from[[i]])),
    numeric(1)
  )
  product_sizes <- vapply(
    seq_len(nrow(edges)),
    \(i) .virtual_glycan_size(glyparse::auto_parse(edges$to[[i]])),
    numeric(1)
  )
  expect_equal(product_sizes - substrate_sizes, rep(1L, nrow(edges)))
  expect_equal(edges$step, product_sizes - 1L)
})

test_that("virtual enzymes desulfate a terminal residue before deleting it", {
  target <- "Gal6S(b1-3)GalNAc(a1-"

  path <- trace_biosynthesis_virtual(target)
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "6SulfoT"))
  expect_equal(edges$from[[1]], "GalNAc(a1-")
  expect_equal(edges$to[[1]], "Gal(b1-3)GalNAc(a1-")
  expect_equal(edges$from[[2]], edges$to[[1]])
  expect_equal(edges$to[[2]], target)
})

test_that("automatically selected virtual starts are desulfated", {
  o_path <- trace_biosynthesis_virtual("GalNAc6S(a1-")
  o_edges <- igraph::as_data_frame(o_path, what = "edges")

  expect_equal(o_edges$from, "GalNAc(a1-")
  expect_equal(o_edges$to, "GalNAc6S(a1-")
  expect_equal(o_edges$enzyme, "6SulfoT")
  expect_equal(o_edges$step, 1L)

  n_target <- paste0(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)",
    "GlcNAc(b1-4)GlcNAc6S(b1-"
  )
  n_path <- trace_biosynthesis_virtual(n_target)
  n_edges <- igraph::as_data_frame(n_path, what = "edges")

  expect_equal(nrow(n_edges), 1L)
  expect_equal(n_edges$from, as.character(glyrepr::n_glycan_core()))
  expect_equal(n_edges$to, n_target)
  expect_equal(n_edges$enzyme, "6SulfoT")
})

test_that("explicit virtual starts preserve sulfate groups", {
  from <- "GalNAc6S(a1-"
  to <- "Gal3S(b1-3)GalNAc6S(a1-"

  path <- path_biosynthesis_virtual(from, to)
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "3SulfoT"))
  expect_equal(edges$from[[1]], from)
  expect_equal(edges$to[[2]], to)
  expect_equal(edges$step, 1:2)
})

test_that("explicit sulfate starts must be a subset of the target", {
  expect_snapshot(
    error = TRUE,
    path_biosynthesis_virtual("GalNAc6S(a1-", "GalNAc(a1-")
  )
  expect_snapshot(
    error = TRUE,
    path_biosynthesis_virtual("GalNAc3S(a1-", "GalNAc6S(a1-")
  )
})

test_that("virtual enzymes include every order for independent branches", {
  target <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"

  path <- trace_biosynthesis_virtual(target)
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

  path <- trace_biosynthesis_virtual(target)
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
    trace_biosynthesis_virtual(topological)
  )
  partial_path <- suppressWarnings(
    trace_biosynthesis_virtual(partial)
  )
  basic_path <- suppressWarnings(
    trace_biosynthesis_virtual(basic)
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

test_that("virtual sulfate labels support reduced structures", {
  intact <- glyparse::auto_parse("Gal6S(b1-3)GalNAc(a1-")
  topological <- glyrepr::reduce_structure_level(intact, "topological")
  basic <- glyrepr::reduce_structure_level(intact, "basic")
  partial <- glyparse::auto_parse("Gal?S(b1-?)GalNAc(a1-")

  paths <- list(
    topological = suppressWarnings(
      trace_biosynthesis_virtual(topological)
    ),
    basic = suppressWarnings(trace_biosynthesis_virtual(basic)),
    partial = suppressWarnings(trace_biosynthesis_virtual(partial))
  )

  expect_equal(
    unname(lapply(paths, \(path) sort(igraph::E(path)$enzyme))),
    list(
      c("6SulfoT", "GalT"),
      c("6SulfoT", "HexT"),
      c("?SulfoT", "GalT")
    )
  )
})

test_that("virtual tracing uses the root for non-N-glycans", {
  target <- "Gal(b1-4)Xyl(b1-"

  path <- trace_biosynthesis_virtual(target)
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

  path <- trace_biosynthesis_virtual(targets)
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
    intact = trace_biosynthesis_virtual(intact),
    partial = suppressWarnings(
      trace_biosynthesis_virtual(partial)
    ),
    topological = suppressWarnings(
      trace_biosynthesis_virtual(topological)
    ),
    basic = suppressWarnings(
      trace_biosynthesis_virtual(basic)
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

  path <- path_biosynthesis_virtual(from, to)
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

  path <- suppressWarnings(path_biosynthesis_virtual(
    from,
    to
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

  path <- suppressWarnings(trace_biosynthesis_virtual(
    targets
  ))
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(igraph::vcount(path), 2L)
  expect_equal(edges$from, "GalNAc(a1-")
  expect_equal(edges$to, as.character(targets[[2]]))
  expect_equal(edges$enzyme, "GalT")
})

test_that("virtual paths handle a trivial starting target", {
  glycan <- "Gal(b1-3)GalNAc(a1-"

  path <- path_biosynthesis_virtual(glycan, glycan)

  expect_equal(igraph::vcount(path), 1L)
  expect_equal(igraph::ecount(path), 0L)
  expect_equal(igraph::V(path)$name, glycan)
})

test_that("enzymatic tracing has no method argument", {
  target <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  enzymes <- "ST6GAL1"

  default <- path_biosynthesis(
    "Gal(b1-4)GlcNAc(b1-",
    target,
    enzymes = enzymes
  )
  expect_equal(igraph::E(default)$enzyme, "ST6GAL1")
})

test_that("virtual tracing only accepts enzymes when annotating", {
  expect_snapshot(
    error = TRUE,
    trace_biosynthesis_virtual(
      "Gal(b1-3)GalNAc(a1-",
      enzymes = "C1GALT1"
    )
  )
})

test_that("virtual tracing reports unreachable paths and has no search controls", {
  target <- "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"

  expect_identical(
    names(formals(trace_biosynthesis_virtual)),
    c("glycans", "enzymes", "annotate_enzymes")
  )
  expect_identical(
    names(formals(path_biosynthesis_virtual)),
    c("from", "to", "enzymes", "annotate_enzymes")
  )
  expect_snapshot(
    error = TRUE,
    path_biosynthesis_virtual(
      "Man(a1-3)GlcNAc(b1-",
      target
    )
  )
})

test_that("virtual tracing annotates exact transitions in candidate order", {
  target <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  candidates <- c("ST6GAL2", "B4GALT2", "ST6GAL1", "B4GALT1")

  path <- trace_biosynthesis_virtual(
    target,
    enzymes = candidates,
    annotate_enzymes = TRUE
  )
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b4GalT", "a6Neu5AcT"))
  expect_equal(
    edges$concrete_enzymes,
    list(c("B4GALT2", "B4GALT1"), c("ST6GAL2", "ST6GAL1"))
  )

  default_path <- path_biosynthesis_virtual(
    "Gal(b1-4)GlcNAc(b1-",
    target,
    annotate_enzymes = TRUE
  )
  expect_equal(
    igraph::E(default_path)$concrete_enzymes,
    list(c("ST6GAL1", "ST6GAL2"))
  )
})

test_that("virtual tracing annotates sulfate transitions", {
  enzyme <- make_enzyme(
    name = "CUSTOM_ST",
    type = "ST",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-3)GalNAc(a1-",
      acceptor_alignment = "whole",
      rejects = NULL,
      product = "Gal3S(b1-3)GalNAc(a1-"
    ))
  )

  path <- path_biosynthesis_virtual(
    "Gal(b1-3)GalNAc(a1-",
    "Gal3S(b1-3)GalNAc(a1-",
    enzymes = list(enzyme),
    annotate_enzymes = TRUE
  )
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(edges$enzyme, "3SulfoT")
  expect_equal(edges$concrete_enzymes, list("CUSTOM_ST"))
})

test_that("virtual tracing retains unsupported annotated transitions", {
  target <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"

  path <- trace_biosynthesis_virtual(
    target,
    enzymes = "C1GALT1",
    annotate_enzymes = TRUE
  )
  edges <- igraph::as_data_frame(path, what = "edges")
  before_branch <- edges$from == "GalNAc(a1-" & edges$enzyme == "b3GalT"
  after_branch <- edges$from == "GlcNAc(b1-6)GalNAc(a1-" &
    edges$enzyme == "b3GalT"

  expect_equal(edges$concrete_enzymes[before_branch], list("C1GALT1"))
  expect_equal(edges$concrete_enzymes[after_branch], list(character()))
})

test_that("virtual tracing annotates reduced-level transitions", {
  intact_from <- glyrepr::o_glycan_core_1()
  intact_to <- glyrepr::o_glycan_core_2()
  topological_from <- glyrepr::reduce_structure_level(
    intact_from,
    "topological"
  )
  topological_to <- glyrepr::reduce_structure_level(
    intact_to,
    "topological"
  )
  basic_from <- glyrepr::reduce_structure_level(intact_from, "basic")
  basic_to <- glyrepr::reduce_structure_level(intact_to, "basic")
  partial_to <- glyparse::auto_parse(
    "Gal(b1-3)[GlcNAc(b1-?)]GalNAc(a1-"
  )

  paths <- list(
    partial = suppressWarnings(path_biosynthesis_virtual(
      intact_from,
      partial_to,
      enzymes = "GCNT1",
      annotate_enzymes = TRUE
    )),
    topological = suppressWarnings(path_biosynthesis_virtual(
      topological_from,
      topological_to,
      enzymes = "GCNT1",
      annotate_enzymes = TRUE
    )),
    basic = suppressWarnings(path_biosynthesis_virtual(
      basic_from,
      basic_to,
      enzymes = "GCNT1",
      annotate_enzymes = TRUE
    ))
  )

  expect_equal(
    unname(purrr::map(paths, ~ igraph::E(.x)$concrete_enzymes)),
    rep(list(list("GCNT1")), 3L)
  )
})

test_that("virtual tracing annotation retains custom enzyme S3 dispatch", {
  calls <- 0L
  action_method <- function(glycans, enzyme, structure_level = "intact") {
    calls <<- calls + 1L
    expect_length(glycans, 1L)
    rep(
      list(glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")),
      length(glycans)
    )
  }
  rlang::local_bindings(
    .apply_enzyme.test_hybrid_enzyme = action_method,
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
  class(custom) <- c("test_hybrid_enzyme", class(custom))

  path <- trace_biosynthesis_virtual(
    "Gal(b1-3)GalNAc(a1-",
    enzymes = list(custom),
    annotate_enzymes = TRUE
  )

  expect_equal(calls, 1L)
  expect_equal(igraph::E(path)$concrete_enzymes, list("CUSTOM"))
})

test_that("virtual annotation preserves topology and trivial paths", {
  targets <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
  )
  virtual <- trace_biosynthesis_virtual(targets)
  hybrid <- trace_biosynthesis_virtual(
    targets,
    enzymes = c("C1GALT1", "GCNT1"),
    annotate_enzymes = TRUE
  )
  virtual_edges <- igraph::as_data_frame(virtual, what = "edges")
  hybrid_edges <- igraph::as_data_frame(hybrid, what = "edges")

  expect_equal(
    hybrid_edges[c("from", "to", "enzyme", "step")],
    virtual_edges
  )
  expect_length(hybrid_edges$concrete_enzymes, nrow(virtual_edges))

  glycan <- "Gal(b1-3)GalNAc(a1-"
  trivial <- path_biosynthesis_virtual(
    glycan,
    glycan,
    enzymes = "C1GALT1",
    annotate_enzymes = TRUE
  )

  expect_equal(igraph::ecount(trivial), 0L)
  expect_equal(igraph::edge_attr(trivial, "concrete_enzymes"), list())
})
