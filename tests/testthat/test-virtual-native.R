test_that("native trimming preserves ordered networks across N-glycan targets", {
  targets <- read.csv(test_path("fixtures", "abstract-n-glycans.csv"))$glycan
  targets <- glyparse::auto_parse(targets)
  from <- .decide_virtual_starting_glycan(targets[1])
  native <- .perform_virtual_synthesis(from, targets)
  local_mocked_bindings(.virtual_native_synthesis = function(...) NULL)
  reference <- .perform_virtual_synthesis(from, targets)

  expect_identical(
    igraph::as_data_frame(native),
    igraph::as_data_frame(reference)
  )
  expect_identical(
    igraph::as_data_frame(native, "vertices"),
    igraph::as_data_frame(reference, "vertices")
  )
})

test_that("native trimming preserves sulfation, alditols and explicit starts", {
  cases <- list(
    c("GalNAc(a1-", "Gal3S6S(b1-3)GalNAc(a1-"),
    c("Gal3S(b1-3)GalNAc(a1-", "Gal3S6S(b1-3)GalNAc(a1-"),
    c("GalNAc-ol(a1-", "Gal(b1-3)[GlcNAc(b1-6)]GalNAc-ol(a1-"),
    c("GlcNAc(b1-", "Gal(b1-4)GlcNAc(b1-")
  )
  native <- lapply(cases, function(x) path_biosynthesis_virtual(x[1], x[2]))
  local_mocked_bindings(.virtual_native_synthesis = function(...) NULL)
  reference <- lapply(cases, function(x) path_biosynthesis_virtual(x[1], x[2]))
  expect_identical(
    lapply(native, igraph::as_data_frame),
    lapply(reference, igraph::as_data_frame)
  )
})

test_that("native annotation preserves default and reordered enzyme matches", {
  targets <- c(
    "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
    "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"
  )
  enzymes <- .process_enzymes_arg(NULL, apply_prefilter = FALSE)
  for (target in targets) {
    graph <- trace_biosynthesis_virtual(target)
    native <- .amplify_virtual_edges(graph, enzymes)
    reference <- .amplify_virtual_edges_r(graph, enzymes)
    expect_identical(
      igraph::as_data_frame(native),
      igraph::as_data_frame(reference)
    )
  }
  enzymes <- .process_enzymes_arg(
    c("ST6GAL2", "B4GALT2", "ST6GAL1", "B4GALT1", "ST6GAL2"),
    apply_prefilter = FALSE
  )
  graph <- trace_biosynthesis_virtual(targets[1])
  expect_identical(
    igraph::E(.amplify_virtual_edges(graph, enzymes))$concrete_enzymes,
    igraph::E(.amplify_virtual_edges_r(graph, enzymes))$concrete_enzymes
  )
})

test_that("native annotation preserves glycan-class restrictions at beta-Man roots", {
  roots <- c("Man(b1-", "Man6S(b1-", "Man3S6S(b1-", "Man-ol(b1-")
  enzymes <- .process_enzymes_arg("MGAT1", apply_prefilter = FALSE)
  for (root in roots) {
    from <- paste0("Man(a1-3)Man(b1-4)", root)
    to <- paste0("GlcNAc(b1-2)", from)
    graph <- path_biosynthesis_virtual(from, to)
    native <- .amplify_virtual_edges(graph, enzymes)
    reference <- .amplify_virtual_edges_r(graph, enzymes)
    expected <- if (identical(root, "Man(b1-")) "MGAT1" else character()

    expect_identical(igraph::E(reference)$concrete_enzymes, list(expected))
    expect_identical(
      igraph::E(native)$concrete_enzymes,
      igraph::E(reference)$concrete_enzymes
    )
  }
})

test_that("modified built-in rules retain their annotation semantics", {
  candidate <- enzyme("C1GALT1")
  candidate$glycan_type <- "N"
  graph <- trace_biosynthesis_virtual("Gal(b1-3)GalNAc(a1-")
  annotated <- .amplify_virtual_edges(graph, list(candidate))
  expect_identical(igraph::E(annotated)$concrete_enzymes, list(character()))
})

test_that("reduced structures retain R matching and resolve default enzymes", {
  target <- glyrepr::remove_linkages(glyrepr::o_glycan_core_2())
  expect_null(.virtual_native_records(target))
  native <- suppressWarnings(trace_biosynthesis_virtual(
    target,
    annotate_enzymes = TRUE
  ))
  local_mocked_bindings(.bfs_default_enzymes = function(...) NULL)
  reference <- suppressWarnings(trace_biosynthesis_virtual(
    target,
    annotate_enzymes = TRUE
  ))
  expect_identical(
    igraph::as_data_frame(native),
    igraph::as_data_frame(reference)
  )
})
