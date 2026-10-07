# The frozen engine is an independent pre-native oracle, not a selectable backend.
bfs_reference <- function() {
  env <- new.env(parent = asNamespace("glyenzy"))
  sys.source(test_path("fixtures", "bfs-reference.R"), envir = env)
  env$bfs_synthesis_search
}

bfs_snapshot <- function(x) {
  for (name in c("parent", "parent_enzyme", "parent_step")) {
    value <- as.list(x[[name]])
    x[[name]] <- if (length(value)) value[order(names(value))] else list()
  }
  x
}

test_that("native BFS matches frozen R traversal for every bundled rule", {
  reference <- bfs_reference()
  enzymes <- c(
    db_enzymes(include_starter_gt = FALSE, include_npre_gt = FALSE),
    abstract_enzymes()
  )
  for (e in enzymes) {
    for (r in e$rules) {
      for (topological in c(FALSE, TRUE)) {
        from <- r$acceptor
        to <- r$product
        if (topological) {
          from <- glyrepr::remove_linkages(from)
          to <- glyrepr::remove_linkages(to)
        }
        args <- list(
          from_g = from,
          to_gs = to,
          enzymes = list(e),
          max_steps = 1L,
          structure_level = if (topological) "topological" else "intact",
          target_match = if (topological) "whole" else "key",
          allow_partial = TRUE
        )
        actual <- do.call(bfs_synthesis_search, args)
        expected <- do.call(reference, args)
        expect_identical(
          bfs_snapshot(actual),
          bfs_snapshot(expected),
          info = paste(e$name, as.character(r$acceptor), topological)
        )
      }
    }
  }
})

test_that("native search preserves generic, sulfate and multiple-target traversal", {
  reference <- bfs_reference()
  from <- glyrepr::as_glycan_structure("GalNAc(a1-")
  targets <- glyrepr::as_glycan_structure(c(
    "Gal3S(b1-3)GalNAc(a1-",
    "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
  ))
  enzymes <- lapply(c("C1GALT1", "GCNT1", "B4GALT1", "GAL3ST4"), enzyme)
  for (generic in c(FALSE, TRUE)) {
    for (topological in c(FALSE, TRUE)) {
      start <- from
      to <- targets
      if (generic) {
        start <- glyrepr::convert_to_generic(start)
        to <- glyrepr::convert_to_generic(to)
      }
      if (topological) {
        start <- glyrepr::remove_linkages(start)
        to <- glyrepr::remove_linkages(to)
      }
      args <- list(
        from_g = start,
        to_gs = to,
        enzymes = enzymes,
        max_steps = 4L,
        structure_level = if (topological) "topological" else "intact",
        target_match = if (topological || generic) "whole" else "key",
        allow_partial = TRUE
      )
      expect_identical(
        bfs_snapshot(do.call(bfs_synthesis_search, args)),
        bfs_snapshot(do.call(reference, args))
      )
    }
  }
})

test_that("native scalar execution preserves filter values and invocation order", {
  reference <- bfs_reference()
  execute <- function(search) {
    seen <- list()
    result <- search(
      glyrepr::as_glycan_structure("GalNAc(a1-"),
      glyrepr::as_glycan_structure(
        "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-"
      ),
      lapply(c("C1GALT1", "GCNT1", "B4GALT1", "B4GALT2"), enzyme),
      max_steps = 4L,
      filter = function(x) {
        seen[[length(seen) + 1L]] <<- as.character(x)
        rep(length(seen) %% 3L != 0L, length(x))
      },
      allow_partial = TRUE
    )
    list(result = bfs_snapshot(result), seen = seen)
  }
  expect_identical(execute(bfs_synthesis_search), execute(reference))
})

test_that("native BFS retains locale ordering and alditol keys", {
  reference <- bfs_reference()
  from <- glyrepr::as_glycan_structure("GalNAc-ol(a1-")
  to <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc-ol(a1-")
  for (locale in unique(c(Sys.getlocale("LC_COLLATE"), "C"))) {
    withr::local_collate(locale)
    args <- list(
      from_g = from,
      to_gs = to,
      enzymes = list(enzyme("C1GALT1")),
      max_steps = 1L,
      allow_partial = TRUE
    )
    expect_identical(
      bfs_snapshot(do.call(bfs_synthesis_search, args)),
      bfs_snapshot(do.call(reference, args))
    )
  }
})

test_that("attribute-bearing graphs use callbacks inside the native engine", {
  reference <- bfs_reference()
  graph <- glyrepr::get_structure_graphs(glyrepr::as_glycan_structure(
    "GalNAc(a1-"
  ))
  igraph::graph_attr(graph, "experiment") <- "retained"
  from <- glyrepr::new_glycan_structure(
    "GalNAc(a1-",
    list("GalNAc(a1-" = graph)
  )
  to <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  engine <- BfsSynthesisSearch$new(
    from,
    to,
    list(enzyme("C1GALT1")),
    max_steps = 1L
  )
  actual <- engine$run()
  expected <- reference(from, to, list(enzyme("C1GALT1")), max_steps = 1L)
  expect_identical(bfs_snapshot(actual), bfs_snapshot(expected))
  expect_identical(engine$native_stats$callback_cells, 1L)
  expect_identical(
    igraph::graph_attr(engine$queue_graphs[[1]], "experiment"),
    "retained"
  )
})

test_that("native BFS handles initial hits, duplicate goals, limits and no enzymes", {
  reference <- bfs_reference()
  from <- glyrepr::as_glycan_structure("GalNAc(a1-")
  to <- glyrepr::as_glycan_structure(c(
    "GalNAc(a1-",
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-3)GalNAc(a1-"
  ))
  for (steps in c(0:2, 1.5, Inf)) {
    for (enzymes in list(list(), list(enzyme("C1GALT1"), enzyme("C1GALT1")))) {
      args <- list(
        from_g = from,
        to_gs = to,
        enzymes = enzymes,
        max_steps = steps,
        allow_partial = TRUE
      )
      expect_identical(
        bfs_snapshot(do.call(bfs_synthesis_search, args)),
        bfs_snapshot(do.call(reference, args))
      )
    }
  }
})

test_that("native multi-level GT and GH searches preserve all first-parent ties", {
  reference <- bfs_reference()
  env <- new.env(parent = asNamespace("glyenzy"))
  sys.source(test_path("fixtures", "bfs-native-cases.R"), envir = env)
  for (case in env$cases) {
    args <- list(
      from_g = glyrepr::as_glycan_structure(case[[1]]),
      to_gs = glyrepr::as_glycan_structure(case[[2]]),
      enzymes = lapply(case[[3]], enzyme),
      max_steps = case[[4]],
      allow_partial = TRUE
    )
    engine <- do.call(BfsSynthesisSearch$new, args)
    actual <- engine$run()
    expect_identical(engine$native_stats$callback_cells, 0L)
    expect_length(actual$missing_target_keys, 0L)
    expect_identical(
      bfs_snapshot(actual),
      bfs_snapshot(do.call(reference, args))
    )
  }
})

test_that("native traversal preserves floating graph callbacks", {
  reference <- bfs_reference()
  from <- glyrepr::as_glycan_structure("{6S}Gal(b1-3)GalNAc(a1-")
  to <- glyrepr::as_glycan_structure("{6S}Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
  args <- list(
    from_g = from,
    to_gs = to,
    enzymes = list(enzyme("GCNT1")),
    max_steps = 1L,
    target_match = "whole",
    allow_partial = TRUE
  )
  expect_identical(
    bfs_snapshot(do.call(bfs_synthesis_search, args)),
    bfs_snapshot(do.call(reference, args))
  )
})

test_that("native filters retain validation errors and custom action conditions", {
  from <- glyrepr::as_glycan_structure("GalNAc(a1-")
  to <- glyrepr::as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  capture <- function(search, filter) {
    tryCatch(
      search(from, to, list(enzyme("C1GALT1")), 1L, filter = filter),
      error = function(e) list(class = class(e), message = conditionMessage(e))
    )
  }
  reference <- bfs_reference()
  for (filter in list(function(x) NA, function(x) logical(), function(x) 1)) {
    expect_identical(
      capture(bfs_synthesis_search, filter),
      capture(reference, filter)
    )
  }
})
