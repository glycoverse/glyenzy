test_that("native prefilter preserves concrete, generic and topological selections", {
  enzymes <- .enzymes_from_arg(NULL)
  targets <- glyparse::auto_parse(c(
    "Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-",
    "Neu5Ac(a2-3)Gal(b1-4)Glc(b1-",
    "Gal6S(b1-4)GlcNAc6S(b1-3)Gal(b1-3)GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  for (generic in c(FALSE, TRUE)) {
    for (topological in c(FALSE, TRUE)) {
      glycans <- targets
      if (generic) {
        glycans <- glyrepr::convert_to_generic(glycans)
      }
      if (topological) {
        glycans <- glyrepr::remove_linkages(glycans)
      }
      for (target in c(
        lapply(seq_along(glycans), \(i) glycans[i]),
        list(glycans)
      )) {
        expected <- vapply(
          enzymes,
          .enzyme_contributes_to_targets,
          logical(1),
          targets = target
        )
        actual <- .bfs_prefilter(enzymes, target)
        expect_identical(actual, expected)
      }
    }
  }
})

test_that("native prefilter preserves all N-glycan trimming decisions", {
  enzymes <- .enzymes_from_arg(c(
    "MOGS",
    "GANAB",
    "MAN1B1",
    "MAN1A1",
    "MAN1A2",
    "MAN1C1",
    "MAN2A1",
    "MAN2A2",
    "MGAT1"
  ))
  targets <- unique(do.call(
    c,
    lapply(unname(enzymes), function(e) {
      do.call(c, lapply(e$rules, function(r) c(r$acceptor, r$product)))
    })
  ))
  for (generic in c(FALSE, TRUE)) {
    glycans <- if (generic) glyrepr::convert_to_generic(targets) else targets
    for (target in lapply(seq_along(glycans), \(i) glycans[i])) {
      expected <- vapply(
        enzymes,
        .enzyme_contributes_to_targets,
        logical(1),
        targets = target
      )
      expect_identical(.bfs_prefilter(enzymes, target), expected)
    }
  }
})

test_that("native prefilter falls back for modified enzymes and graph metadata", {
  target <- glyparse::auto_parse("Gal(b1-3)GalNAc(a1-")
  custom <- enzyme("C1GALT1")
  custom$rules <- enzyme("B3GNT6")$rules
  enzymes <- list(original = enzyme("C1GALT1"), edited = custom)
  expect_identical(
    .bfs_prefilter(enzymes, target),
    c(original = TRUE, edited = NA)
  )
  expect_identical(
    .can_enzymes_contribute(enzymes, target),
    c(original = TRUE, edited = FALSE)
  )
  graph <- igraph::set_graph_attr(
    glyrepr::get_structure_graphs(target),
    "custom",
    "kept"
  )
  target <- glyrepr::new_glycan_structure(
    as.character(target),
    stats::setNames(list(graph), as.character(target))
  )
  expect_identical(
    .bfs_prefilter(enzymes, target),
    c(original = NA, edited = NA)
  )
  expect_identical(
    .can_enzymes_contribute(enzymes, target),
    c(original = TRUE, edited = FALSE)
  )
})

test_that("native prefilter failures retain the original per-enzyme fallback", {
  enzymes <- .enzymes_from_arg(c("C1GALT1", "B3GNT6"))
  target <- glyparse::auto_parse("Gal(b1-3)GalNAc(a1-")
  local_mocked_bindings(cpp_bfs_have_motifs = function(...) stop("unavailable"))
  expect_identical(
    .can_enzymes_contribute(enzymes, target),
    c(C1GALT1 = TRUE, B3GNT6 = FALSE)
  )
})
