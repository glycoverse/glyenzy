test_that("find_synthesis_path finds shortest path for single enzyme", {
  suppressMessages(devtools::load_all())
  
  from <- "Gal(b1-4)GlcNAc(b1-"
  to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  enzymes <- "ST6GAL1"
  
  g <- suppressMessages(find_synthesis_path(from, to, enzymes, max_steps = 3, return = "shortest"))
  
  expect_s3_class(g, "igraph")
  
  # Check edges
  edges <- igraph::as_data_frame(g, what = "edges")
  expect_equal(nrow(edges), 1L)
  expect_equal(edges$enzyme, "ST6GAL1")
  expect_equal(edges$step, 1L)
  expect_equal(edges$from, from)
  expect_equal(edges$to, to)
  
  # Check vertices
  vertices <- igraph::as_data_frame(g, what = "vertices")
  expect_equal(nrow(vertices), 2L)
  expect_true(from %in% vertices$name)
  expect_true(to %in% vertices$name)
})

test_that("find_synthesis_path works with glyrepr_structure input", {
  suppressMessages(devtools::load_all())
  
  from_g <- glyparse::auto_parse("Gal(b1-4)GlcNAc(b1-")
  to_g <- glyparse::auto_parse("Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-")
  enzymes <- "ST6GAL1"
  
  g <- suppressMessages(find_synthesis_path(from_g, to_g, enzymes, max_steps = 3, return = "shortest"))
  
  expect_s3_class(g, "igraph")
  edges <- igraph::as_data_frame(g, what = "edges")
  expect_equal(nrow(edges), 1L)
  expect_equal(edges$enzyme, "ST6GAL1")
})

test_that("find_synthesis_path returns empty graph when from equals to", {
  suppressMessages(devtools::load_all())
  
  from <- "Gal(b1-4)GlcNAc(b1-"
  to <- "Gal(b1-4)GlcNAc(b1-"
  
  g <- suppressMessages(find_synthesis_path(from, to, enzymes = "ST6GAL1", return = "shortest"))
  
  expect_s3_class(g, "igraph")
  edges <- igraph::as_data_frame(g, what = "edges")
  expect_equal(nrow(edges), 0L)
  
  vertices <- igraph::as_data_frame(g, what = "vertices")
  expect_equal(nrow(vertices), 1L)
  expect_equal(vertices$name, from)
})

test_that("find_synthesis_path fails when no path exists", {
  suppressMessages(devtools::load_all())

  from <- "Gal(b1-4)GlcNAc(b1-"
  to <- "Man(a1-3)GlcNAc(b1-"  # Unlikely to be synthesizable from from
  enzymes <- "ST6GAL1"  # This enzyme won't help

  expect_error(
    suppressMessages(find_synthesis_path(from, to, enzymes, max_steps = 2, return = "shortest")),
    "No enzymes are predicted to contribute"
  )
})

test_that("find_synthesis_path works with enzyme objects", {
  suppressMessages(devtools::load_all())
  
  from <- "Gal(b1-4)GlcNAc(b1-"
  to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  enzymes <- list(enzyme("ST6GAL1"))
  
  g <- suppressMessages(find_synthesis_path(from, to, enzymes, max_steps = 3, return = "shortest"))
  
  expect_s3_class(g, "igraph")
  edges <- igraph::as_data_frame(g, what = "edges")
  expect_equal(nrow(edges), 1L)
  expect_equal(edges$enzyme, "ST6GAL1")
})

test_that("find_synthesis_path fails with unknown enzyme", {
  suppressMessages(devtools::load_all())
  
  from <- "Gal(b1-4)GlcNAc(b1-"
  to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  enzymes <- "UNKNOWN_ENZYME"
  
  expect_error(
    find_synthesis_path(from, to, enzymes, max_steps = 3, return = "shortest"),
    "Unknown enzymes"
  )
})

test_that("find_synthesis_path works with NULL enzymes (uses all)", {
  suppressMessages(devtools::load_all())
  
  from <- "Gal(b1-4)GlcNAc(b1-"
  to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  
  g <- suppressMessages(find_synthesis_path(from, to, enzymes = NULL, max_steps = 3, return = "shortest"))
  
  expect_s3_class(g, "igraph")
  edges <- igraph::as_data_frame(g, what = "edges")
  expect_equal(nrow(edges), 1L)
  # Should find ST6GAL1 among all enzymes
  expect_true(edges$enzyme %in% names(glyenzy_enzymes))
})

test_that("find_synthesis_path return='all' includes multiple paths", {
  suppressMessages(devtools::load_all())
  
  # Use a simple case where we know there should be paths
  from <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  to <- "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzymes <- "MGAT2"
  
  g <- suppressMessages(find_synthesis_path(from, to, enzymes, max_steps = 3, return = "all"))
  
  expect_s3_class(g, "igraph")
  edges <- igraph::as_data_frame(g, what = "edges")
  expect_gte(nrow(edges), 1L)  # Should have at least one edge
})

test_that("find_synthesis_path validates input lengths", {
  suppressMessages(devtools::load_all())
  
  from <- c("Gal(b1-4)GlcNAc(b1-", "Man(a1-3)GlcNAc(b1-")  # Length > 1
  to <- "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-"
  
  expect_error(
    find_synthesis_path(from, to, enzymes = "ST6GAL1"),
    "length\\(from_g\\) == 1L && length\\(to_g\\) == 1L"
  )
})

test_that("find_synthesis_path works with filter function", {
  suppressMessages(devtools::load_all())
  
  from <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  to <- "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzymes <- "MGAT2"
  
  # Filter that keeps everything (should work normally)
  filter_fn <- function(glycans) rep(TRUE, length(glycans))
  
  g <- suppressMessages(find_synthesis_path(from, to, enzymes, max_steps = 3, 
                                           filter = filter_fn, return = "shortest"))
  
  expect_s3_class(g, "igraph")
  edges <- igraph::as_data_frame(g, what = "edges")
  expect_gte(nrow(edges), 1L)
})
