test_that("path fallback resumes concrete synthesis after an early virtual step", {
  path <- path_biosynthesis(
    "GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-",
    enzymes = "ST3GAL1",
    max_steps = 2,
    max_virtual_steps = 1
  )
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "ST3GAL1"))
  expect_identical(edges$is_virtual, c(TRUE, FALSE))
  expect_equal(edges$step, 1:2)
})

test_that("an early fallback can precede multiple concrete steps", {
  add_glcnac <- make_enzyme(
    name = "ADD_GLCNAC",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-3)GalNAc(a1-",
      acceptor_alignment = "whole",
      rejects = NULL,
      product = "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
    ))
  )
  add_sia <- make_enzyme(
    name = "ADD_SIA",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
      acceptor_alignment = "whole",
      rejects = NULL,
      product = "Neu5Ac(a2-3)GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-"
    ))
  )

  path <- path_biosynthesis(
    "GalNAc(a1-",
    "Neu5Ac(a2-3)GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
    enzymes = list(add_glcnac, add_sia),
    max_steps = 3,
    max_virtual_steps = 1
  )
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "ADD_GLCNAC", "ADD_SIA"))
  expect_identical(edges$is_virtual, c(TRUE, FALSE, FALSE))
  expect_equal(edges$step, 1:3)
})

test_that("trace fallback resumes concrete synthesis after a virtual step", {
  target <- "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"

  path <- trace_biosynthesis(
    target,
    enzymes = "ST3GAL1",
    max_steps = 2,
    max_virtual_steps = 1
  )
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "ST3GAL1"))
  expect_identical(edges$is_virtual, c(TRUE, FALSE))
  expect_equal(edges$to[[2]], target)
})

test_that("trace fallback combines strict and rescued targets", {
  targets <- c(
    "Gal(b1-3)GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
  )

  path <- trace_biosynthesis(
    targets,
    enzymes = "C1GALT1",
    max_steps = 2,
    max_virtual_steps = 1
  )
  edges <- igraph::as_data_frame(path, what = "edges")

  expect_equal(
    igraph::V(path)$name,
    c(
      "GalNAc(a1-",
      targets[[1]],
      targets[[2]]
    )
  )
  expect_equal(edges$enzyme, c("C1GALT1", "a3Neu5AcT"))
  expect_identical(edges$is_virtual, c(FALSE, TRUE))
})

test_that("fallback minimizes virtual steps before total steps", {
  target <- "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"

  path <- path_biosynthesis(
    "GalNAc(a1-",
    target,
    enzymes = "GCNT1",
    max_steps = 2,
    max_virtual_steps = 2
  )
  edges <- igraph::as_data_frame(path, what = "edges")
  edges <- edges[order(edges$step), ]

  expect_equal(edges$enzyme, c("b3GalT", "GCNT1"))
  expect_identical(edges$is_virtual, c(TRUE, FALSE))
  expect_equal(
    igraph::V(path)$name,
    c(
      "GalNAc(a1-",
      "Gal(b1-3)GalNAc(a1-",
      target
    )
  )
})

test_that("fallback respects the virtual-step limit", {
  call <- function(max_virtual_steps) {
    path_biosynthesis(
      "GalNAc(a1-",
      "GlcNAc(b1-4)Gal(b1-3)GalNAc(a1-",
      enzymes = "ST6GAL1",
      max_steps = 2,
      max_virtual_steps = max_virtual_steps
    )
  }

  expect_snapshot(error = TRUE, call(1L))

  path <- call(2L)
  edges <- igraph::as_data_frame(path, what = "edges")
  expect_equal(edges$enzyme, c("b3GalT", "b4GlcNAcT"))
  expect_identical(edges$is_virtual, c(TRUE, TRUE))
})

test_that("filters apply to virtual products", {
  reject_core_1 <- function(glycans) {
    as.character(glycans) != "Gal(b1-3)GalNAc(a1-"
  }

  expect_snapshot(
    error = TRUE,
    path_biosynthesis(
      "GalNAc(a1-",
      "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-",
      enzymes = "ST3GAL1",
      max_steps = 2,
      filter = reject_core_1,
      max_virtual_steps = 1
    )
  )
})

test_that("fully enzymatic paths are unchanged when fallback is enabled", {
  args <- list(
    from = "Gal(b1-4)GlcNAc(b1-",
    to = "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-",
    enzymes = "ST6GAL1",
    max_steps = 1
  )

  strict <- do.call(path_biosynthesis, args)
  fallback <- do.call(
    path_biosynthesis,
    c(args, max_virtual_steps = 1L)
  )

  expect_equal(
    igraph::as_data_frame(fallback, what = "edges"),
    igraph::as_data_frame(strict, what = "edges")
  )
  expect_null(igraph::edge_attr(fallback, "is_virtual"))
})

test_that("max_virtual_steps is appended to preserve positional filters", {
  expect_identical(
    names(formals(trace_biosynthesis)),
    c("glycans", "enzymes", "max_steps", "filter", "max_virtual_steps")
  )
  expect_identical(
    names(formals(path_biosynthesis)),
    c(
      "from",
      "to",
      "enzymes",
      "max_steps",
      "filter",
      "max_virtual_steps"
    )
  )
})
