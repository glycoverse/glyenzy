test_that("abstract trimming produces exactly the concrete enzyme products", {
  groups <- list(
    GlcH = c("MOGS", "GANAB"),
    ManI = c("MAN1B1", "MAN1A1", "MAN1A2", "MAN1C1"),
    ManII = c("MAN2A1", "MAN2A2")
  )
  reactions <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  glycans <- glyparse::auto_parse(unique(c(
    reactions$substrate,
    reactions$product
  )))
  for (name in names(groups)) {
    actual <- apply_enzyme(glycans, abstract_enzymes()[[name]])
    concrete <- lapply(groups[[name]], function(gene) {
      apply_enzyme(glycans, enzyme(gene))
    })
    for (i in seq_along(glycans)) {
      expected <- unique(unlist(
        lapply(concrete, function(products) {
          as.character(products[[i]])
        }),
        use.names = FALSE
      ))
      expect_setequal(
        as.character(actual[[i]]),
        as.character(expected)
      )
    }
  }
})

test_that("ManI does not trim glucosylated or GlcNAc-substituted precursors", {
  reactions <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  glucosylated <- reactions$substrate[reactions$enzyme == "GlcH"]
  man8 <- paste0(
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  substituted <- paste0("GlcNAc(b1-2)", man8)
  glycans <- c(glucosylated, substituted)
  expect_identical(
    lengths(apply_enzyme(glycans, abstract_enzymes()$ManI)),
    rep(0L, length(glycans))
  )
})

test_that("high-mannose tracing has the concrete vertices and directed edges", {
  reactions <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  targets <- unique(reactions$product[reactions$enzyme %in% c("GlcH", "ManI")])
  targets <- as.character(glyparse::auto_parse(targets))
  abstract <- trace_biosynthesis(
    targets,
    enzymes = abstract_enzymes(),
    max_steps = 7L,
    max_virtual_steps = 0L
  )
  concrete <- trace_biosynthesis(
    targets,
    max_steps = 7L,
    max_virtual_steps = 0L
  )
  expect_equal(igraph::vcount(abstract), 15L)
  expect_equal(igraph::ecount(abstract), 19L)
  expect_setequal(igraph::V(abstract)$name, igraph::V(concrete)$name)
  expect_setequal(igraph::V(abstract)$name[igraph::V(abstract)$target], targets)
  expect_setequal(igraph::V(concrete)$name[igraph::V(concrete)$target], targets)

  transitions <- function(network) {
    edges <- igraph::as_data_frame(network, what = "edges")
    edges <- edges[
      order(edges$from, edges$to),
      c("from", "to", "step", "is_virtual")
    ]
    rownames(edges) <- NULL
    edges
  }
  expect_identical(transitions(abstract), transitions(concrete))
  expect_identical(igraph::E(abstract)$is_virtual, rep(FALSE, 19L))
  expect_equal(igraph::any_multiple(abstract), FALSE)
  expect_setequal(unlist(igraph::E(abstract)$enzymes), c("GlcH", "ManI"))
})

test_that("ManII has the same two-route trimming network as MAN2A1 and MAN2A2", {
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  man5 <- paste0("GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]", core)
  man3 <- paste0("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]", core)
  intermediates <- paste0(
    c(
      "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]",
      "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]"
    ),
    core
  )
  for (fucosylated in c(FALSE, TRUE)) {
    glycans <- c(man5, man3, intermediates)
    if (fucosylated) {
      glycans <- sub(
        "GlcNAc(b1-4)GlcNAc(b1-",
        "GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",
        glycans,
        fixed = TRUE
      )
    }
    glycans <- as.character(glyparse::auto_parse(glycans))
    abstract <- path_biosynthesis(
      glycans[[1]],
      glycans[[2]],
      enzymes = abstract_enzymes()["ManII"],
      max_steps = 2L,
      max_virtual_steps = 0L
    )
    concrete <- path_biosynthesis(
      glycans[[1]],
      glycans[[2]],
      enzymes = c("MAN2A1", "MAN2A2"),
      max_steps = 2L,
      max_virtual_steps = 0L
    )
    expect_setequal(igraph::V(abstract)$name, glycans)
    expect_setequal(igraph::V(abstract)$name, igraph::V(concrete)$name)
    expect_equal(igraph::ecount(abstract), 4L)
    edge_keys <- function(network) {
      edges <- igraph::as_data_frame(network, what = "edges")
      paste(edges$from, edges$to, sep = " -> ")
    }
    expect_setequal(edge_keys(abstract), edge_keys(concrete))
    expect_setequal(
      igraph::V(abstract)$name[igraph::V(abstract)$target],
      glycans[[2]]
    )
    expect_identical(igraph::E(abstract)$is_virtual, rep(FALSE, 4L))
    expect_equal(igraph::any_multiple(abstract), FALSE)
    expect_setequal(unlist(igraph::E(abstract)$enzymes), "ManII")
  }
})

test_that("tracing through ManII matches concrete trimming from the default precursor", {
  target <- paste0(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  abstract <- trace_biosynthesis(
    target,
    enzymes = abstract_enzymes(),
    max_virtual_steps = 0L
  )
  concrete <- trace_biosynthesis(target, max_virtual_steps = 0L)
  expect_equal(igraph::vcount(abstract), 19L)
  expect_equal(igraph::ecount(abstract), 24L)
  expect_setequal(igraph::V(abstract)$name, igraph::V(concrete)$name)
  edges <- function(network) {
    res <- igraph::as_data_frame(network, what = "edges")
    paste(res$from, res$to, res$step, res$is_virtual, sep = "|")
  }
  expect_setequal(edges(abstract), edges(concrete))
  reactions <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  start <- as.character(glyparse::auto_parse(
    reactions$substrate[reactions$enzyme == "GlcH"][[1]]
  ))
  expect_equal(
    as.numeric(igraph::distances(
      abstract,
      v = start,
      to = as.character(glyparse::auto_parse(target)),
      mode = "out"
    )),
    10
  )
})
