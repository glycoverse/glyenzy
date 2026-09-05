test_that("abstract core transferases retain concrete acceptors and rejects", {
  groups <- list(
    GnTI = "MGAT1",
    GnTII = "MGAT2",
    GnTIII = "MGAT3",
    a6FucT = "FUT8",
    GnTIV = c("MGAT4A", "MGAT4B"),
    GnTV = c("MGAT5", "MGAT5B")
  )
  cases <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  arm3 <- "GlcNAc(b1-2)Man(a1-3)"
  arm6 <- "[GlcNAc(b1-2)Man(a1-6)]"
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  base <- paste0(arm3, arm6, core)
  extra <- c(
    paste0("Man(a1-3)[Man(a1-6)]", core),
    paste0(arm3, core),
    paste0("Gal(b1-4)", base),
    paste0(arm3, "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]", core),
    paste0(arm3, arm6, "[GlcNAc(b1-4)]", core),
    paste0(arm3, arm6, "[Gal(b1-4)GlcNAc(b1-4)]", core),
    paste0("GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)", arm6, core),
    paste0("GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)[Man(a1-6)]", core)
  )
  glycans <- glyparse::auto_parse(unique(c(
    cases$substrate,
    cases$product,
    extra
  )))
  for (name in names(groups)) {
    actual <- apply_enzyme(glycans, abstract_enzymes()[[name]])
    concrete <- lapply(groups[[name]], function(gene) {
      apply_enzyme(glycans, enzyme(gene))
    })
    for (i in seq_along(glycans)) {
      expected <- unlist(
        lapply(concrete, function(products) {
          as.character(products[[i]])
        }),
        use.names = FALSE
      )
      expect_setequal(as.character(actual[[i]]), as.character(expected))
    }
  }
})

test_that("FUT8 and MGAT1-2 restrictions are retained despite empty rejects", {
  enzymes <- abstract_enzymes()
  for (name in c("a6FucT", "GnTI", "GnTII")) {
    expect_identical(
      sum(lengths(lapply(enzymes[[name]]$rules, `[[`, "rejects"))),
      0L
    )
  }
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  man3 <- paste0("Man(a1-3)[Man(a1-6)]", core)
  substrate <- paste0("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]", core)
  expect_identical(
    as.character(apply_enzyme(man3, enzymes$GnTI)),
    as.character(glyparse::auto_parse(substrate))
  )
  premature_branch <- paste0(
    "GlcNAc(b1-2)[GlcNAc(b1-4)]Man(a1-3)[Man(a1-6)]",
    core
  )
  galactosylated <- paste0("Gal(b1-4)", substrate)
  expect_identical(
    lengths(apply_enzyme(
      c(substrate, premature_branch, galactosylated),
      enzymes$GnTII
    )),
    c(1L, 0L, 0L)
  )
  missing_arm <- paste0("GlcNAc(b1-2)Man(a1-3)", core)
  expect_identical(
    lengths(apply_enzyme(c(substrate, missing_arm, man3), enzymes$a6FucT)),
    c(1L, 0L, 0L)
  )
})

test_that("MGAT3-5 rejects distinguish bisecting GlcNAc and individual antennae", {
  arm3 <- "GlcNAc(b1-2)Man(a1-3)"
  arm6 <- "[GlcNAc(b1-2)Man(a1-6)]"
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  glycans <- c(
    paste0(arm3, arm6, core),
    paste0("Gal(b1-4)", arm3, arm6, core),
    paste0(arm3, "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]", core),
    paste0(arm3, arm6, "[GlcNAc(b1-4)]", core),
    paste0(arm3, arm6, "[Gal(b1-4)GlcNAc(b1-4)]", core),
    paste0("GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)", arm6, core)
  )
  expected <- list(
    GnTIII = c(1L, 0L, 1L, 0L, 0L, 0L),
    GnTIV = c(1L, 0L, 1L, 0L, 0L, 0L),
    GnTV = c(1L, 1L, 0L, 0L, 0L, 1L)
  )
  enzymes <- abstract_enzymes()
  for (name in names(expected)) {
    expect_identical(
      lengths(apply_enzyme(glycans, enzymes[[name]])),
      expected[[name]]
    )
    enzyme <- enzymes[[name]]
    parsed <- glyparse::auto_parse(glycans)
    graphs <- glyrepr::get_structure_graphs(parsed, return_list = TRUE)
    plan <- .prepare_bfs_rule_plan(list(enzyme))
    products <- .apply_bfs_rule_frontier(
      graphs,
      rep("strict", length(graphs)),
      plan,
      "intact"
    )
    for (i in seq_along(graphs)) {
      direct <- apply_enzyme(parsed[i], enzyme)
      prepared <- .apply_enzyme_prepared(
        graphs[[i]],
        enzyme,
        lapply(enzyme$rules, .prepare_rule_graphs),
        "intact",
        "strict"
      )
      raw <- .collect_bfs_rule_products(products, plan$enzyme_rule_ids[[1]], i)
      batched <- .materialize_product_graphs(raw, enzyme)
      expect_setequal(as.character(prepared), as.character(direct))
      expect_setequal(as.character(batched), as.character(direct))
    }
  }
})

test_that("a3SiaT excludes alpha3-fucosylated acceptors one site at a time", {
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  alpha3 <- "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)Man(a1-3)"
  alpha6 <- "Gal(b1-4)GlcNAc(b1-2)Man(a1-6)"
  substrate <- paste0(alpha3, "[", alpha6, "]", core)
  product <- paste0(alpha3, "[Neu5Ac(a2-3)", alpha6, "]", core)
  expect_identical(
    as.character(apply_enzyme(substrate, abstract_enzymes()$a3SiaT)),
    as.character(glyparse::auto_parse(product))
  )
  for (gene in c("ST3GAL4", "ST3GAL6")) {
    expect_identical(
      as.character(apply_enzyme(substrate, gene)),
      as.character(glyparse::auto_parse(product))
    )
  }
})

test_that("core branching and bisecting networks match the concrete database", {
  cases <- utils::read.csv(test_path("fixtures", "abstract-n-glycans.csv"))
  cases <- cases[cases$name %in% c("A4G0", "G0FB"), ]
  for (i in seq_len(nrow(cases))) {
    target <- cases$glycan[[i]]
    abstract <- trace_biosynthesis(
      target,
      enzymes = abstract_enzymes(),
      max_virtual_steps = 0L
    )
    concrete <- trace_biosynthesis(target, max_virtual_steps = 0L)
    expect_setequal(igraph::V(abstract)$name, igraph::V(concrete)$name)
    transitions <- function(network) {
      edges <- igraph::as_data_frame(network, what = "edges")
      paste(edges$from, edges$to, edges$step, edges$is_virtual, sep = "|")
    }
    expect_setequal(transitions(abstract), transitions(concrete))
    expect_setequal(
      igraph::V(abstract)$name[igraph::V(abstract)$target],
      as.character(glyparse::auto_parse(target))
    )
    expect_identical(
      igraph::E(abstract)$is_virtual,
      rep(FALSE, igraph::ecount(abstract))
    )
    expect_equal(igraph::any_multiple(abstract), FALSE)
    if (cases$name[[i]] == "G0FB") {
      expect_equal("GnTIII" %in% unlist(igraph::E(abstract)$enzymes), TRUE)
    }
  }
})
