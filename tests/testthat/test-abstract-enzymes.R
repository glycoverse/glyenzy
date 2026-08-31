test_that("abstract enzymes have an independent collection and subclass", {
  ordinary <- db_enzymes()
  enzymes <- abstract_enzymes()
  expected <- c(
    "ManI",
    "GnTI",
    "ManII",
    "GnTII",
    "a6FucT",
    "GnTIV",
    "GnTV",
    "b4GalT",
    "iGnT",
    "a3SiaT",
    "GlcH"
  )

  expect_named(enzymes, expected)
  expect_identical(abstract_enzymes(TRUE), expected)
  expect_identical(sum(lengths(lapply(enzymes, `[[`, "rules"))), 28L)
  expect_identical(
    unname(vapply(enzymes, `[[`, character(1), "localization")),
    c("cis", "cis", rep("medial", 4), rep("trans", 4), "ER")
  )
  for (enzyme in enzymes) {
    subtype <- if (enzyme$name %in% c("ManI", "ManII", "GlcH")) "gh" else "gt"
    expect_identical(
      class(enzyme),
      c(
        "glyenzy_abstract_enzyme",
        paste0("glyenzy_", subtype, "_enzyme"),
        "glyenzy_enzyme"
      )
    )
    expect_identical(enzyme$glycan_type, "N")
    expect_identical(enzyme$species, "unspecified")
    expect_invisible(validate_enzyme(enzyme))
  }
  expect_identical(db_enzymes(), ordinary)
  enzymes$ManI$name <- "changed"
  expect_identical(abstract_enzymes()$ManI$name, "ManI")
})

test_that("abstract enzymes remain absent from existing discovery interfaces", {
  abstract_names <- abstract_enzymes(TRUE)
  for (starter in c(FALSE, TRUE)) {
    for (precursor in c(FALSE, TRUE)) {
      expect_length(
        intersect(
          db_enzymes(TRUE, starter, precursor),
          abstract_names
        ),
        0L
      )
    }
  }
  mat <- matrix(
    10,
    nrow = length(abstract_names),
    ncol = 1,
    dimnames = list(abstract_names, "sample")
  )
  expect_length(enzymes_from_rnaseq(mat), 0L)
  expect_length(intersect(names(.enzymes_from_arg(NULL)), abstract_names), 0L)

  target <- paste0(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_length(intersect(find_enzyme(target), abstract_names), 0L)
  network <- trace_biosynthesis(target)
  expect_length(
    intersect(unlist(igraph::E(network)$enzymes), abstract_names),
    0L
  )
  expect_snapshot(error = TRUE, enzyme("ManI"))
  expect_snapshot(error = TRUE, .enzymes_from_arg(abstract_names))
  expect_snapshot(error = TRUE, abstract_enzymes(NA))
})

test_that("printing abstract enzymes includes localization and rejects", {
  without_blank_lines <- function(lines) lines[nzchar(trimws(lines))]
  expect_snapshot(
    print(abstract_enzymes()$ManII),
    transform = without_blank_lines
  )
  expect_snapshot(
    print(abstract_enzymes()$iGnT),
    transform = without_blank_lines
  )
})

test_that("each abstract reaction produces its independently specified product", {
  cases <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  enzymes <- abstract_enzymes()
  expect_equal(nrow(cases), 31L)
  for (i in seq_len(nrow(cases))) {
    expected <- as.character(glyparse::auto_parse(cases$product[[i]]))
    actual <- as.character(apply_enzyme(
      cases$substrate[[i]],
      enzymes[[cases$enzyme[[i]]]]
    ))
    expect_equal(expected %in% actual, TRUE, info = paste("reaction", i))
  }
})

test_that("ManII permits both trimming orders with concrete context restrictions", {
  cases <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  cases <- cases[cases$enzyme == "ManII", ]
  enzyme <- abstract_enzymes()$ManII
  for (substrate in unique(cases$substrate)) {
    expect_setequal(
      as.character(apply_enzyme(substrate, enzyme)),
      as.character(glyparse::auto_parse(cases$product[
        cases$substrate == substrate
      ]))
    )
  }
  substrates <- c(
    sub("GlcNAc(b1-2)", "", cases$substrate[[1]], fixed = TRUE),
    paste0("Gal(b1-4)", cases$substrate[[1]]),
    sub(
      "Man(b1-4)",
      "[GlcNAc(b1-4)]Man(b1-4)",
      cases$substrate[[1]],
      fixed = TRUE
    )
  )
  expect_identical(lengths(apply_enzyme(substrates, enzyme)), rep(0L, 3))
})

test_that("ManII can remove a terminal mannose beside an extended extra branch", {
  substrate <- paste0(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  product <- paste0(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-2)Man(a1-3)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  enzyme <- abstract_enzymes()$ManII
  expect_identical(
    as.character(apply_enzyme(substrate, enzyme)),
    as.character(glyparse::auto_parse(product))
  )
})

test_that("a6FucT accepts an extended prerequisite and cannot fucosylate twice", {
  cases <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  case <- cases[cases$enzyme == "a6FucT", ]
  enzyme <- abstract_enzymes()$a6FucT
  substrate <- paste0("Gal(b1-4)", case$substrate)
  product <- as.character(glyparse::auto_parse(paste0(
    "Gal(b1-4)",
    case$product
  )))
  expect_identical(as.character(apply_enzyme(substrate, enzyme)), product)
  expect_length(apply_enzyme(product, enzyme), 0L)
  missing <- sub("GlcNAc(b1-2)", "", case$substrate, fixed = TRUE)
  expect_length(apply_enzyme(missing, enzyme), 0L)
})

test_that("iGnT blocks the first alpha3 extension and permits repeated alpha6 extension", {
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  enzyme <- abstract_enzymes()$iGnT
  alpha3 <- "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)"
  for (repeats in c(1L, 3L)) {
    arm <- paste0(
      paste(rep("Gal(b1-4)GlcNAc(b1-3)", repeats - 1L), collapse = ""),
      "Gal(b1-4)GlcNAc(b1-2)"
    )
    substrate <- paste0(alpha3, "[", arm, "Man(a1-6)]", core)
    expected <- paste0(alpha3, "[GlcNAc(b1-3)", arm, "Man(a1-6)]", core)
    expect_identical(
      as.character(apply_enzyme(substrate, enzyme)),
      as.character(glyparse::auto_parse(expected))
    )
  }
})

test_that("iGnT rejects both alpha3 antennae without disabling the alpha6 arm", {
  alpha3 <- "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)"
  alpha6 <- "Gal(b1-4)GlcNAc(b1-2)Man(a1-6)"
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  substrate <- paste0(alpha3, "[", alpha6, "]", core)
  product <- paste0(alpha3, "[GlcNAc(b1-3)", alpha6, "]", core)
  expect_identical(
    as.character(apply_enzyme(substrate, abstract_enzymes()$iGnT)),
    as.character(glyparse::auto_parse(product))
  )
})

test_that("iGnT can extend an externally pre-extended alpha3 arm", {
  alpha3 <- "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)"
  alpha6 <- "Gal(b1-4)GlcNAc(b1-2)Man(a1-6)"
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  substrate <- paste0(alpha3, "[", alpha6, "]", core)
  products <- c(
    paste0("GlcNAc(b1-3)", alpha3, "[", alpha6, "]", core),
    paste0(alpha3, "[GlcNAc(b1-3)", alpha6, "]", core)
  )
  expect_setequal(
    as.character(apply_enzyme(substrate, abstract_enzymes()$iGnT)),
    as.character(glyparse::auto_parse(products))
  )
})

test_that("iGnT product statistics retain ordinary motif semantics", {
  arm <- "Gal(b1-4)GlcNAc(b1-2)"
  extended <- paste0("GlcNAc(b1-3)", arm)
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  glycans <- glyparse::auto_parse(c(
    paste0(extended, "Man(a1-3)[", arm, "Man(a1-6)]", core),
    paste0(arm, "Man(a1-3)[", extended, "Man(a1-6)]", core),
    paste0(extended, "Man(a1-3)[", extended, "Man(a1-6)]", core)
  ))
  enzyme <- abstract_enzymes()$iGnT
  expect_identical(have_enzyme(glycans, enzyme), c(TRUE, TRUE, TRUE))
  expect_equal(count_enzyme(glycans, enzyme), c(1L, 1L, 2L))
  expect_identical(lengths(match_enzyme(glycans, enzyme)), c(1L, 1L, 2L))
})

test_that("abstract rejects inherit lenient motif matching for unknown linkages", {
  substrate <- glyparse::auto_parse(paste0(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-?)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  enzyme <- abstract_enzymes()$iGnT
  ordinary <- enzyme
  class(ordinary) <- class(enzyme)[-1]
  expect_identical(
    as.character(suppressWarnings(apply_enzyme(substrate, enzyme))),
    as.character(suppressWarnings(apply_enzyme(substrate, ordinary)))
  )
  graph <- glyrepr::get_structure_graphs(substrate)
  rule <- enzyme$rules[[1]]
  expect_equal(
    .match_rule_graph(
      graph,
      rule,
      .prepare_rule_graphs(rule),
      mode = "lenient"
    ),
    .match_rule(substrate, rule)[[1]]
  )

  topological <- glyrepr::remove_linkages(substrate)
  rule_matches <- .match_rule(topological, rule)[[1]]
  expect_equal(
    .match_rule_graph(
      glyrepr::get_structure_graphs(topological),
      rule,
      .prepare_rule_graphs(rule),
      mode = "lenient"
    ),
    rule_matches
  )
})

test_that("abstract enzymes reject known non-N glycans and occupied sites", {
  enzymes <- abstract_enzymes()
  non_n <- c("GlcNAc(b1-6)GalNAc(a1-", "Gal(b1-4)GlcNAc(b1-6)GalNAc(a1-")
  for (enzyme in enzymes) {
    expect_identical(lengths(apply_enzyme(non_n, enzyme)), c(0L, 0L))
  }
  cases <- utils::read.csv(test_path("fixtures", "abstract-n-glycans.csv"))
  capped <- cases$glycan[cases$name == "G2S2F"]
  for (name in c("b4GalT", "iGnT", "a3SiaT")) {
    expect_length(apply_enzyme(capped, enzymes[[name]]), 0L)
  }
})

test_that("abstract subclasses retain standard actions without trusting custom ones", {
  enzymes <- abstract_enzymes()
  for (enzyme in enzymes) {
    expect_identical(.uses_standard_graph_action(enzyme), TRUE)
    expect_identical(.can_batch_bfs_enzyme(enzyme), TRUE)
    class(enzyme) <- c("custom_action", class(enzyme))
    expect_identical(.uses_standard_graph_action(enzyme), FALSE)
    expect_identical(.can_batch_bfs_enzyme(enzyme), FALSE)
  }
})

test_that("direct, prepared and batched actions agree for all abstract rules", {
  cases <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  enzymes <- abstract_enzymes()
  for (i in seq_len(nrow(cases))) {
    enzyme <- enzymes[[cases$enzyme[[i]]]]
    glycan <- glyparse::auto_parse(cases$substrate[[i]])
    graph <- glyrepr::get_structure_graphs(glycan)
    expected <- apply_enzyme(glycan, enzyme)
    prepared <- .apply_enzyme_prepared(
      graph,
      enzyme,
      lapply(enzyme$rules, .prepare_rule_graphs),
      "intact",
      "strict"
    )
    expect_setequal(as.character(prepared), as.character(expected))
    plan <- .prepare_bfs_rule_plan(list(enzyme))
    products <- .apply_bfs_rule_frontier(list(graph), "strict", plan, "intact")
    raw <- .collect_bfs_rule_products(products, plan$enzyme_rule_ids[[1]], 1L)
    actual <- .materialize_product_graphs(raw, enzyme)
    expect_setequal(as.character(actual), as.character(expected))
  }
})

test_that("rejects agree across direct, prepared and frontier actions", {
  reactions <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  manii <- reactions[reactions$enzyme == "ManII", ]
  core <- "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  alpha6 <- "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]"
  substrates <- list(
    ManII = c(
      unique(manii$substrate),
      paste0(
        "GlcNAc(b1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]",
        core
      )
    ),
    iGnT = c(
      paste0("Gal(b1-4)GlcNAc(b1-2)Man(a1-3)", alpha6, core),
      paste0(
        "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)",
        alpha6,
        core
      ),
      paste0(
        "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
        alpha6,
        core
      )
    )
  )
  expected_counts <- list(ManII = c(2L, 1L, 1L, 1L), iGnT = c(1L, 1L, 2L))
  enzymes <- abstract_enzymes()[names(substrates)]
  for (name in names(substrates)) {
    enzyme <- enzymes[[name]]
    glycans <- glyparse::auto_parse(substrates[[name]])
    graphs <- glyrepr::get_structure_graphs(glycans, return_list = TRUE)
    direct <- apply_enzyme(glycans, enzyme)
    expect_identical(lengths(direct), expected_counts[[name]])
    plan <- .prepare_bfs_rule_plan(list(enzyme))
    products <- .apply_bfs_rule_frontier(
      graphs,
      rep("strict", length(graphs)),
      plan,
      "intact"
    )
    for (i in seq_along(graphs)) {
      prepared <- .apply_enzyme_prepared(
        graphs[[i]],
        enzyme,
        lapply(enzyme$rules, .prepare_rule_graphs),
        "intact",
        "strict"
      )
      raw <- .collect_bfs_rule_products(products, plan$enzyme_rule_ids[[1]], i)
      batched <- .materialize_product_graphs(raw, enzyme)
      expect_setequal(as.character(prepared), as.character(direct[[i]]))
      expect_setequal(as.character(batched), as.character(direct[[i]]))
    }
  }
})

test_that("BFS never shares a rule with rejects with an unrestricted rule", {
  enzyme <- abstract_enzymes()$iGnT
  unrestricted <- enzyme
  unrestricted$name <- "unrestricted"
  unrestricted$rules[[1]]$rejects <- glyrepr::glycan_structure()
  plan <- .prepare_bfs_rule_plan(list(enzyme, unrestricted))
  expect_length(plan$rules, 2L)
  substrate <- glyparse::auto_parse(paste0(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
    "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
  products <- .apply_bfs_rule_frontier(
    list(glyrepr::get_structure_graphs(substrate)),
    "strict",
    plan,
    "intact"
  )
  expect_length(
    .collect_bfs_rule_products(products, plan$enzyme_rule_ids[[1]], 1L),
    1L
  )
  expect_length(
    .collect_bfs_rule_products(products, plan$enzyme_rule_ids[[2]], 1L),
    2L
  )
})

test_that("abstract-only tracing reconstructs typical N-glycans with valid edges", {
  cases <- utils::read.csv(test_path("fixtures", "abstract-n-glycans.csv"))
  reactions <- utils::read.csv(test_path("fixtures", "abstract-reactions.csv"))
  glucoses <- reactions[reactions$enzyme == "GlcH", ]
  enzymes <- abstract_enzymes()
  start <- as.character(glyparse::auto_parse(glucoses$substrate[[1]]))
  intermediates <- as.character(glyparse::auto_parse(glucoses$product))
  cache <- new.env(parent = emptyenv())

  expect_equal(nrow(cases), 7L)
  for (i in seq_len(nrow(cases))) {
    target <- as.character(glyparse::auto_parse(cases$glycan[[i]]))
    path <- trace_biosynthesis(
      target,
      enzymes = enzymes,
      max_virtual_steps = 0L
    )
    edges <- igraph::as_data_frame(path, what = "edges")
    vertices <- igraph::as_data_frame(path, what = "vertices")
    label <- cases$name[[i]]
    expect_s3_class(path, "glyenzy_biosynthesis_network")
    expect_identical(vertices$name[vertices$target], target, info = label)
    expect_identical(
      vertices$name[igraph::degree(path, mode = "in") == 0L],
      start,
      info = label
    )
    expect_equal(
      as.numeric(igraph::distances(path, start, target, mode = "out")),
      cases$steps[[i]],
      info = label
    )
    expect_equal(
      setdiff(intermediates, vertices$name),
      character(),
      info = label
    )
    expect_equal(igraph::any_multiple(path), FALSE, info = label)
    expect_identical(edges$is_virtual, rep(FALSE, nrow(edges)))
    expected_enzymes <- strsplit(cases$enzymes[[i]], "|", fixed = TRUE)[[1]]
    expect_equal(
      sort(unique(unlist(edges$enzymes))),
      sort(expected_enzymes),
      info = label
    )
    expect_length(setdiff(unlist(edges$enzymes), names(enzymes)), 0L)

    for (j in seq_len(nrow(edges))) {
      for (name in edges$enzymes[[j]]) {
        key <- paste(name, edges$from[[j]], sep = "\r")
        if (!exists(key, envir = cache, inherits = FALSE)) {
          assign(
            key,
            as.character(apply_enzyme(edges$from[[j]], enzymes[[name]])),
            envir = cache
          )
        }
        expect_equal(
          edges$to[[j]] %in% get(key, envir = cache),
          TRUE,
          info = paste(label, j, name)
        )
      }
    }
  }
})

test_that("abstract tracing retains precursor enzymes and rejects impossible routes", {
  cases <- utils::read.csv(test_path("fixtures", "abstract-n-glycans.csv"))
  enzymes <- abstract_enzymes()
  target <- glyparse::auto_parse(cases$glycan[cases$name == "G2S2F"])
  expect_identical(
    unname(.can_enzymes_contribute(enzymes, target)),
    rep(TRUE, length(enzymes))
  )
  expect_snapshot(
    error = TRUE,
    trace_biosynthesis(
      cases$glycan[cases$name == "Man5"],
      enzymes = enzymes[names(enzymes) != "GlcH"],
      max_virtual_steps = 0L
    )
  )
  wrong_arm <- paste0(
    "Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
    "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_snapshot(
    error = TRUE,
    trace_biosynthesis(wrong_arm, enzymes = enzymes, max_virtual_steps = 0L)
  )
})
