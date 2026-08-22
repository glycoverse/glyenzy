test_that("glycan types are classified from reducing-end structures", {
  glycans <- glyparse::auto_parse(c(
    n = "GlcNAc(b1-4)GlcNAc(?1-",
    o_galnac = "Gal(b1-3)GalNAc(a1-",
    o_man = "GlcNAc(b1-2)Man(a1-",
    o_fuc = "GlcNAc(b1-3)Fuc(a1-",
    o_glcnac = "Gal(b1-3)GlcNAc(b1-",
    o_xyl = "Gal(b1-4)Xyl(a1-",
    o_glc = "Xyl(??-?)Glc(?1-",
    o_gal = "Glc(a1-2)Gal(?1-",
    lipid_glc = "Gal(b1-4)Glc(b1-",
    lipid_gal = "Glc(b1-4)Gal(b1-",
    incomplete_n = "Man(a1-3)Man(b1-",
    unknown = "Neu5Ac(a2-3)Neu5Ac(?1-"
  ))

  expect_identical(
    .glycan_types(glycans),
    c(
      n = "N",
      o_galnac = "O",
      o_man = "O",
      o_fuc = "O",
      o_glcnac = "O",
      o_xyl = "O",
      o_glc = "O",
      o_gal = "O",
      lipid_glc = "lipid/free",
      lipid_gal = "lipid/free",
      incomplete_n = NA_character_,
      unknown = NA_character_
    )
  )
})

test_that("beta-Man reducing ends retain incomplete N-glycan context", {
  fragments <- glyparse::auto_parse(c(
    mgat1 = "Man(a1-3)[Man(a1-6)]Man(a1-6)[Man(a1-3)]Man(b1-",
    mgat2 = "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-"
  ))

  expect_identical(
    .glycan_types(fragments),
    c(mgat1 = NA_character_, mgat2 = NA_character_)
  )
  expect_true(all(.enzyme_glycan_type_mask(fragments, enzyme("MGAT1"))))
  expect_true(all(.enzyme_glycan_type_mask(fragments, enzyme("MGAT2"))))

  o_mannose <- glyparse::auto_parse("Man(a1-3)Man(a1-")
  expect_identical(.glycan_types(o_mannose), "O")
  expect_false(.enzyme_glycan_type_mask(o_mannose, enzyme("MGAT1")))
})

test_that("enzyme glycan types are validated and normalized", {
  enz <- make_enzyme(
    name = "TYPED_GT",
    rules = list(list(
      acceptor = "Glc(b1-",
      acceptor_alignment = "whole",
      rejects = NULL,
      product = "Gal(b1-3)Glc(b1-"
    )),
    type = "GT",
    species = "human",
    glycan_type = c("free", "lipid", "N", "lipid")
  )

  expect_identical(enz$glycan_type, c("N", "lipid", "free"))
  expect_null(
    make_enzyme(
      name = "ALL_GT",
      rules = list(list(
        acceptor = "Glc(b1-",
        acceptor_alignment = "whole",
        rejects = NULL,
        product = "Gal(b1-3)Glc(b1-"
      )),
      type = "GT",
      species = "human"
    )$glycan_type
  )
  expect_error(
    make_enzyme(
      name = "BAD_GT",
      rules = list(list(
        acceptor = "Glc(b1-",
        acceptor_alignment = "whole",
        rejects = NULL,
        product = "Gal(b1-3)Glc(b1-"
      )),
      type = "GT",
      species = "human",
      glycan_type = "protein"
    ),
    "additional elements"
  )
})

test_that("built-in enzyme metadata comes from broad pathway classes", {
  expect_identical(enzyme("MGAT5B")$glycan_type, c("N", "O"))
  expect_identical(enzyme("ST3GAL3")$glycan_type, c("lipid", "free"))
  expect_identical(enzyme("ST3GAL5")$glycan_type, "lipid")
  expect_null(enzyme("FUT1")$glycan_type)
  expect_identical(enzyme("DPAGT1")$glycan_type, "N")
  expect_true(all(vapply(
    db_enzymes(),
    function(x) {
      is.null(x$glycan_type) ||
        all(x$glycan_type %in% c("N", "O", "lipid", "free"))
    },
    logical(1)
  )))
  expect_equal(
    sum(vapply(db_enzymes(), function(x) !is.null(x$glycan_type), logical(1))),
    132L
  )
  expect_equal(
    sum(vapply(db_enzymes(), function(x) is.null(x$glycan_type), logical(1))),
    31L
  )
})

test_that("typed enzymes handle known and ambiguous glycan classes", {
  n_glycan <- suppressWarnings(glyparse::auto_parse(
    "Gal(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  ))
  o_glycan <- suppressWarnings(glyparse::auto_parse(
    "Gal(b1-4)GlcNAc(b1-"
  ))
  lipid_or_free <- suppressWarnings(glyparse::auto_parse("Gal(b1-4)Glc(b1-"))

  expect_true(suppressWarnings(have_enzyme(n_glycan, "B4GALT1")))
  expect_true(suppressWarnings(have_enzyme(o_glycan, "B4GALT1")))
  expect_true(suppressWarnings(have_enzyme(lipid_or_free, "B4GALT1")))
  expect_equal(suppressWarnings(count_enzyme(lipid_or_free, "B4GALT1")), 1L)
  expect_equal(
    suppressWarnings(match_enzyme(lipid_or_free, "B4GALT1")),
    list(1L)
  )
  expect_true(suppressWarnings("B4GALT1" %in% find_enzyme(lipid_or_free)))
  expect_true(.enzyme_glycan_type_mask(lipid_or_free, enzyme("ST3GAL5")))
  expect_false(.enzyme_glycan_type_mask(lipid_or_free, enzyme("DPAGT1")))

  free_glc <- suppressWarnings(glyparse::auto_parse("Glc(b1-"))
  grown <- suppressWarnings(grow_glycans_step(free_glc, "B4GALT1"))
  expect_length(grown, 1L)
  expect_false(.enzyme_glycan_type_mask(n_glycan, enzyme("ST3GAL5")))
  expect_false(.enzyme_glycan_type_mask(o_glycan, enzyme("DPAGT1")))
  expect_true(.enzyme_glycan_type_mask(
    glyparse::auto_parse("Neu5Ac(a2-3)Neu5Ac(?1-"),
    enzyme("B4GALT1")
  ))
})

test_that("free-compatible enzymes work in known and virtual biosynthesis", {
  from <- "Glc(b1-"
  to <- "Gal(b1-4)Glc(b1-"

  b4galt5 <- path_biosynthesis(from, to, enzymes = "B4GALT5", max_steps = 1)
  expect_identical(igraph::E(b4galt5)$enzyme, "B4GALT5")
  known <- path_biosynthesis(from, to, enzymes = "B4GALT6", max_steps = 1)
  expect_identical(igraph::E(known)$enzyme, "B4GALT6")

  virtual <- trace_biosynthesis_virtual(
    to,
    enzymes = c("B4GALT5", "B4GALT6"),
    annotate_enzymes = TRUE
  )
  expect_identical(
    igraph::E(virtual)$concrete_enzymes[[1]],
    c("B4GALT5", "B4GALT6")
  )
})
