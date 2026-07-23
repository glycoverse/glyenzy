test_that("make_enzyme works with GT", {
  enz <- make_enzyme(
    name = "MySiaT",
    rules = list(
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = NULL,
        product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      )
    ),
    type = "GT",
    species = "human"
  )
  expect_snapshot(enz)
})

test_that("make_enzyme works with GT with rejects", {
  enz <- make_enzyme(
    name = "MySiaT",
    rules = list(
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = c(
          # Reject a1-3 Fuc on GlcNAc
          "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-",
          # Reject a1-2 Fuc on Gal
          "Fuc(a1-2)Gal(b1-4)GlcNAc(b1-"
        ),
        product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      )
    ),
    type = "GT",
    species = "human"
  )
  expect_snapshot(enz)
})

test_that("make_enzyme works with GH", {
  enz <- make_enzyme(
    name = "MyGalH",
    rules = list(
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = NULL,
        product = "GlcNAc(b1-"
      )
    ),
    type = "GH",
    species = "human"
  )
  expect_snapshot(enz)
})

test_that("make_enzyme works with ST requirements", {
  enz <- make_enzyme(
    name = "MySulfoT",
    rules = list(
      list(
        acceptor = "Gal(b1-",
        acceptor_alignment = "substructure",
        rejects = NULL,
        requires = list(
          list(motif = "GalNAc(a1-", alignment = "core"),
          list(
            motif = "GlcNAc(b1-4)GlcNAc(b1-",
            alignment = "core"
          )
        ),
        product = "Gal6S(b1-"
      )
    ),
    type = "ST",
    species = "human"
  )

  expect_s3_class(enz, "glyenzy_st_enzyme")
  expect_length(enz$rules[[1]]$requires, 2L)
  messages <- testthat::capture_messages(print(enz))
  expect_match(paste(messages, collapse = "\n"), "Requires \\(any\\):")
})

test_that("make_enzyme rejects duplicate requirement keys", {
  requirement <- list(
    motif = "GalNAc(a1-",
    motif = "GlcNAc(b1-"
  )

  expect_error(
    make_enzyme(
      name = "MySulfoT",
      rules = list(list(
        acceptor = "Gal(b1-",
        acceptor_alignment = "whole",
        rejects = NULL,
        requires = list(requirement),
        product = "Gal6S(b1-"
      )),
      type = "ST",
      species = "human"
    ),
    "must contain only"
  )
})

test_that("make_enzyme works with GT with multiple rules", {
  enz <- make_enzyme(
    name = "MySiaT",
    rules = list(
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = NULL,
        product = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      ),
      list(
        acceptor = "Gal(b1-4)GlcNAc(b1-",
        acceptor_alignment = "terminal",
        rejects = NULL,
        product = "Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-"
      )
    ),
    type = "GT",
    species = "human"
  )
  expect_snapshot(enz)
})
