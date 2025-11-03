test_that("create_enzyme works with GT", {
  enz <- create_enzyme(
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

test_that("create_enzyme works with GT with rejects", {
  enz <- create_enzyme(
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

test_that("create_enzyme works with GH", {
  enz <- create_enzyme(
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

test_that("create_enzyme works with GT with multiple rules", {
  enz <- create_enzyme(
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