# ===== Input type =====
test_that("is_synthesized_by works for `glyrepr_structure` and `glyenzy_enzyme`", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
  enzyme <- enzyme("ST3GAL3")
  expect_true(is_synthesized_by(glycan, enzyme))
})

test_that("is_synthesized_by works for characters", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  enzyme <- "ST3GAL3"
  expect_true(is_synthesized_by(glycan, enzyme))
})

test_that("is_synthesized_by works vectorizedly", {
  glycans <- c(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    "Gal(b1-3)GlcNAc(b1-"
  )
  expect_equal(is_synthesized_by(glycans, "ST3GAL3"), c(TRUE, FALSE))
})

test_that("is_synthesized_by rejects invalid inputs", {
  expect_error(is_synthesized_by(123, "ST3GAL3"), "`glycans` must be")
  expect_error(is_synthesized_by("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-", 123), "`enzyme` must be")
})

# ===== Special cases for N-glycans =====
test_that("is_synthesized_by works for ALG enzymes and DPAGT1", {
  glycan1 <- "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"  # the basic N-glycan core
  glycan2 <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"  # not an N-glycan
  enzymes <- c("ALG1", "ALG2", "ALG3", "DPAGT1")
  expect_true(all(purrr::map_lgl(enzymes, ~ is_synthesized_by(glycan1, .x))))
  expect_false(any(purrr::map_lgl(enzymes, ~ is_synthesized_by(glycan2, .x))))
})

test_that("is_synthesized_by works for MGAT1", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(is_synthesized_by(glycans, "MGAT1"), c(TRUE, FALSE))
})

test_that("is_synthesized_by works for MOGS", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(is_synthesized_by(glycans, "MOGS"), c(FALSE, TRUE, TRUE))
})

test_that("is_synthesized_by works for GANAB", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(is_synthesized_by(glycans, "GANAB"), c(FALSE, FALSE, TRUE, TRUE))
})

test_that("is_synthesized_by works for MAN2A1 and MAN2A2", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(is_synthesized_by(glycans, "MAN2A1"), c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(is_synthesized_by(glycans, "MAN2A2"), c(FALSE, TRUE, TRUE, TRUE))
})

test_that("is_synthesized_by works correctly for MAN1B1, MAN1A1, MAN1A2, and MAN1C1", {
  # Please check Fig. 115.2 of Handbook of Glycosyltransferases and Related Genes for details.
  # For each composition, structures are assigned from top to bottom.
  glycans <- c(
    # Man(9)GlcNAc(2)
    Man9 = "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # Man(8)GlcNAc(2)
    Man8_1 = "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    Man8_2 = "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    Man8_3 = "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # Man(7)GlcNAc(2)
    Man7_1 = "Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    Man7_2 = "Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    Man7_3 = "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    Man7_4 = "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # Man(6)GlcNAc(2)
    Man6_1 = "Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    Man6_2 = "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    Man6_3 = "Man(a1-3)[Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    # Man(5)GlcNAc(2)
    Man5 = "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )

  MAN1B1_res <- c(FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE)
  MAN1A1_res <- c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  # MAN1A2 and MAN1C1 are the same as MAN1A1
  expect_equal(is_synthesized_by(glycans, "MAN1B1"), MAN1B1_res)
  expect_equal(is_synthesized_by(glycans, "MAN1A1"), MAN1A1_res)
  expect_equal(is_synthesized_by(glycans, "MAN1A2"), MAN1A1_res)
  expect_equal(is_synthesized_by(glycans, "MAN1C1"), MAN1A1_res)
})