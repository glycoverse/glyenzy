# ===== Input Type =====
test_that("count_enzyme_steps works for `glyrepr_structure` and `glyenzy_enzyme`", {
  glycans <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
  enzyme <- enzyme("ST3GAL3")
  expect_equal(count_enzyme_steps(glycans, enzyme), 1L)
})

test_that("count_enzyme_steps works for characters", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  enzyme <- "ST3GAL3"
  expect_equal(count_enzyme_steps(glycan, enzyme), 1L)
})

test_that("count_enzyme_steps works vectorizedly", {
  glycans <- c(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    "Gal(b1-3)GalNAc(a1-"
  )
  expect_equal(count_enzyme_steps(glycans, "ST3GAL3"), c(1L, 0L))
})

# ===== Normal Cases =====
test_that("count_enzyme_steps works for FUT8", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"  # O-glycan
  )
  expect_equal(count_enzyme_steps(glycans, "FUT8"), c(0L, 1L, 0L))
})

test_that("count_enzyme_steps works for B4GALT1", {
  glycans <- c(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"  # O-glycan
  )
  expect_equal(count_enzyme_steps(glycans, "B4GALT1"), c(2L, 1L, 0L, 0L))
})

# ===== Special cases for N-glycans =====
test_that("count_enzyme_steps works for MGAT1", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",
    "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"  # O-glycan
  )
  expect_equal(count_enzyme_steps(glycans, "MGAT1"), c(1L, 1L, 0L, 0L))
})

test_that("count_enzyme_steps works for MOGS", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"  # O-glycan
  )
  expect_equal(count_enzyme_steps(glycans, "MOGS"), c(0L, 1L, 1L, 0L))
})

test_that("count_enzyme_steps works for MAN2A1 and MAN2A2", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"  # O-glycan
  )
  expect_equal(count_enzyme_steps(glycans, "MAN2A1"), c(0L, 1L, 1L, 2L, 0L))
  expect_equal(count_enzyme_steps(glycans, "MAN2A2"), c(0L, 1L, 1L, 2L, 0L))
})

test_that("count_enzyme_steps works for GANAB", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(count_enzyme_steps(glycans, "GANAB"), c(0L, 0L, 1L, 2L))
})

test_that("count_enzyme_steps works for MAN1A1, MAN1A2, and MAN1C1", {
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

  MAN1B1_res <- c(0L, 0L, 0L, 1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L)
  MAN1A1_res <- c(0L, 1L, 1L, 0L, 1L, 2L, 1L, 1L, 3L, 2L, 2L, 3L)
  MAN1C1_res <- c(0L, 1L, 1L, 0L, 2L, 2L, 1L, 1L, 3L, 2L, 2L, 3L)
  expect_equal(count_enzyme_steps(glycans, "MAN1B1"), MAN1B1_res)
  expect_equal(count_enzyme_steps(glycans, "MAN1A1"), MAN1A1_res)
  expect_equal(count_enzyme_steps(glycans, "MAN1A2"), MAN1A1_res)
  expect_equal(count_enzyme_steps(glycans, "MAN1C1"), MAN1C1_res)
})