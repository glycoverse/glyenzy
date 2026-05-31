# ===== Input type =====
test_that("have_enzyme works for `glyrepr_structure` and `glyenzy_enzyme`", {
  glycan <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
  enzyme <- enzyme("ST3GAL3")
  expect_true(have_enzyme(glycan, enzyme))
})

test_that("have_enzyme works for characters", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  enzyme <- "ST3GAL3"
  expect_true(have_enzyme(glycan, enzyme))
})

test_that("have_enzyme works vectorizedly", {
  glycans <- c(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    "Gal(b1-3)GlcNAc(b1-"
  )
  expect_equal(have_enzyme(glycans, "ST3GAL3"), c(TRUE, FALSE))
})

test_that("have_enzyme rejects invalid inputs", {
  expect_error(have_enzyme(123, "ST3GAL3"), "`glycans` must be")
  expect_error(
    have_enzyme("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-", 123),
    "`enzyme` must be"
  )
})

# ===== Special cases for N-glycans =====
test_that("have_enzyme works for MGAT1", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(have_enzyme(glycans, "MGAT1"), c(TRUE, FALSE))
})

test_that("have_enzyme works for MOGS", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(have_enzyme(glycans, "MOGS"), c(FALSE, TRUE, TRUE))
})

test_that("have_enzyme works for GANAB", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(have_enzyme(glycans, "GANAB"), c(FALSE, FALSE, TRUE, TRUE))
})

test_that("have_enzyme works for MAN2A1 and MAN2A2", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(have_enzyme(glycans, "MAN2A1"), c(FALSE, TRUE, TRUE, TRUE))
  expect_equal(have_enzyme(glycans, "MAN2A2"), c(FALSE, TRUE, TRUE, TRUE))
})

test_that("have_enzyme works correctly for MAN1B1, MAN1A1, MAN1A2, and MAN1C1", {
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

  MAN1B1_res <- c(
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    FALSE,
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE
  )
  MAN1A1_res <- c(
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE
  )
  # MAN1A2 and MAN1C1 are the same as MAN1A1
  expect_equal(have_enzyme(glycans, "MAN1B1"), MAN1B1_res)
  expect_equal(have_enzyme(glycans, "MAN1A1"), MAN1A1_res)
  expect_equal(have_enzyme(glycans, "MAN1A2"), MAN1A1_res)
  expect_equal(have_enzyme(glycans, "MAN1C1"), MAN1A1_res)
})

# ===== Edge Cases =====
test_that("have_enzyme handles product alignments", {
  expect_true(have_enzyme("GlcNAc(b1-3)GalNAc(a1-", "B3GNT6"))
  expect_false(have_enzyme("GlcNAc(b1-3)GalNAc(a1-3)GlcNAc(b1-", "B3GNT6"))
})

# ===== Starter Cases =====
test_that("have_enzyme works for DPAGT1", {
  expect_true(have_enzyme(glyrepr::n_glycan_core(), "DPAGT1"))
  expect_false(have_enzyme(glyrepr::o_glycan_core_1(), "DPAGT1"))
})

test_that("have_enzyme works for FUT10", {
  expect_true(have_enzyme("Gal(b1-4)GlcNAc(b1-3)Fuc(a1-", "FUT10"))
  expect_false(have_enzyme(glyrepr::o_glycan_core_1(), "FUT10"))
})

test_that("have_enzyme works for POGLUT2", {
  expect_true(have_enzyme("Xyl(a1-3)Glc(a1-", "POGLUT2"))
  expect_false(have_enzyme(glyrepr::o_glycan_core_1(), "POGLUT2"))
})
