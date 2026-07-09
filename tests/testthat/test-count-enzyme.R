# ===== Input Type =====
test_that("count_enzyme works for `glyrepr_structure` and `glyenzy_enzyme`", {
  glycans <- glyparse::parse_iupac_condensed("Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-")
  enzyme <- enzyme("ST3GAL3")
  expect_equal(count_enzyme(glycans, enzyme), 1L)
})

test_that("count_enzyme works for characters", {
  glycan <- "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
  enzyme <- "ST3GAL3"
  expect_equal(count_enzyme(glycan, enzyme), 1L)
})

test_that("count_enzyme works vectorizedly", {
  glycans <- c(
    "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-",
    "Gal(b1-3)GalNAc(a1-"
  )
  expect_equal(count_enzyme(glycans, "ST3GAL3"), c(1L, 0L))
})

test_that("count_enzyme uses lenient matching for non-intact glycans", {
  glycan <- "Gal(b1-?)GalNAc(a1-"

  expect_warning(
    expect_equal(count_enzyme(glycan, "C1GALT1"), 1L),
    "non-intact glycan structures"
  )
})

test_that("count_enzyme uses lenient matching for non-intact special cases", {
  mgat1_glycan <- "GlcNAc(b1-2)Man(?1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  mogs_glycan <- "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(?1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  man2_glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(?1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

  expect_warning(
    expect_equal(count_enzyme(mgat1_glycan, "MGAT1"), 1L),
    "non-intact glycan structures"
  )
  expect_warning(
    expect_equal(count_enzyme(mogs_glycan, "MOGS"), 1L),
    "non-intact glycan structures"
  )
  expect_warning(
    expect_equal(count_enzyme(man2_glycan, "MAN2A1"), 2L),
    "non-intact glycan structures"
  )
})

# ===== Normal Cases =====
test_that("count_enzyme works for FUT8", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-" # O-glycan
  )
  expect_equal(count_enzyme(glycans, "FUT8"), c(0L, 1L, 0L))
})

test_that("count_enzyme works for B4GALT1", {
  glycans <- c(
    "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-" # O-glycan
  )
  expect_equal(count_enzyme(glycans, "B4GALT1"), c(2L, 1L, 0L, 0L))
})

# ===== Special cases for N-glycans =====
test_that("count_enzyme works for MGAT1", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",
    "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-" # O-glycan
  )
  expect_equal(count_enzyme(glycans, "MGAT1"), c(1L, 1L, 0L, 0L))
})

test_that("count_enzyme works for MOGS", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-" # O-glycan
  )
  expect_equal(count_enzyme(glycans, "MOGS"), c(0L, 1L, 1L, 0L))
})

test_that("count_enzyme works for MAN2A1 and MAN2A2", {
  glycans <- c(
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-3)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-" # O-glycan
  )
  expect_equal(count_enzyme(glycans, "MAN2A1"), c(0L, 1L, 1L, 2L, 0L))
  expect_equal(count_enzyme(glycans, "MAN2A2"), c(0L, 1L, 1L, 2L, 0L))
})

test_that("count_enzyme works for GANAB", {
  glycans <- c(
    "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  )
  expect_equal(count_enzyme(glycans, "GANAB"), c(0L, 0L, 1L, 2L))
})

test_that("count_enzyme works for MAN1A1, MAN1A2, and MAN1C1", {
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
  expect_equal(count_enzyme(glycans, "MAN1B1"), MAN1B1_res)
  expect_equal(count_enzyme(glycans, "MAN1A1"), MAN1A1_res)
  expect_equal(count_enzyme(glycans, "MAN1A2"), MAN1A1_res)
  expect_equal(count_enzyme(glycans, "MAN1C1"), MAN1C1_res)
})

# ===== Starter cases =====
test_that("count_enzyme works for DPAGT1", {
  expect_equal(count_enzyme(glyrepr::n_glycan_core(), "DPAGT1"), 1L)
  expect_equal(count_enzyme(glyrepr::o_glycan_core_1(), "DPAGT1"), 0L)
})

test_that("count_enzyme works for FUT10", {
  expect_equal(count_enzyme("Gal(b1-4)GlcNAc(b1-3)Fuc(a1-", "FUT10"), 1L)
  expect_equal(count_enzyme(glyrepr::o_glycan_core_1(), "FUT10"), 0L)
})

test_that("count_enzyme works for POGLUT2", {
  expect_equal(count_enzyme("Xyl(a1-3)Glc(a1-", "POGLUT2"), 1L)
  expect_equal(count_enzyme(glyrepr::o_glycan_core_1(), "POGLUT2"), 0L)
})

# ===== N-glycan precursor cases =====
test_that("count_enzyme works for ALG2, ALG11, ALG9, and ALG6", {
  glycan <- glyrepr::n_glycan_core()
  expect_equal(count_enzyme(glycan, "ALG2"), 2L)
  expect_equal(count_enzyme(glycan, "ALG11"), 2L)
  expect_equal(count_enzyme(glycan, "ALG9"), 2L)
  expect_equal(count_enzyme(glycan, "ALG6"), 1L)
})

test_that("count_enzyme is vectorized for N-glycan precursor GTs", {
  glycans <- c(glyrepr::n_glycan_core(), glyrepr::o_glycan_core_1())

  expect_equal(count_enzyme(glycans, "ALG2"), c(2L, 0L))
  expect_equal(count_enzyme(glycans, "ALG6"), c(1L, 0L))
})

test_that("count_enzyme path method treats N-glycan precursor GTs specially", {
  glycans <- c(glyrepr::n_glycan_core(), glyrepr::o_glycan_core_1())

  expect_equal(
    suppressMessages(count_enzyme(glycans, "ALG2", method = "path")),
    c(2L, 0L)
  )
  expect_equal(
    suppressMessages(count_enzyme(glycans, "ALG6", method = "path")),
    c(1L, 0L)
  )
})

test_that("count_enzyme path method supports custom enzyme objects", {
  enz <- make_enzyme(
    name = "TEST_ST3GAL",
    type = "GT",
    species = "human",
    rules = list(list(
      acceptor = "Gal(b1-3)GalNAc(a1-",
      acceptor_alignment = "core",
      rejects = NULL,
      product = "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
    ))
  )

  expect_equal(
    suppressMessages(
      count_enzyme("Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-", enz, method = "path")
    ),
    1L
  )
})

test_that("count_enzyme can use trace-derived path enzymes", {
  glycans <- c(
    sialyl_lewis_x = "Neu5Ac(a2-3)Gal(b1-3)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]GalNAc(a1-",
    core1 = "Gal(b1-3)GalNAc(a1-"
  )

  expect_equal(count_enzyme(glycans, "FUT3"), c(1L, 0L))
  expect_equal(
    suppressMessages(count_enzyme(glycans, "FUT3", method = "path")),
    c(0L, 0L)
  )
  expect_equal(
    suppressMessages(count_enzyme(glycans, "FUT7", method = "path")),
    c(2L, 0L)
  )
})
