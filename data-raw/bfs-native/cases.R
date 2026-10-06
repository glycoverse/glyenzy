cases <- list(
  o_branched = list(
    "GalNAc(a1-",
    "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-6)[Neu5Ac(a2-3)Gal(b1-3)]GalNAc(a1-",
    c("C1GALT1", "GCNT1", "B4GALT1", "B4GALT2", "ST3GAL1", "ST3GAL4"),
    6L
  ),
  n_complex = list(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    c("MGAT1", "MGAT2", "B4GALT1", "B4GALT2", "ST6GAL1"),
    6L
  ),
  n_precursor = list(
    as.character(internal(".n_glycan_starting_glycan")("enzymatic")),
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    c(
      "MOGS",
      "GANAB",
      "MAN1A1",
      "MAN1A2",
      "MAN1C1",
      "MGAT1",
      "MAN2A1",
      "MAN2A2",
      "MGAT2"
    ),
    12L
  )
)
cases$n_multi <- list(
  cases$n_complex[[1]],
  c(
    cases$n_complex[[2]],
    "Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)[Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",
    "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)][GlcNAc(b1-4)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ),
  c(
    "MGAT1",
    "MGAT2",
    "MGAT3",
    "MGAT4A",
    "MGAT5",
    "B4GALT1",
    "B4GALT2",
    "FUT8",
    "FUT3",
    "ST6GAL1",
    "ST3GAL4"
  ),
  10L
)
