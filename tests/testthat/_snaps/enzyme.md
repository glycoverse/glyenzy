# print.glyenzy_enzyme works correctly

    Code
      print(enzyme_obj)
    Message
      
      -- Enzyme: ST3GAL3 -------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (3) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-3)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
      > Rule 2: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      > Rule 3: terminal alignment
      Acceptor: "Gal(b1-3)GalNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-3)GalNAc(b1-"

# DPAGT1 enzyme works correctly

    Code
      print(dpagt1)
    Message
      
      -- Enzyme: DPAGT1 --------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (1) --
      
      > Rule 1: de novo synthesis
      Acceptor: none
      Product: "GlcNAc(b1-"

