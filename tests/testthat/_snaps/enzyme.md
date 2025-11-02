# print.glyenzy_enzyme works correctly

    Code
      print(enzyme_obj)
    Message
      
      -- Enzyme: ST3GAL3 -------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (2) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-3)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-3)GlcNAc(b1-"
      > Rule 2: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"

# print.glyenzy_enzyme works with enzyme with no rules

    Code
      print(enzyme_obj)
    Message
      
      -- Enzyme: EMPTY_ENZYME --------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (0) --
      
      No rules defined

# print.glyenzy_enzyme works with enzyme with rejects

    Code
      print(enzyme_obj)
    Message
      
      -- Enzyme: TEST_ENZYME ---------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (1) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-3)GalNAc(a1-"
      Product: "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
      Rejects:
      "Gal(b1-4)GalNAc(a1-"
      "Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-"

