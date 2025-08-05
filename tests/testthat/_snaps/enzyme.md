# print.glyenzy_enzyme works correctly

    Code
      print(enzyme)
    Message
      
      -- Enzyme: ST3GAL2 -------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (1) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-3)GalNAc(a1-"
      Product: "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
      
      -- Rejects (1) --
      
      1: "Gal(b1-3)GalNAc(a1-" (terminal)

# print.glyenzy_enzyme handles empty motif sets

    Code
      print(enzyme)
    Message
      
      -- Enzyme: TEST ----------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (1) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-3)GalNAc(a1-"
      Product: "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
      
      -- Rejects (0) --
      
      No reject motifs defined

# print.glyenzy_enzyme handles multiple rules

    Code
      print(enzyme)
    Message
      
      -- Enzyme: ST3GAL1 -------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (2) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-3)GalNAc(a1-"
      Product: "Neu5Ac(a2-3)Gal(b1-3)GalNAc(a1-"
      > Rule 2: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      
      -- Rejects (0) --
      
      No reject motifs defined

