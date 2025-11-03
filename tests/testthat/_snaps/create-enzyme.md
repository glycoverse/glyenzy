# create_enzyme works with GT

    Code
      enz
    Message
      
      -- Enzyme: MySiaT --------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (1) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"

# create_enzyme works with GT with rejects

    Code
      enz
    Message
      
      -- Enzyme: MySiaT --------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (1) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      Rejects:
      "Fuc(a1-3)[Gal(b1-4)]GlcNAc(b1-"
      "Fuc(a1-2)Gal(b1-4)GlcNAc(b1-"

# create_enzyme works with GH

    Code
      enz
    Message
      
      -- Enzyme: MyGalH --------------------------------------------------------------
      i Type: "GH" (Glycoside hydrolase)
      i Species: "human"
      
      -- Rules (1) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "GlcNAc(b1-"

# create_enzyme works with GT with multiple rules

    Code
      enz
    Message
      
      -- Enzyme: MySiaT --------------------------------------------------------------
      i Type: "GT" (Glycosyltransferase)
      i Species: "human"
      
      -- Rules (2) --
      
      > Rule 1: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-"
      > Rule 2: terminal alignment
      Acceptor: "Gal(b1-4)GlcNAc(b1-"
      Product: "Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-"

