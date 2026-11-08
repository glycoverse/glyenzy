# normalization rejects invalid input types

    Code
      normalize_n_glycan(1)
    Condition
      Error in `normalize_n_glycan()`:
      ! `glycans` must be a glycan structure vector or a character vector.

# normalization rejects generic monosaccharides

    Code
      normalize_n_glycan("Hex(??-?)HexNAc(??-")
    Condition
      Error in `normalize_n_glycan()`:
      ! `glycans` must use concrete monosaccharides.

# normalization rejects intact structures

    Code
      normalize_n_glycan("Gal(b1-4)GlcNAc(b1-")
    Condition
      Error in `normalize_n_glycan()`:
      ! `glycans` must contain only topological structures.
      i Use `glyrepr::remove_linkages()` on intact or partial structures first.

# normalization rejects partial structures

    Code
      normalize_n_glycan("Gal(??-?)GlcNAc(b1-")
    Condition
      Error in `normalize_n_glycan()`:
      ! `glycans` must contain only topological structures.
      i Use `glyrepr::remove_linkages()` on intact or partial structures first.

