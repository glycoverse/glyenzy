# trace_biosynthesis rejects incompatible glycans

    Code
      trace_biosynthesis("Hex(b1-3)GalNAc(a1-")
    Condition
      Error in `trace_biosynthesis()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "mixed".
      x Detected structure level(s): "intact".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.

---

    Code
      trace_biosynthesis(c("Gal(b1-3)GalNAc(a1-", "Hex(b1-3)HexNAc(a1-"))
    Condition
      Error in `trace_biosynthesis()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "concrete" and "generic".
      x Detected structure level(s): "intact".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.

---

    Code
      trace_biosynthesis("Gal(b1-?)GalNAc(a1-")
    Condition
      Error in `trace_biosynthesis()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "concrete".
      x Detected structure level(s): "partial".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.

---

    Code
      trace_biosynthesis(NA_character_)
    Condition
      Error in `trace_biosynthesis()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "NA".
      x Detected structure level(s): "NA".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.

# trace_biosynthesis_virtual rejects mixed structure levels

    Code
      trace_biosynthesis_virtual(c("Gal(b1-3)GalNAc(a1-", "Gal(??-?)GalNAc(??-"))
    Condition
      Error in `trace_biosynthesis_virtual()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "concrete".
      x Detected structure level(s): "intact" and "topological".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.

# path_biosynthesis requires compatible endpoint structures

    Code
      path_biosynthesis("GalNAc(a1-", "HexNAc(a1-")
    Condition
      Error in `path_biosynthesis()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "concrete" and "generic".
      x Detected structure level(s): "intact".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.

---

    Code
      path_biosynthesis("GalNAc(?1-", "Gal(b1-3)GalNAc(?1-")
    Condition
      Error in `path_biosynthesis()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "concrete".
      x Detected structure level(s): "partial".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.

# path_biosynthesis_virtual requires compatible endpoints

    Code
      path_biosynthesis_virtual("GalNAc(a1-", "GalNAc(??-")
    Condition
      Error in `path_biosynthesis_virtual()`:
      ! Biosynthesis glycans must share one monosaccharide type and one structure level.
      x Detected monosaccharide type(s): "concrete".
      x Detected structure level(s): "intact" and "topological".
      i Supported monosaccharide types are "concrete" and "generic"; mixed-residue and missing glycans are not supported.
      i Supported structure levels are "intact" and "topological"; partial structures are not supported.
      i Use `glyrepr::convert_to_generic()` or `glyrepr::remove_linkages()` to standardize inputs.
