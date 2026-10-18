# Fixed structures are generated alongside the enzyme database. No target or
# biosynthesis result is stored here.
.build_bfs_structures <- function() {
  type_motifs <- glyrepr::as_glycan_structure(c(
    "GlcNAc(b1-4)GlcNAc(?1-",
    "Xyl(??-?)Glc(?1-",
    "Glc(a1-2)Gal(?1-",
    "Man(a1-3)Man(?1-",
    "Man(a1-6)Man(?1-"
  ))
  list(
    precursor = glyparse::parse_iupac_condensed(
      "Glc(a1-2)Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-3)[Man(a1-2)Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ),
    core = glyrepr::n_glycan_core(),
    core_topological = glyrepr::n_glycan_core(linkage = FALSE),
    pre_mgat2 = glyrepr::as_glycan_structure(
      "Man(a1-3/6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    ),
    type_motifs = lapply(
      glyrepr::get_structure_graphs(type_motifs, return_list = TRUE),
      .bfs_native_graph
    )
  )
}
