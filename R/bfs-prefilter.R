# Bundled target-independent motif preparation. Custom enzymes and graphs
# with additional metadata retain the R/S3 prefilter.
.bfs_prefilter_plan <- function() {
  glyenzy_bfs_prefilter
}

# Regenerated with the enzyme database by data-raw/glyenzy_enzymes.R.
.build_bfs_prefilter_plan <- function(enzymes) {
  motifs <- list()
  keys <- character()
  alignments <- character()
  add <- function(motif, alignment = "substructure") {
    motif <- glyrepr::as_glycan_structure(motif)
    key <- paste(alignment, as.character(motif), sep = "\r")
    id <- match(key, keys)
    if (is.na(id)) {
      record <- .bfs_native_graph(glyrepr::get_structure_graphs(motif))
      if (is.null(record)) {
        return(NA_integer_)
      }
      id <- length(motifs) + 1L
      motifs[[id]] <<- record
      keys[[id]] <<- key
      alignments[[id]] <<- alignment
    }
    id
  }
  special <- list(
    ncore = add("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"),
    MGAT1 = add("GlcNAc(b1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"),
    MOGS = add(enzymes$MOGS$rules[[1]]$acceptor),
    MAN1B1 = add(
      "Man(a1-2)Man(a1-3)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "core"
    ),
    MAN2A1 = add(
      "Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "core"
    ),
    GANAB = add(
      "Glc(a1-3)Glc(a1-3)Man(a1-2)Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "core"
    ),
    man_a = add("Man(a1-2)Man(a1-3)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-", "core"),
    man_c = add(
      "Man(a1-2)Man(a1-6)Man(a1-6)Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
      "core"
    )
  )
  jobs <- lapply(enzymes, function(enzyme) {
    if (inherits(enzyme, "glyenzy_npre_gt_enzyme")) {
      return("npre")
    }
    if (enzyme$name == "MGAT1" || inherits(enzyme, "glyenzy_gh_enzyme")) {
      return(enzyme$name)
    }
    if (!.uses_standard_graph_action(enzyme)) {
      return(NULL)
    }
    lapply(enzyme$rules, function(rule) {
      list(
        product = add(rule$product, .product_alignment(rule)),
        requires = vapply(
          rule$requires,
          function(requirement) {
            add(requirement$motif, requirement$alignment)
          },
          integer(1)
        )
      )
    })
  })
  list(
    motifs = motifs,
    alignments = alignments,
    jobs = jobs,
    special = special
  )
}

.bfs_native_dictionary <- local({
  dictionary <- NULL
  function() {
    if (is.null(dictionary)) {
      monos <- glyrepr::available_monosaccharides()
      dictionary <<- data.frame(
        concrete = monos,
        generic = glyrepr::convert_to_generic(monos)
      )
    }
    dictionary
  }
})

.bfs_prefilter <- function(enzymes, glycans) {
  out <- rep(NA, length(enzymes))
  names(out) <- names(enzymes)
  builtin <- vapply(
    enzymes,
    function(enzyme) {
      .bfs_is_builtin_enzyme(enzyme)
    },
    logical(1)
  )
  if (!any(builtin)) {
    return(out)
  }
  targets <- lapply(
    glyrepr::get_structure_graphs(glycans, return_list = TRUE),
    .bfs_native_graph
  )
  if (any(vapply(targets, is.null, logical(1)))) {
    return(out)
  }
  plan <- .bfs_prefilter_plan()
  matches <- cpp_bfs_have_motifs(
    targets,
    plan$motifs,
    plan$alignments,
    .glymotif_mode(glycans) == "lenient",
    .bfs_native_dictionary()
  )
  special <- plan$special
  is_n <- matches[, special$ncore]
  n_man <- glyrepr::count_mono(glycans, "Man")
  man123 <- matches[, special$MAN1B1]
  man123[!is.na(n_man) & n_man <= 7L] <- TRUE
  man123[!is.na(n_man) & n_man == 9L] <- FALSE
  man123[
    rowSums(matches[,
      unlist(special[c("man_a", "MAN1B1", "man_c")]),
      drop = FALSE
    ]) ==
      0L
  ] <- TRUE
  for (i in which(builtin)) {
    job <- plan$jobs[[enzymes[[i]]$name]]
    if (is.null(job)) {
      next
    }
    if (is.character(job)) {
      mask <- switch(
        job,
        npre = is_n,
        MGAT1 = is_n & matches[, special$MGAT1],
        MAN1A1 = ,
        MAN1A2 = ,
        MAN1C1 = is_n & man123,
        MAN2A1 = ,
        MAN2A2 = is_n & !matches[, special$MAN2A1],
        MOGS = is_n & !matches[, special$MOGS],
        MAN1B1 = is_n & !matches[, special$MAN1B1],
        GANAB = is_n & !matches[, special$GANAB],
        NULL
      )
      if (!is.null(mask)) out[[i]] <- any(mask)
    } else if (!anyNA(unlist(job))) {
      out[[i]] <- any(vapply(
        job,
        function(rule) {
          mask <- matches[, rule$product]
          if (length(rule$requires)) {
            mask <- mask & rowSums(matches[, rule$requires, drop = FALSE]) > 0L
          }
          any(mask)
        },
        logical(1)
      ))
    }
  }
  out
}
