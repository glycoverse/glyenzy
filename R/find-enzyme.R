#' Identify Potentially Involved Enzymes
#'
#' This function returns all possible isoenzymes associated with the biosynthetic
#' steps of the input glycan.
#' Note that this function ignores the residues in glycans
#' that cannot be matched to any enzyme rules.
#'
#' @inheritSection have_enzyme Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#' @param return_list If `NULL` (default),
#'   return a list of character vectors when `glycans` has length greater than 1,
#'   and a single character vector when `glycans` has length 1.
#'   Set to `TRUE` to always return a list.
#'   This can be useful when you are working programmatically with unknown input length.
#'   Note that when `return_list = FALSE` and `length(glycans) > 1`,
#'   an error will be thrown.
#' @param method Method used to infer enzyme involvement.
#'   `"motif"` checks product motifs directly in each glycan.
#'   `"path"` extracts enzymes from [trace_biosynthesis()] results, which is
#'   more accurate but slower.
#'
#' @return A character vector or a list of character vectors (see `return_list` parameter),
#'   each containing the names of enzymes involved in the biosynthesis of the corresponding glycan.
#'
#' @examples
#' library(glyrepr)
#' library(glyparse)
#'
#' # Use `glycan_structure()`
#' glycans <- auto_parse(c(
#'   "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
#'   "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' ))
#' find_enzyme(glycans)
#'
#' # Or use characters directly
#' find_enzyme("GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-")
#'
#' # Use reconstructed biosynthesis paths
#' find_enzyme(glycans, method = "path")
#'
#' @export
find_enzyme <- function(
  glycans,
  return_list = NULL,
  method = c("motif", "path")
) {
  method <- match.arg(method)
  glycans <- .process_glycans_arg(glycans)
  return_list <- .validate_return_list(return_list, length(glycans))
  res <- switch(
    method,
    motif = .find_enzyme_motif(glycans),
    path = .find_enzyme_path(glycans)
  )

  .format_result(res, return_list)
}

#' Identify enzymes using final-glycan motif matching
#'
#' @param glycans A `glyrepr_structure` vector.
#'
#' @returns A list of character vectors.
#' @noRd
.find_enzyme_motif <- function(glycans) {
  if (length(glycans) == 0L) {
    return(list())
  }

  enzymes <- glyenzy_enzymes
  types <- .glycan_types(glycans)
  compatible <- vapply(
    enzymes,
    .enzyme_type_mask_from_types,
    logical(length(glycans)),
    types = types
  )
  if (length(glycans) == 1L) {
    compatible <- matrix(compatible, nrow = 1L)
  }
  colnames(compatible) <- names(enzymes)

  plan <- .find_enzyme_rule_plan()
  enzyme_mask <- matrix(
    FALSE,
    nrow = length(glycans),
    ncol = length(enzymes),
    dimnames = list(names(glycans), names(enzymes))
  )

  enzyme_mask <- .find_enzyme_special_masks(
    glycans,
    enzymes,
    compatible,
    plan$special_groups,
    enzyme_mask
  )
  enzyme_mask <- .find_enzyme_substituent_masks(
    glycans,
    compatible,
    plan$sub_rows,
    enzyme_mask
  )
  enzyme_mask <- .find_enzyme_batch_masks(
    glycans,
    enzymes,
    compatible,
    plan$batch_rows,
    enzyme_mask
  )

  enzyme_mask <- enzyme_mask & compatible
  purrr::map(
    seq_along(glycans),
    ~ names(enzymes)[enzyme_mask[.x, ]]
  )
}

.find_enzyme_rule_plan <- local({
  plan <- NULL
  function() {
    if (!is.null(plan)) {
      return(plan)
    }

    rows <- list()
    i <- 0L
    for (enzyme_idx in seq_along(glyenzy_enzymes)) {
      enzyme <- glyenzy_enzymes[[enzyme_idx]]
      for (rule_idx in seq_along(enzyme$rules)) {
        rule <- enzyme$rules[[rule_idx]]
        i <- i + 1L
        product_graph <- glyrepr::get_structure_graphs(
          rule$product,
          return_list = FALSE
        )
        rows[[i]] <- list(
          enzyme_idx = enzyme_idx,
          rule = rule,
          product = rule$product,
          product_chr = as.character(rule$product),
          alignment = .product_alignment(rule),
          special = inherits(enzyme, "glyenzy_gh_enzyme") ||
            inherits(enzyme, "glyenzy_npre_gt_enzyme") ||
            identical(enzyme$name, "MGAT1"),
          has_sub = any(igraph::V(product_graph)$sub != "")
        )
      }
    }

    is_special <- vapply(rows, `[[`, logical(1), "special")
    has_sub <- vapply(rows, `[[`, logical(1), "has_sub")
    special_enzyme_idx <- unique(vapply(
      rows[is_special],
      `[[`,
      integer(1),
      "enzyme_idx"
    ))
    special_groups <- split(
      special_enzyme_idx,
      vapply(
        special_enzyme_idx,
        .find_enzyme_special_group,
        character(1),
        enzymes = glyenzy_enzymes
      )
    )
    plan <<- list(
      batch_rows = rows[!is_special & !has_sub],
      sub_rows = rows[!is_special & has_sub],
      special_groups = special_groups
    )
    plan
  }
})

.find_enzyme_special_group <- function(enzyme_idx, enzymes) {
  enzyme <- enzymes[[enzyme_idx]]
  if (inherits(enzyme, "glyenzy_npre_gt_enzyme")) {
    return("npre")
  }
  if (enzyme$name %in% c("MAN1A1", "MAN1A2", "MAN1C1")) {
    return("man123")
  }
  if (enzyme$name %in% c("MAN2A1", "MAN2A2")) {
    return("man2a12")
  }
  enzyme$name
}

.enzyme_type_mask_from_types <- function(types, enzyme) {
  if (is.null(enzyme$glycan_type)) {
    return(rep(TRUE, length(types)))
  }
  vapply(
    types,
    .glycan_type_is_compatible,
    logical(1),
    supported_types = enzyme$glycan_type
  )
}

.find_enzyme_special_masks <- function(
  glycans,
  enzymes,
  compatible,
  special_groups,
  enzyme_mask
) {
  for (enzyme_group in special_groups) {
    enzyme_idx <- enzyme_group[[1]]
    group_compatible <- apply(
      compatible[, enzyme_group, drop = FALSE],
      1L,
      any
    )
    shared_mask <- .safe_have_enzyme_with_mask(
      glycans,
      enzymes[[enzyme_idx]],
      group_compatible
    )
    enzyme_mask[, enzyme_group] <- shared_mask
  }
  enzyme_mask
}

.find_enzyme_substituent_masks <- function(
  glycans,
  compatible,
  sub_rows,
  enzyme_mask
) {
  glycans_have_sub <- .has_substituents(glycans)
  for (row in sub_rows) {
    enzyme_idx <- row$enzyme_idx
    eligible <- compatible[, enzyme_idx] & glycans_have_sub
    if (!any(eligible)) {
      next
    }
    matched <- rep(FALSE, length(glycans))
    matched[eligible] <- tryCatch(
      .have_motif_substituent_subset(
        glycans[eligible],
        row$rule$product,
        alignment = row$alignment
      ) &
        .rule_requirements_met(glycans[eligible], row$rule),
      error = function(e) rep(FALSE, sum(eligible))
    )
    enzyme_mask[, enzyme_idx] <- enzyme_mask[, enzyme_idx] | matched
  }
  enzyme_mask
}

.find_enzyme_batch_masks <- function(
  glycans,
  enzymes,
  compatible,
  batch_rows,
  enzyme_mask
) {
  non_intact <-
    .glycan_structure_levels(glycans) != "intact" |
    glyrepr::get_mono_type(glycans) != "concrete"
  non_intact[is.na(non_intact)] <- TRUE

  batch_meta <- do.call(
    rbind,
    lapply(batch_rows, function(row) {
      enzyme_compatible <- compatible[, row$enzyme_idx]
      data.frame(
        enzyme_idx = row$enzyme_idx,
        product_chr = row$product_chr,
        alignment = row$alignment,
        mode = if (any(enzyme_compatible & non_intact)) "lenient" else "strict",
        stringsAsFactors = FALSE
      )
    })
  )
  batch_meta$key <- paste(
    batch_meta$mode,
    batch_meta$alignment,
    batch_meta$product_chr,
    sep = "\r"
  )

  groups <- split(
    seq_along(batch_rows),
    interaction(
      batch_meta$mode,
      batch_meta$alignment,
      drop = TRUE
    )
  )
  for (group_idx in groups) {
    unique_idx <- group_idx[!duplicated(batch_meta$key[group_idx])]
    motifs <- do.call(c, lapply(batch_rows[unique_idx], `[[`, "product"))
    motif_keys <- batch_meta$key[unique_idx]
    names(motifs) <- motif_keys
    matches <- tryCatch(
      glymotif::have_motifs(
        glycans,
        motifs,
        alignments = batch_meta$alignment[unique_idx],
        strict_sub = FALSE,
        mode = batch_meta$mode[[unique_idx[[1]]]]
      ),
      error = function(e) NULL
    )
    if (is.null(matches)) {
      fallback_idx <- unique(batch_meta$enzyme_idx[group_idx])
      for (enzyme_idx in fallback_idx) {
        enzyme_mask[, enzyme_idx] <-
          enzyme_mask[, enzyme_idx] |
          .safe_have_enzyme_with_mask(
            glycans,
            enzymes[[enzyme_idx]],
            compatible[, enzyme_idx]
          )
      }
      next
    }
    if (length(glycans) == 1L) {
      matches <- matrix(
        matches,
        nrow = 1L,
        dimnames = list(names(glycans), motif_keys)
      )
    }
    for (row_idx in group_idx) {
      enzyme_idx <- batch_meta$enzyme_idx[[row_idx]]
      rule_matches <- matches[, batch_meta$key[[row_idx]]]
      rule <- batch_rows[[row_idx]]$rule
      if (length(rule$requires) > 0L) {
        requirements_met <- rep(FALSE, length(glycans))
        eligible <- compatible[, enzyme_idx]
        if (any(eligible)) {
          requirements_met[eligible] <- .rule_requirements_met(
            glycans[eligible],
            rule
          )
        }
        rule_matches <- rule_matches & requirements_met
      }
      enzyme_mask[, enzyme_idx] <-
        enzyme_mask[, enzyme_idx] |
        rule_matches
    }
  }
  enzyme_mask
}

.safe_have_enzyme_with_mask <- function(glycans, enzyme, compatible) {
  result <- rep(FALSE, length(glycans))
  if (!any(compatible)) {
    return(result)
  }
  tryCatch(
    {
      result[compatible] <- .have_enzyme_motif(glycans[compatible], enzyme)
      result
    },
    error = function(e) result
  )
}

#' Identify trace-derived enzymes for each glycan
#'
#' @param glycans A `glyrepr_structure` vector.
#'
#' @returns A list of character vectors.
#' @noRd
.find_enzyme_path <- function(glycans) {
  purrr::map(glycans, .find_enzyme_path_single)
}

#' Identify trace-derived enzymes for one glycan
#'
#' @param glycan A length-one `glyrepr_structure` vector.
#'
#' @returns A character vector of enzyme names.
#' @noRd
.find_enzyme_path_single <- function(glycan) {
  npre_enzymes <- .find_npre_enzymes(glycan)
  traced_enzymes <- .trace_enzyme_edges_single(glycan)
  unique(c(traced_enzymes, npre_enzymes))
}

#' Identify N-glycan precursor enzymes for one glycan
#'
#' @param glycan A length-one `glyrepr_structure` vector.
#'
#' @returns A character vector of enzyme names.
#' @noRd
.find_npre_enzymes <- function(glycan) {
  npre_enzymes <- purrr::keep(glyenzy_enzymes, .is_npre_gt)
  masks <- purrr::map_lgl(
    npre_enzymes,
    ~ .enzyme_glycan_type_mask(glycan, .x) &&
      .have_enzyme_motif.glyenzy_npre_gt_enzyme(glycan, .x)
  )
  names(npre_enzymes)[masks]
}

# Like `.have_enzyme_motif()`, but returns FALSE instead of throwing error.
.safe_have_enzyme <- function(glycans, enzyme) {
  compatible <- .enzyme_glycan_type_mask(glycans, enzyme)
  result <- rep(FALSE, length(glycans))
  if (!any(compatible)) {
    return(result)
  }
  tryCatch(
    {
      result[compatible] <- .have_enzyme_motif(glycans[compatible], enzyme)
      result
    },
    error = function(e) result
  )
}
