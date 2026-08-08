#' Calculate Product-Substrate Ratios
#'
#' Calculate the ratio between the motif quantification of enzyme products and
#' substrates. For glycoproteomics data, ratios are calculated independently
#' for each glycosite.
#'
#' Each rule contributes its product and acceptor motif once. Quantifications
#' are summed across all rules of an enzyme before division, including when
#' multiple rules contain the same motif. Rule `rejects` and `requires` are
#' ignored. Ratios with zero substrate quantification are returned as `NA`.
#'
#' Motifs are matched strictly for intact glycan structures and leniently for
#' partial or reduced structures, consistent with other glyenzy motif-based
#' functions.
#'
#' @param exp A [glyexp::GlycomicSE()] or [glyexp::GlycoproteomicSE()] object.
#' @param enzymes A character vector of enzyme names or a list of [enzyme()]
#'   objects. Starter glycosyltransferases are not supported because their
#'   substrates are not glycans.
#'
#' @returns A plain [SummarizedExperiment::SummarizedExperiment()] with a
#'   `product_substrate_ratio` assay. For glycomics data, `rowData()` contains
#'   `enzyme`. For glycoproteomics data, it contains `enzyme`, `protein`, and
#'   `protein_site`. `colData()` and metadata are preserved from `exp`.
#'
#' @examples
#' exp <- glyexp::real_experiment2[seq_len(10), seq_len(3)]
#' product_substrate_ratio(exp, "ST6GAL1")
#'
#' @export
product_substrate_ratio <- function(exp, enzymes) {
  exp_type <- .product_substrate_exp_type(exp)
  .check_product_substrate_exp(exp)
  enzymes <- .product_substrate_enzymes(enzymes)
  sites <- .product_substrate_sites(exp, exp_type)
  entries <- .product_substrate_entries(enzymes)

  n_sites <- nrow(sites)
  n_enzymes <- length(enzymes)
  n_samples <- ncol(exp)
  n_rows <- n_sites * n_enzymes
  matrix_dimnames <- list(NULL, colnames(exp))
  product <- matrix(
    0,
    nrow = n_rows,
    ncol = n_samples,
    dimnames = matrix_dimnames
  )
  substrate <- matrix(
    0,
    nrow = n_rows,
    ncol = n_samples,
    dimnames = matrix_dimnames
  )

  for (alignment in unique(entries$alignment)) {
    entry_idx <- which(entries$alignment == alignment)
    motif_keys <- vapply(
      entries$motif[entry_idx],
      as.character,
      character(1)
    )
    unique_keys <- unique(motif_keys)
    unique_idx <- match(unique_keys, motif_keys)
    motifs <- do.call(c, entries$motif[entry_idx[unique_idx]])
    motif_names <- paste0("motif_", seq_along(motifs))
    names(motifs) <- motif_names
    motif_quantifications <- .quantify_product_substrate_motifs(
      exp,
      exp_type,
      sites,
      motifs,
      alignment
    )

    for (i in seq_along(entry_idx)) {
      entry <- entry_idx[[i]]
      motif_name <- motif_names[[match(motif_keys[[i]], unique_keys)]]
      output_rows <- seq(
        from = entries$enzyme_index[[entry]],
        by = n_enzymes,
        length.out = n_sites
      )
      values <- motif_quantifications[[motif_name]]
      if (entries$role[[entry]] == "product") {
        product[output_rows, ] <- product[output_rows, , drop = FALSE] +
          values
      } else {
        substrate[output_rows, ] <- substrate[output_rows, , drop = FALSE] +
          values
      }
    }
  }

  ratio <- product / substrate
  ratio[substrate == 0] <- NA_real_
  .new_product_substrate_result(exp, exp_type, enzymes, sites, ratio)
}

.product_substrate_exp_type <- function(exp) {
  if (glyexp::is_glycomic_se(exp)) {
    return("glycomics")
  }
  if (glyexp::is_glycoproteomic_se(exp)) {
    return("glycoproteomics")
  }
  cli::cli_abort(
    paste0(
      "{.arg exp} must be a {.cls GlycomicSE} or ",
      "{.cls GlycoproteomicSE} object."
    )
  )
}

.check_product_substrate_exp <- function(exp) {
  row_data <- SummarizedExperiment::rowData(exp)
  if (!"glycan_structure" %in% colnames(row_data)) {
    cli::cli_abort(
      "{.arg exp} must contain {.field glycan_structure} in its row data."
    )
  }
  invisible(exp)
}

.product_substrate_enzymes <- function(enzymes) {
  if (is.null(enzymes)) {
    cli::cli_abort(
      "{.arg enzymes} must be a character vector or a list of enzyme objects."
    )
  }
  enzymes <- .enzymes_from_arg(enzymes)
  if (length(enzymes) == 0L) {
    cli::cli_abort("{.arg enzymes} must not be empty.")
  }

  no_rules <- purrr::map_lgl(enzymes, \(enzyme) length(enzyme$rules) == 0L)
  if (any(no_rules)) {
    cli::cli_abort(
      "Enzymes must contain at least one rule: {.val {names(enzymes)[no_rules]}}."
    )
  }

  starter <- purrr::map_lgl(enzymes, .is_starter_gt)
  empty_acceptor <- purrr::map_lgl(
    enzymes,
    function(enzyme) {
      any(
        purrr::map_int(
          enzyme$rules,
          \(rule) length(rule$acceptor)
        ) ==
          0L
      )
    }
  )
  unsupported <- starter | empty_acceptor
  if (any(unsupported)) {
    cli::cli_abort(c(
      "Enzymes with non-glycan substrates are not supported.",
      "x" = paste0(
        "Their substrates cannot be quantified: ",
        "{.val {names(enzymes)[unsupported]}}."
      )
    ))
  }

  unname(enzymes)
}

.product_substrate_sites <- function(exp, exp_type) {
  if (exp_type == "glycomics") {
    return(data.frame(.site = 1L))
  }

  row_data <- SummarizedExperiment::rowData(exp)
  sites <- data.frame(
    protein = as.character(row_data$protein),
    protein_site = as.integer(row_data$protein_site),
    stringsAsFactors = FALSE
  )
  sites[!duplicated(sites), , drop = FALSE]
}

.product_substrate_entries <- function(enzymes) {
  enzyme_index <- integer()
  role <- character()
  alignment <- character()
  motif <- list()

  for (i in seq_along(enzymes)) {
    for (rule in enzymes[[i]]$rules) {
      enzyme_index <- c(enzyme_index, i, i)
      role <- c(role, "product", "substrate")
      alignment <- c(
        alignment,
        .product_alignment(rule),
        rule$acceptor_alignment
      )
      motif <- c(motif, list(rule$product), list(rule$acceptor))
    }
  }

  list(
    enzyme_index = enzyme_index,
    role = role,
    alignment = alignment,
    motif = motif
  )
}

.quantify_product_substrate_motifs <- function(
  exp,
  exp_type,
  sites,
  motifs,
  alignment
) {
  row_data <- SummarizedExperiment::rowData(exp)
  motif_counts <- .count_motifs(
    row_data$glycan_structure,
    motifs,
    alignments = rep(alignment, length(motifs))
  )
  abundance <- as.matrix(SummarizedExperiment::assay(exp, 1))
  abundance[is.na(abundance)] <- 0
  motif_names <- names(motifs)

  if (exp_type == "glycomics") {
    quantified <- crossprod(motif_counts, abundance)
    return(stats::setNames(
      lapply(seq_along(motif_names), \(i) quantified[i, , drop = FALSE]),
      motif_names
    ))
  }

  site_idx <- match(
    .product_substrate_site_key(row_data$protein, row_data$protein_site),
    .product_substrate_site_key(sites$protein, sites$protein_site)
  )
  stats::setNames(
    lapply(seq_along(motif_names), function(i) {
      weighted_abundance <- abundance * motif_counts[, i]
      rowsum(weighted_abundance, site_idx, reorder = FALSE, na.rm = TRUE)
    }),
    motif_names
  )
}

.product_substrate_site_key <- function(protein, protein_site) {
  paste0(
    protein,
    "\r",
    ifelse(is.na(protein_site), "<NA>", as.character(protein_site))
  )
}

.new_product_substrate_result <- function(
  exp,
  exp_type,
  enzymes,
  sites,
  ratio
) {
  enzyme_names <- purrr::map_chr(enzymes, "name")
  if (exp_type == "glycomics") {
    row_ids <- make.unique(enzyme_names)
    row_data <- S4Vectors::DataFrame(
      enzyme = enzyme_names,
      row.names = row_ids
    )
  } else {
    site_idx <- rep(seq_len(nrow(sites)), each = length(enzymes))
    enzyme_idx <- rep(seq_along(enzymes), times = nrow(sites))
    protein_site_label <- ifelse(
      is.na(sites$protein_site[site_idx]),
      "NA",
      as.character(sites$protein_site[site_idx])
    )
    row_ids <- make.unique(paste(
      sites$protein[site_idx],
      protein_site_label,
      enzyme_names[enzyme_idx],
      sep = "-"
    ))
    row_data <- S4Vectors::DataFrame(
      enzyme = enzyme_names[enzyme_idx],
      protein = sites$protein[site_idx],
      protein_site = sites$protein_site[site_idx],
      row.names = row_ids
    )
  }

  rownames(ratio) <- row_ids
  SummarizedExperiment::SummarizedExperiment(
    assays = list(product_substrate_ratio = ratio),
    rowData = row_data,
    colData = SummarizedExperiment::colData(exp),
    metadata = S4Vectors::metadata(exp)
  )
}
