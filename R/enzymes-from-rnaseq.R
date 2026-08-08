#' Filter enzymes using RNA-seq expression
#'
#' `enzymes_from_rnaseq()` selects built-in enzymes whose genes are expressed
#' at or above a TPM threshold. When `mat` contains multiple columns, the mean
#' TPM across columns is used for each gene.
#'
#' @param mat A numeric matrix of TPM values with gene symbols as unique row
#'   names and samples or replicates as columns.
#' @param threshold A non-negative numeric scalar giving the minimum mean TPM
#'   for a gene to be considered expressed. The default follows Huang et al.
#'   (2021): genes with TPM below 1 are treated as not expressed or rarely
#'   expressed, whereas genes with TPM at least 1 are treated as expressed.
#'
#' @returns A named list of [enzyme()] objects. Genes not represented in the
#'   built-in enzyme database are ignored.
#'
#' @references
#' Huang, Y.-F., Aoki, K., Akase, S., et al. (2021). Global mapping of
#' glycosylation pathways in human-derived cells. *Developmental Cell*, 56,
#' 1195-1209. \doi{10.1016/j.devcel.2021.02.023}.
#'
#' @examples
#' tpm <- matrix(
#'   c(0.2, 1, 3),
#'   ncol = 1,
#'   dimnames = list(c("FUT8", "ST3GAL3", "MOGS"), "sample")
#' )
#' enzymes_from_rnaseq(tpm)
#'
#' @export
enzymes_from_rnaseq <- function(mat, threshold = 1) {
  checkmate::assert_matrix(
    mat,
    mode = "numeric",
    any.missing = FALSE,
    min.rows = 1,
    min.cols = 1
  )
  checkmate::assert_names(rownames(mat), type = "unique")
  checkmate::assert_numeric(
    as.vector(mat),
    lower = 0,
    finite = TRUE,
    any.missing = FALSE
  )
  checkmate::assert_number(threshold, lower = 0, finite = TRUE)

  mean_tpm <- rowMeans(mat)
  expressed <- rownames(mat)[mean_tpm >= threshold]
  glyenzy_enzymes[names(glyenzy_enzymes) %in% expressed]
}
