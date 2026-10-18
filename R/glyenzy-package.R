# Import the namespaces used by every native search at package load time.
# No search or target-specific preparation is performed while loading.
#' @useDynLib glyenzy, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @keywords internal
#' @import glyrepr
#' @importFrom igraph graph_from_data_frame
#' @importFrom checkmate assert_class
#' @importFrom fastmap fastmap
#' @importFrom purrr map
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

ignore_unused_imports <- function() {
  dplyr::case_when
  R6::is.R6
  fastmap::fastmap
}
