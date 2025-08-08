#' Make a function that only applies to N-glycans
#'
#' This function is used to wrap a function that only applies to N-glycans.
#' For N-glycans in the input `glycans`, it uses the values in `.f`.
#' Otherwise, it returns FALSE.
#'
#' @param .f A function that takes two arguments: `glycans` and `enzyme`.
#' @noRd
.make_n_glycan_guard <- function(.f, type = "logical") {
  force(.f)
  function(glycans, enzyme, is_n = NULL) {
    if (is.null(is_n)) {
      is_n <- glymotif::is_n_glycan(glycans)
    }
    if (type == "logical") {
      res <- rep(FALSE, length(glycans))
    } else {  # type is integer
      res <- rep(0, length(glycans))
    }
    if (any(is_n)) {
      res[is_n] <- .f(glycans[is_n], enzyme)
    }
    res
  }
}

.process_glycans_arg <- function(x) {
  if (is.character(x)) {
    x <- glyparse::auto_parse(x)
  } else if (!glyrepr::is_glycan_structure(x)) {
    cli::cli_abort(c(
      "{.arg glycans} must be a {.cls glyrepr_structure} vector or a character vector of glycan structure strings.",
      "x" = "Got {.cls {class(x)}}."
    ))
  }
  return(x)
}

.process_enzyme_arg <- function(x) {
  if (is.character(x)) {
    x <- enzyme(x)
  } else if (!inherits(x, "glyenzy_enzyme")) {
    cli::cli_abort(c(
      "{.arg enzyme} must be a {.cls glyenzy_enzyme} object or a character string of gene symbol.",
      "x" = "Got {.cls {class(x)}}."
    ))
  }
  return(x)
}