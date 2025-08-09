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

# Validate and process return_list parameter early
.validate_return_list <- function(return_list, input_length) {
  checkmate::assert_flag(return_list, null.ok = TRUE)
  checkmate::assert_integerish(input_length, len = 1, lower = 1)

  if (is.null(return_list)) {
    return_list <- input_length > 1
  }
  if (!return_list && input_length > 1) {
    cli::cli_abort(c(
      "When {.arg return_list} is FALSE, input must have length 1.",
      "x" = "Input length: {.val {input_length}}."
    ))
  }

  return_list
}

# Format result based on return_list setting
.format_result <- function(result_list, return_list) {
  checkmate::assert_list(result_list)
  checkmate::assert_flag(return_list)

  if (!return_list) {
    result_list[[1]]
  } else {
    result_list
  }
}