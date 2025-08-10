#' Apply Enzymes to Spawn Glycans
#'
#' This function simulates the action of enzymes on glycans.
#' Think of it like a primordial soup where you put in a few glycans and enzymes,
#' and let them interact to generate new glycans.
#'
#' `spawn_glycans_step()` performs one round of enzyme action,
#' while `spawn_glycans()` performs multiple rounds.
#' The only difference between `spawn_glycans_step()` and `spawn_glycans(n_steps = 1)`
#' is that the latter returns the original input glycans as well.
#' For both, a vector of unique glycan structures is returned.
#'
#' @inheritSection is_synthesized_by Important notes
#'
#' @param glycans A [glyrepr::glycan_structure()], or a character vector of
#'   glycan structure strings supported by [glyparse::auto_parse()].
#' @param enzymes A character vector of gene symbols,
#'   or a list of [enzyme()] objects.
#' @param n_steps The maximum number of rounds to perform.
#'   The actual number of rounds may be less if no new glycans can be generated.
#'
#' @returns A [glyrepr::glycan_structure()] vector of all unique glycans generated.
#'
#' @examples
#' # Use `spawn_glycans_step()` to build glycans step by step
#' glycan <- "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' glycan |>
#'   spawn_glycans_step("MGAT2") |>
#'   spawn_glycans_step("B4GALT1") |>
#'   spawn_glycans_step("ST3GAL3")
#'
#' # Use `spawn_glycans()` to simulate a primordial soup
#' glycans <- c(
#'   "GlcNAc(b1-2)Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
#'   "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#' )
#' enzymes <- c("B4GALT1", "ST3GAL3")
#' spawn_glycans(glycans, enzymes, n_steps = 5)
#'
#' @export
spawn_glycans_step <- function(glycans, enzymes) {
  glycans <- .process_glycans_arg(glycans)
  enzymes <- .process_enzymes_arg(enzymes)
  .spawn_glycans_step(glycans, enzymes)
}

#' @rdname spawn_glycans_step
#' @export
spawn_glycans <- function(glycans, enzymes, n_steps = 10) {
  glycans <- .process_glycans_arg(glycans)
  enzymes <- .process_enzymes_arg(enzymes)
  .spawn_glycans(glycans, enzymes, n_steps)
}

.spawn_glycans_step <- function(glycans, enzymes) {
  new_glycans <- purrr::map(enzymes, ~ .apply_enzyme_flatten(glycans, .x))
  unique(do.call(c, new_glycans))
}

.apply_enzyme_flatten <- function(glycans, enzyme) {
  do.call(c, .apply_enzyme(glycans, enzyme))
}

.spawn_glycans <- function(glycans, enzymes, n_steps) {
  pool <- vector("list", n_steps + 1)
  pool[[1]] <- unique(glycans)

  # Initialize progress bar
  cli::cli_progress_bar(
    "Spawning glycans",
    total = n_steps,
    format = "{cli::pb_spin} Step {cli::pb_current}/{cli::pb_total} | {cli::pb_bar} {cli::pb_percent}"
  )

  for (i in 1:n_steps) {
    cli::cli_progress_update()

    new_glycans <- .spawn_glycans_step(glycans, enzymes)
    if (length(new_glycans) == 0L) {
      cli::cli_progress_done()
      break
    }
    glycans <- new_glycans
    pool[[i + 1]] <- new_glycans
  }

  cli::cli_progress_done()
  unique(do.call(c, pool))
}