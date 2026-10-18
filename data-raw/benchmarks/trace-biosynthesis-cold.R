# Usage: Rscript trace-biosynthesis-cold.R BASELINE_LIB OPTIMIZED_LIB OUTPUT [REPEATS]
# Each observation runs in a fresh R process. No pre-GC, warm-up, or prior trace.
args <- commandArgs(trailingOnly = TRUE)
cases <- list(
  high_mannose_N = "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  complex_N = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  O_GalNAc = "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-",
  O_Man = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-",
  example_single = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
)
if (identical(args[[1]], "--worker")) {
  .libPaths(c(normalizePath(args[[2]]), .libPaths()))
  input <- cases[[args[[3]]]]
  start <- proc.time()[["elapsed"]]
  library(glyenzy)
  load_ms <- (proc.time()[["elapsed"]] - start) * 1000
  start <- proc.time()[["elapsed"]]
  network <- trace_biosynthesis(input)
  first_ms <- (proc.time()[["elapsed"]] - start) * 1000
  saveRDS(
    list(
      load_ms = load_ms,
      first_ms = first_ms,
      network = list(
        class = class(network),
        graph = igraph::graph_attr(network),
        data = igraph::as_data_frame(network, what = "both")
      ),
      session = capture.output(sessionInfo())
    ),
    args[[4]]
  )
  quit(status = 0L)
}
stopifnot(length(args) >= 3L)
libraries <- c(
  baseline = normalizePath(args[[1]]),
  optimized = normalizePath(args[[2]])
)
output <- args[[3]]
repeats <- if (length(args) >= 4L) as.integer(args[[4]]) else 20L
stopifnot(!is.na(repeats), repeats > 0L)
dir.create(output, recursive = TRUE, showWarnings = FALSE)
script <- normalizePath(sub(
  "^--file=",
  "",
  commandArgs()[grepl("^--file=", commandArgs())]
))
worker_file <- tempfile(fileext = ".rds")
rows <- list()
reference <- list()
set.seed(20261019)
for (iteration in seq_len(repeats)) {
  for (case in sample(names(cases))) {
    for (version in sample(names(libraries))) {
      status <- system2(
        file.path(R.home("bin"), "Rscript"),
        c(
          "--vanilla",
          shQuote(script),
          "--worker",
          shQuote(libraries[[version]]),
          case,
          shQuote(worker_file)
        )
      )
      stopifnot(status == 0L)
      result <- readRDS(worker_file)
      if (is.null(reference[[case]])) {
        reference[[case]] <- result$network
      }
      stopifnot(identical(reference[[case]], result$network))
      rows[[length(rows) + 1L]] <- data.frame(
        version = version,
        case = case,
        iteration = iteration,
        load_ms = result$load_ms,
        first_ms = result$first_ms
      )
      if (iteration == 1L) {
        writeLines(
          result$session,
          file.path(output, paste0(version, "-session.txt"))
        )
      }
    }
  }
  cat("Completed fresh-process round", iteration, "\n")
}
unlink(worker_file)
raw <- do.call(rbind, rows)
summary <- do.call(
  rbind,
  lapply(names(libraries), function(version) {
    do.call(
      rbind,
      lapply(names(cases), function(case) {
        sample <- raw[raw$version == version & raw$case == case, ]
        data.frame(
          version = version,
          case = case,
          n = nrow(sample),
          load_mean_ms = mean(sample$load_ms),
          first_mean_ms = mean(sample$first_ms),
          first_median_ms = median(sample$first_ms),
          first_p95_ms = unname(quantile(sample$first_ms, 0.95)),
          first_min_ms = min(sample$first_ms),
          first_max_ms = max(sample$first_ms)
        )
      })
    )
  })
)
write.csv(raw, file.path(output, "cold-measurements.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "cold-summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
cat("All networks are identical across versions and repetitions.\n")
