# Usage: Rscript data-raw/benchmarks/trace-biosynthesis.R LIBRARY OUTPUT [REPEATS]
# Run sequentially in otherwise idle R sessions, using normal installed builds.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
.libPaths(c(normalizePath(args[[1]]), .libPaths()))
library(glyenzy)
output <- args[[2]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
repeats <- if (length(args) >= 3L) as.integer(args[[3]]) else 50L
stopifnot(!is.na(repeats), repeats > 0L)
cases <- list(
  high_mannose_N = "Man(a1-2)Man(a1-3)[Man(a1-3)[Man(a1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  complex_N = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
  O_GalNAc = "Gal(b1-4)GlcNAc(b1-6)[Gal(b1-3)]GalNAc(a1-",
  O_Man = "Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-2)Man(a1-",
  example_single = "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-",
  example_multi = c(
    "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-",
    "Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-3)GalNAc(a1-"
  )
)
first_ms <- numeric(length(cases))
networks <- vector("list", length(cases))
names(networks) <- names(cases)
for (i in seq_along(cases)) {
  start <- proc.time()[["elapsed"]]
  network <- trace_biosynthesis(cases[[i]])
  first_ms[[i]] <- (proc.time()[["elapsed"]] - start) * 1000
  stopifnot(sum(igraph::vertex_attr(network, "target")) == length(cases[[i]]))
  networks[[i]] <- list(
    class = class(network),
    graph = igraph::graph_attr(network),
    data = igraph::as_data_frame(network, what = "both")
  )
}
saveRDS(networks, file.path(output, "networks.rds"))
set.seed(20261018)
rows <- vector("list", repeats * length(cases))
k <- 0L
for (iteration in seq_len(repeats)) {
  for (name in sample(names(cases))) {
    start <- proc.time()[["elapsed"]]
    network <- trace_biosynthesis(cases[[name]])
    elapsed_ms <- (proc.time()[["elapsed"]] - start) * 1000
    k <- k + 1L
    rows[[k]] <- data.frame(
      case = name,
      iteration = iteration,
      elapsed_ms = elapsed_ms
    )
  }
  if (iteration %% 10L == 0L) cat("Completed round", iteration, "\n")
}
raw <- do.call(rbind, rows)
summary <- do.call(
  rbind,
  lapply(seq_along(cases), function(i) {
    x <- raw$elapsed_ms[raw$case == names(cases)[[i]]]
    data.frame(
      case = names(cases)[[i]],
      targets = length(cases[[i]]),
      n = length(x),
      mean_ms = mean(x),
      median_ms = median(x),
      sd_ms = sd(x),
      p95_ms = unname(quantile(x, 0.95)),
      min_ms = min(x),
      max_ms = max(x),
      first_ms = first_ms[[i]]
    )
  })
)
write.csv(raw, file.path(output, "raw.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
print(summary, row.names = FALSE)
cat(
  "Single-target equal-weight mean:",
  mean(summary$mean_ms[summary$targets == 1L]),
  "ms\n"
)
