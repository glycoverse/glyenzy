# Run against an installed optimized build, not pkgload's debug compilation.
# BFS_NATIVE_LIB may select an isolated R library containing that build.
lib <- Sys.getenv("BFS_NATIVE_LIB", "")
if (nzchar(lib)) {
  .libPaths(c(lib, .libPaths()))
}
library(glyenzy)
ns <- asNamespace("glyenzy")
reference_env <- new.env(parent = ns)
sys.source("tests/testthat/fixtures/bfs-reference.R", envir = reference_env)
case_env <- new.env(parent = ns)
sys.source("tests/testthat/fixtures/bfs-native-cases.R", envir = case_env)
normalize <- function(x) {
  for (name in c("parent", "parent_enzyme", "parent_step")) {
    values <- as.list(x[[name]])
    x[[name]] <- if (length(values)) values[order(names(values))] else list()
  }
  x
}
rows <- list()
for (label in names(case_env$cases)) {
  case <- case_env$cases[[label]]
  args <- list(
    from_g = glyrepr::as_glycan_structure(case[[1]]),
    to_gs = glyrepr::as_glycan_structure(case[[2]]),
    enzymes = lapply(case[[3]], enzyme),
    max_steps = case[[4]],
    allow_partial = TRUE
  )
  searches <- list(
    R = reference_env$bfs_synthesis_search,
    native = get("bfs_synthesis_search", ns)
  )
  expected <- normalize(do.call(searches$R, args))
  stopifnot(identical(expected, normalize(do.call(searches$native, args))))
  for (i in seq_len(7L)) {
    for (arm in if (i %% 2L) c("R", "native") else c("native", "R")) {
      elapsed <- system.time(result <- do.call(searches[[arm]], args))[[
        "elapsed"
      ]]
      stopifnot(identical(expected, normalize(result)))
      rows[[length(rows) + 1L]] <- data.frame(
        case = label,
        run = i,
        arm = arm,
        elapsed = elapsed,
        edges = length(result$all_edges),
        visited = length(as.list(result$parent)) + 1L
      )
      cat(label, arm, i, elapsed, "\n")
    }
  }
}
timings <- do.call(rbind, rows)
write.csv(
  timings,
  "data-raw/bfs-native/integrated-timings.csv",
  row.names = FALSE
)
summary <- do.call(
  rbind,
  lapply(split(timings, timings$case), function(x) {
    r <- x$elapsed[x$arm == "R"]
    cpp <- x$elapsed[x$arm == "native"]
    data.frame(
      case = x$case[[1]],
      R_median = median(r),
      native_median = median(cpp),
      speedup = median(r) / median(cpp),
      R_min = min(r),
      R_max = max(r),
      native_min = min(cpp),
      native_max = max(cpp),
      paired_wins = sum(r > cpp)
    )
  })
)
write.csv(
  summary,
  "data-raw/bfs-native/integrated-summary.csv",
  row.names = FALSE
)
capture.output(
  sessionInfo(),
  file = "data-raw/bfs-native/integrated-session.txt"
)
print(summary, row.names = FALSE)

manifest <- tools::md5sum(c(
  "DESCRIPTION",
  "R/bfs-native.R",
  "R/bfs-path.R",
  "src/bfs-search.cpp",
  "src/bfs-matcher.h",
  "src/bfs-vf2.h",
  getLoadedDLLs()[["glyenzy"]][["path"]]
))
writeLines(
  c(
    paste(
      "Git revision:",
      system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    ),
    paste(names(manifest), manifest)
  ),
  "data-raw/bfs-native/integrated-build.txt"
)
