# Biosynthesis tracing benchmarks

Install the baseline and optimized revisions into separate R libraries. Run the
scripts sequentially, after compilation and tests have finished, using an output
directory outside the source tree.

## Fresh-process first calls

```sh
Rscript data-raw/benchmarks/trace-biosynthesis-cold.R /path/to/baseline-library /path/to/optimized-library /tmp/trace-cold 20
```

Each sample starts a fresh R process, loads the package, and measures one
`trace_biosynthesis()` call with default arguments. Package loading is timed
separately. The script uses no warm-up or pre-GC and checks network equality
across both versions and every repetition.

## Warm-session calls

```sh
Rscript data-raw/benchmarks/trace-biosynthesis.R /path/to/baseline-library /tmp/trace-baseline 50
Rscript data-raw/benchmarks/trace-biosynthesis.R /path/to/optimized-library /tmp/trace-optimized 50
```

The scripts contain the exact input structures and write timing summaries, raw
measurements, and session information to the requested output directories. The
warm-session script also saves network snapshots for equality comparisons.
Benchmark outputs are not retained in this repository.
