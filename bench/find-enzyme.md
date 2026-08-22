# `find_enzyme()` motif benchmark

Benchmark run on 2026-08-22 with R 4.6.1 on arm64 macOS, using glyenzy
0.8.2.9000 and glymotif 0.18.1.9000.

The workload in `find-enzyme.R` contains five intact N-, O-, sialylated, and
sulfated glycans. The script first requires exact equality between the legacy
enzyme-wise result and the batched result, then records three elapsed-time runs
per implementation.

| Implementation | Median seconds | Speedup |
|---|---:|---:|
| Legacy enzyme-wise loop | 16.082 | 1.0x |
| Batched product matching | 0.446 | 36.1x |

A separate fresh-session run of the batched implementation took 0.499 seconds,
including construction of the cached rule plan.
