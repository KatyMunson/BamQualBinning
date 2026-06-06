# benchmarking

Performance and file-size benchmarks for the QV-binning pipeline. These are
run on demand against real data — they are **not** correctness tests and do not
gate the workflow.

## `benchmark_sizes.sh`

Comparison of compressed size **and processing time** that isolates the levers
independently — quality-score binning **scheme** and kinetics-tag removal:

|                    | keep kinetics       | strip kinetics       |
|--------------------|---------------------|----------------------|
| **raw QVs**        | `raw_keep`          | `raw_strip`          |
| **binned (default)** | `bin_default_keep` | `bin_default_strip` |
| **binned (scheme A)** | `bin_A_keep`     | `bin_A_strip`        |
| **binned (…)**     | …                   | …                    |

The two raw cells use `samtools view` (a neutral htslib baseline, not this
project's own code); the binned cells use `../bin_qv.py`. The default PacBio
Revio scheme is always benchmarked; pass `--bins-file` / `--bins-dir` to add
custom schemes (one binned keep/strip pair per scheme). All cells share thread
count and bgzf compression level, so differences reflect content only.

This is the direct way to answer "how much smaller does scheme X make the file"
— compare each scheme's `ratio_vs_raw_keep` against the others and against
`bin_default`.

**Timing sources:**
- *samtools cells* — wall time from `date +%s%3N`; CPU from `/usr/bin/time` (GNU time)
- *bin_qv.py cells* — wall time, CPU, RSS, and reads/sec read directly from the
  per-run `--metrics` TSV (measured inside the script; more accurate than
  external process timing)

### Usage

```bash
# Whole file (default scheme only)
benchmarking/benchmark_sizes.sh --input sample.hifi_reads.bam --threads 8

# Compare the default scheme against two custom schemes, fast subset
benchmarking/benchmark_sizes.sh --input sample.hifi_reads.bam \
    --subsample 200000 --clean \
    --bins-file bins/coarse.tsv --bins-file bins/binary.tsv

# Benchmark every *.tsv scheme in a directory (plus the always-on default)
benchmarking/benchmark_sizes.sh --input sample.hifi_reads.bam \
    --bins-dir bins/ --clean
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | (required) | Input BAM/UBAM |
| `-o, --outdir` | `./benchmark_results` | Output directory |
| `-s, --subsample N` | `0` (whole file) | Benchmark the first N reads only |
| `-t, --threads N` | `8` | Threads for samtools and bin_qv.py |
| `-b, --bins-file PATH` | (none) | Add a custom bin scheme (repeatable); named by TSV stem |
| `--bins-dir DIR` | (none) | Add every `*.tsv` in DIR as a custom scheme |
| `--clean` | off | Delete generated BAMs after measuring, keep TSV + logs |

The default scheme is always benchmarked. Custom schemes are added on top.
Duplicate scheme names (same TSV stem) are de-duplicated with a numeric suffix.

### Output

`<outdir>/benchmark.tsv`. The file opens with a block of `#INFO` provenance
lines capturing the run parameters (date, input, mode, reads benchmarked,
threads, kinetics tags, schemes compared, samtools/python versions, whether GNU
time was available), followed by one row per cell. Downstream parsers should
skip lines starting with `#`.

Per-cell columns:

| Column | Description |
|---|---|
| `variant` | e.g. `raw_keep`, `bin_default_strip`, `bin_<scheme>_keep` |
| `scheme` | Bin scheme name (`-` for raw cells, `default`, or the custom TSV stem) |
| `kinetics` | `keep` or `strip` |
| `qual` | `raw` or `binned` |
| `bytes` | Compressed BAM size |
| `MB` | Compressed BAM size in MB |
| `bytes_per_read` | Normalized size |
| `ratio_vs_raw_keep` | Size relative to the keep-kinetics/raw-QV baseline |
| `wallclock_sec` | Wall-clock processing time (seconds) |
| `cpu_sec` | CPU time user+sys (seconds); `N/A` for samtools if GNU time unavailable |
| `reads_per_sec` | Read throughput |
| `peak_rss_mb` | Peak RSS (MB); available for `bin_qv.py` cells only |

> **Disk note:** in whole-file mode each output BAM can approach the input size.
> With S schemes that's `2 + 2*S` output BAMs — budget disk accordingly, or use
> `--clean`.
