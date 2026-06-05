# benchmarking

Performance and file-size benchmarks for the QV-binning pipeline. These are
run on demand against real data — they are **not** correctness tests and do not
gate the workflow.

## `benchmark_sizes.sh`

4-way comparison of compressed size **and processing time** that isolates the
two levers independently — quality-score binning and kinetics-tag removal:

|              | keep kinetics | strip kinetics |
|--------------|---------------|----------------|
| **raw QVs**  | `raw_keep`    | `raw_strip`    |
| **binned QVs** | `bin_keep`  | `bin_strip`    |

The two raw cells use `samtools view` (a neutral htslib baseline, not this
project's own code); the two binned cells use `../bin_qv.py`. All four share
thread count and bgzf compression level, so differences reflect content only.

**Timing sources:**
- *samtools cells* — wall time from `date +%s%3N`; CPU from `/usr/bin/time` (GNU time)
- *bin_qv.py cells* — wall time, CPU, RSS, and reads/sec read directly from the
  per-run `--metrics` TSV (measured inside the script; more accurate than
  external process timing)

### Usage

```bash
# Whole file (default)
benchmarking/benchmark_sizes.sh --input sample.hifi_reads.bam --threads 8

# Fast: fixed first-200k-read subset, clean up BAMs afterward
benchmarking/benchmark_sizes.sh --input sample.hifi_reads.bam \
    --subsample 200000 --outdir benchmark_results --clean
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i, --input` | (required) | Input BAM/UBAM |
| `-o, --outdir` | `./benchmark_results` | Output directory |
| `-s, --subsample N` | `0` (whole file) | Benchmark the first N reads only |
| `-t, --threads N` | `8` | Threads for samtools and bin_qv.py |
| `--clean` | off | Delete generated BAMs after measuring, keep TSV + logs |

### Output

`<outdir>/benchmark.tsv` with one row per cell:

| Column | Description |
|---|---|
| `variant` | `raw_keep`, `raw_strip`, `bin_keep`, `bin_strip` |
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

> **Disk note:** in whole-file mode the four output BAMs can each approach the
> input size — budget up to ~4× the input in free space, or use `--clean`.
