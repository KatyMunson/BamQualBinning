# benchmarking

Performance and file-size benchmarks for the QV-binning pipeline. These are
run on demand against real data — they are **not** correctness tests and do not
gate the workflow.

## `benchmark_sizes.sh`

4-way compressed-size comparison that isolates the two size levers
independently — quality-score binning and kinetics-tag removal:

|              | keep kinetics | strip kinetics |
|--------------|---------------|----------------|
| **raw QVs**  | `raw_keep`    | `raw_strip`    |
| **binned QVs** | `bin_keep`  | `bin_strip`    |

The two raw cells are produced with `samtools view` (a neutral htslib baseline,
not this project's own code); the two binned cells come from `../bin_qv.py`.
All four use the same thread count and bgzf compression level, so the only
variable is content.

### Usage

```bash
# Whole file (default)
benchmarking/benchmark_sizes.sh --input sample.hifi_reads.bam --threads 8

# Fast: fixed first-200k-read subset, all four cells see identical input
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

`<outdir>/benchmark_sizes.tsv` with one row per cell: `variant`, `kinetics`,
`qual`, `bytes`, `MB`, `bytes_per_read`, and `ratio_vs_raw_keep`.

> **Disk note:** in whole-file mode the four output BAMs can each approach the
> input size — budget up to ~4× the input in free space, or use `--clean`.
