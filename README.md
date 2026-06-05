# ubam-qv-bin

Remaps per-base quality scores in PacBio UBAM files to the standard
**Revio/CCS 7-bin scheme**, matching the output of `ccs --qv-binning`.

## QV bin table

| Input range | Bin mean | ASCII |
|-------------|----------|-------|
| [0, 6]      | Q3       | `$`   |
| [7, 13]     | Q10      | `+`   |
| [14, 19]    | Q17      | `2`   |
| [20, 24]    | Q22      | `7`   |
| [25, 29]    | Q27      | `<`   |
| [30, 39]    | Q35      | `D`   |
| [40, 93]    | Q40      | `I`   |

---

## Repository layout

```
BamQualBinning/
├── Snakefile
├── config.yaml
├── bin_qv.py          ← core remapping script
├── manifest.tsv       ← you create this (see Quick start)
└── README.md
```

---

## Quick start

### 1. Populate the manifest

Edit `manifest.tsv` — one sample per line, tab-delimited:

```
sample01    /net/eichler/vol28/.../sample01.reads.bam
sample02    /net/eichler/vol28/.../sample02.reads.bam    true
```

| Column | Required | Description |
|--------|----------|-------------|
| 1      | yes      | Sample name |
| 2      | yes      | Path to input BAM/UBAM |
| 3      | no       | Strip kinetics tags (`true`/`yes`/`1` to remove; default: keep) |

Lines beginning with `#` are ignored.

Setting column 3 to `true` removes the PacBio kinetics tags (`ip`, `pw`, `fi`,
`ri`, `fp`, `rp`) from every read, which reduces file size when kinetics
information is no longer needed downstream.

### 2. Edit config.yaml (optional)

Adjust resource allocations under `resources:` if needed.  
`mem` is **GB per slot (per core)**; SGE multiplies by `threads` internally.

### 3. Run

#### Standalone (no cluster)

With `python`, `pysam`, and `samtools` already available on your `PATH`
(e.g. an activated conda/venv environment):

```bash
snakemake -s Snakefile --configfile config.yaml --cores 8
```

> **Note:** The Snakefile's `conda:` directive references
> `envs/ubam_qvbin.yaml`, which is **not** bundled in this repository. To use
> `--use-conda`, first create that file with the dependencies listed under
> [Dependencies](#dependencies).

#### SGE cluster (liger / e002 / e004)

```bash
snakemake -s Snakefile --configfile config.yaml \
    --use-envmodules \
    --drmaa " -l mfree={resources.mem}G -l h_rt={resources.hrs}:00:00 -pe serial {threads}" \
    --jobs 20 \
    --retries 3
```

> **Note:** Verify that `python/3.11` and `pysam/0.22.0` are available on your
> cluster via `module avail`. Update `envmodules:` directives in the Snakefile
> to match what your environment exposes.

### 4. Output

```
results/
  summary_metrics.tsv             ← aggregated run-time table (one row per sample)
  {sample}/
    {sample}.qvbin.bam
    {sample}.qvbin.bam.bai
    {sample}.metrics.tsv          ← per-sample metrics
  logs/
    {sample}/
      bin_qv.log
      index.log
```

---

## Running the script standalone (no Snakemake)

```bash
python bin_qv.py \
    --input   in.reads.bam \
    --output  out.reads.qvbin.bam \
    --threads 8 \
    --log     bin_qv.log \
    --metrics out.metrics.tsv \
    --sample  my_sample \
    [--strip-kinetics]
```

All tags and BAM headers are preserved by default. Pass `--strip-kinetics` to
remove PacBio kinetics tags (`ip`, `pw`, `fi`, `ri`, `fp`, `rp`). The script
streams reads without loading the full file into memory.

`--metrics` and `--sample` are optional; omitting them produces no metrics file
and does not change processing behavior.

---

## Run-time metrics

When `--metrics` is provided (set automatically by the Snakemake workflow),
the script writes a single-row TSV with the following columns:

| Column | Description |
|---|---|
| `sample` | Sample name |
| `n_reads` | Total reads written |
| `n_bases` | Total bases across all reads |
| `n_no_qual` | Reads with no QUAL field (passed through unchanged) |
| `wallclock_sec` | Wall-clock time from open to close (seconds) |
| `cpu_sec` | CPU time (user + system) via `process_time()` |
| `cpu_efficiency` | `cpu_sec / (wallclock_sec × threads)` — near 1.0 = fully parallel |
| `peak_rss_mb` | Peak resident set size (MB) |
| `reads_per_sec` | Read throughput |
| `bases_per_sec` | Base throughput |
| `input_size_bytes` | Input BAM file size |
| `output_size_bytes` | Output BAM file size |
| `threads` | BAM I/O thread count used |
| `strip_kinetics` | Whether kinetics tags were stripped |

The Snakemake workflow aggregates all per-sample files into
`results/summary_metrics.tsv` automatically.

---

## Dependencies

| Tool     | Version  |
|----------|----------|
| Python   | ≥ 3.11   |
| pysam    | ≥ 0.22.0 |
| samtools | ≥ 1.18   |

---

## Notes

- The script passes through reads with no QUAL field unchanged (and warns in the log).
- Phred scores > 93 are clamped to Q40.
- Quality remapping is vectorized via `bytes.translate()` (one C-level call per
  read), so throughput is bound by BAM compression/decompression I/O rather than
  the binning itself — scale `--threads` accordingly.
- Expected **wall-clock time** for a 90 Gbp / 5M-read UBAM: a few minutes on
  cluster scratch, longer on slow NFS, depending on I/O throughput.
- For a ~40% file size reduction matching standard Revio output, combine with
  CRAM conversion after binning:
  ```bash
  samtools view -C -T ref.fa out.qvbin.bam -o out.qvbin.cram
  ```
