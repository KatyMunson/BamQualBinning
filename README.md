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

Bin boundaries and means follow PacBio's authoritative specification:
<https://ccs.how/faq/qv-binning.html>.

---

## Repository layout

```
BamQualBinning/
├── Snakefile
├── config.yaml
├── bin_qv.py                  ← core remapping script
├── manifest.tsv               ← you create this (see Quick start)
├── bins/
│   └── pacbio_revio.tsv       ← reference copy of the default bin scheme
├── envs/
│   └── ubam_qvbin.yaml        ← conda environment (for --use-conda)
├── benchmarking/
│   └── benchmark_sizes.sh     ← 4-way file-size comparison
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

```bash
snakemake -s Snakefile --configfile config.yaml \
    --use-conda --cores 8
```

The `conda:` directive builds the environment from `envs/ubam_qvbin.yaml` on
first run. If `python`, `pysam`, and `samtools` are already on your `PATH`, you
can drop `--use-conda` and just run `snakemake -s Snakefile --configfile
config.yaml --cores 8`.

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
    {sample}.metrics.tsv          ← per-sample metrics
  logs/
    {sample}/
      bin_qv.log
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
    [--strip-kinetics] \
    [--bins-file bins/pacbio_revio.tsv]
```

All tags and BAM headers are preserved by default. Pass `--strip-kinetics` to
remove PacBio kinetics tags (`ip`, `pw`, `fi`, `ri`, `fp`, `rp`). The script
streams reads without loading the full file into memory.

`--metrics` and `--sample` are optional; omitting them produces no metrics file
and does not change processing behavior.

`--bins-file` is optional; omitting it uses the default PacBio Revio 7-bin scheme.

---

## Custom bin schemes

Pass `--bins-file PATH` to use an alternate QV binning scheme. The file is a
tab-delimited TSV with three columns (no header required; `#` lines are comments):

```
# lo	hi	bin_mean
0	6	3
7	13	10
...
```

| Column | Description |
|--------|-------------|
| `lo` | Lower bound of the input Phred range (inclusive) |
| `hi` | Upper bound of the input Phred range (inclusive) |
| `bin_mean` | Output Phred value for all scores in `[lo, hi]` |

**Rules:**
- Bins must be contiguous (no gaps, no overlaps).
- `bin_mean` must be within `[lo, hi]`.
- Phred scores above the highest `hi` are clamped to that bin's `bin_mean`.
- At least one bin must be defined.

`bins/pacbio_revio.tsv` is included as a reference copy of the default scheme —
use it as a starting point for custom schemes.

To use a custom scheme in the Snakemake workflow, set `bins_file` in `config.yaml`:

```yaml
bins_file: "bins/my_custom_scheme.tsv"
```

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

| Tool     | Version  | Used by |
|----------|----------|---------|
| Python   | ≥ 3.11   | workflow + script |
| pysam    | ≥ 0.22.0 | workflow + script |
| samtools | ≥ 1.18   | `benchmarking/` only |

---

## Notes

- The script passes through reads with no QUAL field unchanged (and warns in the log).
- Phred scores > 93 are clamped to the highest defined bin's mean (Q40 for the default scheme; configurable with `--bins-file`).
- Quality remapping is vectorized via `bytes.translate()` (one C-level call per
  read), so throughput is bound by BAM compression/decompression I/O rather than
  the binning itself — scale `--threads` accordingly.
- For a ~40% file size reduction matching standard Revio output, combine with
  CRAM conversion after binning:
  ```bash
  samtools view -C -T ref.fa out.qvbin.bam -o out.qvbin.cram
  ```

---

## Performance & file-size estimates

> **Note:** the figures below are **extrapolated** from a 200k-read benchmark
> subset (linear scaling to 2.5M reads) and will be replaced with measured
> numbers from a full production run. Wall-clock time is I/O-bound and scales
> with your storage subsystem (fast scratch vs. busy NFS), not just read count;
> CPU-time scales more reliably. See `benchmarking/` to reproduce.

### Run time (bin_qv.py, ~2.5M reads, multi-threaded I/O)

| Operation | Wall-clock | CPU-time |
|---|---|---|
| Default binning, keep kinetics | ~37 min | ~93 min |
| Default binning, strip kinetics | ~13 min | ~45 min |

Strip-kinetics runs are ~3× faster than keep, because the output is ~6× smaller
and the cost is dominated by bytes written, not the binning math.

### File size — cost of retaining more quality than the default scheme

The PacBio Revio default binning is the standard. Retaining *more* quality
information (a finer custom scheme, or full unbinned QVs) costs extra disk —
roughly a **fixed absolute amount per file**, independent of kinetics handling.
Per 2.5M-read file:

| What you keep | keep kinetics | strip kinetics |
|---|---|---|
| Default Revio binning (standard) | 161 GB | 24 GB |
| Finer custom scheme (e.g. more high-end bins) | +~4 GB (+2%) | +~4 GB (+16%) |
| Full base quality (no binning) | +~23 GB (+14%) | +~23 GB (+93%) |

The percentage looks small in a kinetics-retained workflow but nearly doubles
the file once kinetics are stripped — so the "is full quality worth the disk?"
trade-off is sharpest when kinetics are already removed.
