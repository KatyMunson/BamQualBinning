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
ubam-qv-bin/
├── Snakefile
├── config.yaml
├── manifest.tsv              ← edit this
└── workflow/
    ├── envs/
    │   └── ubam_qvbin.yaml   ← conda environment
    └── scripts/
        └── bin_qv.py         ← core remapping script
```

---

## Quick start

### 1. Populate the manifest

Edit `manifest.tsv` — one sample per line, tab-delimited:

```
sample01    /net/eichler/vol28/.../sample01.reads.bam
sample02    /net/eichler/vol28/.../sample02.reads.bam
```

Lines beginning with `#` are ignored.

### 2. Edit config.yaml (optional)

Adjust resource allocations under `resources:` if needed.  
`mem` is **GB per slot (per core)**; SGE multiplies by `threads` internally.

### 3. Run

#### Standalone (no cluster)

```bash
snakemake -s Snakefile --configfile config.yaml \
    --use-conda --cores 8
```

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

Results are written to:

```
results/{sample}/{sample}.qvbin.bam
results/{sample}/{sample}.qvbin.bam.bai
results/logs/{sample}/bin_qv.log
results/logs/{sample}/index.log
```

---

## Running the script standalone (no Snakemake)

```bash
python workflow/scripts/bin_qv.py \
    --input  in.reads.bam \
    --output out.reads.qvbin.bam \
    --threads 8 \
    --log    bin_qv.log
```

All tags and BAM headers are preserved. The script streams reads without
loading the full file into memory.

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
- Expected **wall-clock time** for a 90 Gbp / 5M-read UBAM: 3–6 min on cluster
  scratch, 10–20 min on slow NFS, depending on I/O throughput.
- For a ~40% file size reduction matching standard Revio output, combine with
  CRAM conversion after binning:
  ```bash
  samtools view -C -T ref.fa out.qvbin.bam -o out.qvbin.cram
  ```
