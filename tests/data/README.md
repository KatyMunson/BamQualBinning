# tests/data

Holds the (git-ignored) reference BAM used by `tests/test_reference_binned.py`.

Sequencing data is **not** committed — see `.gitignore`. Populate it with:

```bash
tests/data/fetch_reference.sh <bam_url> 2000
```

where `<bam_url>` is an on-instrument-binned HiFi `reads.bam` from
<https://downloads.pacbcloud.com/public/2026Q1/HG002-SPRQ-Nx/>.

Alternatively, point the test at any binned BAM you already have:

```bash
export BINQV_REVIO_REF=/path/to/binned/reads.bam
```
