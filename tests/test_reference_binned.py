"""Validation against real on-instrument-binned Revio HiFi data.

These tests use a local reference BAM that is NOT committed to the repo (the
public PacBio data is large and not ours to redistribute). They skip cleanly
when no reference is available.

Provide the reference one of two ways:
  * run tests/data/fetch_reference.sh <url> to drop revio_binned_ref.bam into
    tests/data/, or
  * set BINQV_REVIO_REF=/path/to/any/binned/reads.bam

What this validates that the synthetic tests cannot:
  * The 7 canonical bin *values* {3,10,17,22,27,35,40} are exactly what real
    Revio output uses for this chemistry (test_only_canonical_bins).
  * Re-binning already-binned Revio reads is a no-op — i.e. this tool's bin
    targets coincide with Revio's on-instrument bins (test_idempotent). If a
    bin mean were even slightly wrong, some value would change here.

Limitation: already-binned data contains only the 7 output values, so these
tests confirm the bin *means* but not the input-range *boundaries* (no
intermediate Phred values exist in the data to exercise them). Boundaries are
covered by the documented-table tests in tests/test_binning.py.
"""

import os
from pathlib import Path

import pytest

pysam = pytest.importorskip("pysam")

import bin_qv

CANONICAL_BINS = {3, 10, 17, 22, 27, 35, 40}


def _ref_path() -> Path:
    env = os.environ.get("BINQV_REVIO_REF")
    if env:
        return Path(env)
    return Path(__file__).resolve().parent / "data" / "revio_binned_ref.bam"


@pytest.fixture(scope="module")
def reference_bam() -> Path:
    p = _ref_path()
    if not p.exists():
        pytest.skip(
            f"No Revio reference BAM at {p} — set BINQV_REVIO_REF or run "
            f"tests/data/fetch_reference.sh <url>"
        )
    return p


def _iter_quals(path: Path, limit: int | None = None):
    n = 0
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as f:
        for read in f:
            q = read.query_qualities
            if q is not None:
                yield q
            n += 1
            if limit is not None and n >= limit:
                return


def test_only_canonical_bins(reference_bam):
    """Real Revio output must contain only the 7 canonical bin values."""
    seen: set[int] = set()
    for q in _iter_quals(reference_bam):
        seen.update(q)
    assert seen, "no QUAL values found in reference BAM"
    extra = seen - CANONICAL_BINS
    assert not extra, f"reference contains non-canonical QV values: {sorted(extra)}"


def test_idempotent_on_reference(reference_bam):
    """Re-binning already-binned Revio reads must change nothing.

    This ties the tool's bin targets to real on-instrument output: a no-op
    here means remap_quals maps each canonical Revio bin to itself.
    """
    checked = 0
    for q in _iter_quals(reference_bam, limit=5000):
        assert bin_qv.remap_quals(q) == q
        checked += 1
    assert checked > 0, "reference BAM had no reads with QUAL to check"
