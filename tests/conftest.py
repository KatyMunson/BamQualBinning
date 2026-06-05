"""Shared pytest fixtures and path setup for the bin_qv test suite."""

import sys
from pathlib import Path

import pytest

# Make the repo root (for `import bin_qv`) and this tests/ dir (for `import
# _helpers`) importable regardless of where pytest is invoked from.
REPO_ROOT = Path(__file__).resolve().parent.parent
TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(TESTS_DIR))

BINQV = REPO_ROOT / "bin_qv.py"


@pytest.fixture
def binqv_path() -> str:
    """Absolute path to bin_qv.py for subprocess invocation."""
    return str(BINQV)


@pytest.fixture
def make_bam():
    """Return a factory that writes a small unaligned BAM.

    Usage:
        make_bam(path, [
            {"name": "r1", "seq": "ACGT", "qual": "IIII", "tags": [("np", 3)]},
            {"name": "r2", "seq": "ACGT", "qual": None},   # no QUAL field
        ])
    """
    pysam = pytest.importorskip("pysam")

    def _make(path, reads):
        header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
        with pysam.AlignmentFile(str(path), "wb", header=header) as out:
            for r in reads:
                a = pysam.AlignedSegment(header)
                a.query_name = r["name"]
                a.flag = 4  # unmapped
                a.query_sequence = r["seq"]
                if r.get("qual") is not None:
                    a.query_qualities = pysam.qualitystring_to_array(r["qual"])
                if r.get("tags"):
                    a.set_tags(r["tags"])
                out.write(a)
        return path

    return _make
