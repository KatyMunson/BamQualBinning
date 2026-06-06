"""Tests for custom bin scheme support: load_bins_file, build_translation_table,
and end-to-end CLI behavior with --bins-file."""

import array
import os
import textwrap

import pytest
import pysam

from bin_qv import (
    load_bins_file,
    build_translation_table,
    remap_quals,
    TRANSLATE_TABLE,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REVIO_BINS_TSV = os.path.join(REPO_ROOT, "bins", "pacbio_revio.tsv")

EXPECTED_REVIO_RULES = [
    (0,  6,  3),
    (7,  13, 10),
    (14, 19, 17),
    (20, 24, 22),
    (25, 29, 27),
    (30, 39, 35),
    (40, 93, 40),
]


def write_bins(tmp_path, content):
    p = tmp_path / "bins.tsv"
    p.write_text(textwrap.dedent(content))
    return str(p)


# ---------------------------------------------------------------------------
# load_bins_file — valid inputs
# ---------------------------------------------------------------------------

def test_load_bins_file_valid():
    rules = load_bins_file(REVIO_BINS_TSV)
    assert rules == EXPECTED_REVIO_RULES


def test_load_bins_file_single_phred_bin(tmp_path):
    """A bin spanning exactly one Phred value (lo == hi) must be accepted."""
    path = write_bins(tmp_path, """\
        0\t0\t0
        1\t93\t40
    """)
    rules = load_bins_file(path)
    assert rules[0] == (0, 0, 0)
    assert rules[1] == (1, 93, 40)


def test_load_bins_file_comments_and_blank_lines(tmp_path):
    """Comment and blank lines must be ignored."""
    path = write_bins(tmp_path, """\
        # header comment
        0\t6\t3

        # another comment
        7\t93\t40
    """)
    rules = load_bins_file(path)
    assert len(rules) == 2
    assert rules[0] == (0, 6, 3)


# ---------------------------------------------------------------------------
# load_bins_file — validation failures
# ---------------------------------------------------------------------------

def test_load_bins_file_gap(tmp_path):
    path = write_bins(tmp_path, "0\t5\t3\n7\t13\t10\n")
    with pytest.raises(ValueError, match="gap"):
        load_bins_file(path)


def test_load_bins_file_overlap(tmp_path):
    path = write_bins(tmp_path, "0\t8\t3\n7\t13\t10\n")
    with pytest.raises(ValueError, match="overlap"):
        load_bins_file(path)


def test_load_bins_file_bad_mean_above(tmp_path):
    path = write_bins(tmp_path, "0\t6\t7\n7\t93\t40\n")
    with pytest.raises(ValueError, match="bin_mean"):
        load_bins_file(path)


def test_load_bins_file_bad_mean_below(tmp_path):
    path = write_bins(tmp_path, "5\t10\t3\n11\t93\t40\n")
    with pytest.raises(ValueError, match="bin_mean"):
        load_bins_file(path)


def test_load_bins_file_lo_gt_hi(tmp_path):
    path = write_bins(tmp_path, "10\t5\t7\n")
    with pytest.raises(ValueError, match="lo.*>.*hi"):
        load_bins_file(path)


def test_load_bins_file_empty(tmp_path):
    path = write_bins(tmp_path, "# only comments\n\n")
    with pytest.raises(ValueError, match="empty"):
        load_bins_file(path)


def test_load_bins_file_must_start_at_zero(tmp_path):
    """A scheme that doesn't start at Phred 0 leaves low values uncovered."""
    path = write_bins(tmp_path, "10\t93\t40\n")
    with pytest.raises(ValueError, match="start at Phred 0"):
        load_bins_file(path)


def test_load_bins_file_must_cover_93(tmp_path):
    """A scheme that doesn't reach Phred 93 leaves high values uncovered."""
    path = write_bins(tmp_path, "0\t50\t25\n")
    with pytest.raises(ValueError, match="cover at least Phred 93"):
        load_bins_file(path)


def test_load_bins_file_top_bin_may_exceed_93(tmp_path):
    """A top bin with hi > 93 is allowed and must not crash (capped at 255)."""
    path = write_bins(tmp_path, "0\t6\t3\n7\t300\t40\n")
    rules = load_bins_file(path)
    table = build_translation_table(rules)
    assert table[93] == 40
    assert table[255] == 40   # clamp, no IndexError


def test_load_bins_file_wrong_columns(tmp_path):
    path = write_bins(tmp_path, "0\t6\n")
    with pytest.raises(ValueError, match="3 tab-delimited"):
        load_bins_file(path)


# ---------------------------------------------------------------------------
# build_translation_table — custom rules
# ---------------------------------------------------------------------------

def test_build_table_custom_three_bins():
    """A simple 3-bin scheme maps values correctly across all 256 entries."""
    rules = [(0, 9, 5), (10, 29, 20), (30, 93, 40)]
    table = build_translation_table(rules)
    assert len(table) == 256
    # spot-check each bin
    for q in range(0, 10):
        assert table[q] == 5,  f"Q{q} should map to 5"
    for q in range(10, 30):
        assert table[q] == 20, f"Q{q} should map to 20"
    for q in range(30, 94):
        assert table[q] == 40, f"Q{q} should map to 40"


def test_custom_bins_clamp_to_highest_mean():
    """Phred values 94–255 must clamp to the highest bin's bin_mean."""
    rules = [(0, 6, 3), (7, 50, 30)]   # highest bin_mean = 30
    table = build_translation_table(rules)
    for q in range(94, 256):
        assert table[q] == 30, f"Q{q} should clamp to 30"


def test_build_table_single_phred_bin():
    """A bin with lo == hi must set exactly one table entry."""
    rules = [(7, 7, 7), (8, 93, 40)]
    table = build_translation_table(rules)
    assert table[7] == 7
    assert table[8] == 40
    assert table[6] == 0   # not covered by any bin → 0 (bytearray default)


def test_default_table_unchanged():
    """build_translation_table() with no args must equal the module constant."""
    assert build_translation_table() == TRANSLATE_TABLE


# ---------------------------------------------------------------------------
# remap_quals with custom table
# ---------------------------------------------------------------------------

def test_remap_quals_custom_table():
    rules = [(0, 20, 10), (21, 93, 40)]
    table = build_translation_table(rules)
    quals = array.array("B", [0, 10, 20, 21, 40, 93])
    result = remap_quals(quals, table)
    assert list(result) == [10, 10, 10, 40, 40, 40]


# ---------------------------------------------------------------------------
# CLI end-to-end with --bins-file
# ---------------------------------------------------------------------------

def _make_bam(path, qualities):
    """Write a minimal UBAM with one read at the given Phred qualities."""
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        read = pysam.AlignedSegment(header)
        read.query_name = "read1"
        read.query_sequence = "A" * len(qualities)
        read.query_qualities = array.array("B", qualities)
        read.flag = 4
        bam.write(read)


def test_cli_bins_file(tmp_path):
    """--bins-file remaps qualities per the custom scheme end-to-end."""
    in_bam  = str(tmp_path / "in.bam")
    out_bam = str(tmp_path / "out.bam")
    bins    = str(tmp_path / "bins.tsv")

    # Custom scheme: Q0–20 → 10, Q21–93 → 40
    with open(bins, "w") as fh:
        fh.write("0\t20\t10\n21\t93\t40\n")

    input_quals = [0, 10, 20, 21, 50, 93]
    _make_bam(in_bam, input_quals)

    import subprocess, sys
    result = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "bin_qv.py"),
         "--input", in_bam, "--output", out_bam,
         "--bins-file", bins, "--threads", "1"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    with pysam.AlignmentFile(out_bam, "rb", check_sq=False) as bam:
        reads = list(bam)
    assert len(reads) == 1
    assert list(reads[0].query_qualities) == [10, 10, 10, 40, 40, 40]


def test_cli_bins_file_reference_equals_default(tmp_path):
    """--bins-file bins/pacbio_revio.tsv must produce identical output to the default."""
    in_bam       = str(tmp_path / "in.bam")
    out_default  = str(tmp_path / "out_default.bam")
    out_ref_file = str(tmp_path / "out_ref.bam")

    quals = list(range(94))   # one read covering Q0–Q93
    _make_bam(in_bam, quals)

    import subprocess, sys
    script = os.path.join(REPO_ROOT, "bin_qv.py")

    r1 = subprocess.run(
        [sys.executable, script, "--input", in_bam, "--output", out_default, "--threads", "1"],
        capture_output=True,
    )
    r2 = subprocess.run(
        [sys.executable, script, "--input", in_bam, "--output", out_ref_file,
         "--bins-file", REVIO_BINS_TSV, "--threads", "1"],
        capture_output=True,
    )
    assert r1.returncode == 0
    assert r2.returncode == 0

    with pysam.AlignmentFile(out_default,  "rb", check_sq=False) as b:
        q_default = list(list(b)[0].query_qualities)
    with pysam.AlignmentFile(out_ref_file, "rb", check_sq=False) as b:
        q_ref = list(list(b)[0].query_qualities)

    assert q_default == q_ref


def test_cli_invalid_bins_file_exits_nonzero(tmp_path):
    """An invalid bins file must cause a non-zero exit with an error message."""
    in_bam  = str(tmp_path / "in.bam")
    out_bam = str(tmp_path / "out.bam")
    bad_bins = str(tmp_path / "bad.tsv")

    _make_bam(in_bam, [10, 20])
    with open(bad_bins, "w") as fh:
        fh.write("0\t5\t3\n10\t93\t40\n")   # gap at 6–9

    import subprocess, sys
    result = subprocess.run(
        [sys.executable, os.path.join(REPO_ROOT, "bin_qv.py"),
         "--input", in_bam, "--output", out_bam,
         "--bins-file", bad_bins, "--threads", "1"],
        capture_output=True, text=True,
    )
    assert result.returncode != 0
    assert "gap" in result.stderr.lower() or "gap" in result.stdout.lower()
