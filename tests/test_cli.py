"""End-to-end tests: run bin_qv.py as a subprocess on a tiny BAM."""

import csv
import subprocess
import sys
from array import array

import pytest

pysam = pytest.importorskip("pysam")

from _helpers import expected_bin


def _read_all(path):
    with pysam.AlignmentFile(str(path), "rb", check_sq=False) as f:
        return list(f)


# A 10-base read whose QUAL spans several bins, a saturated read, and a read
# with no QUAL field (must pass through and be counted in n_no_qual).
QUAL_R1 = "6;:?@ABII~"   # Phred 21,26,25,30,31,32,33,40,40,93
READS = [
    {"name": "r1", "seq": "ACGTACGTAC", "qual": QUAL_R1,
     "tags": [("ip", array("B", list(range(10)))),
              ("pw", array("B", list(range(10)))),
              ("sn", array("f", [1.0, 2.0, 3.0, 4.0])),
              ("np", 7)]},
    {"name": "r2", "seq": "ACGT", "qual": "~~~~"},
    {"name": "r3", "seq": "ACGTA", "qual": None},
]


def _run(binqv_path, inp, out, *extra):
    subprocess.run(
        [sys.executable, binqv_path, "-i", str(inp), "-o", str(out), *extra],
        check=True,
    )


def test_binning_applied_and_passthrough(tmp_path, make_bam, binqv_path):
    inp = make_bam(tmp_path / "in.bam", READS)
    out = tmp_path / "out.bam"
    _run(binqv_path, inp, out)

    recs = _read_all(out)
    assert len(recs) == 3

    # r1: every base binned per the oracle
    expected = [expected_bin(q) for q in pysam.qualitystring_to_array(QUAL_R1)]
    assert list(recs[0].query_qualities) == expected

    # r2: saturated Q93 -> all Q40
    assert list(recs[1].query_qualities) == [40, 40, 40, 40]

    # r3: no QUAL field passes through unchanged
    assert recs[2].query_qualities is None

    # data integrity: names and sequences untouched
    assert [r.query_name for r in recs] == ["r1", "r2", "r3"]
    assert recs[0].query_sequence == "ACGTACGTAC"


def test_metrics_file_contents(tmp_path, make_bam, binqv_path):
    inp = make_bam(tmp_path / "in.bam", READS)
    out = tmp_path / "out.bam"
    metrics = tmp_path / "m.tsv"
    _run(binqv_path, inp, out, "--metrics", str(metrics), "--sample", "s1")

    row = list(csv.DictReader(open(metrics), delimiter="\t"))[0]
    assert row["sample"] == "s1"
    assert int(row["n_reads"]) == 3
    assert int(row["n_no_qual"]) == 1
    assert int(row["n_bases"]) == 10 + 4 + 5
    assert row["strip_kinetics"] == "False"


def test_keeps_kinetics_by_default(tmp_path, make_bam, binqv_path):
    inp = make_bam(tmp_path / "in.bam", READS)
    out = tmp_path / "out.bam"
    _run(binqv_path, inp, out)
    names = {t for t, _ in _read_all(out)[0].get_tags()}
    assert {"ip", "pw", "sn", "np"} <= names


def test_strip_kinetics_flag(tmp_path, make_bam, binqv_path):
    inp = make_bam(tmp_path / "in.bam", READS)
    out = tmp_path / "out.bam"
    metrics = tmp_path / "m.tsv"
    _run(binqv_path, inp, out, "--strip-kinetics", "--metrics", str(metrics))

    rec = _read_all(out)[0]
    names = {t for t, _ in rec.get_tags()}
    assert "ip" not in names and "pw" not in names
    assert {"sn", "np"} <= names           # non-kinetics preserved
    assert list(rec.query_qualities)[0] == 22   # binning still applied (Q21->Q22)

    row = list(csv.DictReader(open(metrics), delimiter="\t"))[0]
    assert row["strip_kinetics"] == "True"
