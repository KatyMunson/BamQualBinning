"""Pure-logic tests for the QV binning core (no BAM I/O required)."""

import csv
from array import array

import pytest

pytest.importorskip("pysam")  # bin_qv imports pysam at module load

import bin_qv
from _helpers import expected_bin


# ---------------------------------------------------------------------------
# Translation table
# ---------------------------------------------------------------------------
def test_table_is_256_bytes():
    assert len(bin_qv.TRANSLATE_TABLE) == 256


def test_table_matches_oracle_for_every_phred():
    """Exhaustive: all 94 valid Phred inputs map per the documented table."""
    for q in range(94):
        assert bin_qv.TRANSLATE_TABLE[q] == expected_bin(q), f"Phred {q} mis-binned"


def test_table_clamps_out_of_range_to_q40():
    """Indices 94–255 (Phred > 93) all clamp to Q40."""
    for q in range(94, 256):
        assert bin_qv.TRANSLATE_TABLE[q] == 40


def test_bin_values_are_self_mapping():
    """Each of the 7 bin means maps to itself (binning is idempotent)."""
    for b in (3, 10, 17, 22, 27, 35, 40):
        assert bin_qv.TRANSLATE_TABLE[b] == b


# ---------------------------------------------------------------------------
# remap_quals
# ---------------------------------------------------------------------------
def test_remap_exhaustive_single_values():
    for q in range(94):
        out = bin_qv.remap_quals(array("B", [q]))
        assert list(out) == [expected_bin(q)]


def test_remap_returns_uint8_array():
    out = bin_qv.remap_quals(array("B", [40, 93, 0]))
    assert isinstance(out, array)
    assert out.typecode == "B"


def test_remap_empty():
    assert bin_qv.remap_quals(array("B", [])) == array("B", [])


def test_remap_does_not_mutate_input():
    inp = array("B", [93, 0, 30])
    _ = bin_qv.remap_quals(inp)
    assert inp == array("B", [93, 0, 30])  # original untouched


def test_remap_is_idempotent():
    """Binning an already-binned string returns it unchanged."""
    raw = array("B", list(range(94)))
    once = bin_qv.remap_quals(raw)
    twice = bin_qv.remap_quals(once)
    assert once == twice


def test_clamp_above_93():
    out = bin_qv.remap_quals(array("B", [94, 100, 200, 255]))
    assert list(out) == [40, 40, 40, 40]


# ---------------------------------------------------------------------------
# Regression: real before/after quality strings that were spot-checked by hand
# (Phred+33 ASCII). These anchor the behavior observed on actual PacBio data.
# ---------------------------------------------------------------------------
SPOT_CHECK_VECTORS = [
    (
        r"uX~~~}~~~{st|T~~xc~~MjZ~xT~~whune~|fh~d~k~r~~\~`~~~~_~C~rq~I~~`~QOZ^~vYqO~wTD~~w_zxc~nY~~~}~1{~]~_~a",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIDIIIIIIIIIIIIIIIIIIIIIDIIIIIIIIIIIIIII2IIIIIII",
    ),
    (
        r"6;^`^J0V\XIU>_TKITKFOF7FG4@``_KVLBF8]]Q<[ENI6?6UI5aaaTG:[H.H'LRH>G9]TP.<@KQKT?bQKDa^?a-Z[QFI.HT)Z]L_",
        "7<IIII2IIIII<IIIIIIDID7DD2DIIIIIIDD7III<IDII7D7II7IIIID<ID+D$IID<D7III+<DIIIIDIIIDIIDI+IIIDI+DI+IIII",
    ),
]


@pytest.mark.parametrize("raw_str,binned_str", SPOT_CHECK_VECTORS)
def test_spot_check_vectors(raw_str, binned_str):
    raw = array("B", [ord(c) - 33 for c in raw_str])
    expected = array("B", [ord(c) - 33 for c in binned_str])
    assert bin_qv.remap_quals(raw) == expected


# ---------------------------------------------------------------------------
# write_metrics
# ---------------------------------------------------------------------------
def test_write_metrics_roundtrip(tmp_path):
    path = tmp_path / "m.tsv"
    metrics = {"sample": "s1", "n_reads": 100, "cpu_sec": 1.5, "strip_kinetics": False}
    bin_qv.write_metrics(str(path), metrics)
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["sample"] == "s1"
    assert rows[0]["n_reads"] == "100"
    assert rows[0]["strip_kinetics"] == "False"
