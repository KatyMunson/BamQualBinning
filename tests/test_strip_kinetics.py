"""Tests for kinetics-tag stripping on synthetic reads."""

from array import array

import pytest

pysam = pytest.importorskip("pysam")

import bin_qv


def _make_read():
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    a = pysam.AlignedSegment(header)
    a.query_name = "r1"
    a.flag = 4
    a.query_sequence = "ACGTACGT"
    a.query_qualities = pysam.qualitystring_to_array("IIIIIIII")
    a.set_tags([
        # kinetics tags (should all be removed)
        ("ip", array("B", [1, 2, 3, 4, 5, 6, 7, 8])),
        ("pw", array("B", [8, 7, 6, 5, 4, 3, 2, 1])),
        ("fi", array("B", [1] * 8)),
        ("ri", array("B", [2] * 8)),
        ("fp", array("B", [3] * 8)),
        ("rp", array("B", [4] * 8)),
        # non-kinetics tags (should be preserved)
        ("np", 5),
        ("rq", 0.999),
        ("sn", array("f", [10.0, 11.0, 12.0, 13.0])),  # float B-array (regression)
        ("zm", 12345),
    ])
    return a


def test_strip_removes_all_kinetics_tags():
    read = _make_read()
    bin_qv.strip_kinetics_tags(read)
    names = {t for t, _ in read.get_tags()}
    assert names == {"np", "rq", "sn", "zm"}


def test_strip_preserves_other_tag_values():
    read = _make_read()
    bin_qv.strip_kinetics_tags(read)
    d = dict(read.get_tags())
    assert d["np"] == 5
    assert abs(d["rq"] - 0.999) < 1e-6
    assert list(d["sn"]) == [10.0, 11.0, 12.0, 13.0]
    assert d["zm"] == 12345


def test_strip_preserves_sequence_and_quals():
    """Stripping tags must not disturb SEQ or QUAL."""
    read = _make_read()
    bin_qv.strip_kinetics_tags(read)
    assert read.query_sequence == "ACGTACGT"
    assert list(read.query_qualities) == [40] * 8


def test_strip_is_safe_when_no_kinetics_present():
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6"}})
    a = pysam.AlignedSegment(header)
    a.query_name = "r1"
    a.flag = 4
    a.query_sequence = "ACGT"
    a.query_qualities = pysam.qualitystring_to_array("IIII")
    a.set_tags([("np", 3), ("zm", 9)])
    bin_qv.strip_kinetics_tags(a)
    assert {t for t, _ in a.get_tags()} == {"np", "zm"}
