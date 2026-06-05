"""Independent test oracle for the QV bin mapping.

This deliberately re-implements the documented 7-bin table from scratch rather
than importing it from bin_qv, so the binning tests check the implementation
against an *independent* source of truth instead of against itself.
"""

# (low, high, binned_value) — straight from the PacBio CCS documentation table
# in bin_qv.py's module docstring.
_BIN_TABLE = [
    (0,  6,  3),
    (7,  13, 10),
    (14, 19, 17),
    (20, 24, 22),
    (25, 29, 27),
    (30, 39, 35),
    (40, 93, 40),
]


def expected_bin(q: int) -> int:
    """Return the binned Phred for a raw Phred `q` (values > 93 clamp to Q40)."""
    for lo, hi, binned in _BIN_TABLE:
        if lo <= q <= hi:
            return binned
    return 40
