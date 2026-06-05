# tests

Correctness tests for `bin_qv.py`, run with [pytest](https://docs.pytest.org/).
Unlike `benchmarking/` (performance, on-demand), these are fast deterministic
checks meant to gate changes.

## Running

```bash
pip install pytest pysam     # if not already available
pytest tests/ -v
```

All tests `importorskip("pysam")`, so they skip cleanly (rather than erroring)
in an environment without pysam.

## Coverage

| File | What it covers |
|------|----------------|
| `test_binning.py` | Core QV mapping — checked against an **independent oracle** (`_helpers.py`), not the code's own table. Exhaustive over all 94 Phred inputs, the >93 clamp, idempotency, no-mutation, empty input, plus hand-verified real-data before/after vectors. Also `write_metrics`. |
| `test_strip_kinetics.py` | `strip_kinetics_tags` removes only `ip/pw/fi/ri/fp/rp`, preserves other tags (incl. the float `sn` B-array that previously crashed `set_tags`), and leaves SEQ/QUAL intact. |
| `test_cli.py` | End-to-end subprocess runs: binning applied per-base, saturated Q93→Q40, no-QUAL reads pass through and count in `n_no_qual`, read names/sequences preserved, `--metrics` contents, and `--strip-kinetics` behavior. |

`_helpers.py` re-implements the documented bin table from scratch so the
binning tests validate the implementation against an independent source of
truth rather than against itself.
