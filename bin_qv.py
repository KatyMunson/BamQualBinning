#!/usr/bin/env python3
"""
bin_qv.py — Remap per-base quality scores in a PacBio UBAM to the standard
Revio/CCS 7-bin scheme.

QV bin table (from PacBio CCS documentation):
    [ 0,  6] → Q3   (ASCII '$', offset 3)
    [ 7, 13] → Q10  (ASCII '+', offset 10)
    [14, 19] → Q17  (ASCII '2', offset 17)
    [20, 24] → Q22  (ASCII '7', offset 22)
    [25, 29] → Q27  (ASCII '<', offset 27)
    [30, 39] → Q35  (ASCII 'D', offset 35)
    [40, 93] → Q40  (ASCII 'I', offset 40)

Usage:
    python bin_qv.py --input in.bam --output out.bam [--threads N] [--log LOG]
                     [--strip-kinetics] [--metrics METRICS] [--sample SAMPLE]

Notes:
    - Works on coordinate-sorted, queryname-sorted, or unsorted BAMs/UBAMs.
    - All read tags and headers are preserved exactly unless --strip-kinetics
      is used, which removes the ip, pw, fi, ri, fp, and rp tags.
    - pysam writes a BGZF-compressed BAM by default; pipe through
      `samtools view -b` if you need a different format.
    - For very large files, use --threads to enable multi-threaded BAM I/O.
    - Pass --metrics to write a single-row TSV of run-time performance metrics
      suitable for aggregation into a per-run summary table.
"""

import argparse
import array
import csv
import logging
import os
import resource
import sys
import time

import pysam


# ---------------------------------------------------------------------------
# QV bin lookup table
# Build a 94-element list (Phred 0–93) mapping each score to its bin mean.
# ---------------------------------------------------------------------------
def build_lookup() -> array.array:
    """Return an array of length 94 mapping raw Phred → binned Phred."""
    BIN_RULES = [
        (0,  6,  3),
        (7,  13, 10),
        (14, 19, 17),
        (20, 24, 22),
        (25, 29, 27),
        (30, 39, 35),
        (40, 93, 40),
    ]
    lut = array.array("B", [0] * 94)
    for lo, hi, binned in BIN_RULES:
        for q in range(lo, hi + 1):
            lut[q] = binned
    return lut


LUT = build_lookup()

# PacBio kinetics tags: IPD and pulse-width arrays (native and per-strand)
KINETICS_TAGS = frozenset({"ip", "pw", "fi", "ri", "fp", "rp"})


def strip_kinetics_tags(read: pysam.AlignedSegment) -> None:
    """Remove PacBio kinetics tags from a read in-place."""
    tags = [(t, v, tp) for t, v, tp in read.get_tags(with_value_type=True)
            if t not in KINETICS_TAGS]
    read.set_tags(tags)


def remap_quals(quals: array.array) -> array.array:
    """Apply LUT to a pysam quality array (array of uint8) in-place."""
    for i in range(len(quals)):
        q = quals[i]
        quals[i] = LUT[q] if q < 94 else 40   # clamp anything > 93 → Q40
    return quals


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Remap PacBio UBAM per-base QVs to the standard 7-bin CCS scheme."
    )
    p.add_argument("-i", "--input",   required=True,  help="Input BAM/UBAM path (or '-' for stdin)")
    p.add_argument("-o", "--output",  required=True,  help="Output BAM path (or '-' for stdout)")
    p.add_argument("-t", "--threads", type=int, default=4,
                   help="Threads for BAM compression I/O (default: 4)")
    p.add_argument("--log", default=None,
                   help="Log file path (default: stderr)")
    p.add_argument("--log-interval", type=int, default=100_000,
                   help="Log progress every N reads (default: 100000)")
    p.add_argument("--strip-kinetics", action="store_true", default=False,
                   help="Remove PacBio kinetics tags (ip, pw, fi, ri, fp, rp) "
                        "from each read (default: keep tags)")
    p.add_argument("--metrics", default=None, metavar="PATH",
                   help="Write per-sample run-time metrics to this TSV file "
                        "(default: no metrics file written)")
    p.add_argument("--sample", default=None, metavar="NAME",
                   help="Sample name written to the metrics TSV "
                        "(default: output filename stem)")
    return p.parse_args()


def setup_logging(log_path: str | None) -> logging.Logger:
    logger = logging.getLogger("bin_qv")
    logger.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s  %(levelname)s  %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    handler = (logging.FileHandler(log_path) if log_path
               else logging.StreamHandler(sys.stderr))
    handler.setFormatter(fmt)
    logger.addHandler(handler)
    return logger


def _gb(n_bytes: int) -> str:
    return f"{n_bytes / 1e9:.2f} GB"


def write_metrics(path: str, metrics: dict) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(metrics.keys()), delimiter="\t")
        w.writeheader()
        w.writerow(metrics)


def main() -> None:
    args = parse_args()
    log  = setup_logging(args.log)

    sample_name = args.sample or os.path.splitext(os.path.basename(args.output))[0]

    log.info("Input          : %s", args.input)
    log.info("Output         : %s", args.output)
    log.info("Sample         : %s", sample_name)
    log.info("Threads        : %d", args.threads)
    log.info("Strip kinetics : %s", args.strip_kinetics)
    if args.metrics:
        log.info("Metrics file   : %s", args.metrics)

    t0     = time.perf_counter()
    cpu_t0 = time.process_time()

    n_reads   = 0
    n_bases   = 0
    n_no_qual = 0

    open_mode_in  = "rb"
    open_mode_out = "wb"

    with pysam.AlignmentFile(
        args.input,  open_mode_in,  check_sq=False, threads=args.threads
    ) as bam_in, \
    pysam.AlignmentFile(
        args.output, open_mode_out, header=bam_in.header, threads=args.threads
    ) as bam_out:

        for read in bam_in:
            n_reads += 1
            n_bases += read.query_length or 0

            quals = read.query_qualities   # array.array("B", ...) or None
            if quals is None:
                n_no_qual += 1
            else:
                read.query_qualities = remap_quals(quals)

            if args.strip_kinetics:
                strip_kinetics_tags(read)

            bam_out.write(read)

            if n_reads % args.log_interval == 0:
                wall = time.perf_counter() - t0
                cpu  = time.process_time() - cpu_t0
                rate = n_reads / wall
                log.info("  %10d reads  %8.1f wall-sec  %8.1f cpu-sec  %.0f reads/sec",
                         n_reads, wall, cpu, rate)

    wall_sec = time.perf_counter() - t0
    cpu_sec  = time.process_time() - cpu_t0

    # ru_maxrss is in KB on Linux
    peak_rss_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    input_bytes  = os.path.getsize(args.input)  if args.input  != "-" else 0
    output_bytes = os.path.getsize(args.output) if args.output != "-" else 0

    reads_per_sec = n_reads / wall_sec if wall_sec > 0 else 0.0
    bases_per_sec = n_bases / wall_sec if wall_sec > 0 else 0.0
    cpu_efficiency = (cpu_sec / (wall_sec * args.threads)
                      if wall_sec > 0 and args.threads > 0 else 0.0)
    size_ratio = output_bytes / input_bytes if input_bytes > 0 else 0.0

    log.info("=== Summary ===")
    log.info("  Reads processed : %d", n_reads)
    log.info("  Bases processed : %d", n_bases)
    log.info("  Reads w/o QUAL  : %d", n_no_qual)
    log.info("  Wallclock time  : %.2f sec", wall_sec)
    log.info("  CPU time        : %.2f sec", cpu_sec)
    log.info("  CPU efficiency  : %.4f  (threads=%d)", cpu_efficiency, args.threads)
    log.info("  Peak RSS        : %.1f MB", peak_rss_mb)
    log.info("  Throughput      : %.0f reads/sec  |  %.0f kb/sec",
             reads_per_sec, bases_per_sec / 1000)
    log.info("  Input size      : %s", _gb(input_bytes))
    log.info("  Output size     : %s", _gb(output_bytes))
    log.info("  Size ratio      : %.3fx", size_ratio)

    if n_no_qual:
        log.warning("%d reads had no QUAL field and were passed through unchanged.", n_no_qual)

    if args.metrics:
        metrics = {
            "sample":            sample_name,
            "n_reads":           n_reads,
            "n_bases":           n_bases,
            "n_no_qual":         n_no_qual,
            "wallclock_sec":     round(wall_sec,  2),
            "cpu_sec":           round(cpu_sec,   2),
            "cpu_efficiency":    round(cpu_efficiency, 4),
            "peak_rss_mb":       round(peak_rss_mb, 1),
            "reads_per_sec":     round(reads_per_sec, 1),
            "bases_per_sec":     round(bases_per_sec, 1),
            "input_size_bytes":  input_bytes,
            "output_size_bytes": output_bytes,
            "threads":           args.threads,
            "strip_kinetics":    args.strip_kinetics,
        }
        write_metrics(args.metrics, metrics)
        log.info("Metrics written : %s", args.metrics)


if __name__ == "__main__":
    main()
