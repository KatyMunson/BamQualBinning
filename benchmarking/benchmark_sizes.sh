#!/usr/bin/env bash
# =============================================================================
# benchmark_sizes.sh — 4-way file-size comparison for QV binning + kinetics
#
# Builds all four combinations from a single input UBAM and reports their
# compressed sizes so the effect of (a) quality-score binning and (b) kinetics
# tag removal can be compared independently:
#
#                  | keep kinetics | strip kinetics
#     -------------+---------------+----------------
#     raw  QVs     |  raw_keep     |  raw_strip      ← samtools (neutral baseline)
#     binned QVs   |  bin_keep     |  bin_strip      ← bin_qv.py
#
# The two raw cells are produced with `samtools view` (htslib, bgzf level 6) so
# the baseline does not depend on this project's own code. The two binned cells
# come from ../bin_qv.py. All four use the same thread count and compression
# level, so size differences reflect *content* (binned-vs-raw quals,
# kinetics-present-vs-absent) rather than tooling.
#
# Usage:
#     benchmarking/benchmark_sizes.sh --input in.bam [options]
#
# Options:
#     -i, --input PATH     Input BAM/UBAM (required)
#     -o, --outdir DIR     Output directory (default: ./benchmark_results)
#     -s, --subsample N    Benchmark the first N reads only (default: 0 = whole
#                          file). A fixed prefix keeps all four cells comparable
#                          and runs in a couple of minutes instead of the full
#                          per-cell processing time.
#     -t, --threads N      Threads for samtools and bin_qv.py (default: 8)
#         --clean          Delete the generated BAMs after measuring, keeping
#                          only the TSV table and logs.
#     -h, --help           Show this help.
#
# Output:
#     <outdir>/raw_keep.bam, raw_strip.bam, bin_keep.bam, bin_strip.bam
#     <outdir>/bin_keep.log, bin_strip.log
#     <outdir>/benchmark_sizes.tsv   ← the comparison table
#
# Requirements: samtools >= 1.18, python with pysam (same as the main workflow).
#
# NOTE: in whole-file mode the four output BAMs can each approach the size of
# the input, so budget up to ~4x the input size in free disk (use --clean to
# reclaim it after the table is printed).
# =============================================================================

set -euo pipefail

# Kinetics tags removed by `samtools --remove-tag` and `bin_qv.py --strip-kinetics`.
KINETICS_TAGS="ip,pw,fi,ri,fp,rp"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINQV="$SCRIPT_DIR/../bin_qv.py"

INPUT=""
OUTDIR="./benchmark_results"
SUBSAMPLE=0
THREADS=8
CLEAN=0

usage() { sed -n '2,/^# ===.*$/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)     INPUT="$2";     shift 2 ;;
        -o|--outdir)    OUTDIR="$2";    shift 2 ;;
        -s|--subsample) SUBSAMPLE="$2"; shift 2 ;;
        -t|--threads)   THREADS="$2";   shift 2 ;;
        --clean)        CLEAN=1;        shift   ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "$INPUT" ]]; then
    echo "ERROR: --input is required." >&2
    usage
    exit 1
fi
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file not found: $INPUT" >&2
    exit 1
fi
if [[ ! -f "$BINQV" ]]; then
    echo "ERROR: could not locate bin_qv.py at $BINQV" >&2
    exit 1
fi

command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH." >&2; exit 1; }

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Optionally build a fixed read-prefix subset so all four cells see identical
# input content. Taking the first N reads is deterministic and avoids a full
# pass over a large file.
# ---------------------------------------------------------------------------
WORK_INPUT="$INPUT"
if [[ "$SUBSAMPLE" -gt 0 ]]; then
    WORK_INPUT="$OUTDIR/subset.bam"
    echo ">> Building subset: first $SUBSAMPLE reads -> $WORK_INPUT" >&2
    # awk exits after N records, sending SIGPIPE upstream to samtools; that is
    # expected, so pipefail is disabled for this one pipeline and the result is
    # validated by a read count below.
    set +o pipefail
    samtools view -h -@ "$THREADS" "$INPUT" \
        | awk -v n="$SUBSAMPLE" 'BEGIN{c=0} /^@/{print; next} {if (c++ < n) print; else exit}' \
        | samtools view -b -@ "$THREADS" -o "$WORK_INPUT" -
    set -o pipefail
    got=$(samtools view -c -@ "$THREADS" "$WORK_INPUT")
    echo ">> Subset written: $got reads" >&2
fi

# ---------------------------------------------------------------------------
# Generate the four variants.
# ---------------------------------------------------------------------------
echo ">> [1/4] raw quals, keep kinetics  (samtools re-encode)" >&2
samtools view -b -@ "$THREADS" -o "$OUTDIR/raw_keep.bam" "$WORK_INPUT"

echo ">> [2/4] raw quals, strip kinetics (samtools --remove-tag)" >&2
samtools view -b -@ "$THREADS" --remove-tag "$KINETICS_TAGS" \
    -o "$OUTDIR/raw_strip.bam" "$WORK_INPUT"

echo ">> [3/4] binned quals, keep kinetics  (bin_qv.py)" >&2
python "$BINQV" --input "$WORK_INPUT" --output "$OUTDIR/bin_keep.bam" \
    --threads "$THREADS" --log "$OUTDIR/bin_keep.log"

echo ">> [4/4] binned quals, strip kinetics (bin_qv.py --strip-kinetics)" >&2
python "$BINQV" --input "$WORK_INPUT" --output "$OUTDIR/bin_strip.bam" \
    --threads "$THREADS" --strip-kinetics --log "$OUTDIR/bin_strip.log"

# ---------------------------------------------------------------------------
# Measure and tabulate. All four cells share the same read count, so count it
# once from the smallest output for speed.
# ---------------------------------------------------------------------------
NREADS=$(samtools view -c -@ "$THREADS" "$OUTDIR/bin_strip.bam")

VARIANTS=(raw_keep raw_strip bin_keep bin_strip)
declare -A KIN=( [raw_keep]=keep [raw_strip]=strip [bin_keep]=keep [bin_strip]=strip )
declare -A QUAL=( [raw_keep]=raw [raw_strip]=raw [bin_keep]=binned [bin_strip]=binned )

TSV="$OUTDIR/benchmark_sizes.tsv"
printf "variant\tkinetics\tqual\tbytes\tMB\tbytes_per_read\tratio_vs_raw_keep\n" > "$TSV"

base=$(stat -c %s "$OUTDIR/raw_keep.bam")
for v in "${VARIANTS[@]}"; do
    sz=$(stat -c %s "$OUTDIR/$v.bam")
    mb=$(awk -v b="$sz" 'BEGIN{printf "%.1f", b/1048576}')
    bpr=$(awk -v b="$sz" -v n="$NREADS" 'BEGIN{printf "%.1f", b/n}')
    ratio=$(awk -v b="$sz" -v a="$base" 'BEGIN{printf "%.3f", b/a}')
    printf "%s\t%s\t%s\t%d\t%s\t%s\t%s\n" \
        "$v" "${KIN[$v]}" "${QUAL[$v]}" "$sz" "$mb" "$bpr" "$ratio" >> "$TSV"
done

echo >&2
echo "===== File-size comparison ($NREADS reads) =====" >&2
column -t -s $'\t' "$TSV"
echo >&2
echo "Table written to: $TSV" >&2

if [[ "$SUBSAMPLE" -eq 0 ]]; then
    orig=$(stat -c %s "$INPUT")
    orig_mb=$(awk -v b="$orig" 'BEGIN{printf "%.1f", b/1048576}')
    echo "Reference — original input as delivered: $orig bytes ($orig_mb MB)" >&2
fi

# ---------------------------------------------------------------------------
# Optional cleanup of the (potentially large) generated BAMs.
# ---------------------------------------------------------------------------
if [[ "$CLEAN" -eq 1 ]]; then
    echo ">> --clean: removing generated BAMs (keeping TSV + logs)" >&2
    rm -f "$OUTDIR"/raw_keep.bam "$OUTDIR"/raw_strip.bam \
          "$OUTDIR"/bin_keep.bam "$OUTDIR"/bin_strip.bam
    [[ "$SUBSAMPLE" -gt 0 ]] && rm -f "$WORK_INPUT"
fi
