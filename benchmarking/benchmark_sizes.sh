#!/usr/bin/env bash
# =============================================================================
# benchmark_sizes.sh — 4-way file-size and timing comparison
#
# Builds all four combinations from a single input UBAM and reports their
# compressed sizes and processing times so the effect of (a) quality-score
# binning and (b) kinetics-tag removal can be compared independently:
#
#                  | keep kinetics | strip kinetics
#     -------------+---------------+----------------
#     raw  QVs     |  raw_keep     |  raw_strip      ← samtools (neutral baseline)
#     binned QVs   |  bin_keep     |  bin_strip      ← bin_qv.py
#
# The two raw cells use `samtools view` (htslib, bgzf level 6) so the baseline
# does not depend on this project's own code. The two binned cells come from
# ../bin_qv.py. All four use the same thread count and compression level, so
# differences reflect *content* only.
#
# Timing sources:
#   samtools cells  — wall time from date +%s%3N; CPU from /usr/bin/time
#   bin_qv.py cells — wall, CPU, RSS, and reads/sec from the metrics TSV
#                     written by --metrics (measured inside the script, more
#                     accurate than external process timing)
#
# Usage:
#     benchmarking/benchmark_sizes.sh --input in.bam [options]
#
# Options:
#     -i, --input PATH     Input BAM/UBAM (required)
#     -o, --outdir DIR     Output directory (default: ./benchmark_results)
#     -s, --subsample N    Benchmark the first N reads only (default: 0 = whole
#                          file). A fixed prefix keeps all four cells comparable
#                          and runs in minutes instead of hours.
#     -t, --threads N      Threads for samtools and bin_qv.py (default: 8)
#         --clean          Delete generated BAMs after measuring; keep TSV + logs.
#     -h, --help           Show this help.
#
# Output:
#     <outdir>/raw_keep.bam, raw_strip.bam, bin_keep.bam, bin_strip.bam
#     <outdir>/bin_keep.bam, bin_keep.log, bin_keep.metrics.tsv
#     <outdir>/bin_strip.bam, bin_strip.log, bin_strip.metrics.tsv
#     <outdir>/benchmark.tsv   ← the combined size + timing table
#
# Requirements: samtools >= 1.18, python with pysam, /usr/bin/time (GNU time,
# for CPU measurement of samtools cells; wall time is always captured).
#
# NOTE: in whole-file mode the four output BAMs can each approach the size of
# the input, so budget up to ~4x the input size in free disk (use --clean to
# reclaim after the table is printed).
# =============================================================================

set -euo pipefail

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
    echo "ERROR: --input is required." >&2; usage; exit 1
fi
if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: input file not found: $INPUT" >&2; exit 1
fi
if [[ ! -f "$BINQV" ]]; then
    echo "ERROR: could not locate bin_qv.py at $BINQV" >&2; exit 1
fi
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH." >&2; exit 1; }

# Check for GNU time (needed for CPU measurement of samtools cells).
HAS_GNU_TIME=0
if /usr/bin/time -f "%e" true 2>/dev/null; then
    HAS_GNU_TIME=1
fi

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Associative arrays to hold per-variant timing, filled during each run.
# ---------------------------------------------------------------------------
declare -A WALL_SEC=()
declare -A CPU_SEC_A=()
declare -A READS_PER_SEC=()
declare -A PEAK_RSS_MB=()

# ---------------------------------------------------------------------------
# Helper: run a samtools command, capture wall and CPU times.
# Usage: run_samtools <variant> <samtools args…>
# ---------------------------------------------------------------------------
run_samtools() {
    local variant="$1"; shift
    local t0 t1 elapsed
    local timefile="$OUTDIR/${variant}.time"

    t0=$(date +%s%3N)
    if [[ "$HAS_GNU_TIME" -eq 1 ]]; then
        # %e = wall sec, %U = user CPU sec, %S = system CPU sec
        /usr/bin/time -f "%e %U %S" -o "$timefile" samtools "$@"
        read -r _wall user sys < "$timefile"
        CPU_SEC_A[$variant]=$(awk -v u="$user" -v s="$sys" \
            'BEGIN{printf "%.2f", u+s}')
    else
        samtools "$@"
        CPU_SEC_A[$variant]="N/A"
    fi
    t1=$(date +%s%3N)

    elapsed=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", (b-a)/1000}')
    WALL_SEC[$variant]="$elapsed"
    PEAK_RSS_MB[$variant]="N/A"
}

# ---------------------------------------------------------------------------
# Helper: parse a bin_qv.py metrics TSV (header + 1 data row) for one field.
# Usage: parse_metric <tsv_file> <column_name>
# ---------------------------------------------------------------------------
parse_metric() {
    local tsv="$1" col="$2"
    awk -F'\t' -v c="$col" \
        'NR==1{for(i=1;i<=NF;i++) if($i==c){idx=i}} NR==2{print $idx}' "$tsv"
}

# ---------------------------------------------------------------------------
# Build subset if requested.
# ---------------------------------------------------------------------------
WORK_INPUT="$INPUT"
if [[ "$SUBSAMPLE" -gt 0 ]]; then
    WORK_INPUT="$OUTDIR/subset.bam"
    echo ">> Building subset: first $SUBSAMPLE reads -> $WORK_INPUT" >&2
    # awk exits early; SIGPIPE to upstream samtools is expected, so pipefail
    # is disabled for this pipeline and the result is validated below.
    set +o pipefail
    samtools view -h -@ "$THREADS" "$INPUT" \
        | awk -v n="$SUBSAMPLE" \
              'BEGIN{c=0} /^@/{print; next} {if (c++ < n) print; else exit}' \
        | samtools view -b -@ "$THREADS" -o "$WORK_INPUT" -
    set -o pipefail
    got=$(samtools view -c -@ "$THREADS" "$WORK_INPUT")
    echo ">> Subset written: $got reads" >&2
fi

# ---------------------------------------------------------------------------
# Generate the four variants.
# ---------------------------------------------------------------------------
echo ">> [1/4] raw quals, keep kinetics  (samtools re-encode)" >&2
run_samtools raw_keep \
    view -b -@ "$THREADS" -o "$OUTDIR/raw_keep.bam" "$WORK_INPUT"

echo ">> [2/4] raw quals, strip kinetics (samtools --remove-tag)" >&2
run_samtools raw_strip \
    view -b -@ "$THREADS" --remove-tag "$KINETICS_TAGS" \
    -o "$OUTDIR/raw_strip.bam" "$WORK_INPUT"

echo ">> [3/4] binned quals, keep kinetics  (bin_qv.py)" >&2
python "$BINQV" \
    --input   "$WORK_INPUT" \
    --output  "$OUTDIR/bin_keep.bam" \
    --threads "$THREADS" \
    --log     "$OUTDIR/bin_keep.log" \
    --metrics "$OUTDIR/bin_keep.metrics.tsv" \
    --sample  bin_keep

echo ">> [4/4] binned quals, strip kinetics (bin_qv.py --strip-kinetics)" >&2
python "$BINQV" \
    --input          "$WORK_INPUT" \
    --output         "$OUTDIR/bin_strip.bam" \
    --threads        "$THREADS" \
    --strip-kinetics \
    --log            "$OUTDIR/bin_strip.log" \
    --metrics        "$OUTDIR/bin_strip.metrics.tsv" \
    --sample         bin_strip

# Pull timing from bin_qv.py metrics TSVs (measured inside the script).
for v in bin_keep bin_strip; do
    m="$OUTDIR/${v}.metrics.tsv"
    WALL_SEC[$v]=$(parse_metric "$m" wallclock_sec)
    CPU_SEC_A[$v]=$(parse_metric "$m" cpu_sec)
    READS_PER_SEC[$v]=$(parse_metric "$m" reads_per_sec)
    PEAK_RSS_MB[$v]=$(parse_metric "$m" peak_rss_mb)
done

# ---------------------------------------------------------------------------
# Count reads once from the smallest output.
# ---------------------------------------------------------------------------
NREADS=$(samtools view -c -@ "$THREADS" "$OUTDIR/bin_strip.bam")

# For samtools cells: derive reads/sec from NREADS / wall time.
for v in raw_keep raw_strip; do
    READS_PER_SEC[$v]=$(awk -v n="$NREADS" -v w="${WALL_SEC[$v]}" \
        'BEGIN{printf "%.1f", (w>0) ? n/w : 0}')
done

# ---------------------------------------------------------------------------
# Build the combined TSV.
# ---------------------------------------------------------------------------
VARIANTS=(raw_keep raw_strip bin_keep bin_strip)
declare -A KIN=(  [raw_keep]=keep  [raw_strip]=strip [bin_keep]=keep  [bin_strip]=strip )
declare -A QUAL=( [raw_keep]=raw   [raw_strip]=raw   [bin_keep]=binned [bin_strip]=binned )

# Run-parameter provenance, written as #INFO comment lines at the top of the
# TSV so the table is self-documenting (downstream parsers should skip lines
# starting with '#').
if [[ "$SUBSAMPLE" -gt 0 ]]; then
    mode="subsample (first $SUBSAMPLE reads)"
else
    mode="whole file"
fi
samtools_ver=$(samtools --version 2>/dev/null | head -1)
python_ver=$(python --version 2>&1)

TSV="$OUTDIR/benchmark.tsv"
{
    printf "#INFO\tdate\t%s\n"            "$(date '+%Y-%m-%d %H:%M:%S')"
    printf "#INFO\tinput\t%s\n"           "$INPUT"
    printf "#INFO\tmode\t%s\n"            "$mode"
    printf "#INFO\treads_benchmarked\t%s\n" "$NREADS"
    printf "#INFO\tthreads\t%s\n"          "$THREADS"
    printf "#INFO\tkinetics_tags\t%s\n"    "$KINETICS_TAGS"
    printf "#INFO\tsamtools\t%s\n"         "$samtools_ver"
    printf "#INFO\tpython\t%s\n"           "$python_ver"
    printf "#INFO\tgnu_time_cpu\t%s\n"     "$([[ "$HAS_GNU_TIME" -eq 1 ]] && echo yes || echo "no (samtools cpu_sec = N/A)")"
} > "$TSV"

printf "variant\tkinetics\tqual\tbytes\tMB\tbytes_per_read\tratio_vs_raw_keep\twallclock_sec\tcpu_sec\treads_per_sec\tpeak_rss_mb\n" \
    >> "$TSV"

base=$(stat -c %s "$OUTDIR/raw_keep.bam")
for v in "${VARIANTS[@]}"; do
    sz=$(stat -c %s "$OUTDIR/$v.bam")
    mb=$(awk    -v b="$sz"   'BEGIN{printf "%.1f", b/1048576}')
    bpr=$(awk   -v b="$sz" -v n="$NREADS" 'BEGIN{printf "%.1f", b/n}')
    ratio=$(awk -v b="$sz" -v a="$base" 'BEGIN{printf "%.3f", b/a}')
    printf "%s\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$v" "${KIN[$v]}" "${QUAL[$v]}" \
        "$sz" "$mb" "$bpr" "$ratio" \
        "${WALL_SEC[$v]}" "${CPU_SEC_A[$v]}" \
        "${READS_PER_SEC[$v]}" "${PEAK_RSS_MB[$v]}" \
        >> "$TSV"
done

echo >&2
echo "===== Benchmark: $NREADS reads =====" >&2
if command -v column >/dev/null; then
    column -t -s $'\t' "$TSV"
else
    cat "$TSV"
fi
echo >&2
echo "Table written to: $TSV" >&2

if [[ "$SUBSAMPLE" -eq 0 ]]; then
    orig=$(stat -c %s "$INPUT")
    orig_mb=$(awk -v b="$orig" 'BEGIN{printf "%.1f", b/1048576}')
    echo "Reference — original input as delivered: $orig bytes ($orig_mb MB)" >&2
fi

# ---------------------------------------------------------------------------
# Optional cleanup.
# ---------------------------------------------------------------------------
if [[ "$CLEAN" -eq 1 ]]; then
    echo ">> --clean: removing generated BAMs (keeping TSV, logs, metrics)" >&2
    rm -f "$OUTDIR"/raw_keep.bam "$OUTDIR"/raw_strip.bam \
          "$OUTDIR"/bin_keep.bam "$OUTDIR"/bin_strip.bam
    [[ "$SUBSAMPLE" -gt 0 ]] && rm -f "$WORK_INPUT"
fi
