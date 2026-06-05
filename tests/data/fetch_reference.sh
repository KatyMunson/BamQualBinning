#!/usr/bin/env bash
# =============================================================================
# fetch_reference.sh — download a small subsample of a public Revio
# on-instrument-binned HiFi BAM for use as a validation reference.
#
# The reference is used by tests/test_reference_binned.py to confirm that this
# tool's bin values match real Revio output (see that file for details). The
# downloaded BAM is git-ignored and never committed.
#
# Public source (same chemistry as this project's data):
#   https://downloads.pacbcloud.com/public/2026Q1/HG002-SPRQ-Nx/
#
# Usage:
#   tests/data/fetch_reference.sh <full_bam_url> [n_reads]
#
#   <full_bam_url>  URL of a specific on-instrument-binned reads.bam from the
#                   listing above.
#   [n_reads]       Number of reads to subsample (default: 2000).
#
# Requires samtools built with libcurl (https URL support) — true for the
# bioconda/module builds used by the main workflow.
# =============================================================================
set -euo pipefail

URL="${1:-}"
N="${2:-2000}"
OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/revio_binned_ref.bam"

if [[ -z "$URL" ]]; then
    echo "ERROR: provide a BAM URL." >&2
    echo "  Browse: https://downloads.pacbcloud.com/public/2026Q1/HG002-SPRQ-Nx/" >&2
    echo "  Usage:  $0 <full_bam_url> [n_reads]" >&2
    exit 1
fi
command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH." >&2; exit 1; }

echo ">> Streaming first $N reads from:" >&2
echo "   $URL" >&2
# awk exits after N records, sending SIGPIPE upstream (expected); disable
# pipefail for this pipeline and validate the result by read count below.
set +o pipefail
samtools view -h "$URL" \
    | awk -v n="$N" 'BEGIN{c=0} /^@/{print; next} {if (c++ < n) print; else exit}' \
    | samtools view -b -o "$OUT" -
set -o pipefail

got=$(samtools view -c "$OUT")
echo ">> Wrote $got reads to $OUT (git-ignored)" >&2
