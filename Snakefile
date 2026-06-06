# =============================================================================
# Snakefile — UBAM QV Binning
#
# Remaps per-base quality scores in PacBio UBAM files to the standard
# Revio/CCS 7-bin scheme using bin_qv.py.
#
# Steps:
#   1. bin_qv            — remap per-base QVs for each sample UBAM
#   2. aggregate_metrics — combine per-sample run-time metrics into one table
#
# Usage: snakemake -s Snakefile --configfile config.yaml --cores <N>
# See README.md for full documentation and cluster submission instructions.
#
# Author:  KM
# Created: 2026-06
# =============================================================================

import os

configfile: "config.yaml"

# Run the target and the lightweight metrics aggregation on the submit host
# rather than the cluster. Snakemake requires `run:` rules to be local under
# cluster execution, and the aggregation is a millisecond pure-Python task that
# should not consume a cluster slot (nor needs the DRMAA resource template).
localrules: all, aggregate_metrics

wildcard_constraints:
    sample="[^/]+",

# ---------------------------------------------------------------------------
# Manifest parsing  (tab-delimited: sample\tinput_bam)
# ---------------------------------------------------------------------------
SAMPLES         = []
BAMS            = {}
STRIP_KINETICS  = {}

_TRUTHY = {"true", "yes", "1"}

with open(config["manifest"]) as fh:
    for lineno, line in enumerate(fh, 1):
        if line.strip() and not line.startswith("#"):
            fields = line.strip().split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Manifest line {lineno} has {len(fields)} field(s), expected ≥ 2 "
                    f"(sample<TAB>bam_path[<TAB>strip_kinetics])"
                )
            sample, bam = fields[0], fields[1]
            strip = fields[2].strip().lower() in _TRUTHY if len(fields) >= 3 else False
            SAMPLES.append(sample)
            BAMS[sample]           = bam
            STRIP_KINETICS[sample] = strip

if not SAMPLES:
    raise ValueError("No samples found in manifest!")


# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        expand("results/{sample}/{sample}.qvbin.bam", sample=SAMPLES),
        "results/summary_metrics.tsv",


# ---------------------------------------------------------------------------
# Rule: bin_qv
# ---------------------------------------------------------------------------
rule bin_qv:
    input:
        bam=lambda wc: BAMS[wc.sample],
    output:
        bam="results/{sample}/{sample}.qvbin.bam",
        metrics="results/{sample}/{sample}.metrics.tsv",
    log:
        "results/logs/{sample}/bin_qv.log",
    threads: config["resources"]["bin_qv"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["bin_qv"]["mem"] * attempt,
        hrs=config["resources"]["bin_qv"]["hrs"],
    conda:
        "envs/ubam_qvbin.yaml"
    envmodules:
        "python/3.11",
        "pysam/0.22.0",
    params:
        strip_kinetics=lambda wc: "--strip-kinetics" if STRIP_KINETICS[wc.sample] else "",
        bins_file=f"--bins-file {config['bins_file']}" if config.get("bins_file") else "",
    shell:
        """
        python {workflow.basedir}/bin_qv.py \
            --input   {input.bam} \
            --output  {output.bam} \
            --threads {threads} \
            --log     {log} \
            --metrics {output.metrics} \
            --sample  {wildcards.sample} \
            {params.strip_kinetics} \
            {params.bins_file} \
            2>> {log}
        """


# ---------------------------------------------------------------------------
# Rule: aggregate_metrics
# ---------------------------------------------------------------------------
rule aggregate_metrics:
    input:
        expand("results/{sample}/{sample}.metrics.tsv", sample=SAMPLES),
    output:
        tsv="results/summary_metrics.tsv",
    run:
        import csv
        rows = []
        for f in input:
            with open(f) as fh:
                rows.extend(list(csv.DictReader(fh, delimiter="\t")))
        if rows:
            with open(output.tsv, "w", newline="") as out:
                w = csv.DictWriter(out, fieldnames=rows[0].keys(), delimiter="\t")
                w.writeheader()
                w.writerows(rows)
