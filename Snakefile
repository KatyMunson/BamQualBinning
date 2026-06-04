# =============================================================================
# Snakefile — UBAM QV Binning
#
# Remaps per-base quality scores in PacBio UBAM files to the standard
# Revio/CCS 7-bin scheme using workflow/scripts/bin_qv.py.
#
# Steps:
#   1. bin_qv  — remap per-base QVs for each sample UBAM
#   2. index   — samtools index on the output BAM
#
# Usage: snakemake -s Snakefile --configfile config.yaml --cores <N>
# See README.md for full documentation and cluster submission instructions.
#
# Author:  KM
# Created: 2026-06
# =============================================================================

import os

configfile: "config.yaml"

wildcard_constraints:
    sample="[^/]+",

# ---------------------------------------------------------------------------
# Manifest parsing  (tab-delimited: sample\tinput_bam)
# ---------------------------------------------------------------------------
SAMPLES = []
BAMS    = {}

with open(config["manifest"]) as fh:
    for lineno, line in enumerate(fh, 1):
        if line.strip() and not line.startswith("#"):
            fields = line.strip().split("\t")
            if len(fields) < 2:
                raise ValueError(
                    f"Manifest line {lineno} has {len(fields)} field(s), expected ≥ 2 "
                    f"(sample<TAB>bam_path)"
                )
            sample, bam = fields[0], fields[1]
            SAMPLES.append(sample)
            BAMS[sample] = bam

if not SAMPLES:
    raise ValueError("No samples found in manifest!")


# ---------------------------------------------------------------------------
# Target rule
# ---------------------------------------------------------------------------
rule all:
    input:
        expand("results/{sample}/{sample}.qvbin.bam.bai", sample=SAMPLES),


# ---------------------------------------------------------------------------
# Rule: bin_qv
# ---------------------------------------------------------------------------
rule bin_qv:
    input:
        bam=lambda wc: BAMS[wc.sample],
    output:
        bam="results/{sample}/{sample}.qvbin.bam",
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
    shell:
        """
        python {workflow.basedir}/scripts/bin_qv.py \
            --input   {input.bam} \
            --output  {output.bam} \
            --threads {threads} \
            --log     {log} \
            2>> {log}
        """


# ---------------------------------------------------------------------------
# Rule: index
# ---------------------------------------------------------------------------
rule index:
    input:
        bam="results/{sample}/{sample}.qvbin.bam",
    output:
        bai="results/{sample}/{sample}.qvbin.bam.bai",
    log:
        "results/logs/{sample}/index.log",
    threads: config["resources"]["index"]["threads"]
    resources:
        mem=lambda wildcards, attempt: config["resources"]["index"]["mem"] * attempt,
        hrs=config["resources"]["index"]["hrs"],
    conda:
        "envs/ubam_qvbin.yaml"
    envmodules:
        "samtools/1.18",
    shell:
        """
        samtools index -@ {threads} {input.bam} 2> {log}
        """
