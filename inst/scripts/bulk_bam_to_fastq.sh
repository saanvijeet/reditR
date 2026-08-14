#!/usr/bin/env bash
# =============================================================================
# bulk_bam_to_fastq.sh  —  reditR pipeline template
# =============================================================================
# Purpose:
#   Convert a bulk RNA-seq BAM file to FASTQ for input to SPRINT.
#   Only needed if you start from a BAM file. If you already have FASTQs,
#   skip this and run SPRINT directly.
#
# Usage (PBS/HPC — submit once per sample):
#   qsub -v BAM=/path/to/sample.bam,FASTQ=/path/to/output.fastq \
#        bulk_bam_to_fastq.sh
#
# Library layout is detected from the first alignment record's FLAG. For a
# paired-end BAM the output is split into <output>_1.fastq and <output>_2.fastq;
# for single-end it is written to <output> unchanged. Override the detection
# with PAIRED=1 or PAIRED=0 if the BAM has an unusual or mixed header.
#
# Passing -fq alone to bamToFastq on a paired-end BAM writes both mates into a
# single file with no mate information, which SPRINT would then process as
# single-end. The split below is what makes paired input usable.
#
# Dependencies:
#   - samtools  (tested with 1.19)
#   - bedtools  (for bamToFastq; tested with 2.31.1)
#
# After this completes:
#   Run SPRINT on the output FASTQ(s) -- "-1 <r1> -2 <r2>" for paired-end,
#   "-1 <fastq>" for single-end -- then pipe results through
#   extracting_read_counts.sh to produce the count table for reditR.
#
# NOTE: this template is not exercised by the reditR test suite, which has no
# BAM fixtures. Verify the output read counts against the BAM before relying
# on it in a pipeline.
# =============================================================================

#PBS -N bulk_bam_to_fastq
#PBS -l select=1:ncpus=4:mem=16gb
#PBS -l walltime=4:00:00
#PBS -j oe

set -euo pipefail

# ─── LOAD MODULES ─────────────────────────────────────────────────────────────
# Adjust module names for your HPC environment.
# module load SAMtools/<VERSION>
# module load BEDTools/<VERSION>
# ──────────────────────────────────────────────────────────────────────────────

# ─── CONFIGURE ────────────────────────────────────────────────────────────────
SAM=samtools    # full path if samtools is not on $PATH
CPU=4           # should match PBS ncpus above
# ──────────────────────────────────────────────────────────────────────────────

# BAM and FASTQ are passed at submission time via qsub -v
: "${BAM:?BAM variable not set — pass via qsub -v BAM=...}"
: "${FASTQ:?FASTQ variable not set — pass via qsub -v FASTQ=...}"

SORTED_BAM="${BAM%.bam}_namesorted.bam"

# Layout detection: bit 0x1 of the FLAG marks a paired read. Reading one
# record is enough -- a BAM is paired or it is not. pipefail is disabled
# inside the subshell only, so samtools taking SIGPIPE from head is not an
# error here.
if [[ -z "${PAIRED:-}" ]]; then
    FIRST_FLAG=$(set +o pipefail; ${SAM} view "${BAM}" | head -n 1 | cut -f2)
    if [[ -z "${FIRST_FLAG}" ]]; then
        echo "ERROR: no alignment records in ${BAM}" >&2
        exit 1
    fi
    PAIRED=$(( FIRST_FLAG & 1 ))
fi

echo "[$(date)] Sorting BAM by read name..."
${SAM} sort -n -@ "${CPU}" "${BAM}" -o "${SORTED_BAM}"

if (( PAIRED )); then
    FASTQ1="${FASTQ%.fastq}_1.fastq"
    FASTQ2="${FASTQ%.fastq}_2.fastq"
    echo "[$(date)] Paired-end BAM: converting to ${FASTQ1} / ${FASTQ2}..."
    bamToFastq -i "${SORTED_BAM}" -fq "${FASTQ1}" -fq2 "${FASTQ2}"
    OUTS="${FASTQ1} ${FASTQ2}"
    NEXT="sprint main -1 ${FASTQ1} -2 ${FASTQ2} ..."
else
    echo "[$(date)] Single-end BAM: converting to ${FASTQ}..."
    bamToFastq -i "${SORTED_BAM}" -fq "${FASTQ}"
    OUTS="${FASTQ}"
    NEXT="sprint main -1 ${FASTQ} ..."
fi

rm "${SORTED_BAM}"

echo "[$(date)] Done. Written: ${OUTS}"
echo "Next step: ${NEXT}"
