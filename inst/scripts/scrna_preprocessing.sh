#!/usr/bin/env bash
# =============================================================================
# scrna_preprocessing.sh  —  reditR pipeline template
# =============================================================================
# Purpose:
#   Single-cell preprocessing pipeline:
#     1. Split a multiplexed BAM by cell-type barcode annotation
#     2. Convert each per-cell-type BAM to FASTQ
#     3. Run SPRINT on each FASTQ to detect RNA editing sites
#
#   Output directories are named {SAMPLE}_{CellType}_ and are ready for
#   extracting_read_counts.sh (also in inst/scripts/).
#
# Usage (PBS/HPC — submit once per sample):
#   qsub -v BAM=/path/to/sample.bam,\
#            ANNOTATION=/path/to/barcodes.tsv,\
#            SAMPLE=case_rep1 \
#        scrna_preprocessing.sh
#
# Annotation file format (tab-separated, no header):
#   ACGTACGT-1    Cardiomyocyte
#   TTGGCCAA-1    Fibroblast
#   ...
#
# Dependencies:
#   - Python 3 + pysam  (for BAM splitting)
#   - samtools           (tested with 1.19)
#   - BWA                (tested with 0.7.18)
#   - bedtools           (for bamToFastq; tested with 2.31.1)
#   - SPRINT             (https://github.com/jumphone/SPRINT)
#
# The barcode-splitting Python script is shipped alongside this template:
#   bam_extract_barcode_reads_commandline_chr_V2.py
# Find its installed path with:
#   Rscript -e "system.file('scripts',
#     'bam_extract_barcode_reads_commandline_chr_V2.py', package='reditR')"
# =============================================================================

#PBS -N scrna_preprocessing
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=12:00:00
#PBS -j oe

set -euo pipefail

# ─── LOAD MODULES ─────────────────────────────────────────────────────────────
# Adjust module names for your HPC environment.
# module load SAMtools/<VERSION>
# module load BWA/<VERSION>
# module load BEDTools/<VERSION>
# ──────────────────────────────────────────────────────────────────────────────

# ─── CONFIGURE THESE PATHS ────────────────────────────────────────────────────
OUT_DIR=<REPLACE: /path/to/your/output/directory>

# Path to the barcode-splitting script (shipped with reditR):
PYTHON_SCRIPT=<REPLACE: /path/to/bam_extract_barcode_reads_commandline_chr_V2.py>
# Tip: Rscript -e "system.file('scripts', 'bam_extract_barcode_reads_commandline_chr_V2.py', package='reditR')"

FASTA=<REPLACE: /path/to/reference/genome.fa>           # GRCh38 recommended
ANN_GTF=<REPLACE: /path/to/transcript/annotation.gtf>  # Ensembl GTF

BWA=bwa         # full path if BWA is not on $PATH
SAM=samtools    # full path if samtools is not on $PATH

# SPRINT's -c: trim the first N bp of each read (default 0). This is read
# trimming, NOT a coverage or cluster-size threshold. Set to 0 to disable.
# Do not confuse it with SPRINT's cluster-size defaults (-csrg / -cshp, both
# 5 when no repeat annotation is supplied via -rp), which coincidentally
# share the value 5.
CUT=5
CPU=8           # should match PBS ncpus above
# ──────────────────────────────────────────────────────────────────────────────

# BAM, ANNOTATION, and SAMPLE are passed at submission time via qsub -v
: "${BAM:?BAM variable not set — pass via qsub -v BAM=...}"
: "${ANNOTATION:?ANNOTATION variable not set — pass via qsub -v ANNOTATION=...}"
: "${SAMPLE:?SAMPLE variable not set — pass via qsub -v SAMPLE=...}"

WORK_DIR="${OUT_DIR}/${SAMPLE}_scrna_tmp"
mkdir -p "${WORK_DIR}"

echo "[$(date)] Step 1: Splitting BAM by cell type..."
python "${PYTHON_SCRIPT}" "${BAM}" "${ANNOTATION}" "${WORK_DIR}/${SAMPLE}"

echo "[$(date)] Step 2: Converting per-cell-type BAMs to FASTQs..."
for bam_file in "${WORK_DIR}/${SAMPLE}"_*.bam; do
    base=$(basename "${bam_file}" .bam)
    ${SAM} sort -n -@ "${CPU}" "${bam_file}" -o "${WORK_DIR}/${base}_sorted.bam"
    bamToFastq -i "${WORK_DIR}/${base}_sorted.bam" -fq "${WORK_DIR}/${base}.fastq"
    rm "${WORK_DIR}/${base}_sorted.bam"
done

echo "[$(date)] Step 3: Running SPRINT prepare (reference index)..."
sprint prepare -t "${ANN_GTF}" "${FASTA}" "${BWA}"

echo "[$(date)] Step 4: Running SPRINT per cell type..."
for fastq_file in "${WORK_DIR}/${SAMPLE}"_*.fastq; do
    cell_type=$(basename "${fastq_file}" .fastq | sed "s/^${SAMPLE}_//")
    sprint_out="${OUT_DIR}/${SAMPLE}_${cell_type}_"
    mkdir -p "${sprint_out}"
    sprint main -c "${CUT}" -p "${CPU}" -1 "${fastq_file}" \
        "${FASTA}" "${sprint_out}" "${BWA}" "${SAM}"
    echo "[$(date)]   Done: ${sprint_out}"
done

echo "[$(date)] All done. Output directories in ${OUT_DIR}/"
echo "Next step: run extracting_read_counts.sh with a samples.txt listing these directories."
