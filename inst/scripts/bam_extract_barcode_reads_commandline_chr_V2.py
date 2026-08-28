# -*- coding: utf-8 -*-
# =============================================================================
# bam_extract_barcode_reads_commandline_chr_V2.py  —  reditR pipeline utility
# =============================================================================
# Provenance:
#   Derived from the barcode read-extraction utility in PPMS
#   (https://github.com/SrivastavaLab-ICL/PPMS), modified here to extract reads
#   genome wide rather than one chromosome at a time: the chromosome argument
#   is optional, and omitting it processes the whole BAM. Included so the
#   single-cell preprocessing steps can be reproduced; the original
#   implementation is not the work of this package's author.
#
# Purpose:
#   Split a multiplexed single-cell BAM file into per-cell-type BAM files
#   using a barcode-to-cell-type annotation table.
#
# Requirements:
#   Python 3 with pysam installed:
#     pip install pysam
#   Input BAM must be indexed (samtools index sample.bam).
#
# Usage:
#   python bam_extract_barcode_reads_commandline_chr_V2.py \
#       <bam> <annotation> <output_prefix> [chromosome]
#
#   bam:           BAM file with CB (cell barcode) tags — must be indexed
#   annotation:    tab-separated file mapping barcode to cell type (no header):
#                    ACGTACGT-1    Cardiomyocyte
#                    TTGGCCAA-1    Fibroblast
#   output_prefix: path prefix for output BAMs
#                  e.g. /path/to/sample1 → /path/to/sample1_Cardiomyocyte.bam
#   chromosome:    optional — restrict to one chromosome (e.g. chr1)
#                  omit to process all chromosomes
#
# Output:
#   One BAM per cell type: <output_prefix>_<CellType>.bam
#   These are used as input to Step 2 of scrna_preprocessing.sh.
# =============================================================================

import pysam
import sys

bamfile = pysam.AlignmentFile(sys.argv[1], "rb")

ANN = {}
cell_types = {}
with open(sys.argv[2]) as f:
    for line in f:
        line = line.rstrip("\n")
        parts = line.split('\t', 2)
        ANN[parts[0]] = parts[1]
        cell_types[parts[1]] = 1

key_cell_types = list(cell_types.keys())

out_files = []
for CT in key_cell_types:
    out = pysam.AlignmentFile(sys.argv[3] + "_" + CT + ".bam", "wb", template=bamfile)
    out_files.append(out)

# If chromosome is given, fetch from that chromosome only; otherwise fetch all reads
fetch_iter = bamfile.fetch(sys.argv[4]) if len(sys.argv) > 4 else bamfile.fetch()

for read in fetch_iter:
    if read.has_tag('CB'):
        CB = read.get_tag('CB')
        ct = ANN.get(CB)
        if ct is not None and ct in key_cell_types:
            out_files[key_cell_types.index(ct)].write(read)

for f in out_files:
    f.close()
bamfile.close()
