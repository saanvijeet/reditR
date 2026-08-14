# Mouse dehydration dataset — barcode-to-cluster annotation builders

The two R scripts in this directory are the ones that produced the cluster
annotations used in the mouse analysis (Dumas et al. 2020 renal endothelial
cells, ArrayExpress **E-MTAB-8145**). They are copied **verbatim**, with no
tidying, genericisation or path changes, so that a citation points at the code
that generated the published numbers rather than at a cleaned-up equivalent.

| Script | Libraries |
|---|---|
| `build_mec_annotations_full.R` | MEC1–MEC5 |
| `build_cec_gec_annotations.R`  | CEC1, CEC2, CEC3, CEC5, GEC1–GEC5 |

Run order: CellRanger &rarr; these scripts &rarr; `../../scrna_preprocessing.sh`
(or the per-dataset PBS equivalent) &rarr; SPRINT per cluster.

## What they do

1. Read the Single Cell Expression Atlas clustering table
   (`E-MTAB-8145.clusters.tsv`), a wide matrix of 9 resolution rows &times;
   36,348 cell columns.
2. Locate the row flagged `sel.K == "TRUE"` and stop if there is not exactly
   one. **`sel.K` marks the clustering SCEA displays by default, which is the
   output of Scanpy's Louvain algorithm at its default resolution of 1.0** — it
   is a software default, not a resolution selected as best-supported for this
   experiment. For E-MTAB-8145 that resolution yields K = 20. Note that K is
   *derived* from `sel.K` here, never hardcoded.
3. Subset to one library's cells by the `SAMEA…-` accession prefix.
4. Reconcile barcode conventions: strip the accession prefix and append
   CellRanger's `-1` GEM-well suffix.
5. Intersect against that library's CellRanger
   `filtered_feature_bc_matrix/barcodes.tsv.gz`, keeping only cells that have
   both an atlas cluster label and passed local cell calling.
6. Write a headerless two-column TSV (`barcode`, `cluster_<n>`) for the
   splitter.

## No minimum cluster-size filter

**These scripts apply no cluster-size threshold.** Every intersected barcode is
written out.

This differs from `../../build_splitter_annotation.sh`, which sits one level up
and carries `MIN_CELLS=30`. That shell script is the **MEC5 single-library
smoke test** (`SAMPLE_PREFIX="SAMEA11354569"`, paths under `mec5_smoketest/`);
it did not produce the analysed annotations. Do not cite it as the method for
the full dataset — doing so would assert a filter that never ran.

The consequence is visible downstream: the analysed dataset contains 96
cluster-level pseudobulk samples, and the number of detected editing sites per
sample ranges from 1 to 1,238 (median 17; 39 of 96 below 10 sites). The
smallest — e.g. `CEC5_cluster_11`, `GEC2_cluster_13`, `MEC2_cluster_16` — carry
a single site on 10–13 reads. Anyone reproducing this analysis should decide
deliberately whether to keep that behaviour or add a threshold; it is recorded
here rather than silently corrected.

## Libraries built vs libraries analysed

`build_cec_gec_annotations.R` builds annotations for **CEC3**, but CEC3 does not
appear in the final analysis, which uses 13 libraries: MEC1–5, CEC1, CEC2, CEC5,
GEC1–5. CEC4 is in neither. The downstream sample list is fixed in
`extract_counts_all_ec.sh`.

## Portability

Both scripts hardcode absolute paths — `$EPHEMERAL/mec_dehydration/…` for
CellRanger output and `$HOME/msc_prj/test_data_mouse/` for the clusters TSV —
and an explicit accession-to-library map. They are a record of what was run, not
a portable utility. To reuse them, edit `WORK_BASE`, `CR_BASE`, `CLUSTERS_TSV`
and the `samea` vector. The only dependency is **data.table**.

## Reference

Dumas SJ, et al. (2020) *Single-cell RNA sequencing reveals renal endothelium
heterogeneity and metabolic adaptation to water deprivation.* J Am Soc Nephrol.
ArrayExpress accession E-MTAB-8145.
