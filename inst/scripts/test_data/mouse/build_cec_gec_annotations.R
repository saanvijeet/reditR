# build_cec_gec_annotations.R
# Build per-sample annotation TSVs for CEC1,2,3,5 and GEC1-5.
# Run after cellranger_cec_gec.pbs finishes all 9 jobs.

library(data.table)

EPHEMERAL <- Sys.getenv("EPHEMERAL")
WORK_BASE <- file.path(EPHEMERAL, "mec_dehydration")
ANN_DIR   <- file.path(WORK_BASE, "annotations_full")
CR_BASE   <- file.path(WORK_BASE, "cellranger_out_full")

CLUSTERS_TSV <- file.path(
  Sys.getenv("HOME"),
  "msc_prj/test_data_mouse/E-MTAB-8145.clusters.tsv"
)

dir.create(ANN_DIR, showWarnings = FALSE, recursive = TRUE)

samea <- c(
  CEC1 = "SAMEA11354556",
  CEC2 = "SAMEA11354557",
  CEC3 = "SAMEA11354558",
  CEC5 = "SAMEA11354559",
  GEC1 = "SAMEA11354560",
  GEC2 = "SAMEA11354561",
  GEC3 = "SAMEA11354562",
  GEC4 = "SAMEA11354563",
  GEC5 = "SAMEA11354564"
)

cat("Reading clusters TSV...\n")
clusters_raw <- fread(CLUSTERS_TSV, header = TRUE, sep = "\t")

sel_row <- which(clusters_raw$sel.K == "TRUE")
if (length(sel_row) != 1) stop("Expected exactly one sel.K=TRUE row")
K_val <- clusters_raw$K[sel_row]
cat("Using K =", K_val, "(sel.K=TRUE, row", sel_row, ")\n\n")

cell_cols   <- setdiff(names(clusters_raw), c("sel.K", "K"))
cluster_vec <- as.integer(unlist(clusters_raw[sel_row, cell_cols, with = FALSE]))
names(cluster_vec) <- cell_cols

for (sample in names(samea)) {
  sa <- samea[sample]
  cat("---------------------------------------\n")
  cat("Sample:", sample, "(", sa, ")\n")

  is_sample <- startsWith(names(cluster_vec), paste0(sa, "-"))
  sample_clusters <- cluster_vec[is_sample]
  cat("  Cells in clusters TSV:", sum(is_sample), "\n")

  barcodes_cr <- sub(paste0("^", sa, "-"), "", names(sample_clusters))
  barcodes_cr <- paste0(barcodes_cr, "-1")

  bc_file <- file.path(CR_BASE, sample, "outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
  if (!file.exists(bc_file)) {
    cat("  WARNING: CellRanger barcodes not found at", bc_file, "\n\n")
    next
  }
  cr_barcodes <- readLines(gzcon(file(bc_file, "rb")))
  cat("  CellRanger filtered barcodes:", length(cr_barcodes), "\n")

  keep <- barcodes_cr %in% cr_barcodes
  cat("  Intersection (annotation rows):", sum(keep), "\n")

  ann <- data.table(
    barcode = barcodes_cr[keep],
    cluster = paste0("cluster_", sample_clusters[keep])
  )

  tbl <- sort(table(ann$cluster), decreasing = TRUE)
  cat("  Cluster distribution (top 5):\n")
  print(head(tbl, 5))

  out_file <- file.path(ANN_DIR, paste0(sample, "_annotation.tsv"))
  fwrite(ann, out_file, sep = "\t", col.names = FALSE)
  cat("  Written:", out_file, "\n\n")
}

cat("Done. Annotation files in:", ANN_DIR, "\n")
cat("Next step: qsub preprocessing_cec_gec.pbs\n")
