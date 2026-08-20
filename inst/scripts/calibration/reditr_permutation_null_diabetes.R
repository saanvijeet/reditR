# reditr_permutation_null_diabetes.R
#
# DIRECT null-calibration test of the reditR package itself.
#
# The existing null_calibration_sim_diabetes.R and the noise-injection
# simulation validated differential_editing_analysis_v3.R -- the hand-rolled
# predecessor pipeline -- not reditR. They establish that the statistics are
# sound, but not that the package implements them faithfully. This script
# closes that gap by running reditR::differential_editing() itself on data
# where the answer is known in advance.
#
# Design: permute the condition labels across the 12 real samples, keeping
# the 6/6 balance, and re-run the full reditR workflow unchanged. Each
# permutation destroys any real diabetic-vs-control signal while preserving
# every other property of the data -- coverage, editing-rate distribution,
# between-sample variance, site correlation structure. A correctly
# calibrated workflow must therefore return ~0 significant sites after FDR
# correction, and p-values indistinguishable from uniform.
#
# This is a stronger test than a simulation because the null data ARE the
# real data; nothing about the noise model has to be assumed.
#
# The true labelling is excluded (a permutation identical to the real one,
# or its exact complement, would not be a null).
#
# Output: reditr_permutation_null_diabetes.txt  (one row per permutation)
#         reditr_permutation_null_pvals.rds     (pooled null p-values)

suppressPackageStartupMessages({
  library(data.table)
  library(reditR)
})

D <- "/rds/general/user/sj1825/home/diabetes_output"
setwd(D)
# EXHAUSTIVE, NOT SAMPLED. This previously drew N_PERM = 10 labellings at
# random from the 462 distinct 6v6 assignments (C(12,6)/2 after collapsing
# complements). With only 9 usable nulls the empirical p floor was 1/10 =
# 0.100 -- a compute choice mistaken for a design limit, and the reason the
# diabetes result read p = 0.200 when the design supports p as low as
# 1/463 = 0.00216. Enumeration is now deterministic and complete; PERM_FROM /
# PERM_TO select a contiguous chunk so the ~34 h serial workload can be run
# as a PBS array and the chunks concatenated afterwards.
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

meta <- fread("sample_metadata.txt")
samples <- meta$sample
truth <- meta[match(samples, sample), condition]
n_dia <- sum(truth == "diabetic")

tmpdir <- file.path(D, "perm_null_tmp")
dir.create(tmpdir, showWarnings = FALSE)

is_trivial <- function(lab) {
  # identical to the truth, or its exact complement (same partition, relabelled)
  all(lab == truth) || all(lab != truth)
}

# Enumerate every distinct labelling once. A labelling and its complement give
# the same two-sample partition, so only one of each pair is kept.
all_labs <- local({
  seen <- character(0); out <- list()
  for (cb in combn(seq_along(samples), n_dia, simplify = FALSE)) {
    lab <- rep("control", length(samples)); lab[cb] <- "diabetic"
    key  <- paste(substr(lab, 1, 1), collapse = "")
    ckey <- paste(substr(ifelse(lab == "control", "diabetic", "control"), 1, 1), collapse = "")
    if (key %in% seen || ckey %in% seen) next
    seen <- c(seen, key); out[[length(out) + 1]] <- lab
  }
  out[!vapply(out, is_trivial, logical(1))]
})

PERM_FROM <- as.integer(Sys.getenv("PERM_FROM", "1"))
PERM_TO   <- as.integer(Sys.getenv("PERM_TO", as.character(length(all_labs))))
PERM_TO   <- min(PERM_TO, length(all_labs))
idx <- seq.int(PERM_FROM, PERM_TO)
cat(sprintf("Exhaustive labelling space: %d nulls. Running %d-%d (%d this job).\n",
            length(all_labs), PERM_FROM, PERM_TO, length(idx)))

rows <- list(); pvals <- list()
for (done in idx) {
  lab <- all_labs[[done]]

  mp <- file.path(tmpdir, sprintf("meta_perm%03d.txt", done))
  op <- file.path(tmpdir, sprintf("DRE_perm%03d.txt", done))
  fwrite(data.table(sample = samples, condition = lab), mp, sep = "\t")

  cat(sprintf("\n[perm %03d/%d] %s\n", done, length(all_labs),
              paste(paste0(samples, "=", substr(lab, 1, 1)), collapse = " ")))
  t0 <- Sys.time()
  res <- tryCatch(differential_editing(
      data_path = "filtered_sites_clustered_t.txt", meta_path = mp,
      test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "control", case_level = "diabetic",
      random_effects = "(1 | sample)", out_path = op,
      summary_path = file.path(tmpdir, sprintf("summary_perm%03d.txt", done)),
      n_cores = N_CORES),
    error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next

  r <- as.data.table(res)
  ks <- function(p) { p <- p[!is.na(p)]
    if (length(p) > 10) suppressWarnings(stats::ks.test(p, "punif")$p.value) else NA_real_ }
  rows[[length(rows) + 1]] <- data.table(
    perm = done,
    labels = paste(substr(lab, 1, 1), collapse = ""),
    n_tested_glmm  = sum(!is.na(r$glmm_pvalue)),
    GLMM_sig   = sum(r$GLMM_sig,   na.rm = TRUE),
    Fisher_sig = sum(r$Fisher_sig, na.rm = TRUE),
    Wilcox_sig = sum(r$Wilcox_sig, na.rm = TRUE),
    med_p_glmm   = median(r$glmm_pvalue,   na.rm = TRUE),
    med_p_fisher = median(r$fisher_pvalue, na.rm = TRUE),
    ks_glmm   = ks(r$glmm_pvalue),
    ks_fisher = ks(r$fisher_pvalue),
    mins = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
  # Keep the site IDs, not just the p-values. Without them the null results
  # cannot be joined to the VEP annotation, which blocks the follow-up
  # question of whether the coding-sequence enrichment seen on the true
  # labelling also appears under permutation -- i.e. whether that finding is
  # biology or a property of which sites the procedure falsely calls.
  pvals[[length(pvals) + 1]] <- data.table(perm = done, site = r$site,
    glmm = r$glmm_pvalue, fisher = r$fisher_pvalue,
    glmm_sig = r$GLMM_sig, fisher_sig = r$Fisher_sig)
  print(rows[[length(rows)]][, .(GLMM_sig, Fisher_sig, Wilcox_sig,
                                 med_p_glmm = round(med_p_glmm, 3), mins)])
  unlink(op)
}

# OUT_PREFIX lets a follow-up batch run without clobbering an earlier one.
# Array chunks additionally carry their range, so concurrent tasks cannot
# overwrite each other; concatenate the chunk files before scoring.
PREFIX <- Sys.getenv("OUT_PREFIX", "reditr_permutation_null")
CHUNK  <- if (PERM_FROM == 1L && PERM_TO == length(all_labs)) {
  ""
} else {
  sprintf("_p%03d-%03d", PERM_FROM, PERM_TO)
}
out <- rbindlist(rows)
fwrite(out, paste0(PREFIX, "_diabetes", CHUNK, ".txt"), sep = "\t")
saveRDS(rbindlist(pvals), paste0(PREFIX, "_pvals", CHUNK, ".rds"))

cat("\n\n=== reditR PERMUTATION NULL — ", nrow(out), " relabellings ===\n", sep = "")
print(out[, .(perm, labels, GLMM_sig, Fisher_sig, Wilcox_sig,
              med_p_glmm = round(med_p_glmm, 3), ks_glmm = signif(ks_glmm, 3))])
cat(sprintf("\nSignificant sites per permutation — GLMM: mean %.1f (max %d) | Fisher: mean %.1f (max %d) | Wilcoxon: mean %.1f (max %d)\n",
            mean(out$GLMM_sig), max(out$GLMM_sig), mean(out$Fisher_sig),
            max(out$Fisher_sig), mean(out$Wilcox_sig), max(out$Wilcox_sig)))
cat(sprintf("TRUE labelling for comparison       — GLMM: 226 | Fisher: 538 | Wilcoxon: 0\n"))
if (nzchar(CHUNK)) {
  cat("\nChunk only -- empirical p is NOT computable until all",
      length(all_labs), "nulls are concatenated.\n")
}
cat("\nWritten: ", PREFIX, "_diabetes", CHUNK, ".txt\n", sep = "")
