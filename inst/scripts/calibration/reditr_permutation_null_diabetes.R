# reditr_permutation_null_diabetes.R
#
# Test the calibration of the reditR package using permutation of the
# diabetes dataset.
#
# The aim is to test whether reditR produces false-positive discoveries
# when the sample labels are randomly reassigned. Unlike a simulation,
# permutation keeps the real data unchanged, including its coverage,
# editing-rate variation and between-sample variability.
#
# The 12 samples are kept balanced at 6 control and 6 diabetic samples.
# Randomly changing the labels removes the original biological comparison
# while preserving the structure of the observed data.
#
# The true labelling is excluded because it contains the real biological
# signal and is therefore not part of the null distribution.
#
# Output:
#   reditr_permutation_null_diabetes.txt  - summary for each permutation
#   reditr_permutation_null_pvals.rds     - site-level p-values for follow-up


suppressPackageStartupMessages({
  library(data.table)
  library(reditR)
})

D <- "/rds/general/user/sj1825/home/diabetes_output"
setwd(D)

# Number of CPU cores used by reditR.
# This can be changed through the N_CORES environment variable when running
# the script on the HPC system.
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

# Read the sample metadata and store the original condition labels.
# These labels are used as the reference for generating the permutations.
meta <- fread("sample_metadata.txt")
samples <- meta$sample
truth <- meta[match(samples, sample), condition]
n_dia <- sum(truth == "diabetic")

# Create a temporary directory for metadata and result files generated
# during the permutation analysis.
tmpdir <- file.path(D, "perm_null_tmp")
dir.create(tmpdir, showWarnings = FALSE)

# Identify labellings that are not useful as null permutations.
# The original labelling is excluded because it contains the real signal.
# Its exact opposite is also excluded because it represents the same
# two-group partition with the labels reversed.
is_trivial <- function(lab) {
  all(lab == truth) || all(lab != truth)
}

# Generate every possible balanced assignment of samples to the diabetic
# group. A labelling and its exact complement represent the same partition,
# so only one of each pair is retained.
#
# Using all possible assignments makes the permutation test exhaustive
# rather than relying on a small random sample of permutations.
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

# Select which part of the permutation list this job should run.
# PERM_FROM and PERM_TO allow the full analysis to be split across several
# HPC jobs instead of running all permutations in one job.
PERM_FROM <- as.integer(Sys.getenv("PERM_FROM", "1"))
PERM_TO   <- as.integer(Sys.getenv("PERM_TO", as.character(length(all_labs))))
PERM_TO   <- min(PERM_TO, length(all_labs))
idx <- seq.int(PERM_FROM, PERM_TO)

cat(sprintf("Exhaustive labelling space: %d nulls. Running %d-%d (%d this job).\n",
            length(all_labs), PERM_FROM, PERM_TO, length(idx)))

# Store summary statistics and site-level results from each permutation.
rows <- list(); pvals <- list()

for (done in idx) {
  lab <- all_labs[[done]]

  # Create a metadata file containing the current permuted condition labels.
  # This file is passed to differential_editing(), while the editing data
  # themselves remain unchanged.
  mp <- file.path(tmpdir, sprintf("meta_perm%03d.txt", done))
  op <- file.path(tmpdir, sprintf("DRE_perm%03d.txt", done))
  fwrite(data.table(sample = samples, condition = lab), mp, sep = "\t")

  cat(sprintf("\n[perm %03d/%d] %s\n", done, length(all_labs),
              paste(paste0(samples, "=", substr(lab, 1, 1)), collapse = " ")))
  t0 <- Sys.time()

  # Run the complete reditR differential-editing analysis using the
  # permuted labels.
  #
  # All three statistical tests are run so that their false-positive
  # behaviour can be compared under the same null permutation.
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

  # Test whether the p-values are consistent with a uniform distribution.
  # Under a valid null model, p-values should be approximately Uniform(0,1).
  # The test is only performed when there are enough non-missing p-values.
  ks <- function(p) { p <- p[!is.na(p)]
    if (length(p) > 10) suppressWarnings(stats::ks.test(p, "punif")$p.value) else NA_real_ }

  # Summarise the results from this permutation.
  # GLMM convergence is measured as the proportion of sites with a
  # non-missing p-value.
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

  # Keep site-level results as well as summary statistics.
  # This allows the same sites to be examined across permutations and
  # allows later analyses of which types of sites are falsely detected.
  pvals[[length(pvals) + 1]] <- data.table(perm = done, site = r$site,
    glmm = r$glmm_pvalue, fisher = r$fisher_pvalue,
    glmm_sig = r$GLMM_sig, fisher_sig = r$Fisher_sig)

  # Print a short progress summary for the current permutation.
  print(rows[[length(rows)]][, .(GLMM_sig, Fisher_sig, Wilcox_sig,
                                 med_p_glmm = round(med_p_glmm, 3), mins)])

  # Remove temporary output before moving to the next permutation.
  unlink(op)
}

# Use a different output name for each HPC chunk so parallel jobs do not
# overwrite one another.
PREFIX <- Sys.getenv("OUT_PREFIX", "reditr_permutation_null")

# Add the permutation range to chunk filenames when only part of the
# complete permutation space was analysed.
CHUNK  <- if (PERM_FROM == 1L && PERM_TO == length(all_labs)) {
  ""
} else {
  sprintf("_p%03d-%03d", PERM_FROM, PERM_TO)
}

# Combine the results from all permutations completed by this job and
# save both the summary results and site-level p-values.
out <- rbindlist(rows)
fwrite(out, paste0(PREFIX, "_diabetes", CHUNK, ".txt"), sep = "\t")
saveRDS(rbindlist(pvals), paste0(PREFIX, "_pvals", CHUNK, ".rds"))

# Display the number of significant sites found in each permutation.
# Because these permutations represent the null distribution, this gives
# a direct indication of how many discoveries reditR makes when there is
# no true condition assignment.
cat("\n\n=== reditR PERMUTATION NULL — ", nrow(out), " relabellings ===\n", sep = "")
print(out[, .(perm, labels, GLMM_sig, Fisher_sig, Wilcox_sig,
              med_p_glmm = round(med_p_glmm, 3), ks_glmm = signif(ks_glmm, 3))])

# Summarise the average and maximum number of significant sites across
# permutations for each statistical test.
cat(sprintf("\nSignificant sites per permutation — GLMM: mean %.1f (max %d) | Fisher: mean %.1f (max %d) | Wilcoxon: mean %.1f (max %d)\n",
            mean(out$GLMM_sig), max(out$GLMM_sig), mean(out$Fisher_sig),
            max(out$Fisher_sig), mean(out$Wilcox_sig), max(out$Wilcox_sig)))

# Show the results from the real, unpermuted condition labels for comparison.
# The real labelling is not part of the null distribution, but provides
# context for how many discoveries are observed in the actual analysis.
cat(sprintf("TRUE labelling for comparison       — GLMM: 226 | Fisher: 538 | Wilcoxon: 0\n"))

# Warn that an empirical p-value cannot be calculated from an individual
# chunk because the complete null distribution is required.
if (nzchar(CHUNK)) {
  cat("\nChunk only -- empirical p is NOT computable until all",
      length(all_labs), "nulls are concatenated.\n")
}

cat("\nWritten: ", PREFIX, "_diabetes", CHUNK, ".txt\n", sep = "")