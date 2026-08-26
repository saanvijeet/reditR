# reditr_permutation_null_endo2_genes.R
#
# Permutation null analysis for the endothelial2 dataset at both
# site and gene level.
#
# The site-level permutation analysis counts significant editing sites,
# but this script also aggregates significant sites into genes. This allows
# gene-level results to be compared with a matching permutation null.
#
# The same gene annotation used for the observed results is used here so
# that the observed and permuted results are calculated in the same way.
#
# A gene is counted as significant if at least one of its sites is
# significant.
#
# The analysis uses all possible balanced label assignments. For this
# dataset there are 6 samples, with 3 assigned to each condition.
# Complementary assignments represent the same partition, so only one
# of each pair is retained. The true labelling is then removed, leaving
# the null permutations.
#
# Output:
#   reditr_permutation_null_endo2_genes_<arm>.txt


suppressPackageStartupMessages({ library(data.table); library(reditR) })

D <- "/rds/general/user/sj1825/home/endothelial2/output"
setwd(D)

# Number of CPU cores used by reditR.
# This can be changed through the N_CORES environment variable when
# running the analysis on the HPC system.
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

# Read the full sample metadata and filtered editing sites.
# These provide the samples and editing observations used in the
# permutation analysis.
full_meta <- fread("sample_metadata.txt")
sites_all <- fread("filtered_sites_clustered.txt")

# Temporary files are created for the metadata and differential-editing
# results generated during each permutation.
tmpdir <- file.path(D, "perm_null_genes_tmp"); dir.create(tmpdir, showWarnings = FALSE)


run_arm <- function(case) {

  # Keep only the reference samples (scr) and the current case condition.
  # The analysis is therefore performed separately for siADAR1 and siADAR2.
  keep <- full_meta[condition %in% c("scr", case)]

  # Keep editing sites belonging to the samples in this comparison.
  s <- sites_all[sample %in% keep$sample]

  # Write the editing data to a temporary file so it can be passed to
  # differential_editing().
  dp <- file.path(tmpdir, sprintf("sites_%s.txt", case))
  fwrite(s, dp, sep = "\t")

  # Read the VEP annotation used for the observed gene-level results.
  # Only sites with a valid gene symbol are retained because sites without
  # a gene annotation cannot be assigned to a gene.
  vm <- fread(sprintf("vep_merged_%s.txt", case))
  gene_map <- unique(vm[SYMBOL != "" & SYMBOL != "-" & !is.na(SYMBOL), .(site, SYMBOL)])

  # Define the gene universe as genes that have at least one site actually
  # tested by differential_editing() in this arm.
  # This gives the correct denominator for the gene-level analysis rather
  # than counting genes present only in the raw annotation.
  universe <- uniqueN(merge(gene_map, fread(sprintf("DRE_%s_v4.txt", case))[, .(site)],
                            by = "site")$SYMBOL)

  samples <- keep$sample
  truth <- keep$condition
  n_case <- sum(truth == case)

  # Generate every possible assignment of the correct number of samples
  # to the case group.
  combs <- combn(seq_along(samples), n_case, simplify = FALSE)

  # Store only unique sample partitions.
  # A partition and its complementary labelling contain the same two groups,
  # so only one needs to be tested.
  seen <- character(0); labsets <- list()
  for (cb in combs) {
    lab <- rep("scr", length(samples)); lab[cb] <- case
    key <- paste(lab, collapse = "")
    ckey <- paste(ifelse(lab == "scr", case, "scr"), collapse = "")
    if (key %in% seen || ckey %in% seen) next
    seen <- c(seen, key); labsets[[length(labsets) + 1]] <- lab
  }

  # Identify and remove the true biological labelling.
  # The remaining assignments form the null distribution used for
  # calculating empirical significance.
  is_true <- vapply(labsets, function(l) all(l == truth) || all(l != truth), logical(1))
  nulls <- labsets[!is_true]

  cat(sprintf("\n#### %s: %d exhaustive nulls | gene universe %d ####\n", case, length(nulls), universe))


  # Test whether the null p-values are consistent with a uniform
  # distribution, as expected when there is no true condition effect.
  # The test is only performed when there are enough non-missing p-values.
  ks_unif <- function(p) {
    p <- p[!is.na(p)]
    if (length(p) > 10) suppressWarnings(stats::ks.test(p, "punif")$p.value) else NA_real_
  }


  score <- function(r, lab_str, idx) {

    # Join the significant site results to their gene annotations.
    # Sites belonging to the same gene are then grouped together.
    g <- merge(gene_map, r[, .(site, GLMM_sig, Fisher_sig)], by = "site")[
      , .(nG = sum(GLMM_sig, na.rm = TRUE), nF = sum(Fisher_sig, na.rm = TRUE)), by = SYMBOL]

    # Summarise both site-level and gene-level results.
    # A gene is counted as significant when at least one of its sites
    # is significant for the relevant statistical test.
    data.table(arm = case, null = idx, labels = lab_str,
               GLMM_sites = sum(r$GLMM_sig, na.rm = TRUE),
               Fisher_sites = sum(r$Fisher_sig, na.rm = TRUE),
               GLMM_genes = sum(g$nG >= 1), Fisher_genes = sum(g$nF >= 1),
               both_genes = sum(g$nG >= 1 & g$nF >= 1), gene_universe = universe,
               n_tested_glmm = sum(!is.na(r$glmm_pvalue)),
               med_p_glmm    = median(r$glmm_pvalue,   na.rm = TRUE),
               med_p_fisher  = median(r$fisher_pvalue, na.rm = TRUE),
               ks_glmm       = ks_unif(r$glmm_pvalue),
               ks_fisher     = ks_unif(r$fisher_pvalue))
  }


  fit <- function(lab, tag) {

    # Create metadata containing the current condition labels.
    # The editing data remain unchanged; only the sample labels are
    # permuted between conditions.
    mp <- file.path(tmpdir, sprintf("meta_%s_%s.txt", case, tag))
    op <- file.path(tmpdir, sprintf("DRE_%s_%s.txt", case, tag))
    fwrite(data.table(sample = samples, condition = lab), mp, sep = "\t")

    # Run the same differential-editing analysis for the current labelling.
    # GLMM and Fisher's exact test are used here because these are the
    # tests being compared at both site and gene level.
    r <- tryCatch(as.data.table(differential_editing(
        data_path = dp, meta_path = mp, test = c("glmm", "fisher"),
        reference_level = "scr", case_level = case,
        random_effects = "(1 | sample)", out_path = op,
        summary_path = file.path(tmpdir, sprintf("sum_%s_%s.txt", case, tag)),
        n_cores = N_CORES, verbose = FALSE)),
      error = function(e) { cat("  FAILED", tag, ":", conditionMessage(e), "\n"); NULL })

    # Remove the temporary differential-editing output after it has been read.
    unlink(op)
    r
  }


  rows <- list()

  # First analyse the real, unpermuted condition labels.
  # This provides the observed result that will later be compared
  # against the permutation null distribution.
  rt <- fit(truth, "TRUE")
  rows[[1]] <- score(rt, paste0("TRUE:", paste(ifelse(truth == "scr", "R", "C"), collapse = "")), 0L)

  cat("  TRUE  "); print(rows[[1]][, .(GLMM_sites, Fisher_sites, GLMM_genes, Fisher_genes)])


  # Run the differential-editing analysis for every null labelling.
  # These results show how many significant sites and genes can be
  # obtained when the condition labels contain no true signal.
  for (i in seq_along(nulls)) {
    lab <- nulls[[i]]
    r <- fit(lab, sprintf("%02d", i)); if (is.null(r)) next

    rows[[length(rows) + 1]] <- score(r, paste(ifelse(lab == "scr", "R", "C"), collapse = ""), i)

    cat(sprintf("  null %d/%d  ", i, length(nulls)))
    print(rows[[length(rows)]][, .(GLMM_sites, Fisher_sites, GLMM_genes, Fisher_genes)])
  }


  # Combine the observed and null results into one table and save them.
  out <- rbindlist(rows)
  fwrite(out, sprintf("reditr_permutation_null_endo2_genes_%s.txt", case), sep = "\t")


  # Separate the true result from the permutation results.
  tv <- out[null == 0]
  nl <- out[null > 0]

  # Calculate the empirical p-value.
  # The +1 correction prevents an empirical p-value of exactly zero and
  # accounts for the observed result when estimating the tail probability.
  ep <- function(o, v) (sum(v >= o) + 1) / (length(v) + 1)


  # Print the observed result, the null distribution, and the empirical
  # p-value for each site- and gene-level measure.
  cat(sprintf("\n=== %s (n = %d exhaustive nulls) ===\n", case, nrow(nl)))

  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "GLMM sites",  tv$GLMM_sites,  mean(nl$GLMM_sites),  min(nl$GLMM_sites),  max(nl$GLMM_sites), ep(tv$GLMM_sites, nl$GLMM_sites)))

  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "Fisher sites", tv$Fisher_sites, mean(nl$Fisher_sites), min(nl$Fisher_sites), max(nl$Fisher_sites), ep(tv$Fisher_sites, nl$Fisher_sites)))

  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "GLMM genes",  tv$GLMM_genes,  mean(nl$GLMM_genes),  min(nl$GLMM_genes),  max(nl$GLMM_genes), ep(tv$GLMM_genes, nl$GLMM_genes)))

  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "Fisher genes", tv$Fisher_genes, mean(nl$Fisher_genes), min(nl$Fisher_genes), max(nl$Fisher_genes), ep(tv$Fisher_genes, nl$Fisher_genes)))

  # Report the median p-value and KS test results as additional checks
  # of whether the null p-values behave as expected.
  cat(sprintf("  %-14s true %.3f | permuted %.3f-%.3f   (0.5 expected under a clean null)\n",
              "median p GLMM", tv$med_p_glmm, min(nl$med_p_glmm), max(nl$med_p_glmm)))

  cat(sprintf("  %-14s true %.3g | permuted %.3g-%.3g  (KS vs uniform; >0.05 = not rejected)\n",
              "KS p GLMM", tv$ks_glmm, min(nl$ks_glmm), max(nl$ks_glmm)))

  out
}


# Run the complete permutation analysis separately for each siADAR arm.
for (a in c("siADAR1", "siADAR2")) run_arm(a)

cat("\nWritten: reditr_permutation_null_endo2_genes_siADAR1.txt / _siADAR2.txt\n")