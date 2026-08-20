# reditr_permutation_null_endo2_genes.R
#
# Permutation null for endothelial2 with GENE-LEVEL aggregation.
#
# WHY. reditr_permutation_null_endo2.R records site counts only, and its
# per-permutation site-level outputs are deleted as it goes, so the gene
# counts reported in the summary figure (192 GLMM-sig, 250 Fisher-sig for
# siADAR1) have no null to be compared against. This script re-runs the same
# exhaustive permutations and records BOTH levels, so the gene counts get an
# empirical p on the same footing as the site counts.
#
# GENE MAP SOURCE. Uses vep_merged_<arm>.txt, which is what build_gene_tables.R
# uses -- NOT vep_gene_body_dedup.txt. The dedup file still carries the
# underscore-join bug documented in VEP_summary.md: patch contigs appear as
# "HG1343:HG173:HG459:PATCH:1259969" instead of "HG1343_HG173_HG459_PATCH:1259969",
# which silently drops 8 genes and shifts the counts to 717 / 190 / 248.
# Aggregating off vep_merged reproduces the published 725 / 192 / 250 and
# 467 / 4 / 0 exactly, so the permuted counts are computed by identical logic
# to the observed ones.
#
# A gene is significant if >= 1 of its sites is significant, matching
# build_gene_tables.R (n_GLMM_sig >= 1). SYMBOL == "-" / "" / NA excluded,
# also matching.
#
# The exhaustive labelling logic is copied unchanged from the site-level
# script: C(6,3) = 20 labellings, halved by complement symmetry to 10, minus
# the true one leaves 9 valid nulls. The empirical p floor is therefore 1/10 =
# 0.100 and cannot be lowered by running more.
#
# Writes to a SEPARATE file so the verified site-level output is not clobbered.
#
# Output: reditr_permutation_null_endo2_genes_<arm>.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })
D <- "/rds/general/user/sj1825/home/endothelial2/output"
setwd(D)
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

full_meta <- fread("sample_metadata.txt")
sites_all <- fread("filtered_sites_clustered.txt")
tmpdir <- file.path(D, "perm_null_genes_tmp"); dir.create(tmpdir, showWarnings = FALSE)

run_arm <- function(case) {
  keep <- full_meta[condition %in% c("scr", case)]
  s <- sites_all[sample %in% keep$sample]
  dp <- file.path(tmpdir, sprintf("sites_%s.txt", case))
  fwrite(s, dp, sep = "\t")

  vm <- fread(sprintf("vep_merged_%s.txt", case))
  gene_map <- unique(vm[SYMBOL != "" & SYMBOL != "-" & !is.na(SYMBOL), .(site, SYMBOL)])
  # The universe is genes with >= 1 site that differential_editing actually
  # TESTED in this arm -- 725 for siADAR1, 467 for siADAR2, matching
  # nrow(gene_table_<arm>.txt). It must be counted against the DRE result
  # table, not the raw site file: the raw sites for both arms map to the same
  # 778 symbols, so intersecting with filtered_sites_clustered.txt does not
  # reduce it. 778 is the union across arms and is not a per-arm denominator.
  universe <- uniqueN(merge(gene_map, fread(sprintf("DRE_%s_v4.txt", case))[, .(site)],
                            by = "site")$SYMBOL)

  samples <- keep$sample; truth <- keep$condition
  n_case <- sum(truth == case)

  combs <- combn(seq_along(samples), n_case, simplify = FALSE)
  seen <- character(0); labsets <- list()
  for (cb in combs) {
    lab <- rep("scr", length(samples)); lab[cb] <- case
    key <- paste(lab, collapse = "")
    ckey <- paste(ifelse(lab == "scr", case, "scr"), collapse = "")
    if (key %in% seen || ckey %in% seen) next
    seen <- c(seen, key); labsets[[length(labsets) + 1]] <- lab
  }
  is_true <- vapply(labsets, function(l) all(l == truth) || all(l != truth), logical(1))
  nulls <- labsets[!is_true]
  cat(sprintf("\n#### %s: %d exhaustive nulls | gene universe %d ####\n", case, length(nulls), universe))

  # Null p-value distribution diagnostics. Counting significant sites answers
  # "how many false positives", but the calibration question is whether the
  # p-values are uniform when the labels carry no signal -- which is what a
  # test claiming to control FDR asserts. med_p should sit near 0.5 and the
  # KS test against Uniform(0,1) should NOT reject. These were recorded for
  # diabetes and mouse but not here, leaving the dataset with the largest
  # true-vs-null gap without its most direct calibration evidence.
  ks_unif <- function(p) {
    p <- p[!is.na(p)]
    if (length(p) > 10) suppressWarnings(stats::ks.test(p, "punif")$p.value) else NA_real_
  }

  score <- function(r, lab_str, idx) {
    g <- merge(gene_map, r[, .(site, GLMM_sig, Fisher_sig)], by = "site")[
      , .(nG = sum(GLMM_sig, na.rm = TRUE), nF = sum(Fisher_sig, na.rm = TRUE)), by = SYMBOL]
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
    mp <- file.path(tmpdir, sprintf("meta_%s_%s.txt", case, tag))
    op <- file.path(tmpdir, sprintf("DRE_%s_%s.txt", case, tag))
    fwrite(data.table(sample = samples, condition = lab), mp, sep = "\t")
    r <- tryCatch(as.data.table(differential_editing(
        data_path = dp, meta_path = mp, test = c("glmm", "fisher"),
        reference_level = "scr", case_level = case,
        random_effects = "(1 | sample)", out_path = op,
        summary_path = file.path(tmpdir, sprintf("sum_%s_%s.txt", case, tag)),
        n_cores = N_CORES, verbose = FALSE)),
      error = function(e) { cat("  FAILED", tag, ":", conditionMessage(e), "\n"); NULL })
    unlink(op)
    r
  }

  rows <- list()
  rt <- fit(truth, "TRUE")
  rows[[1]] <- score(rt, paste0("TRUE:", paste(ifelse(truth == "scr", "R", "C"), collapse = "")), 0L)
  cat("  TRUE  "); print(rows[[1]][, .(GLMM_sites, Fisher_sites, GLMM_genes, Fisher_genes)])

  for (i in seq_along(nulls)) {
    lab <- nulls[[i]]
    r <- fit(lab, sprintf("%02d", i)); if (is.null(r)) next
    rows[[length(rows) + 1]] <- score(r, paste(ifelse(lab == "scr", "R", "C"), collapse = ""), i)
    cat(sprintf("  null %d/%d  ", i, length(nulls)))
    print(rows[[length(rows)]][, .(GLMM_sites, Fisher_sites, GLMM_genes, Fisher_genes)])
  }

  out <- rbindlist(rows)
  fwrite(out, sprintf("reditr_permutation_null_endo2_genes_%s.txt", case), sep = "\t")

  tv <- out[null == 0]; nl <- out[null > 0]
  ep <- function(o, v) (sum(v >= o) + 1) / (length(v) + 1)
  cat(sprintf("\n=== %s (n = %d exhaustive nulls) ===\n", case, nrow(nl)))
  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "GLMM sites",  tv$GLMM_sites,  mean(nl$GLMM_sites),  min(nl$GLMM_sites),  max(nl$GLMM_sites),  ep(tv$GLMM_sites, nl$GLMM_sites)))
  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "Fisher sites", tv$Fisher_sites, mean(nl$Fisher_sites), min(nl$Fisher_sites), max(nl$Fisher_sites), ep(tv$Fisher_sites, nl$Fisher_sites)))
  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "GLMM genes",  tv$GLMM_genes,  mean(nl$GLMM_genes),  min(nl$GLMM_genes),  max(nl$GLMM_genes),  ep(tv$GLMM_genes, nl$GLMM_genes)))
  cat(sprintf("  %-14s true %5d | permuted mean %6.1f range %d-%d | emp p = %.3f\n",
              "Fisher genes", tv$Fisher_genes, mean(nl$Fisher_genes), min(nl$Fisher_genes), max(nl$Fisher_genes), ep(tv$Fisher_genes, nl$Fisher_genes)))
  cat(sprintf("  %-14s true %.3f | permuted %.3f-%.3f   (0.5 expected under a clean null)\n",
              "median p GLMM", tv$med_p_glmm, min(nl$med_p_glmm), max(nl$med_p_glmm)))
  cat(sprintf("  %-14s true %.3g | permuted %.3g-%.3g  (KS vs uniform; >0.05 = not rejected)\n",
              "KS p GLMM", tv$ks_glmm, min(nl$ks_glmm), max(nl$ks_glmm)))
  out
}

for (a in c("siADAR1", "siADAR2")) run_arm(a)
cat("\nWritten: reditr_permutation_null_endo2_genes_siADAR1.txt / _siADAR2.txt\n")
