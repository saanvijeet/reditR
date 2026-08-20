# reditr_permutation_null_endo2.R
#
# Permutation null for the endothelial2 ADAR-knockdown dataset, run through
# reditR — the same test applied to the diabetes dataset, so the two are
# directly comparable.
#
# This dataset is the project's positive control: siADAR1 knocks down the
# editing enzyme itself, so a large reduction in editing MUST occur if the
# workflow works. Site-level paired analysis already shows 78.5% of 7,964
# sites falling in siADAR1 (median delta -0.093, sign test p ~ 5e-324).
# The question here is whether reditR's significance calls exceed what
# label-shuffling produces — i.e. whether the workflow's DIFFERENTIAL TEST,
# not just the effect size, carries real information.
#
# The contrast has only 6 samples (3 scr, 3 siADAR1), so the permutation set
# is EXHAUSTIVE: C(6,3) = 20 labellings, 10 distinct after accounting for the
# arbitrary swap of group names, minus the true one = 9 nulls. No sampling.
#
# siADAR2 is run identically as a second arm. Biology predicts a much
# smaller effect there (ADAR2 contributes little to global Alu editing), so
# it acts as an internal negative control from the same experiment.
#
# Output: reditr_permutation_null_endo2_<arm>.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })

D <- "/rds/general/user/sj1825/home/endothelial2/output"
setwd(D)
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

full_meta <- fread("sample_metadata.txt")
sites_all <- fread("filtered_sites_clustered.txt")
tmpdir <- file.path(D, "perm_null_endo2_tmp"); dir.create(tmpdir, showWarnings = FALSE)

run_arm <- function(case) {
  keep <- full_meta[condition %in% c("scr", case)]
  s <- sites_all[sample %in% keep$sample]
  dp <- file.path(tmpdir, sprintf("sites_%s.txt", case))
  fwrite(s, dp, sep = "\t")

  samples <- keep$sample
  truth <- keep$condition
  n_case <- sum(truth == case)

  # exhaustive balanced labellings, de-duplicated by group-name symmetry
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
  cat("\n#### arm:", case, "-", length(nulls), "exhaustive null labellings ####\n")

  rows <- list()
  for (i in seq_along(nulls)) {
    lab <- nulls[[i]]
    mp <- file.path(tmpdir, sprintf("meta_%s_%02d.txt", case, i))
    op <- file.path(tmpdir, sprintf("DRE_%s_%02d.txt", case, i))
    fwrite(data.table(sample = samples, condition = lab), mp, sep = "\t")
    cat(sprintf("[%s null %d/%d] %s\n", case, i, length(nulls), paste(substr(lab,1,4), collapse=" ")))
    res <- tryCatch(differential_editing(
        data_path = dp, meta_path = mp, test = c("glmm", "fisher", "wilcoxon"),
        reference_level = "scr", case_level = case,
        random_effects = "(1 | sample)", out_path = op,
        summary_path = file.path(tmpdir, sprintf("sum_%s_%02d.txt", case, i)),
        n_cores = N_CORES),
      error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
    if (is.null(res)) next
    r <- as.data.table(res)
    # NB: do not abbreviate with substr(lab, 1, 1) -- "scr" and "siADAR1"/
    # "siADAR2" share their first letter, so every row logged as "ssssss" and
    # the labelling could not be audited from the output. Encode as
    # reference/case instead.
    rows[[length(rows)+1]] <- data.table(arm = case, null = i,
      labels = paste(ifelse(lab == "scr", "R", "C"), collapse = ""),
      GLMM_sig = sum(r$GLMM_sig, na.rm = TRUE),
      Fisher_sig = sum(r$Fisher_sig, na.rm = TRUE),
      Wilcox_sig = sum(r$Wilcox_sig, na.rm = TRUE))
    print(rows[[length(rows)]][, .(GLMM_sig, Fisher_sig, Wilcox_sig)])
    unlink(op)
  }

  # the TRUE labelling, run through the identical code path
  mp <- file.path(tmpdir, sprintf("meta_%s_TRUE.txt", case))
  fwrite(data.table(sample = samples, condition = truth), mp, sep = "\t")
  cat(sprintf("[%s TRUE labelling]\n", case))
  rt <- as.data.table(differential_editing(
      data_path = dp, meta_path = mp, test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "scr", case_level = case, random_effects = "(1 | sample)",
      out_path = file.path(tmpdir, sprintf("DRE_%s_TRUE.txt", case)),
      summary_path = file.path(tmpdir, sprintf("sum_%s_TRUE.txt", case)),
      n_cores = N_CORES))
  truth_row <- data.table(arm = case, null = 0L,
    labels = paste0("TRUE:", paste(ifelse(truth == "scr", "R", "C"), collapse = "")),
    GLMM_sig = sum(rt$GLMM_sig, na.rm = TRUE),
    Fisher_sig = sum(rt$Fisher_sig, na.rm = TRUE),
    Wilcox_sig = sum(rt$Wilcox_sig, na.rm = TRUE))

  out <- rbind(truth_row, rbindlist(rows))
  fwrite(out, sprintf("reditr_permutation_null_endo2_%s.txt", case), sep = "\t")
  nl <- out[null > 0]
  cat(sprintf("\n=== %s ===\nTRUE: GLMM %d | Fisher %d | Wilcox %d\nNULL (n=%d): GLMM mean %.1f range %d-%d | Fisher mean %.1f range %d-%d\nempirical p (GLMM) = %.3f | (Fisher) = %.3f\n",
    case, truth_row$GLMM_sig, truth_row$Fisher_sig, truth_row$Wilcox_sig, nrow(nl),
    mean(nl$GLMM_sig), min(nl$GLMM_sig), max(nl$GLMM_sig),
    mean(nl$Fisher_sig), min(nl$Fisher_sig), max(nl$Fisher_sig),
    (sum(nl$GLMM_sig >= truth_row$GLMM_sig) + 1) / (nrow(nl) + 1),
    (sum(nl$Fisher_sig >= truth_row$Fisher_sig) + 1) / (nrow(nl) + 1)))
  out
}

all_out <- rbindlist(lapply(c("siADAR1", "siADAR2"), run_arm))
fwrite(all_out, "reditr_permutation_null_endo2.txt", sep = "\t")
cat("\nWritten: reditr_permutation_null_endo2.txt\n")
