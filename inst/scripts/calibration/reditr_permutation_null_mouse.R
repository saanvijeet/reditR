# reditr_permutation_null_mouse.R
#
# Permutation null for the mouse pan-EC dehydration dataset, run through
# reditR with the SAME crossed random-effects specification the real analysis
# used — so this measures the calibration of the analysis as performed, not
# of a simplified stand-in.
#
# Third of three: the same test applied to diabetes (10 sampled permutations)
# and endothelial2 (9 exhaustive). Results are directly comparable.
#
# PERMUTATION UNIT. Labels are shuffled at the LIBRARY level, not the
# pseudobulk-unit level. There are 192 pseudobulk units but only 13
# libraries, and the library is the unit of biological replication —
# permuting units would leave each library split across both arms and
# destroy the null's meaning by making the two groups near-identical.
# 3 of 13 libraries are assigned "control", matching the true 3/10 design.
#
# A CAVEAT specific to this dataset, which weakens the test in a way the
# other two do not share. The true control set is exactly one library per
# compartment (MEC1, CEC1, GEC1), all at 0h. Condition is therefore
# completely confounded with timepoint, and partly with compartment balance.
# Any permutation breaks that structure, so this null asks "does the true
# 3-vs-10 split yield more signal than an arbitrary 3-vs-10 split" — which
# is the right question, but a permutation that happens to isolate one
# compartment or one timepoint is not a pure null either. Interpret an
# equivocal result here more cautiously than the same result on
# endothelial2, where the arms are true experimental replicates.
#
# There are C(13,3) = 286 possible labellings; 10 are sampled.
#
# Output: reditr_permutation_null_mouse.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })

EPH <- "/rds/general/user/sj1825/ephemeral/mec_dehydration/sprint_output_full"
OUT <- "/rds/general/user/sj1825/home/msc_prj/test_data_mouse"
setwd(OUT)
# EXHAUSTIVE, NOT SAMPLED. This previously drew N_PERM = 10 control sets at
# random from the C(13,3) = 286 possible, i.e. 3.2% of the space, pinning the
# empirical p at the 1/10 = 0.100 floor. That is a compute choice, not a
# design limit -- the diabetes dataset moved from p = 0.200 to p = 0.0368 once
# its full space was enumerated, so the distinction changes conclusions.
# Arms are unequal (3 control vs 10 dehydrated), so there is no complement
# symmetry to collapse: 286 sets minus the true one leaves 285 nulls and a
# floor of 1/286 = 0.0035.
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

base <- fread(file.path(EPH, "all_ec_clustered_with_condition.txt"))
libs <- unique(base[, .(library, condition)])[order(library)]
true_ctrl <- sort(libs[condition == "control", library])
n_ctrl <- length(true_ctrl)
cat("libraries:", nrow(libs), " true controls:", paste(true_ctrl, collapse = ", "), "\n")

tmpdir <- file.path(EPH, "perm_null_mouse_tmp"); dir.create(tmpdir, showWarnings = FALSE)

all_ctrl <- local({
  cs <- combn(sort(libs$library), n_ctrl, simplify = FALSE)
  cs[!vapply(cs, function(c) identical(sort(c), true_ctrl), logical(1))]
})
PERM_FROM <- as.integer(Sys.getenv("PERM_FROM", "1"))
PERM_TO   <- min(as.integer(Sys.getenv("PERM_TO", as.character(length(all_ctrl)))),
                 length(all_ctrl))
CHUNK <- if (PERM_FROM == 1L && PERM_TO == length(all_ctrl)) "" else
           sprintf("_p%03d-%03d", PERM_FROM, PERM_TO)
cat(sprintf("Exhaustive space: %d nulls. Running %d-%d.\n",
            length(all_ctrl), PERM_FROM, PERM_TO))

# KS against Uniform(0,1): under a clean null the p-values should be uniform,
# so this should NOT reject. Recorded for diabetes and endothelial2 but not
# here, leaving the mouse without its most direct calibration evidence.
ks_unif <- function(p) { p <- p[!is.na(p)]
  if (length(p) > 10) suppressWarnings(stats::ks.test(p, "punif")$p.value) else NA_real_ }

# Written once per task, condition stripped so the per-permutation metadata
# map supplies it instead.
dp <- file.path(tmpdir, sprintf("data_shared_%d.txt", PERM_FROM))
fwrite(base[, setdiff(names(base), "condition"), with = FALSE], dp, sep = "\t")

rows <- list()
for (done in seq.int(PERM_FROM, PERM_TO)) {
  ctrl <- sort(all_ctrl[[done]])

  # I/O FIX. This previously rewrote the ENTIRE dataset for every permutation
  # just to change one column. Under a 15-way array that write dominated wall
  # time -- each task managed 6-13 of its 19 permutations in a 6 h walltime
  # despite the in-script timer (which starts after the write) reporting
  # 2.5 min per fit. The data are now written once per task, and only a small
  # sample -> condition map varies; differential_editing() merges it by
  # `sample` when the data carry no condition column, which is equivalent.
  mp <- file.path(tmpdir, sprintf("meta_perm%03d_%d.txt", done, PERM_FROM))
  op <- file.path(tmpdir, sprintf("DRE_perm%03d_%d.txt", done, PERM_FROM))
  fwrite(unique(base[, .(sample, condition = fifelse(library %in% ctrl,
                                                     "control", "dehydrated"))]),
         mp, sep = "\t")

  cat(sprintf("\n[perm %03d/%d] control = %s\n", done, length(all_ctrl), paste(ctrl, collapse = " ")))
  t0 <- Sys.time()
  res <- tryCatch(differential_editing(
      data_path = dp, meta_path = mp, test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "control", case_level = "dehydrated",
      random_effects = "(1 | library) + (1 | cluster_id)",
      out_path = op, summary_path = file.path(tmpdir, sprintf("sum_perm%03d.txt", done)),
      n_cores = N_CORES),
    error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) { unlink(c(mp, op)); next }

  r <- as.data.table(res)
  rows[[length(rows) + 1]] <- data.table(perm = done,
    control_libs = paste(ctrl, collapse = "|"),
    n_compartments = uniqueN(substr(ctrl, 1, 3)),
    GLMM_sig = sum(r$GLMM_sig, na.rm = TRUE),
    Fisher_sig = sum(r$Fisher_sig, na.rm = TRUE),
    Wilcox_sig = sum(r$Wilcox_sig, na.rm = TRUE),
    n_tested_glmm = sum(!is.na(r$glmm_pvalue)),
    med_p_glmm = median(r$glmm_pvalue, na.rm = TRUE),
    med_p_fisher = median(r$fisher_pvalue, na.rm = TRUE),
    ks_glmm = ks_unif(r$glmm_pvalue),
    ks_fisher = ks_unif(r$fisher_pvalue),
    mins = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
  print(rows[[length(rows)]][, .(GLMM_sig, Fisher_sig, Wilcox_sig, mins)])
  unlink(c(mp, op))
  fwrite(rbindlist(rows), paste0("reditr_permutation_null_mouse", CHUNK, ".txt"), sep = "\t")  # checkpoint
}

out <- rbindlist(rows)
fwrite(out, paste0("reditr_permutation_null_mouse", CHUNK, ".txt"), sep = "\t")

# true values, from the analysis as originally run
TRUE_GLMM <- 851L; TRUE_FISHER <- 3291L; TRUE_WILCOX <- 0L
cat("\n\n=== MOUSE PERMUTATION NULL ===\n")
print(out[, .(perm, control_libs, n_compartments, GLMM_sig, Fisher_sig, Wilcox_sig)])
cat(sprintf("\nTRUE: GLMM %d | Fisher %d | Wilcoxon %d\n", TRUE_GLMM, TRUE_FISHER, TRUE_WILCOX))
cat(sprintf("NULL: GLMM mean %.1f range %d-%d | Fisher mean %.1f range %d-%d\n",
            mean(out$GLMM_sig), min(out$GLMM_sig), max(out$GLMM_sig),
            mean(out$Fisher_sig), min(out$Fisher_sig), max(out$Fisher_sig)))
cat(sprintf("empirical p  GLMM = %.3f | Fisher = %.3f\n",
            (sum(out$GLMM_sig >= TRUE_GLMM) + 1) / (nrow(out) + 1),
            (sum(out$Fisher_sig >= TRUE_FISHER) + 1) / (nrow(out) + 1)))
cat(sprintf("empirical FDR (null mean / true)  GLMM = %.0f%% | Fisher = %.0f%%\n",
            100 * mean(out$GLMM_sig) / TRUE_GLMM, 100 * mean(out$Fisher_sig) / TRUE_FISHER))
cat("\nWritten: reditr_permutation_null_mouse.txt\n")
