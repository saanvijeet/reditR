# reditr_parametric_null.R
#
# Parametric null + power validation of the reditR PACKAGE, using the
# package's own simulator.
#
# WHY THIS EXISTS. The three reditr_permutation_null_*.R scripts are
# permutation-only. Permutation is the harder test — the null data are the
# real data — but it can only probe whatever between-sample variance the real
# dataset happens to carry. It cannot answer "how does reditR's false-positive
# rate respond as that variance grows?", which is the question separating
# "reditR is miscalibrated" from "reditR is fine but this design carries
# variance it cannot absorb". The parametric sweeps that do answer it live in
# the null_calibration_sim_*.R family, which tests the PREDECESSOR pipelines,
# not reditR.
#
# reditR ships the right instruments and they had not been used:
#   simulate_editing_data(n_null, n_effects, n_per_condition, baseline_rate,
#                         mean_coverage, sample_re_sd, seed)
#   validate_against_truth(results, truth, fdr_threshold)
#
# simulate_editing_data() generates true-null sites AND sites with known
# effect sizes in one dataset, so type-I error and power are measured
# together — which is what makes a fair power comparison possible at all.
# Comparing a 33%-FPR test against a 0.3%-FPR test on power is meaningless
# without matched calibration; the truth labels supply that.
#
# SWEEP
#   sample_re_sd    0, 0.1, 0.25, 0.5, 1.0   (spans the mimic-pipeline grid;
#                   avoids the collapsed {0,1} grid of the diabetes sweep)
#   n_per_condition 3, 6, 10                  (endothelial2, diabetes, hypothetical)
#   replicates      N_REPS independent seeds per cell
#
# CAVEAT TO CARRY. This validates reditR against its OWN generative model:
# binomial counts conditional on a per-sample random intercept, no extra
# overdispersion. A clean pass here would not fully exonerate the package if
# real editing data are overdispersed beyond that. Permutation on real data
# remains the harder test; the two are complementary.
#
# Output: reditr_parametric_null.txt          (one row per cell x test)
#         reditr_parametric_null_raw.rds      (per-site results, for re-scoring)

suppressPackageStartupMessages({
  library(data.table); library(reditR)
})

OUT <- "/rds/general/user/sj1825/home/msc_prj"
setwd(OUT)

SIGMAS  <- as.numeric(strsplit(Sys.getenv("SIGMAS", "0,0.1,0.25,0.5,1.0"), ",")[[1]])
NPC     <- as.integer(strsplit(Sys.getenv("NPC", "3,6,10"), ",")[[1]])
N_REPS  <- as.integer(Sys.getenv("N_REPS", "3"))
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
N_NULL  <- as.integer(Sys.getenv("N_NULL", "800"))
EFFECTS <- c(`0.05` = 100, `0.10` = 100, `0.20` = 100)

tmpdir <- file.path(OUT, "reditr_parametric_tmp"); dir.create(tmpdir, showWarnings = FALSE)
grid <- CJ(sigma = SIGMAS, npc = NPC, rep = seq_len(N_REPS))
cat(sprintf("cells: %d (%d sigma x %d designs x %d reps)\n\n", nrow(grid),
            length(SIGMAS), length(NPC), N_REPS))

rows <- list(); raw <- list()

for (i in seq_len(nrow(grid))) {
  sg <- grid$sigma[i]; np <- grid$npc[i]; rp <- grid$rep[i]
  key <- sprintf("s%.2f_n%d_r%d", sg, np, rp)
  seed <- 20260812L + i

  # simulate_editing_data() returns list(editing, metadata, truth):
  #   editing  site, sample, edited, total, edit_ratio   (no condition column)
  #   metadata sample, condition  -- levels are "control" / "diabetic"
  #   truth    site, true_effect  -- 0 for nulls, else the simulated delta
  sim <- simulate_editing_data(n_null = N_NULL, n_effects = EFFECTS,
                               n_per_condition = np, sample_re_sd = sg, seed = seed)
  case_lv <- setdiff(unique(as.character(sim$metadata$condition)), "control")

  dp <- file.path(tmpdir, paste0("dat_", key, ".txt"))
  mp <- file.path(tmpdir, paste0("meta_", key, ".txt"))
  fwrite(sim$editing, dp, sep = "\t"); fwrite(sim$metadata, mp, sep = "\t")

  res <- tryCatch(differential_editing(
      data_path = dp, meta_path = mp, test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "control", case_level = case_lv,
      random_effects = "(1 | sample)",
      out_path = file.path(tmpdir, paste0("res_", key, ".txt")),
      summary_path = file.path(tmpdir, paste0("sum_", key, ".txt")),
      n_cores = N_CORES, verbose = FALSE),
    error = function(e) { cat("  FAILED", key, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) next

  r <- merge(as.data.table(res), as.data.table(sim$truth), by = "site", all.x = TRUE)
  raw[[key]] <- cbind(sigma = sg, npc = np, rep = rp, r)
  is_null_site <- r$true_effect == 0

  for (tst in c("GLMM", "Fisher", "Wilcox")) {
    sig <- r[[paste0(tst, "_sig")]]
    if (is.null(sig)) next
    pw <- function(e) mean(sig[abs(r$true_effect - e) < 1e-9], na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.table(
      sigma = sg, npc = np, rep = rp, test = tst,
      n_sites    = nrow(r),
      converged  = mean(!is.na(r[[c(GLMM="glmm_pvalue",Fisher="fisher_pvalue",Wilcox="wilcox_pvalue")[tst]]])),
      fpr        = mean(sig[is_null_site], na.rm = TRUE),
      power_all  = mean(sig[!is_null_site], na.rm = TRUE),
      power_05   = pw(0.05), power_10 = pw(0.10), power_20 = pw(0.20))
  }
  cat(sprintf("[%2d/%2d] %-16s done\n", i, nrow(grid), key))
  unlink(c(dp, mp))
}

out <- rbindlist(rows, fill = TRUE)
fwrite(out, "reditr_parametric_null.txt", sep = "\t")
saveRDS(rbindlist(raw, fill = TRUE), "reditr_parametric_null_raw.rds")

cat("\n\n=========== reditR parametric sweep: FPR at FDR<0.05 (nominal ~0) ===========\n")
print(dcast(out, sigma + npc ~ test, value.var = "fpr", fun.aggregate = function(x) round(mean(x), 4)))
cat("\n=========== power at true effect = 0.10 ===========\n")
print(dcast(out, sigma + npc ~ test, value.var = "power_10", fun.aggregate = function(x) round(mean(x), 3)))
cat("\nWritten: reditr_parametric_null.txt\n")
