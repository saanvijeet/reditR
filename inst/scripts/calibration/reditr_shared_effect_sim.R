# reditr_shared_effect_sim.R
#
# Does a SHARED sample effect reproduce the inflation that permutation finds
# and reditR's own simulator does not?
#
# BACKGROUND. reditR::simulate_editing_data() draws its random intercept
# inside the per-site loop:
#
#     rows <- lapply(seq_len(nrow(truth)), function(i) {   # per SITE
#         sample_re <- rnorm(length(all_samples), 0, sample_re_sd)
#         p <- plogis(qlogis(mean_rates) + sample_re)
#
# so a sample that runs hot at one site has no tendency to run hot at the
# next. Measured empirically, that produces a per-sample coherent offset of
# 0.0038 at sample_re_sd = 1.0, against 0.0136 in the real diabetes data and
# 0.1946 in the real mouse data. Real samples carry genome-wide offsets the
# simulator cannot generate at any parameter value — and a shared sample
# effect is precisely what the GLMM's (1 | sample) term exists to absorb.
#
# THIS SCRIPT reproduces reditR's generator exactly, changing one thing: the
# random intercept is drawn ONCE per sample and reused at every site. Counts
# are otherwise identical in construction:
#     mean_rates = baseline_rate (+ effect in the case arm), clipped to [.001,.999]
#     p          = plogis(qlogis(mean_rates) + sample_re)
#     total      = pmax(10, rnbinom(mu = mean_coverage, size = 5))
#     edited     = rbinom(n, total, p)
#
# Both modes are generated at each setting so the comparison is like-for-like,
# and the coherent-offset diagnostic is reported alongside the false-positive
# rate to confirm the intended structure was actually created.
#
# PREDICTION. If the shared effect is the missing ingredient, the shared mode
# should inflate the false-positive rate at sigma values where the per-site
# mode stays near nominal — and should do so at the low sigma that matches the
# real diabetes data, which is where the parametric sweep currently and
# wrongly reports good calibration.
#
# Output: reditr_shared_effect_sim.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })
OUT <- "/rds/general/user/sj1825/home/msc_prj"
setwd(OUT)

SIGMAS  <- as.numeric(strsplit(Sys.getenv("SIGMAS", "0,0.1,0.25,0.5"), ",")[[1]])
NPC     <- as.integer(Sys.getenv("NPC", "6"))
N_REPS  <- as.integer(Sys.getenv("N_REPS", "3"))
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
N_NULL  <- as.integer(Sys.getenv("N_NULL", "800"))
EFFECTS <- c(`0.05` = 100, `0.10` = 100, `0.20` = 100)
BASE    <- as.numeric(Sys.getenv("BASELINE", "0.1"))
COV     <- as.integer(Sys.getenv("COVERAGE", "30"))

tmpdir <- file.path(OUT, "shared_effect_tmp"); dir.create(tmpdir, showWarnings = FALSE)

# reditR's generator, with `shared` controlling where the rnorm is drawn.
sim_data <- function(npc, sigma, shared, seed) {
  set.seed(seed)
  samples <- c(paste0("ctrl_", seq_len(npc)), paste0("case_", seq_len(npc)))
  conds   <- rep(c("control", "diabetic"), each = npc)
  meta    <- data.table(sample = samples, condition = conds)
  truth   <- rbind(data.table(site = paste0("site_null_", seq_len(N_NULL)), true_effect = 0),
                   rbindlist(Map(function(nm, n)
                     data.table(site = paste0("site_eff", nm, "_", seq_len(n)),
                                true_effect = as.numeric(nm)), names(EFFECTS), EFFECTS)))
  re_fixed <- rnorm(length(samples), 0, sigma)     # drawn once; used only if shared
  rows <- rbindlist(lapply(seq_len(nrow(truth)), function(i) {
    eff <- truth$true_effect[i]
    re  <- if (shared) re_fixed else rnorm(length(samples), 0, sigma)
    mr  <- pmax(0.001, pmin(0.999, ifelse(conds == "diabetic", BASE + eff, BASE)))
    p   <- plogis(qlogis(mr) + re)
    tot <- pmax(10L, rnbinom(length(samples), mu = COV, size = 5))
    data.table(site = truth$site[i], sample = samples,
               edited = rbinom(length(samples), tot, p), total = tot)
  }))
  rows[, edit_ratio := edited / total]
  list(editing = rows, metadata = meta, truth = truth)
}

# the diagnostic from the QC: does a sample carry a genome-wide offset?
coherent_offset <- function(e) {
  x <- copy(e)[total >= 10]
  x[, r := edited / total][, dev := r - mean(r), by = .(site, condition)]
  sd(x[, .(off = mean(dev)), by = sample]$off)
}

grid <- CJ(sigma = SIGMAS, shared = c(FALSE, TRUE), rep = seq_len(N_REPS))
cat(sprintf("cells: %d  (npc=%d, baseline=%.2f, coverage=%d)\n\n", nrow(grid), NPC, BASE, COV))
rows <- list()

for (i in seq_len(nrow(grid))) {
  sg <- grid$sigma[i]; sh <- grid$shared[i]; rp <- grid$rep[i]
  key <- sprintf("s%.2f_%s_r%d", sg, ifelse(sh, "shared", "persite"), rp)
  s <- sim_data(NPC, sg, sh, seed = 5000L + i)
  e <- merge(s$editing, s$metadata, by = "sample")
  off <- coherent_offset(e)

  dp <- file.path(tmpdir, paste0("d_", key, ".txt")); mp <- file.path(tmpdir, paste0("m_", key, ".txt"))
  fwrite(s$editing, dp, sep = "\t"); fwrite(s$metadata, mp, sep = "\t")
  res <- tryCatch(differential_editing(
      data_path = dp, meta_path = mp, test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "control", case_level = "diabetic",
      random_effects = "(1 | sample)", out_path = file.path(tmpdir, paste0("r_", key, ".txt")),
      summary_path = file.path(tmpdir, paste0("s_", key, ".txt")),
      n_cores = N_CORES, verbose = FALSE),
    error = function(err) { cat("  FAILED", key, ":", conditionMessage(err), "\n"); NULL })
  unlink(c(dp, mp))
  if (is.null(res)) next

  r <- merge(as.data.table(res), s$truth, by = "site", all.x = TRUE)
  isnull <- r$true_effect == 0
  for (tst in c("GLMM", "Fisher", "Wilcox")) {
    sig <- r[[paste0(tst, "_sig")]]; if (is.null(sig)) next
    rows[[length(rows) + 1]] <- data.table(
      sigma = sg, mode = ifelse(sh, "shared", "per-site"), rep = rp, test = tst,
      coherent_offset = off,
      fpr = mean(sig[isnull], na.rm = TRUE),
      power_10 = mean(sig[abs(r$true_effect - 0.10) < 1e-9], na.rm = TRUE))
  }
  cat(sprintf("[%2d/%2d] %-22s offset=%.4f\n", i, nrow(grid), key, off))
}

out <- rbindlist(rows)
fwrite(out, "reditr_shared_effect_sim.txt", sep = "\t")

cat("\n\n===== coherent per-sample offset (structure check) =====\n")
print(dcast(out[test == "GLMM"], sigma ~ mode, value.var = "coherent_offset",
            fun.aggregate = function(x) round(mean(x), 4)))
cat("  real diabetes = 0.0136 | real mouse = 0.1946\n")
cat("\n===== false-positive rate at FDR<0.05 =====\n")
for (tst in c("GLMM", "Fisher")) {
  cat("--", tst, "\n")
  print(dcast(out[test == tst], sigma ~ mode, value.var = "fpr",
              fun.aggregate = function(x) round(mean(x), 4)))
}
cat("\n===== power at true effect 0.10 =====\n")
print(dcast(out[test == "GLMM"], sigma ~ mode, value.var = "power_10",
            fun.aggregate = function(x) round(mean(x), 3)))
cat("\nWritten: reditr_shared_effect_sim.txt\n")
