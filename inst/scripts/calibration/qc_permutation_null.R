# qc_permutation_null.R
#
# QC on the PERMUTATION arm of the reditR validation.
#
# The existing checks (TRUE row reproduces the published analysis; labellings
# distinct, balanced and non-mirrored; permuted counts vary with the labels)
# verify the plumbing. None of them tests whether the permutation PROCEDURE is
# calibrated -- whether the empirical p it returns means what it claims.
#
# Two controls on data with known truth, using the package's own simulator so
# the answer is not in question:
#
#   NEGATIVE  no planted effect. The true labelling is then just one draw from
#             the same null as the other nine, so the empirical p must be
#             UNIFORM over its 10 attainable values {0.1,...,1.0}. p <= 0.1
#             should occur ~10% of the time. Materially more than that means
#             the procedure treats the real labelling as special when it is
#             not -- i.e. the permutation null itself is anti-conservative and
#             the dataset results built on it are overstated.
#
#   POSITIVE  strong planted effect in every site. The true labelling should
#             rank first essentially always (p = 0.100, the floor). If it does
#             not, the procedure has no power at this design and a null result
#             from it carries no information.
#
# Design mirrors endothelial2 exactly: 3 v 3, so C(6,3) = 20 labellings, 10
# after collapsing complements, 9 usable nulls, floor 1/10 = 0.100.
#
# Output: qc_permutation_null.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })
setwd("/rds/general/user/sj1825/home/msc_prj")

N_REP   <- as.integer(Sys.getenv("N_REP",   "20"))
N_SITES <- as.integer(Sys.getenv("N_SITES", "300"))
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))
SIGMA   <- as.numeric(Sys.getenv("SIGMA",   "0.25"))
tmp <- file.path(getwd(), "qc_perm_tmp"); dir.create(tmp, showWarnings = FALSE)

# the 9 non-trivial labellings for a 3v3 design, complements collapsed
labsets <- local({
  truth <- c(rep("control", 3), rep("diabetic", 3))
  seen <- character(0); out <- list()
  for (cb in combn(6, 3, simplify = FALSE)) {
    lab <- rep("control", 6); lab[cb] <- "diabetic"
    k <- paste(lab, collapse = ""); ck <- paste(ifelse(lab == "control", "diabetic", "control"), collapse = "")
    if (k %in% seen || ck %in% seen) next
    seen <- c(seen, k); out[[length(out) + 1]] <- lab
  }
  out[!vapply(out, function(l) all(l == truth) || all(l != truth), logical(1))]
})
stopifnot(length(labsets) == 9L)

one_rep <- function(arm, rep_i) {
  seed <- 4200L + rep_i + if (arm == "positive") 1000L else 0L
  eff  <- if (arm == "positive") c(`0.20` = N_SITES) else c(`0.20` = 0L)
  sim  <- simulate_editing_data(
    n_null = if (arm == "positive") 0L else N_SITES,
    n_effects = eff, n_per_condition = 3L, sample_re_sd = SIGMA, seed = seed)
  dp <- file.path(tmp, sprintf("d_%s_%02d.txt", arm, rep_i))
  fwrite(sim$editing, dp, sep = "\t")

  count <- function(lab, tag) {
    mp <- file.path(tmp, sprintf("m_%s_%02d_%s.txt", arm, rep_i, tag))
    fwrite(data.table(sample = sim$metadata$sample, condition = lab), mp, sep = "\t")
    r <- tryCatch(as.data.table(differential_editing(
           data_path = dp, meta_path = mp, test = "glmm",
           reference_level = "control", case_level = "case",
           random_effects = "(1 | sample)", n_cores = N_CORES, verbose = FALSE)),
         error = function(e) NULL)
    if (is.null(r)) NA_integer_ else sum(r$GLMM_sig, na.rm = TRUE)
  }

  truth_lab <- as.character(sim$metadata$condition)
  obs  <- count(truth_lab, "TRUE")
  nul  <- vapply(seq_along(labsets), function(i) count(labsets[[i]], sprintf("%02d", i)), integer(1))
  nul  <- nul[!is.na(nul)]
  data.table(arm = arm, rep = rep_i, true_n = obs,
             null_mean = mean(nul), null_max = max(nul), n_nulls = length(nul),
             emp_p = (sum(nul >= obs) + 1) / (length(nul) + 1))
}

res <- rbindlist(lapply(c("negative", "positive"), function(a)
         rbindlist(lapply(seq_len(N_REP), function(i) {
           x <- one_rep(a, i); cat(sprintf("%-8s rep %2d  true %4d  null mean %6.1f  p %.3f\n",
                                           a, i, x$true_n, x$null_mean, x$emp_p)); x }))))
fwrite(res, "qc_permutation_null.txt", sep = "\t")

cat("\n================ PERMUTATION QC ================\n")
for (a in c("negative", "positive")) {
  d <- res[arm == a]
  cat(sprintf("\n-- %s control (%d reps, %d sites, sigma %.2f) --\n", a, nrow(d), N_SITES, SIGMA))
  cat(sprintf("   true count      : median %.0f (range %d-%d)\n", median(d$true_n), min(d$true_n), max(d$true_n)))
  cat(sprintf("   null mean count : median %.1f\n", median(d$null_mean)))
  cat(sprintf("   empirical p     : median %.3f | mean %.3f\n", median(d$emp_p), mean(d$emp_p)))
  cat(sprintf("   frac p <= 0.100 : %.1f%%   (expect ~10%% negative / ~100%% positive)\n",
              100 * mean(d$emp_p <= 0.1000001)))
}
cat("\nWritten: qc_permutation_null.txt\n")
