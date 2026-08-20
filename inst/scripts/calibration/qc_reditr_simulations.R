# qc_reditr_simulations.R
#
# Quality control on reditR::simulate_editing_data(), the generator behind the
# parametric arm of the reditR validation.
#
# WHY. The parametric sweep tests reditR against its OWN generative model. If
# that model does not resemble real editing data, a clean pass there means
# little — and, more importantly, the false-positive rates it reports would be
# LOWER BOUNDS. This script asks two separate questions:
#
#   A. PARAMETER RECOVERY — does the simulator produce what it was asked for?
#      baseline_rate, mean_coverage, sample_re_sd, n_per_condition, and the
#      declared true_effect must all be recoverable from the output. This is a
#      check on the simulator itself.
#
#   B. REALISM — does simulated data resemble the real datasets? The decisive
#      quantity is OVERDISPERSION. The simulator draws counts as binomial
#      conditional on a per-sample random intercept, so within a condition,
#      once the random effect is accounted for, residual variance should be
#      binomial (dispersion phi ~ 1). Real editing data may be more variable
#      than that. If real phi >> simulated phi, the simulation is an easier
#      problem than reality and its FPRs understate the true ones.
#
# Dispersion is measured per site x condition as the Pearson chi-square
# statistic over its degrees of freedom:
#     p_hat = sum(edited)/sum(total)
#     phi   = sum( (edited_j - total_j*p_hat)^2 / (total_j*p_hat*(1-p_hat)) ) / (n-1)
# phi ~ 1 is binomial; phi > 1 is overdispersed.
#
# Output: qc_reditr_simulations.txt (tables), console report

suppressPackageStartupMessages({ library(data.table); library(reditR) })
H <- "/rds/general/user/sj1825/home"
setwd(file.path(H, "msc_prj"))
set.seed(20260812)

dispersion <- function(d) {
  # d: site, sample, edited, total, condition
  d <- d[total > 0]
  d[, {
    n <- .N
    if (n < 2L) NA_real_ else {
      ph <- sum(edited) / sum(total)
      if (ph <= 0 || ph >= 1) NA_real_ else
        sum((edited - total * ph)^2 / (total * ph * (1 - ph))) / (n - 1)
    }
  }, by = .(site, condition)]$V1
}

report <- list()

# =====================================================================
# A. PARAMETER RECOVERY
# =====================================================================
cat("=========== A. PARAMETER RECOVERY ===========\n")
grid <- CJ(sigma = c(0, 0.25, 1.0), npc = c(3L, 6L))
rec <- rbindlist(lapply(seq_len(nrow(grid)), function(i) {
  sg <- grid$sigma[i]; np <- grid$npc[i]
  s <- simulate_editing_data(n_null = 600L, n_effects = c(`0.05` = 100, `0.10` = 100, `0.20` = 100),
                             n_per_condition = np, baseline_rate = 0.1,
                             mean_coverage = 30L, sample_re_sd = sg, seed = 100L + i)
  e <- merge(as.data.table(s$editing), as.data.table(s$metadata), by = "sample")
  e <- merge(e, as.data.table(s$truth), by = "site")
  nulls <- e[true_effect == 0]

  # between-sample SD on the logit scale, from control samples of null sites
  ctrl <- nulls[condition == "control" & edited > 0 & edited < total]
  lg <- ctrl[, .(l = qlogis(sum(edited) / sum(total))), by = .(site, sample)]
  sd_obs <- lg[, .(s = sd(l)), by = site][, median(s, na.rm = TRUE)]

  # realised effect: case minus control mean ratio, per declared effect size
  eff <- e[true_effect > 0, .(r = sum(edited) / sum(total)), by = .(site, condition, true_effect)]
  eff <- dcast(eff, site + true_effect ~ condition, value.var = "r")
  cl <- setdiff(names(eff), c("site", "true_effect", "control"))
  eff[, delta := get(cl) - control]

  data.table(sigma_set = sg, npc_set = np,
             samples_obs = uniqueN(e$sample),
             sites_obs = uniqueN(e$site),
             baseline_obs = round(nulls[, sum(edited) / sum(total)], 4),
             coverage_obs = round(mean(e$total), 1),
             sd_logit_obs = round(sd_obs, 3),
             delta_05 = round(eff[abs(true_effect - 0.05) < 1e-9, median(delta)], 4),
             delta_10 = round(eff[abs(true_effect - 0.10) < 1e-9, median(delta)], 4),
             delta_20 = round(eff[abs(true_effect - 0.20) < 1e-9, median(delta)], 4))
}))
print(rec)
report$recovery <- rec
cat("\n  expected: samples = 2*npc | sites = 900 | baseline ~ 0.100 | coverage ~ 30\n")
cat("  expected: sd_logit_obs tracks sigma_set | delta_XX ~ its nominal effect size\n")

# =====================================================================
# B. REALISM — dispersion, coverage, editing rate vs the real datasets
# =====================================================================
cat("\n\n=========== B. REALISM vs REAL DATA ===========\n")

sim_disp <- rbindlist(lapply(c(0, 0.25, 1.0), function(sg) {
  s <- simulate_editing_data(n_null = 800L, n_effects = c(`0.10` = 100),
                             n_per_condition = 6L, sample_re_sd = sg, seed = 7L)
  e <- merge(as.data.table(s$editing), as.data.table(s$metadata), by = "sample")
  e <- merge(e, as.data.table(s$truth), by = "site")[true_effect == 0]
  data.table(source = sprintf("simulated, sigma=%.2f", sg),
             phi_median = median(dispersion(e), na.rm = TRUE),
             cov_median = median(e$total), ratio_median = median(e$edited / e$total),
             n_site_cond = uniqueN(e[, .(site, condition)]))
}))

real <- list()
dia <- merge(fread(file.path(H, "diabetes_output/filtered_sites_clustered_t.txt")),
             fread(file.path(H, "diabetes_output/sample_metadata.txt")), by = "sample")
real[["diabetes cardiomyocytes (6v6)"]] <- dia[, .(site, sample, edited, total, condition)]

mo <- fread("/rds/general/user/sj1825/ephemeral/mec_dehydration/sprint_output_full/all_ec_clustered_with_condition.txt")
real[["mouse pseudobulk"]] <- mo[, .(site, sample, edited, total, condition)]

real_disp <- rbindlist(lapply(names(real), function(nm) {
  e <- real[[nm]]
  data.table(source = nm, phi_median = median(dispersion(e), na.rm = TRUE),
             cov_median = as.numeric(median(e$total)),
             ratio_median = median(e$edited / e$total),
             n_site_cond = uniqueN(e[, .(site, condition)]))
}))

disp <- rbind(sim_disp, real_disp)
disp[, `:=`(phi_median = round(phi_median, 2), ratio_median = round(ratio_median, 3))]
print(disp)
report$dispersion <- disp

cat("\n  phi ~ 1 means residual variance is binomial, as the simulator assumes.\n")
cat("  phi >> simulated means real data is harder than the simulation, and the\n")
cat("  sweep's false-positive rates are LOWER BOUNDS on the real ones.\n")

fwrite(rec,  "qc_reditr_simulations.txt", sep = "\t")
fwrite(disp, "qc_reditr_simulations.txt", sep = "\t", append = TRUE, col.names = TRUE)
cat("\nWritten: qc_reditr_simulations.txt\n")
