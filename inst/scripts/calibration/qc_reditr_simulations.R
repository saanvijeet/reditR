# qc_reditr_simulations.R
#
# Quality control for reditR::simulate_editing_data(), which is used for the
# parametric validation of reditR.
#
# This script checks two aspects of the simulator:
#
#   A. PARAMETER RECOVERY — whether the simulated data match the parameters
#      used to generate them, including baseline editing rate, coverage,
#      sample-level variation, sample number and effect size.
#
#   B. REALISM — whether simulated data have similar properties to the real
#      datasets. Dispersion is compared because the simulator assumes
#      binomial variation after accounting for sample-level random effects.
#
# Dispersion is calculated for each site x condition group using the Pearson
# chi-square statistic:
#     p_hat = sum(edited) / sum(total)
#     phi   = Pearson statistic / (n - 1)
#
# A phi value near 1 is consistent with binomial variation, while values above
# 1 indicate additional variation.
#
# Dispersion cannot be estimated for groups with only one observation, so the
# proportion of estimable groups and median group size are also reported.
# endothelial2 is reported both using all sites and using sites detected in
# all 9 samples because many site x condition groups contain only one sample.
#
# Output: qc_reditr_simulations.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })
H <- "/rds/general/user/sj1825/home"
setwd(file.path(H, "msc_prj"))
set.seed(20260812)

dispersion <- function(d) {
  # Calculate the dispersion of editing counts for each site and condition.
  # Groups with fewer than two observations cannot have a variance estimate.
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

  # Estimate sample-level variation on the logit scale using null control sites.
  ctrl <- nulls[condition == "control" & edited > 0 & edited < total]
  lg <- ctrl[, .(l = qlogis(sum(edited) / sum(total))), by = .(site, sample)]
  sd_obs <- lg[, .(s = sd(l)), by = site][, median(s, na.rm = TRUE)]

  # Estimate the realised effect size as the difference between case and
  # control editing rates.
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

# Check that the main simulated parameters are recovered.
cat("\n  expected: samples = 2*npc | sites = 900 | baseline ~ 0.100 | coverage ~ 30\n")
cat("  expected: sd_logit_obs tracks sigma_set | delta_XX ~ its nominal effect size\n")

# =====================================================================
# B. REALISM — dispersion, coverage, editing rate vs the real datasets
# =====================================================================
cat("\n\n=========== B. REALISM vs REAL DATA ===========\n")

# Summarise dispersion and basic coverage/editing-rate properties for each
# dataset. Only groups with at least two observations contribute to phi.
disp_row <- function(nm, e) {
  e <- e[total > 0]
  grp <- e[, .N, by = .(site, condition)]
  data.table(source        = nm,
             phi_median    = median(dispersion(e), na.rm = TRUE),
             pct_estimable = round(100 * mean(grp$N >= 2L), 1),
             n_median      = as.numeric(median(grp$N)),
             cov_median    = as.numeric(median(e$total)),
             ratio_median  = median(e$edited / e$total),
             n_site_cond   = nrow(grp))
}

sim_disp <- rbindlist(lapply(c(0, 0.25, 1.0), function(sg) {
  s <- simulate_editing_data(n_null = 800L, n_effects = c(`0.10` = 100),
                             n_per_condition = 6L, sample_re_sd = sg, seed = 7L)
  e <- merge(as.data.table(s$editing), as.data.table(s$metadata), by = "sample")
  e <- merge(e, as.data.table(s$truth), by = "site")[true_effect == 0]
  disp_row(sprintf("simulated, sigma=%.2f", sg), e)
}))

real <- list()

# Load the real diabetes dataset and its sample metadata.
dia <- merge(fread(file.path(H, "diabetes_output/filtered_sites_clustered_t.txt")),
             fread(file.path(H, "diabetes_output/sample_metadata.txt")), by = "sample")
real[["diabetes cardiomyocytes (6v6)"]] <- dia[, .(site, sample, edited, total, condition)]

# Load the real mouse pseudobulk dataset.
mo <- fread("/rds/general/user/sj1825/ephemeral/mec_dehydration/sprint_output_full/all_ec_clustered_with_condition.txt")
real[["mouse pseudobulk"]] <- mo[, .(site, sample, edited, total, condition)]

# Load endothelial2. It is also analysed using a subset containing sites
# detected in all 9 samples, because dispersion requires repeated observations.
en <- merge(fread(file.path(H, "endothelial2/output/filtered_sites_clustered.txt")),
            fread(file.path(H, "endothelial2/output/sample_metadata.txt")), by = "sample")
en <- en[, .(site, sample, edited, total, condition)]
real[["endothelial2 (3x3), all sites"]] <- en
complete <- en[total > 0, uniqueN(sample), by = site][V1 == uniqueN(en$sample), site]
real[["endothelial2 (3x3), 9-sample sites"]] <- en[site %in% complete]

# Combine the simulated and real-data summaries.
real_disp <- rbindlist(lapply(names(real), function(nm) disp_row(nm, real[[nm]])))

disp <- rbind(sim_disp, real_disp)
disp[, `:=`(phi_median = round(phi_median, 2), ratio_median = round(ratio_median, 3))]
print(disp)
report$dispersion <- disp

# Explain how to interpret the dispersion and why the estimable proportion
# is reported alongside it.
cat("\n  phi ~ 1 means residual variance is binomial, as the simulator assumes.\n")
cat("  phi >> simulated means real data is harder than the simulation, and the\n")
cat("  sweep's false-positive rates are LOWER BOUNDS on the real ones.\n")
cat("\n  Read phi together with pct_estimable. endothelial2 has 3 replicates per\n")
cat("  arm and sites called per sample, so most site x condition groups hold a\n")
cat("  single observation and phi does not exist for them. Restricting to sites\n")
cat("  seen in all 9 samples makes phi estimable throughout but selects the\n")
cat("  best-covered sites, and the two subsets straddle phi = 1. Neither is\n")
cat("  precise enough to place against the simulated range on its own.\n")

# Save the parameter-recovery and dispersion summaries.
fwrite(rec,  "qc_reditr_simulations.txt", sep = "\t")
fwrite(disp, "qc_reditr_simulations.txt", sep = "\t", append = TRUE, col.names = TRUE)
cat("\nWritten: qc_reditr_simulations.txt\n")