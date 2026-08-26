# reditr_parametric_null.R
#
# Parametric simulation to assess the false-positive rate and power of reditR.
#
# This complements the permutation analyses. Permutations use the real data
# and therefore preserve the variability present in each dataset, whereas
# this simulation allows sample-level variability to be deliberately changed.
# This lets us test how reditR performs as between-sample variability increases.
#
# The simulation uses reditR's own simulator, which creates both:
#   - null sites with no true condition effect
#   - sites with known effect sizes
#
# This allows false-positive rate and power to be measured against a known
# ground truth.
#
# The simulation varies:
#   sample_re_sd     = between-sample variability
#   n_per_condition = number of samples in each condition
#   replicates       = independent simulation repeats
#
# Important limitation:
# The simulator uses a binomial model with a sample-level random effect.
# Therefore, this assesses reditR under the assumptions of its own simulator
# and does not reproduce every source of variability in real RNA-seq data.
# The parametric and permutation analyses therefore provide complementary
# assessments of performance.
#
# Output:
#   reditr_parametric_null.txt      Summary results
#   reditr_parametric_null_raw.rds  Per-site results for further analysis


suppressPackageStartupMessages({
  library(data.table)
  library(reditR)
})


# Set the directory containing the simulation outputs.
OUT <- "/rds/general/user/sj1825/home/msc_prj"
setwd(OUT)


# Define the simulation settings.
#
# SIGMAS controls the amount of sample-to-sample variability.
# NPC controls the number of samples per condition.
# N_REPS repeats each combination using independent random seeds.
# More repeats give a more stable estimate of the average performance.
SIGMAS  <- as.numeric(strsplit(
  Sys.getenv("SIGMAS", "0,0.1,0.25,0.5,1.0"), ",")[[1]])

NPC     <- as.integer(strsplit(
  Sys.getenv("NPC", "3,6,10"), ",")[[1]])

N_REPS  <- as.integer(Sys.getenv("N_REPS", "3"))

# Number of CPU cores used by differential_editing().
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

# Number of sites with no true condition effect.
N_NULL  <- as.integer(Sys.getenv("N_NULL", "800"))

# Number of sites simulated at each known effect size.
EFFECTS <- c(`0.05` = 100, `0.10` = 100, `0.20` = 100)


# Create a temporary directory for the input and output files generated
# during each simulation.
tmpdir <- file.path(OUT, "reditr_parametric_tmp")
dir.create(tmpdir, showWarnings = FALSE)


# Create every combination of variability, sample size and replicate.
# Each row represents one independent simulation to be performed.
grid <- CJ(
  sigma = SIGMAS,
  npc   = NPC,
  rep   = seq_len(N_REPS)
)

cat(sprintf(
  "cells: %d (%d sigma x %d designs x %d reps)\n\n",
  nrow(grid), length(SIGMAS), length(NPC), N_REPS
))


# Store summary results and raw per-site results separately.
rows <- list()
raw  <- list()


# Run each simulation setting in the grid.
for (i in seq_len(nrow(grid))) {

  # Extract the settings for the current simulation.
  sg <- grid$sigma[i]
  np <- grid$npc[i]
  rp <- grid$rep[i]

  # Create a unique name for the current simulation.
  key <- sprintf("s%.2f_n%d_r%d", sg, np, rp)

  # Give each simulation a different seed so that replicates are independent
  # while remaining reproducible.
  seed <- 20260812L + i


  # Generate synthetic editing data with known ground truth.
  #
  # The simulator returns:
  #   editing  = editing counts for each site and sample
  #   metadata = sample condition information
  #   truth    = the known effect for each simulated site
  #
  # sample_re_sd controls how much the editing rate can vary between samples.
  sim <- simulate_editing_data(
    n_null = N_NULL,
    n_effects = EFFECTS,
    n_per_condition = np,
    sample_re_sd = sg,
    seed = seed
  )


  # Identify the simulated case condition.
  # "control" is used as the reference condition, while the other condition
  # is treated as the case condition.
  case_lv <- setdiff(
    unique(as.character(sim$metadata$condition)),
    "control"
  )


  # Save the simulated editing data and metadata as temporary files.
  # differential_editing() accepts file paths as its input.
  dp <- file.path(tmpdir, paste0("dat_", key, ".txt"))
  mp <- file.path(tmpdir, paste0("meta_", key, ".txt"))

  fwrite(sim$editing, dp, sep = "\t")
  fwrite(sim$metadata, mp, sep = "\t")


  # Run all three differential-editing tests on the simulated data.
  #
  # The GLMM uses sample as a random effect so that sample-level variation
  # is accounted for in the statistical model.
  #
  # tryCatch() prevents one failed simulation from stopping the entire sweep.
  res <- tryCatch(
    differential_editing(
      data_path = dp,
      meta_path = mp,
      test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "control",
      case_level = case_lv,
      random_effects = "(1 | sample)",
      out_path = file.path(tmpdir, paste0("res_", key, ".txt")),
      summary_path = file.path(tmpdir, paste0("sum_", key, ".txt")),
      n_cores = N_CORES,
      verbose = FALSE
    ),
    error = function(e) {
      cat(
        "  FAILED", key, ":",
        conditionMessage(e), "\n"
      )
      NULL
    }
  )

  # Skip this simulation if the analysis failed.
  if (is.null(res)) next


  # Combine the statistical results with the known simulated effect
  # for each site. This allows detected sites to be compared with the
  # true underlying effect.
  r <- merge(
    as.data.table(res),
    as.data.table(sim$truth),
    by = "site",
    all.x = TRUE
  )


  # Keep the full per-site results so that they can be inspected later.
  raw[[key]] <- cbind(
    sigma = sg,
    npc = np,
    rep = rp,
    r
  )


  # Identify sites that were deliberately simulated with no condition effect.
  # These sites are used to calculate the false-positive rate.
  is_null_site <- r$true_effect == 0


  # Calculate performance separately for each statistical test.
  for (tst in c("GLMM", "Fisher", "Wilcox")) {

    # Extract the significance result for the current test.
    sig <- r[[paste0(tst, "_sig")]]

    if (is.null(sig)) next


    # Calculate power for a particular known effect size.
    # Power is the proportion of sites with that true effect that were
    # correctly identified as significant.
    pw <- function(e) {
      mean(
        sig[abs(r$true_effect - e) < 1e-9],
        na.rm = TRUE
      )
    }


    # Store the main performance measures for this simulation.
    rows[[length(rows) + 1]] <- data.table(
      sigma = sg,
      npc = np,
      rep = rp,
      test = tst,

      # Total number of sites analysed.
      n_sites = nrow(r),

      # Proportion of sites for which a p-value was successfully produced.
      # For GLMM this represents the convergence rate.
      converged = mean(
        !is.na(
          r[[c(
            GLMM = "glmm_pvalue",
            Fisher = "fisher_pvalue",
            Wilcox = "wilcox_pvalue"
          )[tst]]]
        )
      ),

      # False-positive rate among sites where the true effect was zero.
      fpr = mean(
        sig[is_null_site],
        na.rm = TRUE
      ),

      # Overall power across all sites with a true non-zero effect.
      power_all = mean(
        sig[!is_null_site],
        na.rm = TRUE
      ),

      # Power at each of the predefined effect sizes.
      power_05 = pw(0.05),
      power_10 = pw(0.10),
      power_20 = pw(0.20)
    )
  }


  # Report progress so that long simulations can be monitored.
  cat(sprintf(
    "[%2d/%2d] %-16s done\n",
    i, nrow(grid), key
  ))


  # Remove the temporary input files before moving to the next simulation.
  unlink(c(dp, mp))
}


# Combine the results from all simulations into one table.
out <- rbindlist(rows, fill = TRUE)


# Save the summary results as a tab-separated file for later analysis.
fwrite(
  out,
  "reditr_parametric_null.txt",
  sep = "\t"
)


# Save the complete per-site results so that the simulations can be
# re-analysed without having to repeat the computationally expensive runs.
saveRDS(
  rbindlist(raw, fill = TRUE),
  "reditr_parametric_null_raw.rds"
)


# Print the average false-positive rate for each combination of
# sample-level variability, sample size and statistical test.
#
# This is the main calibration measure: a well-calibrated test should have
# a false-positive rate close to the expected level under the null.
cat("\n\n=========== reditR parametric sweep: FPR at FDR<0.05 (nominal ~0) ===========\n")

print(
  dcast(
    out,
    sigma + npc ~ test,
    value.var = "fpr",
    fun.aggregate = function(x) round(mean(x), 4)
  )
)


# Print average power for sites with a true effect of 0.10.
# This provides a simple comparison of how well each test detects
# a moderate simulated effect under different conditions.
cat("\n=========== power at true effect = 0.10 ===========\n")

print(
  dcast(
    out,
    sigma + npc ~ test,
    value.var = "power_10",
    fun.aggregate = function(x) round(mean(x), 3)
  )
)


# Report where the main summary results were saved.
cat("\nWritten: reditr_parametric_null.txt\n")