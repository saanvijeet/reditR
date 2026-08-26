# reditr_permutation_null_mouse.R
#
# Permutation null analysis for the mouse pan-endothelial-cell dehydration
# dataset using the reditR differential editing workflow.
#
# The aim is to test whether the number of significant sites in the real
# analysis is greater than expected from random assignment of conditions.
#
# The analysis uses the same crossed random-effects model as the real analysis
# so that the permutation test evaluates the workflow that was actually used.
#
# PERMUTATION UNIT:
# Conditions are shuffled at the library level rather than the pseudobulk
# level. There are 192 pseudobulk samples but only 13 biological libraries,
# so the library is the appropriate unit of replication. Shuffling individual
# pseudobulk samples would allow one library to appear in both groups and
# would not provide a valid null comparison.
#
# Three of the 13 libraries are assigned to the control group, matching the
# original 3-control / 10-dehydrated design.
#
# IMPORTANT DATASET LIMITATION:
# Condition is completely confounded with timepoint because all control
# libraries are from 0 h. The control libraries also represent different
# endothelial compartments. Therefore, permuting the labels does not produce
# a completely clean biological null. The result should be interpreted as
# asking whether the observed 3-vs-10 grouping produces more signal than
# arbitrary 3-vs-10 groupings, rather than as definitive evidence of correct
# statistical calibration.
#
# There are C(13,3) = 286 possible control-library assignments. Because the
# groups are unequal, a labelling and its complement are not equivalent.
# Removing the true labelling leaves 285 null permutations.
#
# Output:
#   reditr_permutation_null_mouse.txt

suppressPackageStartupMessages({ library(data.table); library(reditR) })

EPH <- "/rds/general/user/sj1825/ephemeral/mec_dehydration/sprint_output_full"
OUT <- "/rds/general/user/sj1825/home/msc_prj/test_data_mouse"
setwd(OUT)

# Set the number of CPU cores used by differential_editing().
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

# Load the pseudobulk editing data and identify the unique biological
# libraries and their original conditions.
base <- fread(file.path(EPH, "all_ec_clustered_with_condition.txt"))
libs <- unique(base[, .(library, condition)])[order(library)]
true_ctrl <- sort(libs[condition == "control", library])
n_ctrl <- length(true_ctrl)

cat("libraries:", nrow(libs), " true controls:",
    paste(true_ctrl, collapse = ", "), "\n")

# Create a temporary directory for the metadata and result files generated
# for each permutation.
tmpdir <- file.path(EPH, "perm_null_mouse_tmp")
dir.create(tmpdir, showWarnings = FALSE)

# Generate every possible set of control libraries.
# The true control set is removed so that only null permutations remain.
all_ctrl <- local({
  cs <- combn(sort(libs$library), n_ctrl, simplify = FALSE)
  cs[!vapply(
    cs,
    function(c) identical(sort(c), true_ctrl),
    logical(1)
  )]
})

# PERM_FROM and PERM_TO allow the complete set of permutations to be split
# across multiple jobs, which is useful when running the analysis on HPC.
PERM_FROM <- as.integer(Sys.getenv("PERM_FROM", "1"))
PERM_TO   <- min(
  as.integer(Sys.getenv("PERM_TO", as.character(length(all_ctrl)))),
  length(all_ctrl)
)

CHUNK <- if (PERM_FROM == 1L && PERM_TO == length(all_ctrl)) "" else
           sprintf("_p%03d-%03d", PERM_FROM, PERM_TO)

cat(sprintf(
  "Exhaustive space: %d nulls. Running %d-%d.\n",
  length(all_ctrl), PERM_FROM, PERM_TO
))

# Test whether the p-values are consistent with a Uniform(0,1) distribution.
# Under a well-calibrated null, p-values should be approximately uniform,
# so a significant KS test would indicate a departure from this expectation.
ks_unif <- function(p) {
  p <- p[!is.na(p)]
  if (length(p) > 10)
    suppressWarnings(stats::ks.test(p, "punif")$p.value)
  else
    NA_real_
}

# Write the editing data once per job without the condition column.
# The condition assignment will instead be supplied through the small
# permutation-specific metadata files below.
dp <- file.path(tmpdir, sprintf("data_shared_%d.txt", PERM_FROM))
fwrite(
  base[, setdiff(names(base), "condition"), with = FALSE],
  dp,
  sep = "\t"
)

rows <- list()

for (done in seq.int(PERM_FROM, PERM_TO)) {

  # Select the libraries that will be labelled as controls in this permutation.
  ctrl <- sort(all_ctrl[[done]])

  # Create metadata containing the new condition assignment for each sample.
  # Libraries in the selected control set are labelled "control"; all others
  # are labelled "dehydrated".
  mp <- file.path(
    tmpdir,
    sprintf("meta_perm%03d_%d.txt", done, PERM_FROM)
  )

  op <- file.path(
    tmpdir,
    sprintf("DRE_perm%03d_%d.txt", done, PERM_FROM)
  )

  fwrite(
    unique(
      base[, .(
        sample,
        condition = fifelse(
          library %in% ctrl,
          "control",
          "dehydrated"
        )
      )]
    ),
    mp,
    sep = "\t"
  )

  cat(sprintf(
    "\n[perm %03d/%d] control = %s\n",
    done,
    length(all_ctrl),
    paste(ctrl, collapse = " ")
  ))

  t0 <- Sys.time()

  # Run the complete differential editing workflow using the permuted labels.
  # The same statistical tests and random-effects structure are used for every
  # permutation so that only the condition labels differ between analyses.
  res <- tryCatch(
    differential_editing(
      data_path = dp,
      meta_path = mp,
      test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "control",
      case_level = "dehydrated",
      random_effects = "(1 | library) + (1 | cluster_id)",
      out_path = op,
      summary_path = file.path(
        tmpdir,
        sprintf("sum_perm%03d.txt", done)
      ),
      n_cores = N_CORES
    ),
    error = function(e) {
      cat("  FAILED:", conditionMessage(e), "\n")
      NULL
    }
  )

  if (is.null(res)) {
    unlink(c(mp, op))
    next
  }

  r <- as.data.table(res)

  # Record the number of significant sites and basic p-value diagnostics for
  # this permutation. These values form the empirical null distribution.
  rows[[length(rows) + 1]] <- data.table(
    perm = done,
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
    mins = round(
      as.numeric(difftime(Sys.time(), t0, units = "mins")),
      2
    )
  )

  print(
    rows[[length(rows)]][,
      .(GLMM_sig, Fisher_sig, Wilcox_sig, mins)
    ]
  )

  # Remove temporary files before moving to the next permutation.
  unlink(c(mp, op))

  # Save a checkpoint after each permutation so completed results are not
  # lost if the HPC job stops before the full analysis finishes.
  fwrite(
    rbindlist(rows),
    paste0("reditr_permutation_null_mouse", CHUNK, ".txt"),
    sep = "\t"
  )
}

# Combine all permutation results and write the final output file.
out <- rbindlist(rows)

fwrite(
  out,
  paste0("reditr_permutation_null_mouse", CHUNK, ".txt"),
  sep = "\t"
)

# These are the numbers of significant sites obtained from the original
# analysis using the true condition labels.
TRUE_GLMM <- 851L
TRUE_FISHER <- 3291L
TRUE_WILCOX <- 0L

cat("\n\n=== MOUSE PERMUTATION NULL ===\n")

print(
  out[, .(
    perm,
    control_libs,
    n_compartments,
    GLMM_sig,
    Fisher_sig,
    Wilcox_sig
  )]
)

# Compare the observed number of significant sites with the null distribution.
cat(sprintf(
  "\nTRUE: GLMM %d | Fisher %d | Wilcoxon %d\n",
  TRUE_GLMM,
  TRUE_FISHER,
  TRUE_WILCOX
))

cat(sprintf(
  "NULL: GLMM mean %.1f range %d-%d | Fisher mean %.1f range %d-%d\n",
  mean(out$GLMM_sig),
  min(out$GLMM_sig),
  max(out$GLMM_sig),
  mean(out$Fisher_sig),
  min(out$Fisher_sig),
  max(out$Fisher_sig)
))

# Calculate an empirical p-value: the proportion of null permutations that
# produce at least as many significant sites as the true labelling.
# The +1 correction prevents a p-value of exactly zero.
cat(sprintf(
  "empirical p  GLMM = %.3f | Fisher = %.3f\n",
  (sum(out$GLMM_sig >= TRUE_GLMM) + 1) / (nrow(out) + 1),
  (sum(out$Fisher_sig >= TRUE_FISHER) + 1) / (nrow(out) + 1)
))

# Compare the average number of significant sites under the null with the
# observed number. This gives an estimate of how large the null signal is
# relative to the observed result.
cat(sprintf(
  "empirical FDR (null mean / true)  GLMM = %.0f%% | Fisher = %.0f%%\n",
  100 * mean(out$GLMM_sig) / TRUE_GLMM,
  100 * mean(out$Fisher_sig) / TRUE_FISHER
))

cat("\nWritten: reditr_permutation_null_mouse.txt\n")