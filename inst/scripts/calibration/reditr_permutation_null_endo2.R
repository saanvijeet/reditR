# reditr_permutation_null_endo2.R
#
# Permutation null analysis for the endothelial2 ADAR-knockdown dataset.
# This tests whether the significant sites identified by reditR are greater
# than expected when the sample labels are randomly reassigned.
#
# siADAR1 is expected to produce a strong reduction in RNA editing because
# ADAR1 is directly involved in A-to-I RNA editing. siADAR2 provides a second
# condition where a much smaller effect is expected.
#
# Because each arm contains only 6 samples (3 per condition), all possible
# balanced label assignments can be tested rather than using random sampling.
# This provides an exact permutation null for these data.
#
# The true labelling is analysed separately so that the number of significant
# sites under the real biological conditions can be compared with the null
# distribution produced by permuting the labels.
#
# Output:
#   reditr_permutation_null_endo2_<arm>.txt
#   reditr_permutation_null_endo2.txt

suppressPackageStartupMessages({
  library(data.table)
  library(reditR)
})

D <- "/rds/general/user/sj1825/home/endothelial2/output"
setwd(D)
N_CORES <- as.integer(Sys.getenv("N_CORES", "8"))

# Load the sample metadata and filtered editing sites used for the analysis.
# These provide the sample conditions and the editing counts tested by reditR.
full_meta <- fread("sample_metadata.txt")
sites_all <- fread("filtered_sites_clustered.txt")

# Temporary files are used to provide each permutation with its own metadata
# file and output file without changing the original analysis data.
tmpdir <- file.path(D, "perm_null_endo2_tmp")
dir.create(tmpdir, showWarnings = FALSE)

run_arm <- function(case) {

  # Keep only the reference samples (scr) and the condition being tested.
  # This allows the same function to analyse either siADAR1 or siADAR2.
  keep <- full_meta[condition %in% c("scr", case)]
  s <- sites_all[sample %in% keep$sample]

  # Write the editing data to a temporary file for differential_editing().
  dp <- file.path(tmpdir, sprintf("sites_%s.txt", case))
  fwrite(s, dp, sep = "\t")

  samples <- keep$sample
  truth <- keep$condition
  n_case <- sum(truth == case)

  # Generate every possible balanced assignment of samples to the two
  # conditions. With 3 samples per group, this gives C(6,3) = 20 assignments.
  combs <- combn(seq_along(samples), n_case, simplify = FALSE)

  # A labelling and its exact reverse represent the same two-group partition,
  # so keep only one version of each pair to avoid testing duplicates.
  seen <- character(0)
  labsets <- list()

  for (cb in combs) {
    lab <- rep("scr", length(samples))
    lab[cb] <- case

    key <- paste(lab, collapse = "")
    ckey <- paste(ifelse(lab == "scr", case, "scr"), collapse = "")

    if (key %in% seen || ckey %in% seen)
      next

    seen <- c(seen, key)
    labsets[[length(labsets) + 1]] <- lab
  }

  # Remove the real biological labelling from the permutation null.
  # The remaining assignments represent situations where the labels contain
  # no information about the true condition.
  is_true <- vapply(
    labsets,
    function(l) all(l == truth) || all(l != truth),
    logical(1)
  )

  nulls <- labsets[!is_true]

  cat("\n#### arm:", case, "-", length(nulls),
      "exhaustive null labellings ####\n")

  rows <- list()

  # Run the complete differential editing analysis for every null labelling.
  # Keeping the analysis unchanged means that differences in the number of
  # significant sites are due to the labels rather than changes to the method.
  for (i in seq_along(nulls)) {

    lab <- nulls[[i]]

    mp <- file.path(
      tmpdir,
      sprintf("meta_%s_%02d.txt", case, i)
    )

    op <- file.path(
      tmpdir,
      sprintf("DRE_%s_%02d.txt", case, i)
    )

    # Create metadata containing the permuted condition labels.
    fwrite(
      data.table(sample = samples, condition = lab),
      mp,
      sep = "\t"
    )

    cat(sprintf(
      "[%s null %d/%d] %s\n",
      case,
      i,
      length(nulls),
      paste(ifelse(lab == "scr", "R", "C"), collapse = " ")
    ))

    # Run all three statistical tests using the same settings as the
    # original analysis.
    res <- tryCatch(
      differential_editing(
        data_path = dp,
        meta_path = mp,
        test = c("glmm", "fisher", "wilcoxon"),
        reference_level = "scr",
        case_level = case,
        random_effects = "(1 | sample)",
        out_path = op,
        summary_path = file.path(
          tmpdir,
          sprintf("sum_%s_%02d.txt", case, i)
        ),
        n_cores = N_CORES
      ),
      error = function(e) {
        cat("  FAILED:", conditionMessage(e), "\n")
        NULL
      }
    )

    if (is.null(res))
      next

    r <- as.data.table(res)

    # Count how many sites are significant under this null labelling.
    # These counts form the empirical null distribution against which the
    # results from the true labels can be compared.
    rows[[length(rows) + 1]] <- data.table(
      arm = case,
      null = i,
      labels = paste(ifelse(lab == "scr", "R", "C"), collapse = ""),
      GLMM_sig = sum(r$GLMM_sig, na.rm = TRUE),
      Fisher_sig = sum(r$Fisher_sig, na.rm = TRUE),
      Wilcox_sig = sum(r$Wilcox_sig, na.rm = TRUE)
    )

    print(
      rows[[length(rows)]][,
        .(GLMM_sig, Fisher_sig, Wilcox_sig)
      ]
    )

    # Remove the temporary output before moving to the next permutation.
    unlink(op)
  }

  # Analyse the real condition labels using exactly the same reditR workflow.
  # This provides the observed number of significant sites for comparison
  # with the permutation-based null distribution.
  mp <- file.path(
    tmpdir,
    sprintf("meta_%s_TRUE.txt", case)
  )

  fwrite(
    data.table(sample = samples, condition = truth),
    mp,
    sep = "\t"
  )

  cat(sprintf("[%s TRUE labelling]\n", case))

  rt <- as.data.table(
    differential_editing(
      data_path = dp,
      meta_path = mp,
      test = c("glmm", "fisher", "wilcoxon"),
      reference_level = "scr",
      case_level = case,
      random_effects = "(1 | sample)",
      out_path = file.path(
        tmpdir,
        sprintf("DRE_%s_TRUE.txt", case)
      ),
      summary_path = file.path(
        tmpdir,
        sprintf("sum_%s_TRUE.txt", case)
      ),
      n_cores = N_CORES
    )
  )

  # Store the observed result as null = 0 so it can be separated from the
  # permuted results when calculating the empirical p-value.
  truth_row <- data.table(
    arm = case,
    null = 0L,
    labels = paste0(
      "TRUE:",
      paste(ifelse(truth == "scr", "R", "C"), collapse = "")
    ),
    GLMM_sig = sum(rt$GLMM_sig, na.rm = TRUE),
    Fisher_sig = sum(rt$Fisher_sig, na.rm = TRUE),
    Wilcox_sig = sum(rt$Wilcox_sig, na.rm = TRUE)
  )

  # Combine the true result with all permutation results and save them.
  out <- rbind(truth_row, rbindlist(rows))

  fwrite(
    out,
    sprintf("reditr_permutation_null_endo2_%s.txt", case),
    sep = "\t"
  )

  # Keep only the null permutations for calculating the empirical
  # distribution and comparing it with the true result.
  nl <- out[null > 0]

  cat(sprintf(
    "\n=== %s ===\nTRUE: GLMM %d | Fisher %d | Wilcox %d\nNULL (n=%d): GLMM mean %.1f range %d-%d | Fisher mean %.1f range %d-%d\nempirical p (GLMM) = %.3f | (Fisher) = %.3f\n",
    case,
    truth_row$GLMM_sig,
    truth_row$Fisher_sig,
    truth_row$Wilcox_sig,
    nrow(nl),
    mean(nl$GLMM_sig),
    min(nl$GLMM_sig),
    max(nl$GLMM_sig),
    mean(nl$Fisher_sig),
    min(nl$Fisher_sig),
    max(nl$Fisher_sig),
    (sum(nl$GLMM_sig >= truth_row$GLMM_sig) + 1) /
      (nrow(nl) + 1),
    (sum(nl$Fisher_sig >= truth_row$Fisher_sig) + 1) /
      (nrow(nl) + 1)
  ))

  out
}

# Run the same permutation analysis separately for siADAR1 and siADAR2.
all_out <- rbindlist(
  lapply(c("siADAR1", "siADAR2"), run_arm)
)

# Save the combined results from both experimental arms.
fwrite(
  all_out,
  "reditr_permutation_null_endo2.txt",
  sep = "\t"
)

cat("\nWritten: reditr_permutation_null_endo2.txt\n")