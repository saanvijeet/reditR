#' Run differential RNA editing analysis
#'
#' Identifies differentially RNA-edited (DRE) sites between two conditions
#' using one or more of three statistical tests: GLMM, Fisher's exact test,
#' and Wilcoxon rank-sum test. Each test is run independently on all sites
#' and receives its own Benjamini-Hochberg (BH) correction.
#'
#' @param data_path Path to the filtered and clustered editing-site count file.
#'   Required columns are \code{site}, \code{sample}, \code{edited}, and
#'   \code{total}. Editing ratios are calculated if not already provided.
#'
#' @param meta_path Path to the sample metadata file. Required columns are
#'   \code{sample} and \code{condition}. Not required if \code{condition}
#'   is already included in the input data.
#'
#' @param out_path Path for the output results file. Set to \code{NULL}
#'   to avoid writing a file.
#'
#' @param summary_path Path for a summary file containing the number of
#'   tested and significant sites. Set to \code{NULL} to avoid writing a file.
#'
#' @param test Statistical test(s) to run. Choose from \code{"glmm"},
#'   \code{"fisher"}, and \code{"wilcoxon"}. Multiple tests can be selected.
#'   The default runs all three.
#'
#' @param reference_level Name of the reference/control condition.
#'   Default is \code{"control"}.
#'
#' @param case_level Name of the case/treated condition.
#'   Default is \code{"diabetic"}.
#'
#' @param random_effects Random-effect term used in the GLMM.
#'   Defaults to \code{"(1 | sample)"}. For pseudobulk data, this can be
#'   changed to account for shared libraries or clusters.
#'
#' @param min_obs Minimum number of observations required to fit a GLMM
#'   for a site. Default is 4.
#'
#' @param fdr_threshold BH-adjusted p-value threshold used to identify
#'   significant sites. Default is 0.05.
#'
#' @param n_cores Number of CPU cores used for parallel GLMM fitting.
#'   Default is 1. Windows runs are always performed using one core.
#'
#' @param verbose Logical value indicating whether progress and summary
#'   information should be printed. Default is \code{TRUE}.
#'
#' @return A \code{data.table} containing the site-level p-values,
#'   BH-adjusted p-values, and significance results for each selected test.
#'   Results from each test are reported independently; no combined
#'   significance result is produced.
#'
#' @details
#' The three tests use different aspects of the editing data. GLMM models
#' edited and unedited read counts while accounting for between-sample
#' variability. Fisher's exact test compares pooled edited and unedited
#' counts between conditions. Wilcoxon compares editing ratios between
#' individual samples.
#'
#' The appropriate test depends on the experimental design. GLMM can be
#' useful when there are sufficient independent biological replicates and
#' can account for between-sample variability. Fisher's exact test is
#' simple and can be used for sparse sites, but pooling reads can ignore
#' biological variability. Wilcoxon does not depend on read-count
#' distributions but has lower power at small sample sizes.
#'
#' Because test performance can depend on sample size, replication, and
#' between-sample variability, users should assess calibration before
#' interpreting significant sites. reditR therefore provides separate
#' statistical tests rather than combining them into a single result.
#'
#' @examples
#' \dontrun{
#' res <- differential_editing(
#'   data_path = "filtered_sites_clustered.txt",
#'   meta_path = "sample_metadata.txt",
#'   test = c("glmm", "wilcoxon"),
#'   out_path = "DRE_results.txt"
#' )
#'
#' # View sites significant by GLMM
#' res[GLMM_sig == TRUE]
#'
#' # Find sites significant by both tests
#' res[GLMM_sig == TRUE & Wilcox_sig == TRUE]
#' }
#'
#' @export
differential_editing <- function(data_path,
                                  meta_path       = NULL,
                                  out_path        = NULL,
                                  summary_path    = NULL,
                                  test            = c("glmm", "fisher", "wilcoxon"),
                                  reference_level = "control",
                                  case_level      = "diabetic",
                                  random_effects  = "(1 | sample)",
                                  min_obs         = 4L,
                                  fdr_threshold   = 0.05,
                                  n_cores         = 1L,
                                  verbose         = TRUE) {

# Check that the requested tests are valid.
  # This prevents spelling errors or unsupported test names.
  test <- match.arg(test, choices = c("glmm", "fisher", "wilcoxon"), several.ok = TRUE)
 # Read the editing-site count data.
  # Each row represents an editing site measured in one sample.
  data <- fread(data_path)
   # Add condition information if it is not already in the input data.
  # Keeping condition in one table makes the later per-site analyses simpler.
  if (!"condition" %in% names(data)) {
    if (is.null(meta_path))
      stop("'condition' column not found in data file and no meta_path supplied.")
    meta <- fread(meta_path)
    data <- merge(data, meta, by = "sample")
  }

  # Keep only the two conditions being compared.
  # This prevents unrelated experimental groups from entering the analysis.
  data <- data[condition %in% c(reference_level, case_level)]

# Check that both conditions are present.
  # A comparison cannot be performed if one group is missing.
  absent <- setdiff(c(reference_level, case_level), unique(as.character(data$condition)))
  if (length(absent) > 0)
    stop("condition level(s) not present in data: ", paste(absent, collapse = ", "),
         ". Check reference_level / case_level against the condition column.")

  # Explicitly set the order of the condition factor.
  # This ensures that the GLMM uses the requested reference condition and
  # that the estimated coefficient represents the case-versus-reference
  # comparison.
  data[, condition := factor(as.character(condition),
                             levels = c(reference_level, case_level))]
  
  # Calculate unedited reads if they are not already provided.
  # GLMM and Fisher's exact test require both edited and unedited counts.
  if (!"unedited" %in% names(data))   data[, unedited   := total - edited]
  # Calculate the editing ratio if it is not already provided.
  # Wilcoxon uses this per-sample editing proportion.
  if (!"edit_ratio" %in% names(data)) data[, edit_ratio := edited / total]

  # Print basic information about the analysis.
  # This provides a simple record of the samples, conditions, sites, and tests
  # used in each run.
  if (verbose) {
    cat("Samples:   ", paste(sort(unique(data$sample)), collapse = " "), "\n")
    cat("Conditions:", reference_level, "=", sum(data$condition == reference_level),
        "obs |", case_level, "=", sum(data$condition == case_level), "obs\n")
    cat("Sites:     ", uniqueN(data$site), "\n")
    cat("Tests:     ", paste(test, collapse = ", "), "\n\n")
  }

  # Name of the GLMM coefficient representing the case condition.
  # This is the coefficient whose p-value is used to test for differential
  # editing between the two conditions.
  coef_name <- paste0("condition", case_level)

  # Build the GLMM formula only when GLMM has been requested.
  # Edited and unedited reads are modelled as a binomial response, with
  # condition as the fixed effect and the user-defined random effects
  # accounting for repeated or clustered observations.
  #User also has the option to add a clinical covariate, an example of this can 
  #found in the developer's notes.
  glmm_formula <- NULL
  if ("glmm" %in% test) {
    glmm_formula <- stats::as.formula(
      paste("cbind(edited, unedited) ~ condition +", random_effects)
    )
    # Check that every variable used in the random-effects term exists.
    # This gives a clear error before site-level model fitting begins.
    re_vars <- setdiff(all.vars(glmm_formula), c("edited", "unedited", "condition"))
    missing_re_vars <- setdiff(re_vars, names(data))
    if (length(missing_re_vars) > 0)
      stop("random_effects references column(s) not found in data: ",
           paste(missing_re_vars, collapse = ", "))
  }

  # ── Functions for the three statistical tests ──
   # GLMM:
  # Models the number of edited and unedited reads at each site.
  # The random-effects term accounts for variation between the biological
  # units specified by the user.
  run_glmm <- function(d) {
    # Skip sites with too few observations or only one condition.
    if (nrow(d) < min_obs || uniqueN(d$condition) < 2L) return(NULL)
    p <- tryCatch({
      # Fit the binomial GLMM.
      # BOBYQA is used as the optimiser and maxfun is increased to allow
      # more iterations for difficult site-level models.
      m <- withCallingHandlers(
        lme4::glmer(glmm_formula,
                    data = d, family = binomial,
                    control = lme4::glmerControl(optimizer = "bobyqa",
                                                 optCtrl = list(maxfun = 2e5))),
        # Suppress model-fitting warnings so that one difficult site does
        # not interrupt the full analysis.
        warning = function(w) invokeRestart("muffleWarning"))
        # Extract the p-value for the condition effect.
      summary(m)$coefficients[coef_name, "Pr(>|z|)"]
    }, error = function(e) NA_real_)
    # If a model cannot be fitted, return NA rather than stopping the
      # analysis for all remaining sites.
    data.table(glmm_pvalue = p)
  }
   # Fisher's exact test:
  # Pools edited and unedited reads within each condition and tests whether
  # the overall editing counts differ between conditions.
  run_fisher <- function(d) {
    # Both conditions are required for a comparison.
    if (uniqueN(d$condition) < 2L) return(NULL)
    # Create a 2 x 2 table:
    #             edited     unedited
    # reference   ...        ...
    # case        ...        ...
    #
    # Counts are summed across samples within each condition.
    mat <- matrix(c(
      sum(d[condition == reference_level, edited]),
      sum(d[condition == case_level,      edited]),
      sum(d[condition == reference_level, unedited]),
      sum(d[condition == case_level,      unedited])
    ), nrow = 2, dimnames = list(c(reference_level, case_level), c("edited", "unedited")))
    # Fisher's exact test is used because it provides an exact p-value from
    # the observed count table without relying on large-sample assumptions.
    p <- tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
    data.table(fisher_pvalue = p)
  }
  # Wilcoxon rank-sum test:
  # Compares editing ratios between individual samples rather than pooling reads. 
  # This provides a non-parametric test based on biological replicates.
  run_wilcox <- function(d) {
    if (uniqueN(d$condition) < 2L) return(NULL)
    p <- tryCatch(
      wilcox.test(edit_ratio ~ condition, data = d, exact = FALSE)$p.value,
      error = function(e) NA_real_)
    data.table(wilcox_pvalue = p)
  }
  # -------------------------------------------------------------------------
  # Run each test separately for each editing site
  # -------------------------------------------------------------------------

  # Split the dataset into one table per editing site.
  # This allows each statistical test to be performed independently at
  # every site.
  data_list <- split(data, by = "site")
  # Windows does not support fork-based parallelisation used by mclapply,
# so run serially; otherwise use the requested number of cores.
  cores <- if (.Platform$OS.type == "windows") 1L else n_cores
  # Apply one statistical test to every site.
  # NULL results from sites that cannot be tested are removed before combining
  # the results into one table.
  run_one_test <- function(fn) {
    rbindlist(Filter(Negate(is.null), parallel::mclapply(data_list, fn, mc.cores = cores)),
              idcol = "site")
  }
  # Start the results table with all tested sites.
  # Individual test results are added to this table as each analysis finishes.
  results <- data.table(site = names(data_list))

  #GLMM analysis
  if ("glmm" %in% test) {
    if (verbose) cat("[", format(Sys.time()), "] Running GLMM on", length(data_list), "sites...\n")
    t0 <- Sys.time()
    stage <- run_one_test(run_glmm)
    # Correct all GLMM p-values using Benjamini-Hochberg.
    # This controls the expected false discovery rate across the tested sites.
    stage[, GLMM_FDR := p.adjust(glmm_pvalue, method = "BH")]
    # Mark sites passing the chosen FDR threshold.
    stage[, GLMM_sig := !is.na(GLMM_FDR) & GLMM_FDR < fdr_threshold]
    # Add GLMM results to the main results table.
    results <- merge(results, stage, by = "site", all.x = TRUE)
    if (verbose) cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 2),
                     "min. FDR<", fdr_threshold, ":", sum(stage$GLMM_sig, na.rm = TRUE), "\n\n")
  }
   # Fisher's exact test
  if ("fisher" %in% test) {
    if (verbose) cat("[", format(Sys.time()), "] Running Fisher's exact on", length(data_list), "sites...\n")
    t0 <- Sys.time()
    stage <- run_one_test(run_fisher)
    # Apply BH correction independently of the other tests.
    stage[, Fisher_FDR := p.adjust(fisher_pvalue, method = "BH")]
    # Identify sites below the selected FDR threshold.
    stage[, Fisher_sig := !is.na(Fisher_FDR) & Fisher_FDR < fdr_threshold]
    results <- merge(results, stage, by = "site", all.x = TRUE)
    if (verbose) cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 2),
                     "min. FDR<", fdr_threshold, ":", sum(stage$Fisher_sig, na.rm = TRUE), "\n\n")
  }
    # Wilcoxon rank-sum test
  if ("wilcoxon" %in% test) {
    if (verbose) cat("[", format(Sys.time()), "] Running Wilcoxon on", length(data_list), "sites...\n")
    t0 <- Sys.time()
    stage <- run_one_test(run_wilcox)
    stage[, Wilcox_FDR := p.adjust(wilcox_pvalue, method = "BH")]
    stage[, Wilcox_sig := !is.na(Wilcox_FDR) & Wilcox_FDR < fdr_threshold]
    # Add Wilcoxon results to the main results table.
    results <- merge(results, stage, by = "site", all.x = TRUE)
    if (verbose) cat("  Done in", round(difftime(Sys.time(), t0, units = "mins"), 2),
                     "min. FDR<", fdr_threshold, ":", sum(stage$Wilcox_sig, na.rm = TRUE), "\n\n")
  }
   # Save results
   # Write the complete site-level results if an output path was provided.
  if (!is.null(out_path)) fwrite(results, out_path, sep = "\t")
  # Create a simple summary showing the total number of sites tested and
  # the number significant for each requested test.
  sig_cols <- c(glmm = "GLMM_sig", fisher = "Fisher_sig", wilcoxon = "Wilcox_sig")[test]
  summary_dt <- data.table(
    metric = c("sites_tested", unname(sig_cols)),
    value  = c(nrow(results), vapply(sig_cols, function(cn) sum(results[[cn]], na.rm = TRUE), numeric(1)))
  )
  # Save the summary if requested.
  if (!is.null(summary_path)) fwrite(summary_dt, summary_path, sep = "\t")
  if (verbose) { cat("Summary:\n"); print(summary_dt) }

  results
}
